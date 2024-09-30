
import os
import argparse
from scipy.constants import epsilon_0
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kv as K0, iv as I0
from scipy.linalg import solve, pinv, svd
from EMGinv_fns import *
from geometry_utilities import *
from plotter_utility_functions import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull, Delaunay
from scipy.optimize import minimize
import pandas as pd
import scipy.io
import mat73
from pathlib import Path
import shutil
import datetime
import mne
from tqdm import tqdm
import gc
import pydicom
import pyvistaqt

# Set the folder for batch processing
BATCH_FOLDER = 'Data/Batch'
NUM_ELECTRODES = 256
MIN_DIPOLE_DISTANCE = 2.0 # mm
SAMPLE_RATE = 2000. # Hz
FIXED_TIME_VALUE = 4 # index
ANIMATION_EXPORT_ROOT = 'C:/Data/MetaWB'
EXPORT_FOLDER = 'Data/Batch/Exported'
SUMMARY_FOLDER = 'Data/Summaries'
ANIMATED_SLICE_DEPTHS = [50] # Which depths to export animation frames for.
ARM = 'right'
EXPORT_FRAMES = False

# Process each file in the folder
def run_emg_pipeline_batch(folder=BATCH_FOLDER, export_folder=EXPORT_FOLDER, summary_folder=SUMMARY_FOLDER, t_value=FIXED_TIME_VALUE, arm=ARM, export_frames=EXPORT_FRAMES):
    # Prepare the summary folder
    os.makedirs(summary_folder, exist_ok=True)
    # Initialize an empty list to store summary data
    summary_data = []

    MNE_FWD_SOLUTIONS_FILE = f'Data/simp_arm_{arm}_{MIN_DIPOLE_DISTANCE}mm_{NUM_ELECTRODES}channel_{SAMPLE_RATE}Hz-fwd.fif'
    ELECTRODES_FILE = f'Data/{NUM_ELECTRODES}simparm-{arm}_electrode_pos.npy'

    # Load the necessary data
    electrode_pos = np.load(ELECTRODES_FILE)
    mne_fwd = mne.read_forward_solution(MNE_FWD_SOLUTIONS_FILE, verbose='WARNING')
    fwd = mne_fwd['sol']['data']
    pos = mne_fwd['source_rr']

    # Process forward model and bone removal
    xscaling, yscaling, zscaling = (MIN_DIPOLE_DISTANCE*1e-3, MIN_DIPOLE_DISTANCE*1e-3, MIN_DIPOLE_DISTANCE*1e-3)
    pos, fwd = bone_remover(pos, fwd, -9e-3, 0, 5e-3)  # Ulnar
    pos, fwd = bone_remover(pos, fwd, 5e-3, -9e-3, 5e-3)  # Radius
    # Condense the forward model
    dipole_ori = [0, 0, 1]
    fwd = fwd_convertfixed(fwd, dipole_ori)

    # Get a list of all files and directories
    entries = os.listdir(folder)
    # Filter out only the files
    files = [entry for entry in entries if os.path.isfile(os.path.join(folder, entry))]
    for file in tqdm(files, desc="Localizing MUAPs", unit="File"):
        if file.endswith('.mat') and ('muaps_template-' in file) and (arm in file):
            file_path = os.path.join(folder, file)
            
            # Set the file naming convention based on the file
            file_stem = Path(file_path).stem
            parts = file_stem.split('_')

            # Parse TANK (everything up to the 4th underscore)
            TANK = '_'.join(parts[:4])  # Example: 'MCP01_2024_04_12'

            # Parse BLOCK (the value just before 'muaps')
            BLOCK = int(parts[-5])  # Example: '15'

            # Parse MUAP_ID (the integer right after 'muaps_template-')
            MUAP_ID = int(file_stem.split('muaps_template-')[1].split('_')[0])  # Example: '02'

            EXPERIMENT = f'{TANK}_{BLOCK}'
            TEMPLATES_FILE = file_path
            COVARIANCE_FILE = f'Data/Batch/{EXPERIMENT}_muaps_covariance.mat'                

            waveform = scipy.io.loadmat(TEMPLATES_FILE)['muaps']
            try:
                noise_cov = scipy.io.loadmat(COVARIANCE_FILE)['noise_cov']
                data_cov = scipy.io.loadmat(COVARIANCE_FILE)['data_cov']
                if np.isnan(data_cov).any() or np.isinf(data_cov).any():
                    print(f"Invalid values found in data_cov for {EXPERIMENT}, skipping this iteration.")
                    continue  # Skip this iteration if data_cov contains NaNs or infs
            except FileNotFoundError:
                print(f"Covariance file not found for {EXPERIMENT}, skipping this iteration.")
                continue
            except Exception as e:
                print(f"Unexpected error while loading covariance file for {EXPERIMENT}: {str(e)}")
                continue  # Skip if there's any unexpected error

            w_lcmv = lcmv_beamformer_constructor(fwd, data_cov=data_cov.data, noise_cov= noise_cov.data, pos=pos, arr_gain=True, max_power=False)
            source_activity_time = np.dot(w_lcmv, waveform)
            save_dir = f'{ANIMATION_EXPORT_ROOT}/{TANK}/Animations/BEM/{EXPERIMENT}/MU{MUAP_ID:02d}'
            os.makedirs(save_dir, exist_ok=True)
            fig = plot_waveform_grid(waveform.T, fs=SAMPLE_RATE, cols=NUM_ELECTRODES/8, row_spacing=50.0, filename=f'{save_dir}/{EXPERIMENT}_Surface-Timeseries.png')
            plt.close(fig)

            # Threshold
            thresh = 0.5
            source_activity = source_activity_time[:, t_value]

            # Reshape source activity to condense N source orientations into 1 per voxel - Confirmed works for 3 orientations, should work for more.
            reshape_by = source_activity.shape[0] // pos.shape[0]
            reshaped_act = np.array(source_activity.reshape((reshape_by, -1), order='F'))
            source_activity = np.linalg.norm(reshaped_act, axis=0)

            ind = np.abs(source_activity) > thresh*np.max(np.abs(source_activity))
            source_activity = source_activity[ind]
            pos_t = pos[ind]

            # Plot the convex hull and the moved points
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            # Plot electrode positions
            ax.scatter(electrode_pos[:, 0], electrode_pos[:, 1], electrode_pos[:, 2], c=waveform[:,t_value], marker='o', alpha=0.2, cmap='turbo')
            # Plot the source space
            ax.scatter(pos_t[:, 0], pos_t[:, 1], pos_t[:, 2], c=source_activity, marker='s', alpha=0.8, cmap='viridis')
            # Set labels
            ax.set_xlabel('X Axis (m)')
            ax.set_ylabel('Y Axis (m)')
            ax.set_zlabel('Z Axis (m)')
            n_channels, n_sources = fwd.shape
            ax.set_title(f'{n_channels}-Channel {n_sources}-Dipole LCMV Beamformer')

            # Save the figure to a file in the created folder
            fig.savefig(f'{save_dir}/{EXPERIMENT}_3D-T{t_value}.png', dpi=300)
            plt.close(fig)

            x_min, x_max = np.min(pos[:, 0]), np.max(pos[:, 0])
            y_min, y_max = np.min(pos[:, 1]), np.max(pos[:, 1])
            z_min, z_max = np.min(pos[:, 2]), np.max(pos[:, 2])
            grid = pos_to_3Dgrid_converter(pos, source_activity_time[:,t_value], (xscaling, yscaling, zscaling))
            # # # # Determine the nearest muscle based on highest activity # # # #
            options = {
                'LandMarksFile': "C:/Data/Anatomy/Human Arm/Sections/Landmarks.xlsx",
                'LandMarkSheetExpression': "R_Forearm_Section_%d",
                'GridRows': grid.shape[0],
                'GridCols': grid.shape[1],
                'GridRowTopOffset': 100,
                'GridRowBottomOffset': 500,
                'GridColumnLeftOffset': 330,
                'GridColumnRightOffset': 670,
                'MillimetersPerSection': 8.5,
                'StartingSection': 105,
                'XOffset': -0.020,
                'YOffset': -0.020,
                'XScale': 0.003,
                'YScale': -0.002
            }
            
            nearest_muscle, nearest_muscle_distance = find_nearest_muscle(source_activity_time, pos, t=t_value, options=options)
            summary_data.append({
                'TANK': TANK,
                'BLOCK': BLOCK,
                'MUAP_ID': MUAP_ID,
                'Nearest_Muscle': nearest_muscle,
                'Distance': nearest_muscle_distance,
                'Time_Point': t_value
            })

            # Create a unique folder based on the current runtime
            if export_frames:
                for z in ANIMATED_SLICE_DEPTHS:
                    save_dir = f'{ANIMATION_EXPORT_ROOT}/{TANK}/Animations/BEM/{EXPERIMENT}/MU{MUAP_ID:02d}/z{z}mm'
                    os.makedirs(save_dir, exist_ok=True)  # Create the folder if it doesn't exist
                    for t in tqdm(range(source_activity_time.shape[1]),desc=f'Exporting frames for z={z}-mm',unit='Frame'):
                        source_activity = source_activity_time[:, t]

                        # Faster to reconfigure the scatter points to be in a grid, and then use imshow to plot the activity.
                        grid = pos_to_3Dgrid_converter(pos, source_activity, (xscaling, yscaling, zscaling))

                        # Create a figure
                        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

                        # MRI visualization on the left subplot
                        # img = axs[0].imshow(arm_image[z * arm_image.shape[0] // grid.shape[2], 100:500, 330:670, 2], 
                                            # extent=[y_min, y_max, x_min, x_max], cmap='gray', origin='upper')
                        scatter_left = axs[0].scatter(electrode_pos[:, 0], electrode_pos[:, 1], c=waveform[:, t], cmap='turbo', marker='o')
                        axs[0].set_xticks([])
                        axs[0].set_yticks([])

                        # Turn off the box (spines) around the plot
                        for spine in axs[0].spines.values():
                            spine.set_visible(False)

                        # Source estimate on the right subplot
                        z_idx = z_to_grid(z*1e-3, z_min, zscaling)
                        img_source = axs[1].imshow(grid[:, :, z_idx], origin='upper', cmap='viridis', extent=[y_min, y_max, x_min, x_max],
                                                vmin=0, vmax=150)  # Clamping the colorbar scale
                        scatter_right = axs[1].scatter(electrode_pos[:, 0], electrode_pos[:, 1], c=waveform[:, t], cmap='turbo', marker='o')
                        axs[1].set_title('Source Estimate')
                        axs[1].set_xlabel('X Axis (m)')
                        axs[1].set_ylabel('Y Axis (m)')
                        plt.colorbar(img_source, ax=axs[1], label='Source Activity')  # Associate the colorbar with the imshow (source estimate)

                        # Example call to annotate the slice with landmarks
                        options = {
                            'FileExpression': "C:/Data/Anatomy/Human Arm/Sections/R_Forearm_Section_%d.png",
                            'LandMarksFile': "C:/Data/Anatomy/Human Arm/Sections/Landmarks.xlsx",
                            'LandMarkSheetExpression': "R_Forearm_Section_%d",
                            'LandmarksToAdd': ["ALL"],
                            'GridColumns': grid.shape[1],
                            'GridRows': grid.shape[0],
                            'GridColumnLeftOffset': 100,
                            'GridColumnRightOffset': 500,
                            'GridRowTopOffset': 330,
                            'GridRowBottomOffset': 670,
                            'XScale': 0.003,
                            'YScale': -0.002,
                            'XOffset': -0.0575,
                            'YOffset': -0.005,
                            'MillimetersPerSection': 8.5, 
                            'AddLabels': True,
                            'AddImage': True
                        }

                        plot_annotated_template_slice(axs[0], z, options)
                        options['AddImage'] = False
                        # Plot the annotated landmarks on the right subplot (source estimate plot)
                        plot_annotated_template_slice(axs[1], z, options)
                        # Adjust the layout
                        plt.tight_layout()
                        # Save the figure to a file in the created folder
                        plt.savefig(f'{save_dir}/frame_{t:03d}.png', dpi=300)
                        # Close the figure to save memory
                        plt.close(fig)
            new_file_path = os.path.join(export_folder, file)
            shutil.move(file_path, new_file_path)

    # Convert the summary data to a DataFrame and save it as a CSV file
    summary_df = pd.DataFrame(summary_data)
    summary_filename = f'{summary_folder}/summary_t{t_value}.csv'
    summary_df.to_csv(summary_filename, index=False)
    print("Batch processing complete.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run EMG pipeline batch processing.')
    parser.add_argument('--folder', '-f', type=str, default=BATCH_FOLDER, help='The folder containing the files to process')
    parser.add_argument('--export-folder', '-e', type=str, default=EXPORT_FOLDER, help='The folder where processed files will be moved')
    parser.add_argument('--summary-folder', '-s', type=str, default='Data/Summaries', help='The folder to save the summary table')
    parser.add_argument('--t-value', '-t', type=int, default=4, help='Time point (t) to analyze')
    parser.add_argument('--arm', '-a', type=str, default='right', help='Which arm to analyze')
    parser.add_argument('--animate', '-i', type=int, default=0, help='Export animation frames?')

    args = parser.parse_args()

    # Call the processing function with the folder from the command line
    run_emg_pipeline_batch(args.folder, args.export_folder, args.summary_folder, args.t_value, args.arm, args.animate > 0)