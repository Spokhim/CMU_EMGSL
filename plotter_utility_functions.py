import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

def plot_annotated_template_slice(ax, z, options):
    """
    Plot annotated slice on a provided axis with labels and an image.

    Parameters:
    - ax (matplotlib axis): The axis on which to plot.
    - z (int): The slice index to plot.
    - options (dict): A dictionary containing various options:
        - FileExpression (str): File path pattern for the slice image (default provided).
        - LandMarksFile (str): File path to the Excel file with landmarks (default provided).
        - LandMarkSheetExpression (str): Excel sheet name pattern for landmarks (default provided).
        - LandmarksToAdd (list): List of landmarks to annotate. Use ["ALL"] for all landmarks.
        - AddLabels (bool): Flag to add text labels (default True).
        - AddImage (bool): Flag to add the image (default True).
    """
    # Set default options if not provided
    options.setdefault('FileExpression', "C:/Data/Anatomy/Human Arm/Sections/R_Forearm_Section_%d.png")
    options.setdefault('LandMarksFile', "C:/Data/Anatomy/Human Arm/Sections/Landmarks.xlsx")
    options.setdefault('LandMarkSheetExpression', "R_Forearm_Section_%d")
    options.setdefault('LandmarksToAdd', ["ALL"])
    options.setdefault('GridRows', 19)
    options.setdefault('GridCols', 19)
    options.setdefault('GridRowTopOffset', 100)
    options.setdefault('GridRowBottomOffset', 500)
    options.setdefault('GridColumnLeftOffset', 330)
    options.setdefault('GridColumnRightOffset', 670)
    options.setdefault('MillimetersPerSection', 8.5)
    options.setdefault('StartingSection', 105)
    options.setdefault('XOffset', -0.020)
    options.setdefault('YOffset', -0.020)
    options.setdefault('AddLabels', True)
    options.setdefault('AddImage', True)

    z_index = round(z / options['MillimetersPerSection']) + options['StartingSection']
    # Load the image
    img_file = options['FileExpression'] % z_index
    img = plt.imread(img_file)

    # Load the landmarks from the Excel file
    sheet_name = options['LandMarkSheetExpression'] % z_index
    landmarks_df = pd.read_excel(options['LandMarksFile'], sheet_name=sheet_name)
    grid_row_scale = options['GridRows'] / (options['GridRowBottomOffset']-options['GridRowTopOffset'])
    grid_col_scale = options['GridCols'] / (options['GridColumnRightOffset']-options['GridColumnLeftOffset'])

    # Plot the image if specified
    if options['AddImage']:
        ax.imshow(img, extent=[0, img.shape[1], img.shape[0], 0], cmap='gray')

    # Plot the landmarks if specified
    if options['AddLabels']:
        if "ALL" in options['LandmarksToAdd']:
            for i, row in landmarks_df.iterrows():
                x = options['XScale']*((row['X']-options['GridColumnLeftOffset']) * grid_col_scale) + options['XOffset']
                y = options['YScale']*((row['Y']-options['GridRowTopOffset']) * grid_row_scale) + options['YOffset']
                # print(f"x: {x} | y: {y}")
                ax.text(x, y, row['Landmark'], fontsize=9, color='white',
                        ha='center', va='center', fontweight='bold', fontname='Consolas')
        else:
            for landmark in options['LandmarksToAdd']:
                row = landmarks_df[landmarks_df['Landmark'].str.upper() == landmark.upper()]
                if not row.empty:
                    x = options['XScale']*((row.iloc[0]['X']-options['GridColumnLeftOffset'])) * grid_col_scale + options['XOffset']
                    y = options['YScale']*((row.iloc[0]['Y']-options['GridRowTopOffset'])) * grid_row_scale + options['YOffset']
                    ax.text(x, y, row.iloc[0]['Landmark'], 
                            fontsize=9, color='white', ha='center', va='center', 
                            fontweight='bold', fontname='Consolas')


def plot_waveform_grid(waveform, rows=8, cols=16, fs=4000, row_spacing=50, vertical_scale_units='μV', colors=["#EF3A47", "#008F91", "#FDB515", "#043673"], filename=None):
    """
    Plot a set of time-series snippets arranged on the same axis but as a grid of signals (8 rows, 16 columns).
    Each group of 4 adjacent columns has the same color, in the order specified by 'colors'.
    
    Parameters:
    - waveform: 2D numpy array of shape (n_samples, 128), where n_samples is the number of time samples, and 128 is the number of channels.
    - rows: Number of rows in the grid (default 8).
    - cols: Number of columns in the grid (default 16).
    - fs: Sampling frequency (default 4000 Hz).
    - row_spacing: Vertical distance between rows (default 50).
    - vertical_scale_units: Units for vertical scale (default 'μV').
    - colors: List of 4 colors to use for groups of 4 columns (default provided).
    """
    n_samples, n_channels = waveform.shape
    assert n_channels == rows * cols, "Waveform should have 128 channels for 8x16 grid."

    # Create a figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Set up vertical and horizontal offsets to arrange the grid
    v_offset = row_spacing  # Vertical distance between rows
    h_offset = 1.05 * n_samples / fs  # Horizontal distance between columns

    for channel_idx in range(n_channels):
        col = math.floor(channel_idx // rows)
        row = rows - (channel_idx % rows) - 1

        # Determine the color based on the group of 4 columns
        color = colors[(col // 4) % len(colors)]
        
        # Plot each snippet with offsets for vertical and horizontal arrangement
        ax.plot(np.arange(n_samples) / fs + col * h_offset,  # Time scaling based on fs
                waveform[:, channel_idx] + row * v_offset,  # Vertical offset for each row
                color=color, lw=1.5)

    # Remove axis ticks for a cleaner look
    ax.set_xticks([])
    ax.set_yticks([])

    # Turn off the box (spines) around the plot
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    # Add scalebar
    total_time_sec = n_samples / fs  # Total time in seconds
    time_scalebar_width = h_offset
    time_scalebar_label = f"{time_scalebar_width:.1f}s" if h_offset >= 0.1 else f"{time_scalebar_width * 1000:.1f}ms"
    
    # Draw the horizontal time scalebar (bottom left)
    ax.plot([-0.5 * time_scalebar_width, 0.5 * time_scalebar_width], [-row_spacing, -row_spacing], color='black', lw=2)
    ax.text(time_scalebar_width / 2 -0.5 * h_offset, -row_spacing * 1.25, time_scalebar_label, ha='center', va='center', fontsize=10, color='black')

    # Vertical scale for 2 * row_spacing
    vertical_scalebar_height = 2 * row_spacing
    ax.plot([-0.5 * h_offset, -0.5 * h_offset], [-row_spacing, vertical_scalebar_height - row_spacing], color='black', lw=2)
    ax.text(-0.85 * total_time_sec, vertical_scalebar_height / 2 - row_spacing, f"{vertical_scalebar_height:.1f} {vertical_scale_units}",
            ha='center', va='center', rotation=90, fontsize=10, color='black')

    # Set the limits to adjust for the scalebar
    ax.set_xlim([-1.25 * h_offset, total_time_sec + (cols - 1) * h_offset])
    ax.set_ylim([-row_spacing*1.5, rows * row_spacing])

    if filename is not None:
        fig.savefig(filename)
        print(f'Figure saved in {filename}')