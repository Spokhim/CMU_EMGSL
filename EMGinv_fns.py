from scipy.constants import epsilon_0
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kv as K0, iv as I0
from scipy.linalg import solve, pinv, svd
from joblib import Parallel, delayed

import scipy.io
import numpy as np
import os
import pydicom
import cv2 # opencv
import sklearn

from tmsi_python_interface.TMSiFileFormats.file_readers import Poly5Reader, Xdf_Reader, Edf_Reader

import sys
from os.path import join, dirname, realpath
import tkinter as tk
from tkinter import filedialog
import mne 

def dipole_potential(r_vec, p, r0):
    """If in cylindrical co-ordinates, make sure r_vec[1] = 0.

    Parameters:
    - r_vec: position where potential is calculated
    - p: dipole moment vector
    - r0: position of the dipole
    
    Returns:
    - potential: potential at the position r_vec due to the dipole moment p at position
    """
    R = r_vec - r0
    R_mag = np.linalg.norm(R)
    return (1 / (4 * np.pi * epsilon_0)) * (np.dot(p, R)) / (R_mag**3)

def point_charge_potential(r_vec, q, r0):
    """If in cylindrical co-ordinates, make sure r_vec[1] = 0.
    
    Parameters:
    - r_vec: position where potential is calculated
    - q: charge
    - r0: position of the dipole
    
    Returns:
    - potential: potential at the position r_vec due to the dipole moment p at position
    """
    R = r_vec - r0
    R_mag = np.linalg.norm(R)
    return (1 / (4 * np.pi * epsilon_0)) * q / R_mag

def fwd_generator(potential_func, pos, electrode_positions):
    """ This function generates a forward model projecting from the source space in the arm to electrodes.
    This will also work for generating the leadfield matrix for a subset of sources to the electrodes.

    The source space is a 3D grid of voxels with 3 orthogonal dipoles (charge configurations) at each voxel.  The electrodes are placed on the surface of the arm.  
    The forward model is a matrix that maps the dipole moments in the source space to the potentials at the electrodes.  
    The forward model is generated by calculating the potential at each electrode due to a unit dipole moment at each voxel in the source space.  
    For now, the potential is calculated using the formula for the potential due to a dipole moment in free space.  The forward model is then the matrix of potentials at the electrodes due to unit dipole moments at each voxel in the source space.  

    Parameters:
    - potential_func (function): Analytical function that calculates the potential at a point due to a charge configuration at another point.
    - pos (array): Position of all the source space voxels in 3D space (SI units).
    - electrode_positions (array): Position of all the electrodes in 3D space (SI units).
    
    Returns:
    - fwd (array): Forward model matrix of size N x 3V where N is the number of electrodes and V is the number of voxels.  The 3 is for the x, y, z components of the dipole moment.
    """

    # Number of electrodes
    N = len(electrode_positions)

    # Number of voxels in the source space
    V = pos.shape[0]

    # Initialize the forward model matrix
    fwd = np.zeros((N, 3*V))

    # Define unit vectors
    e_x = np.array([1, 0, 0])
    e_y = np.array([0, 1, 0])
    e_z = np.array([0, 0, 1])

    # Function to compute potentials for a single electrode
    def compute_potentials(electrode_pos):
        potentials = np.zeros(3 * V)
        for j in range(V):
            r0 = pos[j, :]
            potentials[3*j] = potential_func(electrode_pos, e_x, r0)
            potentials[3*j+1] = potential_func(electrode_pos, e_y, r0)
            potentials[3*j+2] = potential_func(electrode_pos, e_z, r0)
        return potentials

    # Parallel processing for each electrode
    results = Parallel(n_jobs=-1)(delayed(compute_potentials)(electrode_pos) for electrode_pos in electrode_positions)

    # Combine results into the forward model matrix
    fwd = np.array(results)

    return fwd

def load_src_template(filename=None, con='muscle', flip_dim=None, xscaling=1.5e-4, yscaling=1.5e-4, zscaling=0.5e-2):
    """Load the source space template from the MAT file.  This will be the positions of the voxels in the source space.
    
    Parameters: 
    - filename (str): the path to the file to be loaded.
    - con (str): the condition to load the source space under.  Options are 'muscle' or 'arm'.
    - flip_dim (int): the dimension to flip the source space along if source orientation is incorrect. 
                    Flipping dimension 1 (y) is useful when flipping between left and right arm orientation.
    - xscaling, yscaling, zscaling (float): the scaling factors for the x, y, z axes respectively to turn into SI units (m).

    Returns:
    - pos (array): Positions of the source space voxels (n_voxels x 3).
    """

    if filename is None:
        filename = 'Data/R_Forearm.mat'
    
    labels = scipy.io.loadmat(filename)['labels']
    print(labels.shape)

    if con == 'muscle':
        # Obtain muscle positions
        cond = labels==5
    elif con == 'arm':
        # Obtain all positions within arm
        cond = labels != 0

    # Generate the positions of each of the voxels which remain in the restricted forward matrix.  Transpose so that the columns are the x, y, z positions
    pos = np.array(np.where(cond), dtype='float').T
    print(pos.shape)
    print(pos.max(axis=1))
    print(pos.min(axis=1))

    # # Right now pos is just indexes... not SI units. Do a conversion to SI units (i.e. m)
    pos[:, 0] = pos[:, 0] * xscaling
    pos[:, 1] = pos[:, 1] * yscaling
    # # For z axis, different scaling - Let's say it's ~13cm
    pos[:, 2] = pos[:, 2] * zscaling

    # Flip dimensions as appropriate.
    if flip_dim is not None:
        pos[:, flip_dim] = 2*pos.mean(axis=0)[flip_dim] - pos[:, flip_dim]

    return pos

def pos_to_3Dgrid_converter(pos, source_activity, scaling):
    """ This function converts the positions of the source space to indices of a 3D grid and uses the source_activity as the value at each voxel.
    This is useful for visualising the source space. 
    
    Parameters:
    - pos (array): Positions of the source space voxels (n_voxels x 3).
    - source_activity (array): Source activity (n_voxels*3 x 1) for each dipole orientation.
    - scaling (tuple): Tuple containing the (x, y, z) scaling factors to convert the positions to indices of the 3D grid.

    Returns:
    - grid (array): 3D grid of the source space with the source activity values at each voxel.
    """
    # Define grid dimensions
    x_min, x_max = np.min(pos[:, 0]), np.max(pos[:, 0])
    y_min, y_max = np.min(pos[:, 1]), np.max(pos[:, 1])
    z_min, z_max = np.min(pos[:, 2]), np.max(pos[:, 2])

    xscaling, yscaling, zscaling = scaling

    # Define the resolution of the grid
    x_res, y_res, z_res = int((x_max-x_min)/xscaling), int((y_max-y_min)/yscaling), int((z_max-z_min)/zscaling)

    # Create an empty grid
    grid = np.zeros((x_res, y_res, z_res))
    grid.fill(np.nan)

    # Map positions to grid indices
    x_indices = ((pos[:, 0] - x_min) / (x_max - x_min) * (x_res - 1)).astype(int)
    y_indices = ((pos[:, 1] - y_min) / (y_max - y_min) * (y_res - 1)).astype(int)
    z_indices = ((pos[:, 2] - z_min) / (z_max - z_min) * (z_res - 1)).astype(int)

    # Assign activity values to the grid
    # Needed because 3 dipole moments per voxel 
    reshaped_act = np.array(source_activity.reshape((3, -1), order='F'))
    slice_act = np.linalg.norm(reshaped_act, axis=0)
    for i in range(len(slice_act)):
        grid[x_indices[i], y_indices[i], z_indices[i]] = slice_act[i]

    return grid

def load_tmsi_data(filename=None, return_mne_object=False):
    """Loading tmsi EMG data from a file.  Taken from 'tmsi_python_interface/examples_reading_data/example_file_reader.py'.

    Parameters:
    - filename (str): the path to the file to be loaded.
    - return_mne_object (bool): whether to return the mne object or not

    Returns: 
    - mne_object 
    /or/
    - samples, ch_names, sample_rate, num_channels
    """


    if filename is None:
        filename = 'Data/Pok_2024_08_21_A_PROX_8.poly5'

    try:
        if filename.lower().endswith('poly5'):
            reader = Poly5Reader(filename)
            mne_object = reader.read_data_MNE()
            # The mne_object reader incorrectly labels the channels as EEG so we can fix that
            # ch_types = np.array(mne_object.get_channel_types())
            # ch_types[ch_types == 'eeg'] = 'emg'
            # mne_object.set_channel_types(dict(zip(mne_object.ch_names, ch_types )))
            # However, lots of MNE-python code is gated for EEG channel types so we will just leave it as is.

            # Extract the samples and channel names from the Poly5Reader object
            samples = reader.samples
            ch_names = reader.ch_names
            sample_rate = reader.sample_rate
            num_channels = reader.num_channels

        elif filename.lower().endswith('xdf'):
            reader = Xdf_Reader(filename)
            data = reader.data[0]
            
            samples = data.get_data(units = {'eeg':'uV'})
            ch_names = data.ch_names
            sample_rate = data.info['sfreq']
            num_channels = len(ch_names)

        elif filename.lower().endswith('edf'):
            reader = Edf_Reader(filename)
            mne_object = reader.read_data_MNE()
            
            samples = data.get_data()
            ch_names = data.ch_names
            sample_rate = data.info['sfreq']
            num_channels = len(ch_names)

        elif not filename:
            tk.messagebox.showerror(title='No file selected', message = 'No data file selected.')

        else:
            tk.messagebox.showerror(title='Could not open file', message = 'File format not supported. Could not open file.')
        
        # Print retrieved data from the files
        print('Sample rate: ', sample_rate, ' Hz')
        print('Channel names: ', ch_names)
        print('Shape samples: ', np.shape(samples))

    except:
        tk.messagebox.showerror(title='Could not open file', message = 'Something went wrong. Could not open file.')        

    if return_mne_object:
        return mne_object
    else:
        return samples, ch_names, sample_rate, num_channels

def load_tmsitomne_combine(f_prox=None, f_dist=None):
    """ Combines the distal and proximal datasets from TMSI together into one MNE object.
    
    Parameters: 
    - f_prox (str): the path to the proximal file to be loaded.
    - f_dist (str): the path to the distal file to be loaded.

    Returns:
    - MNE_raw (MNE raw object): the combined MNE object.
    """
    if f_prox is None:
        f_prox = 'Data/Pok_2024_08_21_A_PROX_8.poly5'
    if f_dist is None:
        f_dist = 'Data/Pok_2024_08_21_B_DIST_8.poly5' 

    prox = load_tmsi_data(filename=f_prox)
    dist = load_tmsi_data(filename=f_dist)      
    ch_names = ["Prox - " + s for s in prox[1]] + ["Dist - " + s for s in dist[1]]
    fs = prox[2]
    ch_types = ['eeg']*64 + ['misc']*(len(prox[1])-64) + ['eeg']*64 + ['misc']*(len(dist[1])-64)
    num_channels = len(ch_names)
    dist_sec = np.where(dist[0][-3,:]==254)[0]
    prox_sec = np.where(prox[0][-3,:]==254)[0]

    dist_sample = dist[0][:, dist_sec[0]:dist_sec[-1]]
    prox_sample = prox[0][:, prox_sec[0]:prox_sec[-1]]

    del dist, prox
    # Create MNE raw object
    info = mne.create_info(ch_names, fs, ch_types)
    # The lengths of the two datasets are stil off, so remove some samples, still misaligned slightly
    MNE_raw = mne.io.RawArray(np.concatenate((prox_sample[:,1:-1], dist_sample), axis=0), info)

    return MNE_raw

def tmsi_eventextractor(channel_data):
    """Purpose of this function is to extract when the events occur.  It will generate an event array in the form (sample index, 0, event value).
    This is the same form as in MNE-python.
    
    Parameters:
    - channel_data (array): the data from the TMSI file / MNE object.

    Returns:
    - events (array): the event array in the form (sample index, 0, event value).
    """

    # If the data is in MNE format, run the 2 lines below before calling function.
    # Extract the specific channel data - if want to zoom in 100000:150000 is a good spot
    # channel_data = MNE_raw.get_data()[-3, :]-252

    # Find where the data goes from positive to negative 
    transitions = np.where(channel_data[:-1]*channel_data[1:] < 0)[0]
    transitions = transitions[channel_data[transitions] > 0]

    # Construct event array in form (sample index, 0, event value)
    events = np.zeros((len(transitions), 3))
    events[:,0] = transitions
    events[:,2] = channel_data[transitions+1]
    events = events.astype(int)

    return events

# Implement Beamformer
def lcmv_beamformer_constructor(fwd, data_cov, noise_cov=None):
    """
    Constructs a Linearly Constrained Minimum Variance (LCMV) beamformer for each source in the forward model.

    Parameters:
    - fwd (array): Leadfield matrix (n_channels x n_sources).
    - data_cov (array): Data covariance matrix of EMG data (n_channels x n_channels) estimated during signal of interested.  
                        Can be calculated using np.cov(data, rowvar=True). Or use MNE Python's mne.cov.compute_covariance().
    - noise_cov (array): Noise covariance matrix (n_channels x n_channels), optional. If None, don't whiten based on noise covariance.

    Returns:
    - weights (array): Beamformer weights for each source (n_sources x n_channels).
    """

    n_channels, n_sources = fwd.shape

    # If noise covariance matrix provided, whiten the data
    if noise_cov is not None:
        # Compute the whitening matrix using the noise covariance - svd should be same as pca since high-pass filtered. Can consider using ZCA.
        U_noise, s_noise, _ = svd(noise_cov)
        noise_whitening = np.diag(1.0 / np.sqrt(s_noise)) @ U_noise.T
        # Whiten the data covariance matrix
        data_cov = noise_whitening @ data_cov @ noise_whitening.T.conj()

    # Compute the inverse of the data covariance matrix
    data_cov_inv = pinv(data_cov)  # Should be able to just use the inverse?

    # Calculate the beamformer weights for each source
    weights = np.zeros((n_sources, n_channels))
    for i in range(n_sources):
        L_i = fwd[:, i]  # Leadfield for the i-th source
        numerator = np.dot(data_cov_inv, L_i)
        denominator = np.dot(L_i.T, np.dot(data_cov_inv, L_i))
        weights[i, :] = numerator / denominator

    # Leave applying the beamformer to another function
    return weights

# Implement Minumum Norm Estimate (MNE) 
def minimum_norm_estimate(fwd, data, noise_cov=None, reg=0.1):
    """
    Compute the Minimum-Norm Estimate (MNE) for source localization.

    Parameters:
    - fwd (array): Leadfield matrix (n_sensors x n_sources).
    - data (array): Recorded data (n_sensors x n_timepoints).
    - noise_cov (array): Noise covariance matrix (n_sensors x n_sensors). If None, identity is used.
    - reg (float): Regularization parameter (lambda).

    Returns:
    - source_estimates (array): Estimated source activities (n_sources x n_timepoints).
    """
    n_sensors, n_sources = fwd.shape

    # If no noise covariance matrix is provided, use the identity matrix
    if noise_cov is None:
        noise_cov = np.eye(n_sensors)

    # Compute the SVD of the forward matrix
    U, s, Vt = svd(fwd, full_matrices=False)
    
    # Regularization: regularize the singular values
    s_inv = s / (s ** 2 + reg ** 2)
    
    # Compute the inverse operator
    inv_op = Vt.T @ np.diag(s_inv) @ U.T
    
    # Apply the noise covariance matrix
    inv_op = inv_op @ np.linalg.inv(noise_cov)
    
    # Compute the source estimates
    source_estimates = inv_op @ data

    return source_estimates

def sloreta(fwd, data, noise_cov=None, reg=0.1):
    """
    Warning: Taken from ChatGPT, will need to be checked with a paper. 
    Compute sLORETA for source localization.

    Parameters:
    - fwd (array): Leadfield matrix (n_sensors x n_sources).
    - data (array): Recorded data (n_sensors x n_timepoints).
    - noise_cov (array): Noise covariance matrix (n_sensors x n_sensors). If None, it's estimated from the data.
    - reg (float): Regularization parameter (lambda).

    Returns:
    - source_estimates (array): Estimated source activities (n_sources x n_timepoints).
    """
    n_sensors, n_sources = fwd.shape

    # If no noise covariance matrix is provided, estimate it from the data
    if noise_cov is None:
        noise_cov = np.cov(data)

    # Compute the whitening matrix using the noise covariance
    U_noise, s_noise, _ = svd(noise_cov)
    noise_whitening = U_noise @ np.diag(1.0 / np.sqrt(s_noise)) @ U_noise.T

    # Whiten the forward matrix
    fwd_whitened = noise_whitening @ fwd

    # Compute the SVD of the whitened forward matrix
    U, s, Vt = svd(fwd_whitened, full_matrices=False)

    # Regularization: regularize the singular values
    s_inv = s / (s ** 2 + reg ** 2)

    # Compute the inverse operator for sLORETA
    inv_op = Vt.T @ np.diag(s_inv) @ U.T @ noise_whitening

    # Calculate the variance normalization factor (sLORETA)
    source_cov = np.sum(inv_op**2, axis=1)
    normalization = np.sqrt(source_cov)

    # Apply the inverse operator to the data
    source_estimates = inv_op @ data

    # Standardize the source estimates
    source_estimates /= normalization[:, np.newaxis]

    return source_estimates

def find_weighted_centroid(pos, source_activity, fixedorient=True):
    """
    Calculate the weighted centroid of the source activity.
    
    Parameters:
    - pos (array): Positions of the source space voxels (n_voxels x 3).
    - source_activity (array): Source activity (n_voxels x 1).
    - fixedorient (bool): Whether the source activity is fixed orientation (True) or free orientation (False).  
        If free orientation, assume that the source activity is a 3*n_voxels x 1 array .
    
    Returns:
    - weighted_centroid (array): Weighted centroid of the source activity (3 x 1).
    """
    
    if fixedorient is False:
        # Reshape the source activity to a n_voxels x 1 array through normalisation across the 3 orientations
        reshaped_act = source_activity.reshape((3, -1), order='F')
        source_activity = np.linalg.norm(reshaped_act, axis=0)

    # Calculate the weighted sum of positions
    weighted_sum = np.sum(pos * source_activity[:, np.newaxis], axis=0)
    # Calculate the total weight
    total_weight = np.sum(source_activity)
    # Calculate the weighted centroid
    weighted_centroid = weighted_sum / total_weight
    
    return weighted_centroid

def find_matching_rows_floats(pos, dpos, tol=1e-8):
    """
    Find indices of rows in `pos` that match rows in `dpos` within a tolerance.

    Parameters:
    - pos (array): Matrix of positions (n_rows, n_cols).
    - dpos (array): Matrix of target rows to find matches for (m_rows, n_cols).
    - tol (float): Tolerance for floating point comparison (default 1e-8).

    Returns:
    - matching_indices (array): Array of indices of rows in `pos` that match rows in `dpos`.
    """
    # Use broadcasting and np.isclose to compare each row in pos with dpos
    matches = np.all(np.isclose(pos, dpos[:, None], atol=tol), axis=-1)
    
    # Get indices of matching rows
    matching_indices = np.where(matches.any(axis=0))[0]
    
    return matching_indices

def find_closest_rows(pos, dpos):
    """
    Find the indices of the rows in `pos` that are closest to each row in `dpos`.

    Parameters:
    - pos (array): Matrix of positions (n_rows, n_cols).
    - dpos (array): Matrix of target rows (m_rows, n_cols).

    Returns:
    - closest_indices (array): Array of indices in `pos` corresponding to the closest rows for each row in `dpos`.
    """
    closest_indices = []
    
    # Iterate over each row in dpos
    for d_row in dpos:
        # Compute the Euclidean distance between d_row and each row in pos
        distances = np.linalg.norm(pos - d_row, axis=1)
        
        # Find the index of the closest row in pos
        closest_index = np.argmin(distances)
        closest_indices.append(closest_index)
    
    return np.array(closest_indices)

def best_dipole_infwd(fwd, data):
    """ Find the best dipole position in the forward model that would match the data when scaled.  Obtains the strengths of each orientation in said location.
    
    Parameters:
    - fwd (array): Forward model matrix (n_sensors x n_dipoles*3).
    - data (array): Recorded data (n_sensors x n_timepoints).

    Returns:
    - best_index (array): Index of the best dipole in the forward model for each timepoint.
    - best_weights (array): Weights of the best dipole in the forward model for each timepoint.
    - save_arr (array): Array containing the residuals and weights for each dipole (and its three orientations) in the forward model.
    """

    # Reshape fwd - so that the dipole orientation is the third dimension
    fwd = fwd.reshape((fwd.shape[0],-1,3), order='C')

    # Sanity check that the reshaping is correct
    # arr = np.arange(12).reshape(2,6)
    # print(arr)
    # arr.reshape(2,-1,3, order='C')[:,0,:]
    # print(arr)

    # Initialise array based on shape of the input data, whether it is 1D (n_electrodes) or 2D (n_electrodes x n_timepoints)
    if len(data.shape) == 1:
        save_arr_dim2 = 1
    else:
        save_arr_dim2 = data.shape[1]
    save_arr = np.zeros((fwd.shape[1], 4, save_arr_dim2))

    # For each dipole, solve for the weights that would match the data with np.linalg.listsq (Xw = y)
    for i in np.arange(fwd.shape[1]):
        w, residuals, _, _ = np.linalg.lstsq(fwd[:,i,:], data, rcond=None) 
        # Save the weights and residuals
        try:
            save_arr[i,0,:] = residuals
            save_arr[i,1:,:] = w.reshape((3,save_arr_dim2))
        except:
            save_arr[i,0,:] = np.nan
            save_arr[i,1:,:] = np.nan

    best_index = np.nanargmin(save_arr[:,0,:], axis=0)
    best_weights = save_arr[best_index,1:,:].diagonal(axis1=0, axis2=2)

    return best_index, best_weights, save_arr

###########################################################################################################################################

# Optimisers tried which don't work well:

# Equivalent Current Dipole (ECD) fitting using analytical function instead of pre-generated leadfield matrix
# Reason for not using: Slow, unable to fit global minima in a sensible time frame, doesn't converge well
def ECD_fit_dipoles_analytical(data, electrode_pos, n_dipoles=1, initial_guess=None, local=True):
    """
    Fit multiple dipoles to the data using the ECD method and analytical dipole equation.

    Parameters:
    - data: Recorded data (n_sensors x timepoints).
    - electrode_pos: Position of all the electrodes (n_electrodes x 3).
    - n_dipoles: Number of dipoles to fit.  Technically we are fitting 3*n_dipoles dipoles due to 3 orientations per position.
    - initial_guess: Initial guess for the dipole parameters (optional).
    - local: (bool): If True, perform a local optimization. If False, perform a global optimization.
    
    Returns:
    - Optimal dipole parameters (positions, orientations, strengths).
    """

    # Define the bounds for the dipole positions
    bound_max = electrode_pos.max(axis=0)
    bound_min = electrode_pos.min(axis=0)
    bounds = [(bound_min[0], bound_max[0]), (bound_min[1], bound_max[1]), (bound_min[2], bound_max[2]), 
            (None, None), (None, None), (None, None)] * n_dipoles

    # Define the initial conditions
    if initial_guess is None:
        # Generate a random initial guess for the dipole parameters
        # Place in middle of source space if only 1 dipole
        if n_dipoles == 1:
            middle_pos = np.mean(electrode_pos, axis=0)
            initial_guess = np.ones(n_dipoles * 6)
            initial_guess[:3] = middle_pos
        # Otherwise, randomly place the dipoles within the source space
        else:
            initial_guess = np.random.rand(n_dipoles * 6)  # 3 for position, 3 for strength of each orientation
            # Scale the random values to the range of the source space
            index = np.arange(len(initial_guess))%6
            initial_guess[index==0] = initial_guess[index==0] * bound_max[0]
            initial_guess[index==1] = initial_guess[index==1] * bound_max[1]
            initial_guess[index==2] = initial_guess[index==2] * bound_max[2]

    # Nested objective function for the optimization
    def objective_function(d_params, electrode_pos, data):
        """
        Objective function to minimize: the difference between measured data and the forward model prediction.

        Parameters:
        - d_params: Dipole parameters (n_dipoles*6). Each successive 6 indices contain the position (x, y, z) of the nth dipole and strength of each of the 3 dipole orientations.
        - electrode_pos: Position of all the electrodes (n_electrodes x 3).
        - data: Recorded data (n_sensors x timepoints).

        Returns:
        - Error (L2-norm) between predicted data and observed data.
        """
        # Extract some useful variables
        n_dipoles = len(d_params) // 6
        # Reshape
        d_params = d_params.reshape((n_dipoles, 6), order='C')
        d_pos = d_params[:, :3]
        d_strength = d_params[:, 3:]
        d_strength = d_strength.flatten(order='C')

        # Calculate the leadfield matrix for the given dipole parameters
        d_fwd = fwd_generator(dipole_potential, d_pos, electrode_pos)

        # Compute the predicted data
        predicted_data =  d_fwd @ d_strength.T
        
        # Calculate the L2 norm of the error
        error = np.linalg.norm(predicted_data - data)
        return error

    # Use an optimization routine to minimize the objective function
    if local:
        # Perform a local optimization
        result = minimize(objective_function, initial_guess, args=(electrode_pos, data), bounds=bounds,)
    else:
        # Perform a global optimization
        result = scipy.optimize.basinhopping(objective_function, initial_guess, minimizer_kwargs={'args': (electrode_pos, data), 'bounds': bounds})

    print(result)
    d_x = result.x.reshape((n_dipoles, 6), order='C')
    d_pos = d_x[:, :3]
    d_strength = d_x[:, 3:]
    
    # Return the optimized dipole parameters
    return d_pos, d_strength

# Equivalent Current Dipole (ECD) fitting - for offline processing 
# Reason for not using: Does not converge well for global.  Local converges okay for one dipole, but not multiple, and depends on initial guess.
# Additionally, it's a bit silly since we are not taking advantage of the fwd, and the discretised source space (which would need a different package)
def ECD_fit_dipoles(data, fwd, pos, n_dipoles=1, initial_guess=None, local=True):
    """
    Fit multiple dipoles to the data using the ECD method.

    Parameters:
    - data (array): Recorded data (n_sensors x timepoints).
    - fwd (array): Leadfield matrix (n_sensors x n_voxels*3).
    - pos (array): Position of all the source space voxels (n_voxels x 3).
    - n_dipoles (int): Number of dipoles to fit.  Technically we are fitting 3*n_dipoles dipoles due to 3 orientations per position.
    - initial_guess (array): Initial guess for the dipole parameters (optional).  You should use one for a local optimization.
    - local: (bool): If True, perform a local optimization. If False, perform a global optimization.
    
    Returns:
    - Optimal dipole parameters (positions, orientations, strengths).
    """

    # Nested objective function for the optimization
    def objective_function(d_params, pos, fwd, data):
        """
        Objective function to minimize: the difference between measured data and the forward model prediction.

        Parameters:
        - d_params: Dipole parameters (n_dipoles*6). Each successive 6 indices contain the position (x, y, z) of the nth dipole and strength of each of the 3 dipole orientations.
        - pos: Position of all the source space voxels (n_voxels x 3).
        - fwd: Leadfield matrix for all of source space (n_sensors x n_voxels).
        - data: Recorded data (n_sensors x timepoints).

        Returns:
        - Error (L2-norm) between predicted data and observed data.
        """
        # Extract some useful variables
        n_dipoles = len(d_params) // 6
        # Reshape
        d_params = d_params.reshape((n_dipoles, 6), order='C')
        d_pos = d_params[:, :3]
        d_strength = d_params[:, 3:]
        d_strength = d_strength.flatten(order='C')

        # Extract the relevant parts of the leadfield matrix based on the estimated positions
        matching_indices = find_closest_rows(pos, d_pos)
        # Check if the number of matching indices is equal to the number of dipoles
        if len(matching_indices) != n_dipoles:
            print(matching_indices)
            print(d_params)
            raise ValueError("Number of matching indices does not match the number of dipoles.")
        
        # For each dipole, the entries in the leadfield matrix are the 3*matching_indices, 3*matching_indices+1, 3*matching_indices+2 for each orientation.
        matching_indices = np.array([ [3*x, 3*x+1, 3*x+2] for x in matching_indices]).flatten()
        d_fwd = fwd[:, matching_indices]
        
        # Compute the predicted data
        predicted_data =  d_fwd @ d_strength.T
        
        # Calculate the L2 norm of the error
        error = np.linalg.norm(predicted_data - data)
        return error

    # Define the bounds for the dipole positions
    bound_max = pos.max(axis=0)
    bound_min = pos.min(axis=0)
    bounds = [(bound_min[0], bound_max[0]), (bound_min[1], bound_max[1]), (bound_min[2], bound_max[2]), 
            (None, None), (None, None), (None, None)] * n_dipoles

    # Define the initial conditions
    if initial_guess is None:
        # Generate a random initial guess for the dipole parameters
        # Place in middle of source space if only 1 dipole
        if n_dipoles == 1:
            middle_pos = np.mean(pos, axis=0)
            initial_guess = np.ones(n_dipoles * 6)
            initial_guess[:3] = middle_pos
        # Otherwise, randomly place the dipoles within the source space
        else:
            initial_guess = np.random.rand(n_dipoles * 6)  # 3 for position, 3 for strength of each orientation
            # Scale the random values to the range of the source space
            index = np.arange(len(initial_guess))%6
            initial_guess[index==0] = initial_guess[index==0] * (bound_max[0] - bound_min[0]) + bound_min[0]
            initial_guess[index==1] = initial_guess[index==1] * (bound_max[1] - bound_min[1]) + bound_min[1]
            initial_guess[index==2] = initial_guess[index==2] * (bound_max[2] - bound_min[2]) + bound_min[2]

    # Use an optimization routine to minimize the objective function

    if local:
        # Perform a local optimization
        result = minimize(objective_function, initial_guess, args=(pos, fwd, data), bounds=bounds, method='Nelder-Mead', )
    else:
        # Perform a global optimization
        result = scipy.optimize.basinhopping(objective_function, initial_guess, minimizer_kwargs={'args': (pos, fwd, data), 'bounds': bounds, 'method': 'Nelder-Mead'})
    print(result)
    # Adjust the position vector since we take the closest one
    d_x = result.x.reshape((n_dipoles, 6), order='C')
    d_pos = d_x[:, :3]
    d_pos = pos[find_closest_rows(pos, d_pos)]
    d_strength = d_x[:, 3:]
    
    # Return the optimized dipole parameters
    return d_pos, d_strength