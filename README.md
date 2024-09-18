# CMU_EMGSL

This repository contains mostly Python code for the purpose of performing inverse source localisation for EMG data.  
`EMG_pipeline.ipynb` is the main file for analysis and demonstrates how to use all the necessary functions contained in `EMGinv_fns.py`. 

## Features

- **Data Preprocessing**: Load and preprocess EMG data from tMSI devices.
- **Inverse Source Localisation**: Perform inverse source localisation on EMG using a few basic algorithms.
- **Visualisation**: Plot the results of the inverse source localisation.

## File Structure

- **EMGinv_fns.py**: Contains all the functions necessary for performing inverse source localisation.
- **EMG_pipeline.ipynb**: The main Jupyter notebook for analysis.
- **README.md**: This documentation file.

## Inputs

The inputs for the EMG inverse source localisation are as follows:

1. **Electrode positions**: The electrode positions should be stored in a file with the following columns: 'Channel_Name', 'X', 'Y', 'Z'. The 'Channel_Name' column should contain the channel nuname, and the 'X', 'Y', 'Z' columns should contain the corresponding coordinates.  
    * Alternatively, positions can be created following the steps `Electrode_pos_mapper.py`.  This involves an initial selection of electrode positions defining the two ends of the grid such as using https://github.com/Neuro-Mechatronics-Interfaces/Annotation_Tools on MRI images and then using the `Electrode_pos_mapper.py` to interpolate the positions of the remaining electrodes.
2. **EMG data**: The EMG data can be loaded with `load_tmsi_data` or `load_tmsi_combine`.  The data should end in the form of an array with dimensions (n_channels, n_samples).
3. **Forward model**: The forward model can be created with `create_forward_model`.  The forward model is a matrix that maps the source activity to the electrode data.  The forward model can be created using the electrode positions and the head model.

## Processing Steps

The processing steps in the EMG source inversion pipeline are as follows:

1. **Load Data**: Load the EMG data and electrode positions.
2. **Preprocess Data**: Preprocess the EMG data by filtering and normalising it.
3. **Construct Noise and Data Covariance Matrices**: Construct the noise and data covariance matrices from the EMG data.
4. **Solver**: Solve the inverse problem using one of the available solvers.
5. **Visualise Results**: Visualise the results of the inverse source localisation on image slices.