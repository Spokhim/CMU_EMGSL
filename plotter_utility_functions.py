import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    options.setdefault('XOffset', -0.020)
    options.setdefault('YOffset', -0.020)
    options.setdefault('AddLabels', True)
    options.setdefault('AddImage', True)

    # Load the image
    img_file = options['FileExpression'] % z
    img = plt.imread(img_file)

    # Load the landmarks from the Excel file
    sheet_name = options['LandMarkSheetExpression'] % z
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
                print(f"x: {x} | y: {y}")
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
