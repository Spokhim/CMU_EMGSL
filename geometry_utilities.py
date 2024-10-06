import numpy as np
import pandas as pd

def closest_point_on_triangle(triangle, point):
    # Using Barycentric coordinates to find the closest point on the triangle
    A, B, C = triangle
    AB = B - A
    AC = C - A
    AP = point - A

    d1 = np.dot(AB, AP)
    d2 = np.dot(AC, AP)
    d3 = np.dot(AB, AB)
    d4 = np.dot(AC, AC)
    d5 = np.dot(AB, AC)

    denom = d3 * d4 - d5 * d5
    v = (d4 * d1 - d5 * d2) / denom
    w = (d3 * d2 - d5 * d1) / denom
    u = 1 - v - w

    if u >= 0 and v >= 0 and w >= 0:
        return u * A + v * B + w * C
    else:
        # If the point is outside the triangle, project it onto the closest edge
        closest_points = [
            project_point_on_segment(A, B, point),
            project_point_on_segment(B, C, point),
            project_point_on_segment(C, A, point)
        ]
        distances = [np.linalg.norm(point - cp) for cp in closest_points]
        return closest_points[np.argmin(distances)]
    
def project_point_on_segment(P, Q, R):
    PQ = Q - P
    t = np.dot(R - P, PQ) / np.dot(PQ, PQ)
    t = np.clip(t, 0, 1)
    return P + t * PQ

# Function to find the closest point on the convex hull surface
def closest_point_on_hull(hull, point):
    closest_point = None
    min_distance = float('inf')
    for simplex in hull.simplices:
        triangle = hull.points[simplex]
        cp = closest_point_on_triangle(triangle, point)
        distance = np.linalg.norm(cp - point)
        if distance < min_distance:
            min_distance = distance
            closest_point = cp
    return closest_point

def find_nearest_muscle(source_activity, pos, t, options):
    """
    Find the nearest muscle based on the highest-intensity source activity.

    Parameters:
    - source_activity (numpy array): Source activity array with shape (n_sources, n_timepoints).
    - pos (numpy array): Source position array with shape (n_sources, 3).
    - t (int): Time point to analyze.
    - options (dict): Dictionary of options similar to the plotting function.

    Returns:
    - nearest_muscle (str): The name of the nearest muscle.
    """
    # Step 1: Find the highest-intensity source in the source_activity at time t
    source_at_t = source_activity[:, t]
    max_index = np.argmax(np.abs(source_at_t))  # Index of the highest activity voxel
    max_pos = pos[max_index]  # 3D position of the highest intensity source

    # Step 2: Determine the nearest transverse slice (z-index)
    z_index = round(max_pos[2]*1000 / options['MillimetersPerSection']) + options['StartingSection']

    # Step 3: Load the landmarks for the nearest slice
    sheet_name = options['LandMarkSheetExpression'] % z_index
    landmarks_df = pd.read_excel(options['LandMarksFile'], sheet_name=sheet_name)

    # Step 4: Compute the nearest landmark in the 2D plane of the nearest slice
    grid_row_scale = options['GridRows'] / (options['GridRowBottomOffset'] - options['GridRowTopOffset'])
    grid_col_scale = options['GridCols'] / (options['GridColumnRightOffset'] - options['GridColumnLeftOffset'])

    # Convert 3D coordinates to 2D by projecting max_pos onto the slice
    max_x_2d = (max_pos[0] - options['XOffset']) / options['XScale']
    max_y_2d = (max_pos[1] - options['YOffset']) / options['YScale']

    # Iterate through landmarks to find the closest one
    min_distance = np.inf
    nearest_muscle = None

    for i, row in landmarks_df.iterrows():
        x_landmark = (row['X'] - options['GridColumnLeftOffset']) * grid_col_scale
        y_landmark = (row['Y'] - options['GridRowTopOffset']) * grid_row_scale
        
        # Compute Euclidean distance between max_pos and the landmark
        distance = np.sqrt((max_x_2d - x_landmark)**2 + (max_y_2d - y_landmark)**2)
        
        if distance < min_distance:
            min_distance = distance
            nearest_muscle = row['Landmark']

    return nearest_muscle, min_distance, z_index