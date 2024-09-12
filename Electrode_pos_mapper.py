"""Script for mapping electrode positions to the voxel positions for use in the forward model."""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull, Delaunay
import pandas as pd
from EMGinv_fns import load_src_template

# Scaling factors
z_scaling = 6.3e-3
x_scaling = 1.7e-4
y_scaling = x_scaling

pos = load_src_template(filename=None, con='arm', flip_dim=1, xscaling=x_scaling, yscaling=y_scaling, zscaling=z_scaling)

# Electrode positions - from excel file
filename = '/Users/pokhims/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/Coding/CMU_EMGSL/Data/annotations_l_ends.csv'
electrode_pos = pd.read_csv(filename, header=0)
# electrode_pos = pd.read_excel(filename, header=0)
print(electrode_pos)
# Order it correctly to be same as source space
electrode_pos = electrode_pos[["Y", "X", "Slice"]].to_numpy(dtype='float')
# Scale the values appropriately
electrode_pos[:, :2] = electrode_pos[:,:2] * x_scaling
electrode_pos[:, 2] = (electrode_pos[:, 2]-105) * z_scaling


# Assuming pos is a 2D numpy array with shape (3, V)

# Compute the convex hull
hull = ConvexHull(pos)

# Create a Delaunay triangulation of the convex hull vertices
tri = Delaunay(pos[hull.vertices])
# Check if points in electrode_pos lie within the convex hull
inside = tri.find_simplex(electrode_pos) >= 0
# Print the indices of the points that lie within the convex hull
inside_indices = np.where(inside)[0]
print("Indices of points in electrode_pos that lie within the convex hull:", inside_indices)
# Optionally, print the points themselves
print("Points in electrode_pos that lie within the convex hull:")
print(electrode_pos[inside])


# Plot the convex hull
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
electrode_pos = electrode_pos[:32]

# Plot some points on the hull
sc = ax.scatter(electrode_pos[:, 0], electrode_pos[:, 1], electrode_pos[:, 2], c=np.arange(32), marker='o')
plt.colorbar(sc)

# # Plot the convex hull
for simplex in hull.simplices:
    ax.plot(pos[simplex, 0], pos[simplex, 1], pos[simplex, 2], 'r-')

# Set labels
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

fig.show()

ends_pos = electrode_pos
# Place the electrodes between points on the end slices. 
# 0->7 maps to 16->23, 8->15 maps to 24->31
# Put 19 evenly spaced points between the two points
electrode_pos = np.zeros((256,3)) 

for i in np.arange(16):
    temp = np.linspace(ends_pos[i], ends_pos[i+16], 19)
    # Drop the points not on any patches, i.e. the 5th, 10th, 15th points
    temp = np.delete(temp, [5, 10, 15], axis=0)     
    electrode_pos[i*16:(i+1)*16] = temp


# Shift electrode positions to the closest point on the convex hull surface

# Extract the vertices of the convex hull
hull_points = pos[hull.vertices]

# Function to find the closest point on a triangle
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
        def project_point_on_segment(P, Q, R):
            PQ = Q - P
            t = np.dot(R - P, PQ) / np.dot(PQ, PQ)
            t = np.clip(t, 0, 1)
            return P + t * PQ

        closest_points = [
            project_point_on_segment(A, B, point),
            project_point_on_segment(B, C, point),
            project_point_on_segment(C, A, point)
        ]
        distances = [np.linalg.norm(point - cp) for cp in closest_points]
        return closest_points[np.argmin(distances)]

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

# Move all points in electrode_pos to the closest point on the convex hull surface
electrode_pos = np.array([closest_point_on_hull(hull, point) for point in electrode_pos])

# Check if points in electrode_pos lie within the convex hull
inside = tri.find_simplex(electrode_pos) >= 0
# Print the indices of the points that lie within the convex hull
inside_indices = np.where(inside)[0]
print("Indices of points in electrode_pos that lie within the convex hull:", inside_indices)
# Optionally, print the points themselves
print("Points in electrode_pos that lie within the convex hull:")
print(electrode_pos[inside])

# Plot the convex hull and the moved points
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Plot the convex hull
for simplex in hull.simplices:
    ax.plot(pos[simplex, 0], pos[simplex, 1], pos[simplex, 2], 'r-')
# Plot the moved points on the hull
ax.scatter(electrode_pos[:, 0], electrode_pos[:, 1], electrode_pos[:, 2], c=np.arange(len(electrode_pos)), marker='o')
# Set labels
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
fig.show()

# Load the electrode number mapping from excel file
filename = '/Users/pokhims/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/Coding/CMU_EMGSL/Data/Electrode_pos.xlsx'
lookup = pd.read_excel(filename, sheet_name="Lookup", header=0)
# Obtain the electrode number from the "Channel" column
lookup["Ch_num"] = lookup["Channel"].str.extract('(\d+)', expand=False).astype(int)

# Will need to map twice. 
# First, map each electrode according to its row and column number
elec_pos_df = pd.DataFrame(electrode_pos, columns=['X', 'Y', 'Z'])
RC_array = np.zeros((256, 2))
# Column 
RC_array[:, 1] = np.tile(np.repeat(np.arange(1, 9), 16), 2)
# Row repeats every 4 electrodes for the first 128 electrodes
RC_array[:128, 0] = np.tile(np.arange(1, 5), 32)
# For the last 128 electrodes, the row repeats every 4, but is descending
RC_array[128:, 0] = np.tile(np.arange(8, 4, -1), 32)  
elec_pos_df["RowColumn"] = [f"R{int(row)}C{int(col)}" for row, col in RC_array]
# Now map each electrode according to its row and column to the electrode number
elec_pos_df = elec_pos_df.merge(lookup, on=["RowColumn"], how="left")
# Side on top of arm
side_1 = np.tile(np.repeat(np.arange(1, 8, 2), 4), 8)
# Side on underside of arm
side_2 = np.tile(np.repeat(np.arange(2, 9, 2), 4), 8)
elec_pos_df["Patch"] = np.concatenate((side_1, side_2))

elec_pos_df = elec_pos_df.sort_values(by=["Patch","Ch_num"]).reset_index(drop=True)

# Save the electrode positions as csv
elec_pos_df.to_csv("256leftarm_electrode_pos.csv", index=False)

# Plot the electrode positions
electrode_pos = elec_pos_df[["X", "Y", "Z"]].to_numpy()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Plot the convex hull
for simplex in hull.simplices:
    ax.plot(pos[simplex, 0], pos[simplex, 1], pos[simplex, 2], 'r-')
p = ax.scatter(electrode_pos[:, 0], electrode_pos[:, 1], electrode_pos[:, 2], c=elec_pos_df["Ch_num"], marker='o')
plt.colorbar(p)
# Set labels
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
fig.show()
plt.show()