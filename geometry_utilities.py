import numpy as np

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