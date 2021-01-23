# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_separate_centroids.py: this module is responsible for separating the elements based on their x direction.
# This is done to speed up the backtracing process (when locating in which element we are currently residing to compute velocity). 
# Depending on the geometry, this can be also re-written for y and z axes (for y-axis and z-axis dominant geometries).
# If modified, the find_velocity -ordered method has to be also adjusted accordingly. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *
from HPM_mathematics import *

def SeparateCentroidsByDx(centroids, xmin, xmax, N):
    centroids_ordered = []
    centroids_ordered_dxs = []
    dx = abs(-xmin + xmax)/N
    for i in range(0, N+1):
        centroids_ordered_dxs.append(xmax - (i*dx))
        centroids_ordered.append([])

    for i in range(0, len(centroids)):
        min_x_dist = 1e9
        dx_region = -1
        for j in range(0, N+1):
            x_dist = abs(centroids_ordered_dxs[j] - centroids[i][0])
            if x_dist < min_x_dist:
                min_x_dist = x_dist
                dx_region = j

        centroids_ordered[dx_region].append([i, centroids[i]])

    centroids_ordered_final = []
    centroids_ordered_dxs_final = []
    for i in range(0, len(centroids_ordered)):
        if len(centroids_ordered[i]) > 0.:
            centroids_ordered_final.append(centroids_ordered[i])
            centroids_ordered_dxs_final.append(centroids_ordered_dxs[i])

    return centroids_ordered_final, centroids_ordered_dxs_final

