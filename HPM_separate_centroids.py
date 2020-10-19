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

