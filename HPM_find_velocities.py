# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_find_velocities.py: this module contains function to determine the local velocity while backtracing.
# Either the centroid-ordering has been performed - then use is made of FindVelVectorOrdered. 
# Otherwise, the slower FindVelVectorNotOrdered should be implemented.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *
from HPM_mathematics import *



def FindVelVectorNotOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids):
    dmin = 1e10
    for e in range(0, len(centroids)):
        centroid = centroids[e]
        d = ((coordinates[0] - centroid[0])**2 + (coordinates[1] - centroid[1])**2 + (coordinates[2] - centroid[2])**2)**0.5
        if d < dmin:
            dmin = d
            c = e
            found_element = e

    local_velocity = [us_mesh[c], vs_mesh[c], ws_mesh[c]]
    return local_velocity, found_element



###############################################################################



def FindVelVectorOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids_ordered_local, centroids):
    dmin = 1e10
    for e in range(0, len(centroids_ordered_local)):
        centroid = centroids_ordered_local[e][1]
        d = ((coordinates[0] - centroid[0])**2 + (coordinates[1] - centroid[1])**2 + (coordinates[2] - centroid[2])**2)**0.5
        if d < dmin:
            dmin = d
            index = centroids_ordered_local[e][0]
            found_element = index

    local_velocity = [us_mesh[index], vs_mesh[index], ws_mesh[index]]
    return local_velocity, found_element


