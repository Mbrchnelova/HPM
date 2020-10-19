# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_find_stagnation: this module contains functions which help determine which node is the stagnation point. 
# Based on user's input, this is currently done based on either tracing the front-most point (geometrically) or from pressure.
# Note that due to the fact that pressure is element-resolved, calculation from pressure cannot be completely accurate.
# For coarse geometries and geometries with one stagnation region only, geometry-based calculation is thus preferred.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be


from HPM_import import *



def FindStagnationPointFromPressure(nodes, pressures):
    pmax = 0.
    for n in range(0, len(nodes)):
        p = pressures[n]
        if p > pmax:
            stag_point = nodes[n]
            stag_p = p
            pmax = p
            nstag = n

    return [stag_point], [stag_p], [nstag]



###############################################################################



def FindStagnationPointFromGeometry(nodes, pressures):
    xmax = 1e9
    stag_points = []
    stag_ps = []
    nstags = []
    for n in range(0, len(nodes)):
        x = nodes[n][0]
        if x < xmax:
            stag_point = nodes[n]
            xmax = x

    for n in range(0, len(nodes)):
        x = nodes[n][0]
        if x == xmax:
            stag_point = nodes[n]
            stag_points.append(stag_point)
            stag_ps.append(pressures[n])
            nstags.append(n)
    return stag_points, stag_ps, nstags



###############################################################################



def FindStagnationPointGeneral(centroids, u_velocity, v_velocity, w_velocity):
    return



