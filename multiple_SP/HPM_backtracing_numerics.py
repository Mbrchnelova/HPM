# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_backtracing_numerics: all functions which are used for backtracing specifically go here. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *



def CrossedStagnationLine(coordinates, stag_node, epsilon):
    dist_to_stag_point = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
    if dist_to_stag_point <= epsilon:
        return [True, dist_to_stag_point]
    else:
        return [False]


###############################################################################



def FindLength(nodes):
    xmax = -1e9
    for n in range(0, len(nodes)):
        if nodes[n][0] > xmax:
            xmax = nodes[n][0]
    return xmax



###############################################################################



def FindBounds(nodes):
    xmax = -1e9
    xmin = 1e9
    ymax = -1e9
    ymin = 1e9
    zmax = -1e9
    zmin = 1e9
    for n in range(0, len(nodes)):
        if nodes[n][0] > xmax:
            xmax = nodes[n][0]
        if nodes[n][0] < xmin:
            xmin = nodes[n][0]

        if nodes[n][1] > ymax:
            ymax = nodes[n][1]
        if nodes[n][1] < ymin:
            ymin = nodes[n][1]

        if nodes[n][2] > zmax:
            zmax = nodes[n][2]
        if nodes[n][2] < zmin:
            zmin = nodes[n][2]

    return xmin, xmax, ymin, ymax, zmin, zmax

