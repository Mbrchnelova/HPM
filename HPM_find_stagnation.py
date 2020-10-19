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



