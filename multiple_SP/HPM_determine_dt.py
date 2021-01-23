# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_determine_dt.py: this function provides an estimate for the desired dt during backtracing. 
# For some geometries or conditions, this function might need to be updated to give better convergence.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *
from HPM_mathematics import *



def DetermineDt(current_coord, u_vector, element, connectivity):
    element_nodes = connectivity[element]
    max_x = -1e9
    for i in range(0, len(element_nodes)):
        node = nodes[element_nodes[i]]
        if node[0] > max_x:
            max_x = node[0]
    dist_to_max_x = abs(max_x - current_coord[0])
    u_vec = abs(u_vector[0])
    if round(u_vec, 2) != 0.:
        dt = dist_to_max_x / u_vec
    else:
        dt = 1e-6
    return dt


