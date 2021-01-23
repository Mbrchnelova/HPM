# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_assign_stagnation: this module (currently) assigns the correct behind-NSW velocity on the stagnation line. 
# The rest of the stagnation properties (apart from the velocities) are computed in HPM_calculate_stagnation.py. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be

from HPM_import import *

def AssgnStagnationPoint(node, nodes, Minf, pinf, rhoinf, Tinf, p_nodes, rho_nodes, T_nodes, mu_nodes, u_nodes, v_nodes, w_nodes, VERIFY):
    gamma = 1.4
    R = 287.05
    for n in range(0, len(nodes)):
        if nodes[n] == node:
            Mach_behindNSW = math.sqrt(((gamma - 1.0)*Minf**2 + 2.0)/(2. * gamma * Minf**2 - gamma + 1.))
            u_nodes[n] = Mach_behindNSW * math.sqrt(gamma * R * T_nodes[n])
            v_nodes[n] = 0.
            w_nodes[n] = 0.
    return p_nodes, T_nodes, mu_nodes, rho_nodes, u_nodes, v_nodes, w_nodes
