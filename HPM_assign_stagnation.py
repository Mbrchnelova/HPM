from HPM_import import *

def AssgnStagnationPoint(node, nodes, Minf, pinf, rhoinf, Tinf, p_nodes, rho_nodes, T_nodes, mu_nodes, u_nodes, v_nodes, w_nodes, VERIFY):
    gamma = 1.4
    R = 287.05
    for n in range(0, len(nodes)):
        if nodes[n] == node:
            #pstag = pinf * (1. + (gamma - 1.)/2. * Minf**2)**(gamma/(gamma - 1.))   
            #pstag_behindNSW = pstag * (((gamma + 1.)*Minf**2)/ (2. + (gamma - 1.)*Minf**2))**(gamma/(gamma - 1.)) * ((gamma + 1.)/(2. * gamma * Minf**2 - (gamma - 1.)))**(1./(gamma - 1.))
            #Tstag = Tinf * (1. + (gamma - 1.)/2. * Minf**2)
            #Tstag_behindNSW = Tstag
            #rhostag_behindNSW = pstag_behindNSW / (R * Tstag_behindNSW)
            #if rhostag_behindNSW > rhoinf * 6.:
            #    rhostag_behindNSW = rhoinf * 6.
   
            #S = 110.
            #mustag_behindNSW = muinf* (Tstag_behindNSW/Tinf)**(3./2.) * (Tinf + S)/(Tstag_behindNSW + S)
   
            #p_nodes[n] = pstag_behindNSW
            #T_nodes[n] = Tstag_behindNSW
            #mu_nodes[n] = mustag_behindNSW
            #rho_nodes[n] = rhostag_behindNSW
            Mach_behindNSW = math.sqrt(((gamma - 1.0)*Minf**2 + 2.0)/(2. * gamma * Minf**2 - gamma + 1.))
            u_nodes[n] = Mach_behindNSW * math.sqrt(gamma * R * T_nodes[n])
            v_nodes[n] = 0.
            w_nodes[n] = 0.
    '''            
    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node[0], node[1], node[2], s=10, c = 'r')
    '''
    return p_nodes, T_nodes, mu_nodes, rho_nodes, u_nodes, v_nodes, w_nodes
