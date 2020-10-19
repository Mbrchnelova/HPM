# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_get_thermodynamics.py: this function computes the "inviscid thermodynamics" from the pressure and velocity field.
# The thermodynamic variables (temperature, density, viscosity) are computed at the node locations.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *
from HPM_mathematics import *
from HPM_SETUP import *



def GetThermodynamics(ps, Minf, nodes, gamma, R, muinf, Tinf, rhoinf, VERIFY):
    if gamma == 1.:
        print "Gamma cannot be 1! Adjusted to prevent singularities."
        gamma += 0.001

    if gamma == 0.:
        print "Gamma cannot be 0! Adjusted to prevent singularities. Results will be unphysical until you get your shit together with the input data."
        gamma += 0.001

    pstag = pinf * (1. + (gamma - 1.)/2. * Minf**2)**(gamma/(gamma - 1.))
    pstag_behindNSW = pstag * (((gamma + 1.)*Minf**2)/ (2. + (gamma - 1.)*Minf**2))**(gamma/(gamma - 1.)) * ((gamma + 1.)/(2. * gamma * Minf**2 - (gamma - 1.)))**(1./(gamma - 1.))
    Tstag = Tinf * (1. + (gamma - 1.)/2. * Minf**2)
    Tstag_behindNSW = Tstag
    rhostag_behindNSW = pstag_behindNSW / (R * Tstag_behindNSW)
    if rhostag_behindNSW > rhoinf * 6.:
        print "Accuracy compromised, rho behind NSW had to be adjusted to the maximum value (6rho_inf). The resulting data are unreliable."
        rhostag_behindNSW = rhoinf * 6.
    T_nodes = []
    mu_nodes = []
    rho_nodes = []
    M_nodes = []
    for n in range(0, len(nodes)):
        p = ps[n]
        if p == 0.:
            "Zero pressure detected on the windward side! Check the suitability of the mesh and input data. Adjusted, but unphysical."
            p += 0.001

        T = Tstag_behindNSW * (p/pstag)**((gamma - 1.)/gamma)
        T = Tstag_behindNSW * (p/pstag_behindNSW)**((gamma - 1.)/gamma)

        T_nodes.append(T)
        if p > pstag_behindNSW:
            p = pstag_behindNSW
        M = math.sqrt( 2. / (gamma - 1.) * ((p/pstag_behindNSW)**((gamma - 1.)/(-gamma)) - 1.))
        M_nodes.append(M)
        S = 110.
        mu = muinf* (T/Tinf)**(3./2.) * (Tinf + S)/(T + S)
        mu_nodes.append(mu)
        #rho = p/(R * T)
        rho = (p/pstag_behindNSW)**(1./gamma) * rhostag_behindNSW
        rho_nodes.append(rho)


    if VERIFY:
        x_coords_ele = []
        y_coords_ele = []
        z_coords_ele = []
        for n in range(0, len(nodes)):
            x_coords_ele.append(nodes[n][0])
            y_coords_ele.append(nodes[n][1])
            z_coords_ele.append(nodes[n][2])

        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(T_nodes)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('T of nodes on a cone using modified Newton technique')


        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(rho_nodes)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Rho of nodes on a cone using modified Newton technique')

        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(mu_nodes)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('mu of nodes on a cone using modified Newton technique')

    return T_nodes, mu_nodes, rho_nodes, M_nodes
