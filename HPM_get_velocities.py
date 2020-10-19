# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_get_velocities.py: this function computes the surface (inviscid) velocities. 
# They are determined from the pressure field, freestream conditions and projection onto the surface.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be


from HPM_import import *
from HPM_mathematics import *
from HPM_SETUP import *




def GetVelocities(V_vector, normals, centroids, Cps, pinf, Minf, Tinf, gamma, R, VERIFY):
    if gamma == 1.:
        print "Gamma cannot be 1! Adjusted to prevent singularities."
        gamma += 0.001

    if gamma == 0.:
        print "Gamma cannot be 0! Adjusted to prevent singularities. Results will be unphysical until you get your shit together with the input data."
        gamma += 0.001

    velocities = []
    pressures = []
    pstag = pinf * (1. + (gamma - 1.)/2. * Minf**2)**(gamma/(gamma - 1.))
    pstag_behindNSW = pstag * (((gamma + 1.)*Minf**2)/ (2. + (gamma - 1.)*Minf**2))**(gamma/(gamma - 1.)) * ((gamma + 1.)/(2. * gamma * Minf**2 - (gamma - 1.)))**(1./(gamma - 1.))
    Tstag = Tinf * (1. + (gamma - 1.)/2. * Minf**2)
    Tstag_behindNSW = Tstag
    rhostag_behindNSW = pstag_behindNSW / (R * Tstag_behindNSW)
    if rhostag_behindNSW > rhoinf * 6.:
        rhostag_behindNSW = rhoinf * 6.

    for c in range(0, len(centroids)):
        if Cps[c] == 0.:
            Cps[c] = - 1./Minf**2
        n = normals[c]
        p = pinf * (Cps[c] * gamma * Minf**2 / 2. + 1.)
        if p == 0.:
            "Zero pressure detected on the windward side! Check the suitability of the mesh and input data. Adjusted, but unphysical."
            p += 0.001
        T = Tstag_behindNSW * (p/pstag_behindNSW)**((gamma - 1.)/gamma)
        pressures.append(p)
        #print 2. / (gamma - 1.) * ((p/pstag_behindNSW)**((gamma - 1.)/(-gamma)) - 1.), p/pstag_behindNSW
        if p/pstag_behindNSW > 1.:
            p = pstag_behindNSW
        M = math.sqrt( 2. / (gamma - 1.) * ((p/pstag_behindNSW)**((gamma - 1.)/(-gamma)) - 1.))
        #M = math.sqrt(2./(gamma * Cps[c]) * ((p/pinf) - 1.))
        if (gamma * R * T ) < 0.:
            print "The product of gamma * R * T is negative. Calculation will be terminated."
        V_mag = math.sqrt(gamma * R * T ) * M
        #print M, (gamma * Cps[c]),  ((p/pinf) - 1.), p, V_mag
        if Mag(n) == 0.:
            print "Normal of a zero magnitude detected. Adjusted to prevent singularities. Please fix your mesh for physical results."
            n[0] += 0.001
        projection = Dot(V_vector, n) / Mag(n) * 1. / Mag(n)
        v_proj = [n[0] * projection, n[1] * projection, n[2] * projection]
        v_rej = [V_vector[0] - v_proj[0], V_vector[1] - v_proj[1], V_vector[2] - v_proj[2] ]
        v_rej_mag = (v_rej[0]**2 + v_rej[1]**2 + v_rej[2]**2)**0.5
        if v_rej_mag > 0.:
            v_rej = [v_rej[0]/ v_rej_mag, v_rej[1]/ v_rej_mag, v_rej[2]/ v_rej_mag]
            V = [V_mag * v_rej[0], V_mag * v_rej[1], V_mag * v_rej[2]]
        else:
            V = [0., 0., 0.]
        velocities.append(V)


    if VERIFY:
        x_coords_ele = []
        y_coords_ele = []
        z_coords_ele = []
        u_vel_mag = []
        v_vel_mag = []
        w_vel_mag = []
        for c in range(0, len(centroids)):
            if centroids[c][0] > -2.:
                x_coords_ele.append(centroids[c][0])
                y_coords_ele.append(centroids[c][1])
                z_coords_ele.append(centroids[c][2])
                u_vel_mag.append(velocities[c][0])
                v_vel_mag.append(velocities[c][1])
                w_vel_mag.append(velocities[c][2])

        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(u_vel_mag)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('U velocity on a cone using modified Newton technique')

        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(v_vel_mag)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('V velocity on a cone using modified Newton technique')

        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(w_vel_mag)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('W velocity on a cone using modified Newton technique')


    return velocities, pressures
