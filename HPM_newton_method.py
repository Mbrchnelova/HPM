from HPM_import import *
from HPM_mathematics import *

def Newton(nodes, connectivity, normals, centroids, areas, gamma, V_vector, Mach, VERIFY):
    x_coords = []
    y_coords = []
    z_coords = []

    for node in nodes:
        x_coords.append(node[0])
        y_coords.append(node[1])
        z_coords.append(node[2])

    x_coords_ele = []
    y_coords_ele = []
    z_coords_ele = []

    x_coords_ele_cp = []
    y_coords_ele_cp = []
    z_coords_ele_cp = []

    nonzero_cp = []
    shadowed = []

    for centroid in centroids:
        x_coords_ele.append(centroid[0])
        y_coords_ele.append(centroid[1])
        z_coords_ele.append(centroid[2])


    V_mag = (V_vector[0]**2 + V_vector[1]**2 + V_vector[2]**2)**0.5
    Cps = []

    if gamma == 1.:
        print "Gamma cannot be 1! Adjusted to prevent singularities in calculations."
        gamma += 0.001

    for i in range(0, len(connectivity)):
        normal = normals[i]
        dot_product = -(normal[0] * V_vector[0] + normal[1] * V_vector[1] + normal[2] * V_vector[2])
        normal_mag = (normal[0]**2 + normal[1]**2 + normal[2]**2)**0.5



        #Cp_max = (1. - gamma + 2. * Mach**2 * gamma)/(gamma + 1.) * ((gamma + 1.)**2 * Mach**2/(4. * gamma * Mach**2 - 2.*( gamma - 1.)))**(-gamma / (gamma - 1.))
        Cp_max = ((gamma + 1.)**2 / (4. * gamma))**(gamma / (gamma - 1.))  * (4. / (gamma + 1.))
        if dot_product < 0.:
            theta = math.acos(dot_product/(V_mag * normal_mag))
            Cp = Cp_max * math.sin(math.pi/2. - theta)**2
            shadowed.append('False')
            #print theta
            #print dot_product, normal_mag
        else:
            Cp = 0.
            shadowed.append('True')

        if Cp == 0:
            centroid = centroids[i]
            x_coords_ele_cp.append(centroid[0])
            y_coords_ele_cp.append(centroid[1])
            z_coords_ele_cp.append(centroid[2])
            nonzero_cp.append(Cp)

        Cps.append(Cp)
    #print Cp_max

    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(Cps)

        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Cp on a cone using modified Newton technique')
        ax.axes.set_aspect('equal')

    return Cps, shadowed

