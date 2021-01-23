# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_process_mesh.py: here, the main mesh process for .ply files is defined.
# Additional one may be added to handle other than ply formats. 
# If an additional mesh processer is added, the outputs must be the same. 
# The inputs must be at least the node coordinates, connectivity and the directions of the outward normals.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *

def process_mesh(filename, VERIFY):
    with open(filename, "r") as ins:
        array = []
        for line in ins:
            array.append(line)

    alllines = []
    COUNT = False
    for line in array:
        if line[-1] == '\n':
            line = line[0:-1]

        line = line.strip(' ')

        if COUNT:
            alllines.append(line.split(' '))

        if line == 'end_header':
            COUNT = True

    nodes = []
    normals_nodes = []
    connectivity = []
    centroids = []
    normals = []
    test_dot_products = []
    areas = []

    for line in alllines:
        if len(line) >= 6:
            nodes.append([float(line[0]), float(line[1]), float(line[2])])
            normals_nodes.append([float(line[3]), float(line[4]), float(line[5])])

        if len(line) == 4:
            connectivity.append([int(line[1]), int(line[2]), int(line[3])])

    elements_num = len(connectivity)

    if len(connectivity[0]) > 3:
        print "Incorrect mesh formar! Triangular elements required, please re-format and try again."

    for i in range(0, elements_num):
        first_point = connectivity[i][0]
        second_point = connectivity[i][1]
        third_point = connectivity[i][2]

        ctr_x = (nodes[first_point][0] + nodes[second_point][0] + nodes[third_point][0])/3.
        ctr_y = (nodes[first_point][1] + nodes[second_point][1] + nodes[third_point][1])/3.
        ctr_z = (nodes[first_point][2] + nodes[second_point][2] + nodes[third_point][2])/3.

        centroids.append([ctr_x, ctr_y, ctr_z])

    for i in range(0, elements_num):
        first_point = nodes[connectivity[i][0]]
        second_point = nodes[connectivity[i][1]]
        third_point = nodes[connectivity[i][2]]

        vector1 = [second_point[0] - first_point[0], second_point[1] - first_point[1], second_point[2] - first_point[2]]
        vector2 = [third_point[0] - first_point[0], third_point[1] - first_point[1], third_point[2] - first_point[2]]

        n_x = vector1[1] * vector2[2] - vector1[2] * vector2[1]
        n_y = vector1[2] * vector2[0] - vector1[0] * vector2[2]
        n_z = vector1[0] * vector2[1] - vector1[1] * vector2[0]

        normal = [n_x, n_y, n_z]

        dot_product = n_x * centroids[i][0] +  n_y * centroids[i][1] +  n_z * centroids[i][2]

        if dot_product < 0.:
            normal = [-n_x, -n_y, -n_z]

        test_normal_1 = normals_nodes[connectivity[i][0]]
        test_normal_2 = normals_nodes[connectivity[i][1]]
        test_normal_3 = normals_nodes[connectivity[i][2]]

        ave_x = (test_normal_1[0] + test_normal_2[0] + test_normal_3[0])/3.
        ave_y = (test_normal_1[1] + test_normal_2[1] + test_normal_3[1])/3.
        ave_z = (test_normal_1[2] + test_normal_2[2] + test_normal_3[2])/3.

        test_normal_ave = [ave_x, ave_y, ave_z]
        #print test_normal_1, test_normal_2, test_normal_3, test_normal_4
        test_dot_product = normal[0] * test_normal_ave[0] +  normal[1] * test_normal_ave[1] +  normal[2] * test_normal_ave[2]

        test_dot_products.append(test_dot_product)

        if test_dot_product > 0:
            normal = [-n_x, -n_y, -n_z]
            #print 'problem', test_dot_product

        normals.append(normal)

    for i in range(0, elements_num):
        first_point = nodes[connectivity[i][0]]
        second_point = nodes[connectivity[i][1]]
        third_point = nodes[connectivity[i][2]]

        vector1 = [second_point[0] - first_point[0], second_point[1] - first_point[1], second_point[2] - first_point[2]]
        vector2 = [third_point[0] - first_point[0], third_point[1] - first_point[1], third_point[2] - first_point[2]]

        a_x = vector1[1] * vector2[2] - vector1[2] * vector2[1]
        a_y = vector1[2] * vector2[0] - vector1[0] * vector2[2]
        a_z = vector1[0] * vector2[1] - vector1[1] * vector2[0]

        area = 0.5 * (a_x **2 + a_y**2 + a_z**2)*0.5
        areas.append(area)

    #RETURN node coordinates, node connectivity, element normals, element centroids, element areas

    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)

        x_coords_ele = []
        y_coords_ele = []
        z_coords_ele = []
        #u_vel_mag = []
        for c in range(0, len(centroids)):
            x_coords_ele.append(centroids[c][0])
            y_coords_ele.append(centroids[c][1])
            z_coords_ele.append(centroids[c][2])
            #u_vel_mag.append(velocities[c][2])
        c = np.array(areas)
        ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Nearest points resolved for parameter fitting')


        # Here you have to insert the node indices which were resolved as being used:
        t_check = [96, 45, 234, 97, 232, 800]
        x_coords_neigh_check = [nodes[t_check[0]][0], nodes[t_check[1]][0], nodes[t_check[2]][0], nodes[t_check[3]][0], nodes[t_check[4]][0], nodes[t_check[5]][0]]
        y_coords_neigh_check = [nodes[t_check[0]][1], nodes[t_check[1]][1], nodes[t_check[2]][1], nodes[t_check[3]][1], nodes[t_check[4]][1], nodes[t_check[5]][1]]
        z_coords_neigh_check = [nodes[t_check[0]][2], nodes[t_check[1]][2], nodes[t_check[2]][2], nodes[t_check[3]][2], nodes[t_check[4]][2], nodes[t_check[5]][2]]

        ax.scatter(x_coords_neigh_check, y_coords_neigh_check, z_coords_neigh_check, s=15, c = 'r')

    return[nodes, connectivity, normals, centroids, areas]


