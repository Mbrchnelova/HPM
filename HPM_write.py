from HPM_import import *




def WriteRunInfo(filename, data):

    pinf = data[0]
    Tinf = data[1]
    Mach = data[2]
    rhoinf = data[3]
    muinf = data[4]
    R = data[5]
    cp = data[6]
    gamma = data[7]
    name = data[8]
    USE_IDS = data[9]
    IDWP = data[10]
    COMPUTE_STAGNATION_POINT_FROM_GEOMETRY = data[11]
    epsilon = data[12]
    T_w = data[13]
    PERZHIKAR = data[14]
    INTERPOLATE_N = data[15]
    beta_tol = data[16]
    OWN_TRANSITION = data[17]
    max_x_user  = data[18]
    dist_yz_min = data[19]

    f = open(filename, "w")
   
    f.write("Simulation at p(inf), T(inf), M(inf), rho(inf), mu(inf), R, cp and gamma: \n" + str(pinf) + ', ' + str(Tinf) + ', ' + str(Mach) + ', ' + str(rhoinf) + ', ' + str(muinf) + ', ' + str(R) + ', ' + str(cp) + ', ' + str(gamma) + '\n')
    f.write("Mesh used: " + str(name) + '\n \n')
   
    f.write("####### Inviscid computation data " + '\n') 
    f.write("IDW used: " + str(USE_IDS) + " with IDWP: " + str(IDWP) + '\n' )
    f.write("SP computed from geometry: " + str(COMPUTE_STAGNATION_POINT_FROM_GEOMETRY) + '\n\n' )

    f.write("####### Viscous computation data " + '\n')
    f.write("Epsilon: " + str(epsilon) + '\n')
    f.write("T wall: " + str(T_w) + '\n')
    f.write("Perzhikar used: " + str(PERZHIKAR) + " " + '\n')
    f.write("Turbulent N was interpolated rather than computed: " + str(INTERPOLATE_N) + '\n')
    f.write("Beta tolerance while projecting (Hamilton): " + str(beta_tol) + '\n')
    f.write("Own transition vector: " + str(OWN_TRANSITION) + '\n')
    f.write("Data ignored beyond: " + str(max_x_user) + '\n')
    f.write("Min. distance between streamlines (Perzhikar): " + str(dist_yz_min) + '\n')
   
    f.close()

    return True



###############################################################################



def WriteDataForSurfaceFit(filenamevel, filenamecoord, connectivity, centroids, u_nodes, v_nodes, w_nodes, nodes):
    f = open(filenamevel, "w")
   
    for c in range(0, len(centroids)):
        conn = connectivity[c]
        vertex_1 = conn[0]
        vertex_2 = conn[1]
        vertex_3 = conn[2]

        u = [u_nodes[vertex_1], u_nodes[vertex_2], u_nodes[vertex_3]]
        v = [v_nodes[vertex_1], v_nodes[vertex_2], v_nodes[vertex_3]]
        w = [w_nodes[vertex_1], w_nodes[vertex_2], w_nodes[vertex_3]]

        line = str(u[0]) + ' ' +  str(u[1]) + ' ' +  str(u[2]) + ' ' +  str(v[0]) + ' ' +  str(v[1]) + ' ' +  str(v[2]) + ' ' +  str(w[0]) + ' ' +  str(w[1]) + ' ' +  str(w[2]) + ' ' + '\n'
        f.write(line)
    f.close()
   
    f = open(filenamecoord, "w")

    for c in range(0, len(centroids)):
        conn = connectivity[c]
        vertex_1 = conn[0]
        vertex_2 = conn[1]
        vertex_3 = conn[2]

        coord_1 = nodes[vertex_1]
        coord_2 = nodes[vertex_2]
        coord_3 = nodes[vertex_3]

        line = str(coord_1[0]) + ' ' + str(coord_1[1]) + ' ' + str(coord_1[2]) + ' ' + str(coord_2[0]) + ' '+ str(coord_2[1]) + ' '+ str(coord_2[2]) + ' ' + str(coord_3[0]) + ' '+ str(coord_3[1]) + ' ' + str(coord_3[2]) + '\n'
        f.write(line)
    f.close()
    return True



###############################################################################



def WriteDataNodes(filename, data, nodes):
    f = open(filename, "w")

    for n in range(0, len(nodes)):
        node = nodes[n]
        data_n = data[n]
        line = str(node[0]) + ' ' +  str(node[1]) + ' ' +  str(node[2]) + ' ' + str(data_n) + '\n'
        f.write(line)
    f.close()

    return True



###############################################################################



def WriteDataCentroids(filename, filenamecon, centroids, connectivity):
    f = open(filename, "w")

    for n in range(0, len(centroids)):
        centroid = centroids[n]
        line = str(centroid[0]) + ' ' +  str(centroid[1]) + ' ' +  str(centroid[2]) + '\n'
        f.write(line)
    f.close()
   
    f = open(filenamecon, "w")

    for n in range(0, len(connectivity)):
        conn = connectivity[n]
        line = str(conn[0]) + ' ' +  str(conn[1]) + ' ' +  str(conn[2]) + '\n'
        f.write(line)
    f.close()

    return True



###############################################################################



def WriteStagPointData(filename, stag_points, stag_idxs, epsilon):
    f = open(filename, "w")

    i = 0
    for stag_point in stag_points:
        stag_idx = stag_idxs[i]
        line = str(stag_point[0]) + ' ' +  str(stag_point[1]) + ' ' +  str(stag_point[2]) + ' ' + str(stag_idx) + ' ' + str(epsilon) + '\n'
        f.write(line)
        i += 1
    f.close()
   
    return True



###############################################################################



def WriteNormals(filename, normals):
    f = open(filename, "w")

    for n in range(0, len(normals)):
        normal = normals[n]
        line = str(normal[0]) + ' ' +  str(normal[1]) + ' ' +  str(normal[2]) + '\n'
        f.write(line)
    f.close()
   


###############################################################################



def WriteBackTracing(filename_intersections, filename_nodepaths, filename_nodecoords, filename_nodesresolved, filename_epsilonnodes, intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs):
    f = open(filename_intersections, "w")
    for n in range(0, len(intersections)):
        intersection = intersections[n]
        line = str(intersection[0]) + ' ' +  str(intersection[1]) + ' ' +  str(intersection[2]) + '\n'
        f.write(line)
    f.close()
   
    f = open(filename_nodesresolved, "w")
    for n in range(0, len(nodes_resolved_idxs)):
        line = str(nodes_resolved_idxs[n]) + '\n'
        f.write(line)
    f.close()
   
    f = open(filename_epsilonnodes, "w")
    for n in range(0, len(epsilon_nodes)):
        epsilon_node = epsilon_nodes[n]
        line = str(epsilon_node[0]) + ' ' +  str(epsilon_node[1]) + ' ' +  str(epsilon_node[2]) + '\n'
        f.write(line)
    f.close()
   
    f = open(filename_nodepaths, "w")
    for n in range(0, len(node_paths_elem)):
        node_path = node_paths_elem[n]
        line = ""
        for m in range(0, len(node_path)):
            line += str(node_path[m]) + ' '
   
        line += '\n'
        f.write(line)
    f.close()

    f = open(filename_nodecoords, "w")
    for n in range(0, len(node_paths_elem)):
        node_path_coord = node_paths_coord[n]
        line = ""
        for m in range(0, len(node_path_coord)):
            coord = node_path_coord[m]
            for l in range(0, len(coord)):
                line += str(coord[l]) + ' '

            line += ","

        line += '\n'
        f.write(line)
    f.close()

    return True
