# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_read.py: here, all functions for reading of the input data and save data are defined. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be




from HPM_import import *






def ReadNullSpace(filename):
    j = 0.
    null_spaces = []
    with open(filename, "r") as ins:
        nullspace_array = []
        for line in ins:
            line = line.strip()
            line = line.split()
            if len(line) == 1:
                nullspace_array.append(float(line[0]))
            if len(line) > 1:
                null_spaces.append(nullspace_array)
                nullspace_array = []

            j += 1
    null_spaces.append(nullspace_array)

    null_spaces_formatted = []
    for i in range(0, len(null_spaces)):
        null_spaces_formatted.append(null_spaces[i])

    return null_spaces_formatted[1:]



###############################################################################



def ReadInviscidData(filename):
    data = []
    nodes = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])

            nodes.append([x, y, z])
            data.append(float(line[3]))

    return data, nodes



###############################################################################



def ReadCentroids(filename):
    centroids = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])

            centroids.append([x, y, z])

    return centroids



###############################################################################



def ReadNormals(filename):
    normals = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])

            normals.append([x, y, z])

    return normals



###############################################################################



def ReadStagPointData(filename):
    stag_points = []
    stg_idxs = []
    epss = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            stag_point = [float(line[0]), float(line[1]), float(line[2])]
            stg_idx = int(line[3])
            eps = float(line[4])

            stag_points.append(stag_point)
            stg_idxs.append(stg_idx)
            epss.append(eps)

    return stag_points, stg_idxs, epss


def ReadLePointData(filename):
    stag_points = []
    epss = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            stag_point = [float(line[0]), float(line[1]), float(line[2])]
            eps = float(line[3])
            stag_points.append(stag_point)
            epss.append(eps)

    return stag_points, epss


###############################################################################


def ReadLeThermalData(filename):
    thermal = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            th = [float(line[0]), float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6]), float(line[7])]
            thermal.append(th)

    return thermal



###############################################################################



def ReadConnectivity(filename):
    connectivity = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = int(line[0])
            y = int(line[1])
            z = int(line[2])

            connectivity.append([x, y, z])

    return connectivity



###############################################################################



def ReadSolvedFitting(filename):
    j = 0.
    vectors = []
    with open(filename, "r") as ins:
        vector_array = []
        for line in ins:
            line = line.strip()
            line = line.split()
            if len(line) == 1:
                vector_array.append(float(line[0]))
            if len(line) > 1:
                vectors.append(vector_array)
                vector_array = []

            j += 1
    vectors.append(vector_array)

    solutions_formatted = []
    for i in range(0, len(vectors)):
        solutions_formatted.append(vectors[i])

    return solutions_formatted[1:]



###############################################################################



def ReadBacktracingData(filename_intersections, filename_nodepaths, filename_nodecoords, filename_nodesresolved, filename_epsilonnodes, filename_whereended):
    intersections = []
    node_paths_elem = []
    node_paths_coord = []
    epsilon_nodes = []
    nodes_resolved_idxs = []
    where_ended = []

    with open(filename_epsilonnodes, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])

            epsilon_nodes.append([x, y, z])

    with open(filename_intersections, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])

            intersections.append([x, y, z])

    with open(filename_nodesresolved, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            node_idx = int(line[0])

            nodes_resolved_idxs.append(node_idx)


    with open(filename_nodepaths, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            node_path = []
            for i in range(0, len(line)):
                node_path.append(int(line[i]))

            node_paths_elem.append(node_path)


    with open(filename_nodecoords, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            node_path_coord = []
            coord = []
            for i in range(0, len(line)):
                #print line
                if line[i][0] != ',':
                    if line[i] != '' and line[i] != ' ' and line[i] != '\n':
                        coord_cur = float(line[i])
                        coord.append(coord_cur)
                else:
                    if line[i] != '' and line[i] != ' '  and line[i] != '\n':
                        node_path_coord.append(coord)
                        if len(line[i][1:]) > 1:
                            new_coord = line[i][1:]
                            coord = [float(new_coord)]

            node_paths_coord.append(node_path_coord)



    with open(filename_whereended, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            if line[0] == 'SP':
                where_ended.append('SP')
            else:
                end_idx = int(line[0])
                where_ended.append(end_idx)

    return intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs, where_ended


