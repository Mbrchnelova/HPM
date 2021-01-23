# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_compute_parzhikar.py: this function computes the metric coefficients for the Parzhikar method. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be


from HPM_import import *
from HPM_mathematics import *
from HPM_numerics import *






def ComputeMetricPerzhikar(nodes, node_resolved_idxs, node_path_coordinates_orig, stag_point, dist_yz_min, where_ends):
    node_path_coordinates = copy.deepcopy(node_path_coordinates_orig)
    #for n in range(0, len(node_path_coordinates)):
    #    first_point = nodes[node_resolved_idxs[n]]
    #    node_path_coordinates[n] = [first_point] + node_path_coordinates[n][:]

    total_ss_all = []
    #For each node that is traced
    for n in range(0, len(node_resolved_idxs)):
        if n%10 == 0:
            print "Node... ", n
        total_ss = []
        node_path_coordinates_curr = node_path_coordinates[n]

        #For each step made during tracing
        for m in range(1, len(node_path_coordinates_curr)):
            total_s = 0.

            #For each coordinate it has passed
            for k in range(m, len(node_path_coordinates_curr)):
                x_i = node_path_coordinates_curr[k][0]
                y_i = node_path_coordinates_curr[k][1]
                z_i = node_path_coordinates_curr[k][2]
                x_imin1 = node_path_coordinates_curr[k-1][0]
                y_imin1 = node_path_coordinates_curr[k-1][1]
                z_imin1 = node_path_coordinates_curr[k-1][2]
                total_s += ((x_i - x_imin1)**2 + (y_i - y_imin1)**2 + (z_i - z_imin1)**2)**0.5

            total_ss.append(total_s)
        #print "For node ", nodes[node_resolved_idxs[n]], " with ", len(node_path_coordinates_curr), " path elements, the ss are: ", total_ss
        #print '\n'
        if len(total_ss) == 0.:
            coords_zero_path = nodes[node_resolved_idxs[n]]
            radius = ((coords_zero_path[1]-stag_point[1])**2 + (coords_zero_path[2]-stag_point[2])**2)**0.5
            new_ss_extrapolated = radius

        elif len(total_ss) == 1:
             new_ss_extrapolated = total_ss[0]

        else:
            new_ss_extrapolated = (total_ss[0] - total_ss[1] ) + total_ss[0]

        total_ss = [new_ss_extrapolated] + total_ss[:]
        total_ss_all.append(total_ss)

    hs_perzhikar_all = []
    for n in range(0, len(node_resolved_idxs)):
        if where_ends[n] != 'SP':
            hs_perzhikar = []

            for m in range(0, len(node_path_coordinates[n])):
                hs_perzhikar.append(1.0)

            hs_perzhikar.append(1.0)
            hs_perzhikar_all.append(hs_perzhikar)

        else:
            #if n%100 == 0:
            #    print "Node... cycle b ", n
            node_path_coordinates_curr = node_path_coordinates[n]
            hs_perzhikar = []
            streamline_dist_min = 1e9
            nearest_streamline_order = 0
            nearest_streamline_point_idx = 0

            #if n == 6522:
            #print "We have a streamline starting from node: ", node_path_coordinates_curr[0], "with index: ", n

            for l in range(0, len(node_resolved_idxs)):
                node_coord = nodes[node_resolved_idxs[l]]
                node_coord_curr = node_path_coordinates_curr[0]

                last_point_node_coor = node_path_coordinates[l][-1]
                last_point_node_coord_curr = node_path_coordinates_curr[-1]

                dist_yz = (0.*(node_coord[0] - node_coord_curr[0])**2 + (node_coord[1] - node_coord_curr[1])**2 + (node_coord[2] - node_coord_curr[2])**2)**0.5
                dist_yz += 0.3*((last_point_node_coor[0] - last_point_node_coord_curr[0])**2 + (last_point_node_coor[1] - last_point_node_coord_curr[1])**2 + (last_point_node_coor[2] - last_point_node_coord_curr[2])**2)**0.5

                if dist_yz < streamline_dist_min and dist_yz > dist_yz_min and node_coord[0] >= node_coord_curr[0]:
                    streamline_dist_min = dist_yz
                    nearest_streamline_point_idx = node_resolved_idxs[l]
                    nearest_streamline_order = l

            #At this point, the node index from which the nearest y-z streamline originates is known

            node_path_coordinates_nearest = node_path_coordinates[nearest_streamline_order]
            total_s_nearest = total_ss_all[nearest_streamline_order]
            total_s_current = total_ss_all[n]

            #if n == 6522:
            #print "We found the nearest streamline to originate from: ",  node_path_coordinates_nearest[0], " with index ", l, " and with a y-z distance of ", streamline_dist_min


            #March through the points on the streamline
            for m in range(1, len(node_path_coordinates_curr)):
                current_s = total_s_current[m]
                s_dist_min = 1e9
                nearest_point = 0

                for k in range(0, len(node_path_coordinates_nearest)):
                    s_nearest_at_k = total_s_nearest[k]
                    #print current_s, s_nearest_at_k
                    s_dist  = abs(current_s - s_nearest_at_k)
                    if s_dist < s_dist_min:
                        s_dist_min = s_dist
                        nearest_point = k

                second_nearest_point = 0
                second_s_dist_min = 1e9

                for l in range(0, len(node_path_coordinates_nearest)):
                    s_nearest_at_l = total_s_nearest[l]
                    second_s_dist  = abs(current_s - s_nearest_at_l)
                    if second_s_dist < second_s_dist_min and second_s_dist != s_dist_min:
                        second_s_dist_min = second_s_dist
                        second_nearest_point = l

                #Figure out linearly interpolated point between these two
                p_n_nearest = node_path_coordinates_nearest[nearest_point]
                p_n_secondnearest = node_path_coordinates_nearest[second_nearest_point]

                x_p_n1 = p_n_nearest[0]
                y_p_n1 = p_n_nearest[1]
                z_p_n1 = p_n_nearest[2]
                x_p_n2 = p_n_secondnearest[0]
                y_p_n2 = p_n_secondnearest[1]
                z_p_n2 = p_n_secondnearest[2]

                if s_dist_min == 0.0:
                    s_dist_min = 1.0e-9
                if second_s_dist_min == 0.0:
                    second_s_dist_min = 1.0e-9
                x_p_n = ((1./s_dist_min) * x_p_n1 + (1./ second_s_dist_min) * x_p_n2) / ((1./s_dist_min) + (1./ second_s_dist_min))
                y_p_n = ((1./s_dist_min) * y_p_n1 + (1./ second_s_dist_min) * y_p_n2) / ((1./s_dist_min) + (1./ second_s_dist_min))
                z_p_n = ((1./s_dist_min) * z_p_n1 + (1./ second_s_dist_min) * z_p_n2) / ((1./s_dist_min) + (1./ second_s_dist_min))

                #At this point, we know which of the points, k, on the nearest streamline is the closest one to the point of interest, m


                #p_n = node_path_coordinates_nearest[nearest_point]
                p_n = [x_p_n, y_p_n, z_p_n]
                p_i = node_path_coordinates_curr[m]
                p_imin1 = node_path_coordinates_curr[m-1]

                #if n == 6522:
                #    print "Staring with the second step ", p_i, " after ", p_imin1, " and taking ", p_n, "from the second streamline, with an error of ", s_dist_min

                hbetadbeta = ((p_n[0] - p_i[0])**2 + (p_n[1] - p_i[1])**2 + (p_n[2] - p_i[2])**2)**0.5
                xpipimin1 = ((p_imin1[0] - p_i[0])**2 + (p_imin1[1] - p_i[1])**2 + (p_imin1[2] - p_i[2])**2)**0.5
                xpimin1pn = ((p_imin1[0] - p_n[0])**2 + (p_imin1[1] - p_n[1])**2 + (p_imin1[2] - p_n[2])**2)**0.5

                #print xpipimin1, hbetadbeta

                term = -(xpimin1pn**2) + xpipimin1**2 + hbetadbeta**2
                if hbetadbeta == 0. or xpipimin1 == 0.:
                    print "Error for a node ", node_path_coordinates_curr[0], " using streamline from the node ", nodes[nearest_streamline_point_idx], " at coordinates pi, pimin1 and pn ", p_i, '\t', p_imin1, '\t', p_n
                    print node_path_coordinates_curr
                    print node_path_coordinates_nearest


                if (0.5*term / hbetadbeta / xpipimin1 ) > 1.0:
                    gamma = math.acos(1.0)

                elif (0.5*term / hbetadbeta / xpipimin1) < -1.0:
                    gamma = math.acos(-1.0)
                else:
                    gamma = math.acos(0.5*term / hbetadbeta / xpipimin1 )

                hdbeta = hbetadbeta * math.sin(gamma)
                if m == 1:
                    factor = hbetadbeta / ((p_i[1]-stag_point[1])**2 + (p_i[2]-stag_point[2])**2)**0.5

                if n == 6522:
                    print "Comparison: ", hbetadbeta / factor, '\t' ,((p_i[1]-stag_point[1])**2 + (p_i[2]-stag_point[2])**2)**0.5

                #print "The found angle is ", gamma *180./3.14 , " and the metric ", round(hdbeta, 5), " while the actual metric is ", ((p_i[1]-stag_point[1])**2 + (p_i[2]-stag_point[2])**2)**0.5, "\n \n"

                hs_perzhikar.append(hdbeta)
            hs_perzhikar.append(hdbeta)
            hs_perzhikar_all.append(hs_perzhikar)
    print "Pre-computing of metric coefficients has finished."
    return hs_perzhikar_all




