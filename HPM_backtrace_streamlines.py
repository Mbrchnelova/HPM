# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_backtrace_streamlines.py: this module uses the inviscid data (velocity field) to backtrace the streamlines from all nodes back to the epsilon curve.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be




from HPM_import import *
from HPM_mathematics import *
from HPM_numerics import *




def TraceStreamlinesBack(max_x_user, nodes, nodes_relevance, us, vs, ws, dt_end, stag_node, epsilon, us_mesh, vs_mesh, ws_mesh, centroids, centroids_ordered, centroids_ordered_dxs, normals, connectivity, xmin, xmax, ymin, ymax, zmin, zmax, VERIFY, PRINT):
    node_paths_coord = []
    node_paths_elem = []
    nodes_resolved = []
    nodes_resolved_idxs = []
    intersections = []
    streambacked_nodes = []
    #iter_number = int(length / (min(us)*dt))
    threshold_vel = 0.
    #nodes = nodes[:33]
    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Backtracing streamlines to SP from nodes')
        ax.scatter(stag_node[0], stag_node[1], stag_node[2], s=10, c = 'g')

    #nodes = nodes[:17]
    for n in range(0, len(nodes)):
        s_to_travel = epsilon/2.
        #if VERIFY:                                                              
        #    ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=1, c = 'b')
        RESOLVED = False
        FAILED = False
        if PRINT:
            print 'Node: ', n
        elif n % 100 == 0:
            print 'Node: ', n
        local_u_vector = [us[n], vs[n], ws[n]]

        coordinates = copy.deepcopy(nodes[n])
        dist_to_stag_point = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
        nit_crit = int(dist_to_stag_point*20000)
        nit_crit = dist_to_stag_point / s_to_travel * 2
        if PRINT:
            print dist_to_stag_point, coordinates, stag_node, local_u_vector


        if local_u_vector[0] < threshold_vel and local_u_vector[1] < threshold_vel and local_u_vector[2] < threshold_vel:
            if PRINT:
                print 'Node in shadow: ', n                                         #THESE ARE ALL FAILURE OF TRACING POSSIBILITIES
            #if VERIFY:
            #    ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=4, c = 'r')
        elif nodes_relevance[n] == 0:
            if PRINT:
                print 'Node irrelevant: ', n                                        #THESE ARE ALL FAILURE OF TRACING POSSIBILITIES
            #if VERIFY:
            #    ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=4, c = 'r')



        elif nodes[n][0] > max_x_user:
            if PRINT:
                print 'Node irrelevant: ', n                                        #THESE ARE ALL FAILURE OF TRACING POSSIBILITIES
            #if VERIFY:
            #    ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=4, c = 'r')


        else:
            crossed = CrossedStagnationLine(coordinates, stag_node, epsilon)[0]
            if dist_to_stag_point < epsilon:
                path_coord = [[stag_node[0], stag_node[1], stag_node[2]]]
                local_u_vector, element = FindVelVectorNotOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids)
                path_element = [element]                                        #THESE ALL INDICATE THAT THE POINT ALREADY IS IN A STAGNATION REGION

            elif crossed == True:
                path_coord = [[coordinates[0], coordinates[1], coordinates[1]]]
                local_u_vector, element = FindVelVectorNotOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids)
                path_element = [element]                                        #THESE ALL INDICATE THAT THE POINT ALREADY IS IN A STAGNATION REGION                                          #THESE ALL INDICATE THAT THE POINT ALREADY IS IN A STAGNATION REGION
            else:
                local_u_vector, element = FindVelVectorNotOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids)
                path_coord = [[coordinates[0], coordinates[1], coordinates[2]]]
                path_element = [element]

                crossed = CrossedStagnationLine(coordinates, stag_node, epsilon)[0]
                it = 0.
                dt = 1.*dt_end
                dist_to_stag_point_previous = dist_to_stag_point
                while crossed == False:
                    if dist_to_stag_point_previous > 2.*epsilon:
                        s_to_travel = epsilon
                    else:
                        s_to_travel = epsilon/2.

                    '''
                    if FIRST_LOOP:
                        dist_to_stag_point_previous = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
                        dist_to_stag_point_previous_x = (coordinates[0] - stag_node[0])
                        dist_to_stag_point_previous_y = (coordinates[1] - stag_node[1])
                        dist_to_stag_point_previous_z = (coordinates[2] - stag_node[2])
                        #dt = dist_to_stag_point_previous / us[n] 
                        FIRST_LOOP = False
                    else:
                        element = path_element[-1]
                        dist_to_stag_point_previous = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
                        dist_to_stag_point_previous_x = (coordinates[0] - stag_node[0])
                        dist_to_stag_point_previous_y = (coordinates[1] - stag_node[1])
                        dist_to_stag_point_previous_z = (coordinates[2] - stag_node[2])
                        
                        #dt = DetermineDt(coordinates, local_u_vector, element, connectivity, centroids, normals)
                        #dt = DetermineDt(coordinates, local_u_vector, element, connectivity)
                        #print dt
                    '''

                    #print  (local_u_vector)
                    '''
                    if max(abs_vector) != 0.:
                        if max(abs_vector) == abs_vector[0]:
                            dt = abs(dist_to_stag_point_previous_x) / (abs_vector[0])*0.25
                        if max(abs_vector) == abs_vector[1]:
                            dt = abs(dist_to_stag_point_previous_y) / (abs_vector[1])*0.25
                        if max(abs_vector) == abs_vector[2]:
                            dt = abs(dist_to_stag_point_previous_z) / (abs_vector[2])*0.25
                        
                    else:
                        dt = 0.001* dt
                    '''
                    #dt = dt_end*20.
                    total_u = (local_u_vector[0]**2  + local_u_vector[1]**2  + local_u_vector[2]**2)**0.5
                    if total_u == 0.:
                        dt = 0.00001
                    else:
                        dt = s_to_travel / total_u

                    '''
                    if dist_to_stag_point_previous < 2.0 * dt * max(local_u_vector):
                        dt = 1.*dt_end
                        #dt = dist_to_stag_point_previous / max(local_u_vector)
                        
                    if dist_to_stag_point_previous < 1.0 * dt * max(local_u_vector):
                        dt = 0.25*dt_end
                        #dt = dist_to_stag_point_previous / max(local_u_vector)
                        
                    if dist_to_stag_point_previous < 0.1 * dt * max(local_u_vector):
                        dt = 0.1*dt_end
                        #dt = dist_to_stag_point_previous / max(local_u_vector)
                        
                    if dist_to_stag_point_previous < 0.01 * dt * max(local_u_vector):
                        dt = 0.01*dt_end
                        #dt = dist_to_stag_point_previous / max(local_u_vector)
                
                    '''
                    coordinates[0] = coordinates[0] - local_u_vector[0]*dt
                    coordinates[1] = coordinates[1] - local_u_vector[1]*dt
                    coordinates[2] = coordinates[2] - local_u_vector[2]*dt



                    LOC = 0

                    min_loc = len(centroids_ordered_dxs) - 2
                    max_loc = len(centroids_ordered_dxs) - 1

                    for i in range(0, len(centroids_ordered_dxs)):
                        if coordinates[0] > centroids_ordered_dxs[i] and not LOC:
                            LOC = 1
                            min_loc = i-1
                            max_loc = i

                    centroids_ordered_local = centroids_ordered[min_loc] + centroids_ordered[max_loc]

                    local_u_vector, element = FindVelVectorOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids_ordered_local, centroids)

                    normal = normals[element]
                    centroid = centroids[element]

                    normal = [normal[0]/Mag(normal), normal[1]/Mag(normal), normal[2]/Mag(normal)]
                    diff = [coordinates[0] - centroid[0], coordinates[1] - centroid[1], coordinates[2] - centroid[2]]
                    dot_scalar = Dot(diff, normal)
                    element_proj = [coordinates[0] - dot_scalar * normal[0], coordinates[1] - dot_scalar * normal[1], coordinates[2] - dot_scalar * normal[2] ]
                    coordinates = element_proj

                    if VERIFY:
                        if nodes[n][0] < 1.0:                                  #Only these streamlines will be shown for clarity
                            ax.scatter(coordinates[0], coordinates[1], coordinates[2], s=2, c = 'k')


                    crossed = CrossedStagnationLine(coordinates, stag_node, epsilon)[0]
                    if crossed:
                        dist_to_stag_point_current = CrossedStagnationLine(coordinates, stag_node, epsilon)[1]
                        #coordinates = path_coord[-1]
                        dist_to_stag_point_previous = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
                        dt = dt * dist_to_stag_point_previous / (dist_to_stag_point_previous + dist_to_stag_point_current)

                        normal = normals[element]
                        centroid = centroids[element]

                        normal = [normal[0]/Mag(normal), normal[1]/Mag(normal), normal[2]/Mag(normal)]
                        diff = [coordinates[0] - centroid[0], coordinates[1] - centroid[1], coordinates[2] - centroid[2]]
                        dot_scalar = Dot(diff, normal)
                        element_proj = [coordinates[0] - dot_scalar * normal[0], coordinates[1] - dot_scalar * normal[1], coordinates[2] - dot_scalar * normal[2] ]
                        coordinates = element_proj

                        if VERIFY:
                            ax.scatter(coordinates[0], coordinates[1], coordinates[2], s=4, c = 'g')

                        nodes_resolved.append(coordinates)
                        nodes_resolved_idxs.append(n)
                        path_element.append(element)
                        path_coord.append([coordinates[0], coordinates[1], coordinates[2]])
                        RESOLVED = True
                        break


                    path_element.append(element)
                    path_coord.append([coordinates[0], coordinates[1], coordinates[2]])
                    it += 1.

                    if coordinates[0] < xmin or coordinates[0] > xmax or coordinates[1] < ymin or coordinates[1] > ymax or coordinates[2] < zmin or coordinates[2] > zmax:
                        FAILED = True
                        if PRINT:
                            print 'Out of bounds: ', coordinates,
                        break
                    if it > nit_crit:
                        FAILED = True
                        if PRINT:
                            print 'Too many iterations: ',
                        break
                if FAILED:
                    if PRINT:
                        print 'Could not find tracing for node:', n
                    if VERIFY:
                        ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=4, c = 'r')
            if RESOLVED:
                node_paths_coord.append(path_coord)
                node_paths_elem.append(path_element)
                intersections.append([coordinates[0], coordinates[1], coordinates[2]])                    #THIS WILL BE WRONG, YOU HAVE TO REDO THIS AFTER YOU TEST THE REST
                streambacked_nodes.append(nodes[n])

    return intersections, node_paths_elem, node_paths_coord, nodes_resolved, nodes_resolved_idxs


