# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_INVISCID.py: this module is responsible for processing the mesh and computing the inviscid pressure distribution and thermodynamics.  
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be

from HPM_import import *
from HPM_process_mesh import *
from HPM_transform_mesh import *
from HPM_get_thermodynamics import *
from HPM_vertex_values import *
from HPM_get_velocities import *
from HPM_newton_method import *
from HPM_stagnation import *
from HPM_write import *
from HPM_read import *

from HPM_SETUP import *
from HPM_NAMECONFIG import *




if INVSCID:
    print 'Processing mesh...'

    print 'Mesh Processor start time:', time.time()
    nodes, connectivity, normals, centroids, areas = process_mesh(name, CHECK_MESH)
    print 'Mesh Processor stop time:', time.time()

    print 'Transforming based on AOA, Beta and Scale...'
    centroids, nodes, normals = TransformByAOA(centroids, nodes, normals, alpha)
    centroids, nodes, normals = TransformByBeta(centroids, nodes, normals, beta)
    centroids, nodes, normals = TransformByScaling(centroids, nodes, normals, scale)

    print 'Inviscid Solver start time:', time.time()

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

    for centroid in centroids:
        x_coords_ele.append(centroid[0])
        y_coords_ele.append(centroid[1])
        z_coords_ele.append(centroid[2])


    print 'Computing pressure coefficients... '
    Cps, shadowed = Newton(nodes, connectivity, normals, centroids, areas, gamma, V_vector, Mach, VERIFY_NEWTON)

    print 'Computing velocities and pressures at the mesh elements...'
    V, p = GetVelocities(V_vector, normals, centroids, Cps, pinf, Mach, Tinf, gamma, R, VERIFY_NEWTON)


    u_velocities = []
    v_velocities = []
    w_velocities = []

    for v in V:
        u_velocities.append(v[0])
        v_velocities.append(v[1])
        w_velocities.append(v[2])


    if USE_IDS:
        print 'Computing U velocities at the vertices by IDW interpolation with p of' , IDWP, '...'
        u_nodes = GetVertexValues(nodes, centroids, u_velocities, IDWP, VERIFY_ISD_NODES)
        print 'Computing V velocities at the vertices...'
        v_nodes = GetVertexValues(nodes, centroids, v_velocities, IDWP, VERIFY_ISD_NODES)
        print 'Computing W velocities at the vertices...'
        w_nodes = GetVertexValues(nodes, centroids, w_velocities, IDWP, VERIFY_ISD_NODES)
        print 'Computing pressures at the vertices...'
        p_nodes = GetVertexValues(nodes, centroids, p, IDWP, VERIFY_ISD_NODES)


    else:
        print 'Computing node velocities and pressures by cell averaging... '
        u_nodes, v_nodes, w_nodes, p_nodes = GetVertexValuesFaster(nodes, centroids, connectivity, u_velocities, v_velocities, w_velocities, p)


    print 'Computing therodynamics at nodes...'
    T_nodes, mu_nodes, rho_nodes, M_nodes = GetThermodynamics(p_nodes, Mach, nodes, gamma, R, muinf, Tinf, rhoinf, VERIFY_THERMODYNAMICS)

    print 'Computing therodynamics at elements...'
    T_elements, mu_elements, rho_elements, M_elements = GetThermodynamics(p, Mach, centroids, gamma, R, muinf, Tinf, rhoinf, VERIFY_THERMODYNAMICS)


    p_elements = p
    u_elements = u_velocities
    v_elements = v_velocities
    w_elements = w_velocities

    print 'Searching for stagnation point...'
    if COMPUTE_STAGNATION_POINT_FROM_GEOMETRY:
        stag_points, stag_ps, stag_idxs = FindStagnationPointFromGeometry(nodes, p_nodes)
    else:
        stag_points, stag_ps, stag_idxs = FindStagnationPointFromPressure(nodes, p_nodes)                 #CURRENTLY ONLY WORKS FOR ONE STAGNATION POINT, MUST BE REDONE IN CASE OF A WINGED GEOMETRY

    if VERIFY_STAG_POINT:
        '''
        fig = plt.figure()
        ax = Axes3D(fig)
        ''' 
    for stag_point in stag_points:
        '''
        if VERIFY_STAG_POINT:
            ax.scatter(stag_point[0], stag_point[1], stag_point[2], s=10, c = 'r')
        '''

        p_nodes, T_nodes, mu_nodes, rho_nodes, u_nodes, v_nodes, w_nodes = AssgnStagnationPoint(stag_point, nodes, Mach, pinf, rhoinf, Tinf, p_nodes, rho_nodes, T_nodes, mu_nodes, u_nodes, v_nodes, w_nodes, VERIFY_STAG_POINT)


    print 'Stagnation point found at,', stag_points, ' with stagnation pressure of', stag_ps
    #print 'Backtracing the streamlines...'

    #intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs = TraceStreamlinesBack(nodes, u_nodes, v_nodes, w_nodes, dt, stag_point, stag_idx, epsilon, u_velocities, v_velocities, w_velocities, centroids, normals, VERIFY_BACKTRACING)

    print 'Inviscid Solver stop time:', time.time()

    print 'Creating surface fits... '




if WRITE_FIELDS:

    WRITTEN = WriteDataForSurfaceFit(str(path + 'velocities.txt'), str(path + 'coordinates.txt'), connectivity, centroids, u_nodes, v_nodes, w_nodes, nodes)

    WRITTEN = WriteDataNodes(filename_p_out, p_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_u_out, u_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_v_out, v_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_w_out, w_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_T_out, T_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_rho_out, rho_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_M_out, M_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_mu_out, mu_nodes, nodes)

    WRITTEN = WriteDataNodes(filename_u_ele_out, u_elements, centroids)
    WRITTEN = WriteDataNodes(filename_v_ele_out, v_elements, centroids)
    WRITTEN = WriteDataNodes(filename_w_ele_out, w_elements, centroids)
    WRITTEN = WriteDataNodes(filename_T_ele_out, T_elements, centroids)
    WRITTEN = WriteDataNodes(filename_rho_ele_out, rho_elements, centroids)
    WRITTEN = WriteDataNodes(filename_mu_ele_out, mu_elements, centroids)
    WRITTEN = WriteDataNodes(filename_M_ele_out, M_elements, centroids)

    WRITTEN = WriteStagPointData(filename_stag_point, stag_points, stag_idxs, epsilon)
    WRITTEN = WriteNormals(filename_normals, normals)
    WRITTEN = WriteDataCentroids(filename_centroids, filename_connectivity, centroids, connectivity)

