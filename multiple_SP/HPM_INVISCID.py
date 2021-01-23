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
from HPM_lep_functions import *

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

    print 'Inviscid Solver stop time:', time.time()

    stag_point_x = 0
    stag_point_y = 0
    stag_point_z = 0
   
    print 'Inviscid Solver LEP start time:', time.time()

    for i in range(0, len(stag_points)):
        stag_point_x += stag_points[i][0]
        stag_point_y += stag_points[i][1]
        stag_point_z += stag_points[i][2]
   
    stag_point = [stag_point_x, stag_point_y, stag_point_z]
    
    SPs, xLEPs, yLEPs, zLEPs, SPthermo, xLEPThermo, yLEPThermo, zLEPThermo = LocalizeSPandLEP(centroids, normals, crit_angle, connectivity,
                                                                                              u_velocities, v_velocities, w_velocities, nodes,
                                                                                              T_nodes, p_nodes, mu_nodes, rho_nodes, u_nodes,
                                                                                              v_nodes, w_nodes, M_nodes, nosetype, Mach, beta, R_N)
   
    filtered_SPs, SPthermo = FilterOutSPs(SPs, stag_point, epsilon, SPthermo)
    filtered_xLEPs, xLEPThermo = FilterOutSPs(xLEPs, stag_point, epsilon, xLEPThermo)
    filtered_yLEPs, yLEPThermo = FilterOutSPs(yLEPs, stag_point, epsilon, yLEPThermo)
    filtered_zLEPs, zLEPThermo = FilterOutSPs(zLEPs, stag_point, epsilon, zLEPThermo)
   
    LEPs_all = filtered_xLEPs + filtered_yLEPs + filtered_zLEPs
    Thermo_all = xLEPThermo + yLEPThermo + zLEPThermo
   
    if VERIFY_LEPS:
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(Cps)
    
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=2, c = c, marker = 'x')
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Detected LEPs and SPs.')
        #ax.axes.set_aspect('equal')

        for i in range(0, len(filtered_xLEPs)):
            ax.scatter(filtered_xLEPs[i][0], filtered_xLEPs[i][1], filtered_xLEPs[i][2], s=8, c = 'r')

        for i in range(0, len(filtered_zLEPs)):
            ax.scatter(filtered_zLEPs[i][0], filtered_zLEPs[i][1], filtered_zLEPs[i][2], s=8, c = 'r')

        for i in range(0, len(filtered_yLEPs)):
            ax.scatter(filtered_yLEPs[i][0], filtered_yLEPs[i][1], filtered_yLEPs[i][2], s=8, c = 'r')

        for i in range(0, len(filtered_SPs)):
            ax.scatter(filtered_SPs[i][0], filtered_SPs[i][1], filtered_SPs[i][2], s=20, c = 'k')

    print len(LEPs_all), "LEP found. Check on the Figure whether they are relevant. Computing their thermodynamics..."



   
    print 'Stagnation point found at,', stag_points, ' with stagnation pressure of', stag_ps
    #print 'Backtracing the streamlines...'

    #intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs = TraceStreamlinesBack(nodes, u_nodes, v_nodes, w_nodes, dt, stag_point, stag_idx, epsilon, u_velocities, v_velocities, w_velocities, centroids, normals, VERIFY_BACKTRACING)
    print 'Inviscid Solver LEP stop time:', time.time()

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
    WRITTEN = WriteLEPointData(filename_le_point, LEPs_all, epsilon_LEP)
    WRITTEN = WriteLEPThermoData(filename_le_thermal, LEPs_all, epsilon_LEP, Thermo_all)
    WRITTEN = WriteNormals(filename_normals, normals)
    WRITTEN = WriteDataCentroids(filename_centroids, filename_connectivity, centroids, connectivity)
