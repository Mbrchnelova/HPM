from HPM_import import *
from HPM_process_mesh import *
from HPM_get_thermodynamics import *
from HPM_vertex_values import *
from HPM_get_velocities import *
from HPM_newton_method import *
from HPM_stagnation import *
from HPM_write import *
from HPM_read import *
from HPM_backtrace_streamlines import *
from HPM_detect_irrelevant import *
from HPM_compute_parzhikar import *
from HPM_marchdown_integrate import *

from HPM_SETUP import *
from HPM_NAMECONFIG import *








if VISCOUS:
    print 'Viscous Solver backtrace start time:', time.time()

    connectivity = ReadConnectivity(filename_connectivity)
    centroids = ReadCentroids(filename_centroids)
    normals =  ReadNormals(filename_normals)

    p_nodes, nodes =   ReadInviscidData(filename_p_out)
    u_nodes, nodes =   ReadInviscidData(filename_u_out)
    v_nodes, nodes =   ReadInviscidData(filename_v_out)
    w_nodes, nodes =   ReadInviscidData(filename_w_out)
    T_nodes, nodes =   ReadInviscidData(filename_T_out)
    rho_nodes, nodes = ReadInviscidData(filename_rho_out)
    mu_nodes, nodes =  ReadInviscidData(filename_mu_out)
    M_nodes, nodes =   ReadInviscidData(filename_M_out)

    u_velocities, centroids = ReadInviscidData(filename_u_ele_out)
    v_velocities, centroids = ReadInviscidData(filename_v_ele_out)
    w_velocities, centroids = ReadInviscidData(filename_w_ele_out)
    T_elements, centroids =   ReadInviscidData(filename_T_ele_out)
    mu_elements, centroids =  ReadInviscidData(filename_mu_ele_out)
    rho_elements, centroids = ReadInviscidData(filename_rho_ele_out)
    M_elements, centroids =   ReadInviscidData(filename_M_ele_out)
    u_elements = u_velocities
    v_elements = v_velocities
    w_elements = w_velocities
   
    stag_points, stg_idxs, epss = ReadStagPointData(filename_stag_point)
    if len(stag_points) == 1:
        stag_point = stag_points[0]
        stag_idx = stg_idxs[0]
        epsilon = epss[0]
    if len(stag_points) > 1:
        stag_point_x = 0.
        stag_point_y = 0.
        stag_point_z = 0.
        for i in range(0, len(stag_points)):
            stag_point_x += stag_points[i][0]
            stag_point_y += stag_points[i][1]
            stag_point_z += stag_points[i][2]
        stag_point = [stag_point_x/len(stag_points), stag_point_y/len(stag_points), stag_point_z/len(stag_points)]
        stag_idx = stg_idxs[0]
        epsilon = epss[0]


  
    if BACKTRACE:

        nodes_relevance = DetectIrrelevantNodes(nodes, p_nodes, pinf, Mach, gamma, VERIFY_NODE_RELEVANCE)
        print "Finding intersections with epsilon line..."
        #epsilon = 2. * epsilon
        xmin, xmax, ymin, ymax, zmin, zmax = FindBounds(nodes)
    
        centroids_ordered, centroids_ordered_dxs = SeparateCentroidsByDx(centroids, xmin, xmax, 2)
        PRINT_BTR = True
        intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs = TraceStreamlinesBack(max_x_user, nodes, nodes_relevance, u_nodes, v_nodes, w_nodes, dt, stag_point, epsilon, u_velocities, v_velocities, w_velocities, centroids, centroids_ordered, centroids_ordered_dxs, normals, connectivity, xmin, xmax, ymin, ymax, zmin, zmax, VERIFY_BACKTRACING, PRINT_BTR)
        print "Writing intersection data..."
    
        print 'Viscous Solver backtrace stop time:', time.time()

    if WRITE_BACKTRACING:
        WRITTEN = WriteBackTracing(filename_intersections, filename_nodepaths, filename_nodecoords, filename_nodesresolved, filename_epsilonnodes, intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs)

    else:
        print "Reading intersection data..."
        intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs = ReadBacktracingData(filename_intersections, filename_nodepaths, filename_nodecoords, filename_nodesresolved, filename_epsilonnodes)




    if not PERZHIKAR:
        filename_null = path + "nullspace.txt"
        #filename_p = path + "p_nodes_fitted.txt"
        filename_u = path + "u_nodes_fitted.txt"
        filename_v = path + "v_nodes_fitted.txt"
        filename_w = path + "w_nodes_fitted.txt"
        #filename_T = path + "T_nodes_fitted.txt"
        #filename_rho = path + "rho_nodes_fitted.txt"
        #filename_M = path + "M_nodes_fitted.txt"
        #filename_mu = path + "mu_nodes_fitted.txt"


        null_spaces = ReadNullSpace(filename_null)
        #p_fitted = ReadSolvedFitting(filename_p)
        u_fitted = ReadSolvedFitting(filename_u)
        v_fitted = ReadSolvedFitting(filename_v)
        w_fitted = ReadSolvedFitting(filename_w)
        #T_fitted = ReadSolvedFitting(filename_T)
        #rho_fitted = ReadSolvedFitting(filename_rho)
        #M_fitted = ReadSolvedFitting(filename_M)
        #mu_fitted = ReadSolvedFitting(filename_mu)

    else:
        u_fitted = []
        v_fitted = []
        w_fitted = []
        for i in range(0, len(nodes)):
            u_fitted.append([0.,0.,0.,0.,0.,0.])
            v_fitted.append([0.,0.,0.,0.,0.,0.])
            w_fitted.append([0.,0.,0.,0.,0.,0.])

        null_spaces = []
        for i in range(0, len(centroids)):
            null_spaces.append([0.,0.,0.,0.,0.,0.])


    wd1s = []
    wd2s = []
    wd3s = []
    wd4s = []
    wd5s = []
    wd6s = []
    for i in range(0, len(w_fitted)):
        wd1s.append(w_fitted[i][0])
        wd2s.append(w_fitted[i][1])
        wd3s.append(w_fitted[i][2])
        wd4s.append(w_fitted[i][3])
        wd5s.append(w_fitted[i][4])
        wd6s.append(w_fitted[i][5])

    wds = [wd1s, wd2s, wd3s, wd4s, wd5s, wd6s]

    Fxs = []
    GradFs = []
   
  
    if not PERZHIKAR:
        for c in range(0, len(centroids)):

            a_vector = null_spaces[c]
            centroid = centroids[c]
            x = centroid[0]
            y = centroid[1]
            z = centroid[2]
            Fx = 2. * a_vector[4] * x + a_vector[5]
            Fy = 2. * a_vector[0] * y + a_vector[1]
            Fz = 2. * a_vector[2] * z + a_vector[3]

            GradF = (Fx**2 + Fy**2 + Fz**2)**0.5

            Fxs.append(Fx)

            GradFs.append(GradF)

        if VERIFY_A_COMPUTATION:
            errors = []
            xs_err = []
            ys_err = []
            zs_err = []


            for c in range(0, len(centroids)):
                a_vector = null_spaces[c]
                x = centroids[c][0]
                y = centroids[c][1]
                z = centroids[c][2]

                u = u_elements[c]
                v = v_elements[c]
                w = w_elements[c]
                centroid = centroids[c]
                Fx = 2. * a_vector[4] * centroid[0] + a_vector[5]
                Fy = 2. * a_vector[0] * centroid[1] + a_vector[1]
                Fz = 2. * a_vector[2] * centroid[2] + a_vector[3]
                if abs(Fx) > 0.:
                    u_th = (v * Fy + w * Fz )/  (-Fx)
                else:
                    u_th = 0.
                #print "Error in a vector: ", (u - u_th)/u_elements[c]]
                if u_elements[c] > 0.:
                    err = (u - u_th)/u_elements[c]
                else:
                    err = 0.
                if -1. < err < 1.:
                    xs_err.append(x)
                    ys_err.append(y)
                    zs_err.append(z)
                    if u_elements[c] > 0.:
                        errors.append((u - u_th)/u_elements[c])
                    else:
                        errors.append(0.)
            fig = plt.figure()
            ax = Axes3D(fig)
            c = np.array(errors)
            cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
            p = ax.scatter(xs_err, ys_err, zs_err, s=15, c = c)
            fig.colorbar(p, cax = cax, orientation = 'vertical')
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            ax.set_title('Errors in the matrix inversion')



        x_coords_elefourth = []
        y_coords_elefourth = []
        z_coords_elefourth = []
        GradFsfourth = []
        Fxsfourth = []
        for i in range(0, len(centroids)):
            x = centroids[i][0]
            y = centroids[i][1]
            z = centroids[i][2]

            #if x > -0.5:
            #    if y > 0.0:
            #        if z > 0.0:
            x_coords_elefourth.append(x)
            y_coords_elefourth.append(y)
            z_coords_elefourth.append(z)

            if GradFs[i] > 30:
                GradFs[i] = 30.
            if GradFs[i] < -30:
                GradFs[i] = -30.
            if Fxs[i] > 1:
                Fxs[i] = 1.
            if Fxs[i] < -1:
                Fxs[i] = -1.

            GradFsfourth.append(GradFs[i])
            Fxsfourth.append(Fxs[i])

        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(Fxsfourth)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_elefourth, y_coords_elefourth, z_coords_elefourth, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Fxs on a cone using modified Newton technique')

        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(GradFsfourth)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_elefourth, y_coords_elefourth, z_coords_elefourth, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('GradFs on a cone using modified Newton technique')


    else:
        for i in range(0, len(centroids)):
            GradFs.append(0.0)
            Fxs.append(0.0)




    if not PERZHIKAR:
        print 'Performing epsilon projecton and calculating matric coefficients at the stagnation region...'

        rn_s, betas = EpsilonCircTransformation(epsilon_nodes)
        hs, ss = FindInitialMetricCoefficientsCircle(betas, epsilon, epsilon_nodes, beta_tol, VERIFY_HS_CLUTTERING)
        wdydbetazs = FindInitial_wdydbetazStable(epsilon_nodes, nodes, T_w, T_nodes, mu_nodes, rho_nodes, M_nodes, centroids, u_elements, v_elements, w_elements , hs, ss, cp, k, R, gamma, stag_idx)


    else:
        hs = []
        ss = []
        rn_s = []
        betas = []
        wdydbetazs = []
        for i in range(0, len(epsilon_nodes)):
            hs.append(epsilon)
            ss.append(0.0)
            wdydbetazs.append(0.0)
            rn_s.append(0.0)
            betas.append(0.0)



    final_xs = []
    final_ys = []
    final_zs = []
    final_values_thetas = []
    final_values_dotqws = []
    final_values_hs = []
    final_values_Ves = []
    final_values_rhoes = []
    final_values_rhostars = []
    final_values_mustars = []
    final_values_integs = []
    final_values_skincfs = []

    print "Computing SP heat flux..."
    Tstag_behindNSW, rhostag_behindNSW, pstag_behindNSW, mustag_behindNSW, Vstag_behindNSW = GetStagPointConditions(Mach, pinf, rhoinf, Tinf)
    stg_cnds = [Tstag_behindNSW, rhostag_behindNSW, pstag_behindNSW, mustag_behindNSW, Vstag_behindNSW]
    qstag = GetStagnationHeating(stg_cnds, T_w, 1.0, cp, pinf, epsilon)



    if PERZHIKAR and MARCH:
        print 'Viscous Solver Parzhikar start time:', time.time()

        print "Pre-computing Perzhikar's metric coefficients..."
        metriccoeffs = ComputeMetricPerzhikar(nodes, nodes_resolved_idxs, node_paths_coord, stag_point, dist_yz_min)

        print 'Viscous Solver Parzhikar stop time:', time.time()

    else:
        metriccoeffs = []
        for i in range(0, len(epsilon_nodes)):
            metriccoeffs.append(0.0)



    if MARCH:
        print 'Marching forward and computing solution...'

        print 'Viscous Solver marchdown start time:', time.time()

        for s in range(0, len(epsilon_nodes)):
            metric_coeffs_current = copy.deepcopy(list(reversed(metriccoeffs[s])))
            h_ini = hs[s]
            wdydbetaz_ini = wdydbetazs[s]
            list_of_crossed_elements = copy.deepcopy(list(reversed(node_paths_elem[s])))
            list_of_coordinates = copy.deepcopy(list(reversed(node_paths_coord[s])))

            hs_res, hs_coords, thetas, dotq_ws, deltas, transition_points, V_es, rho_es, mu_stars, rho_stars, integs, dss, skin_cfs = MarchDownstream(stag_point, stg_cnds, wdydbetaz_ini, list_of_crossed_elements, list_of_coordinates, h_ini, GradFs, Fxs, u_elements, v_elements, w_elements, T_elements, mu_elements, rho_elements, M_elements, w_fitted, T_w, cp, k, wds, metric_coeffs_current, INTERPOLATE_N, OWN_TRANSITION, PERZHIKAR)

            if len(hs_coords) > 0.:
                y_hs = []
                x_hs = []
                z_hs = []

                if hs_coords[-1][0] > -1.:
                    for i in range(0, len(hs_coords)):
                        x_hs.append(hs_coords[i][0])
                        y_hs.append(hs_coords[i][1])
                        z_hs.append(hs_coords[i][2])

                    #c = np.array(dss)  


                    final_xs.append(x_hs[-1])
                    final_ys.append(y_hs[-1])
                    final_zs.append(z_hs[-1])
                    final_values_thetas.append(thetas[-1])
                    final_values_dotqws.append(dotq_ws[-1])
                    final_values_hs.append(hs_res[-1])
                    final_values_Ves.append(V_es[-1])
                    final_values_rhoes.append(rho_es[-1])
                    final_values_rhostars.append(rho_stars[-1])
                    final_values_mustars.append(mu_stars[-1])
                    final_values_integs.append(integs[-1])
                    final_values_skincfs.append(skin_cfs[-1])

                    final_xs.append(x_hs[0])
                    final_ys.append(y_hs[0])
                    final_zs.append(z_hs[0])
                    final_values_thetas.append(thetas[0])
                    final_values_dotqws.append(dotq_ws[0])
                    final_values_hs.append(hs_res[0])
                    final_values_Ves.append(V_es[0])
                    final_values_rhoes.append(rho_es[0])
                    final_values_rhostars.append(rho_stars[0])
                    final_values_mustars.append(mu_stars[0])
                    final_values_integs.append(integs[0])
                    final_values_skincfs.append(skin_cfs[0])

                    #p = ax.scatter(x_hs, y_hs, z_hs, s=5, c = c)

    print 'Viscous Solver marchdown stop time:', time.time()

    final_values_dotqws.append(qstag/10000.)
    final_xs_q = copy.deepcopy(final_xs)
    final_ys_q = copy.deepcopy(final_ys)
    final_zs_q = copy.deepcopy(final_zs)
    final_xs_q.append(stag_point[0])
    final_ys_q.append(stag_point[1])
    final_zs_q.append(stag_point[2])


    fig = plt.figure()
    ax = Axes3D(fig)
    cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])

    #for theta in range(0, len(final_values_thetas)):
    #    if final_values_thetas[theta] > 0.005:
    #        final_values_thetas[theta] = 0.00
    c = np.array(final_values_thetas)
    p = ax.scatter(final_xs, final_ys, final_zs, s=5, c = c)
    fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Thetas along a cone, final values only, m')


    fig = plt.figure()
    ax = Axes3D(fig)
    cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])

    for k in range(0, len(final_xs_q)-1):
        dist_to_stag_point = ((final_xs_q[k] - final_xs_q[-1])**2 + (final_ys_q[k] - final_ys_q[-1])**2 + (final_zs_q[k] - final_zs_q[-1])**2)**0.5
        if final_values_dotqws[k] > qstag:
            final_values_dotqws[k] = qstag
        if dist_to_stag_point < 1. * epsilon:
            final_values_dotqws[k] = qstag / 10000.
        else:
            final_values_dotqws[k] /= 10000.


    c = np.array(final_values_dotqws)
    p = ax.scatter(final_xs_q, final_ys_q, final_zs_q, s=5, c = c)
    ax.scatter(final_xs_q[-1], final_ys_q[-1], final_zs_q[-1], s=100, c = 'r')
    fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Dotqws along a cone, final values only, W/cm2')

    fig = plt.figure()
    ax = Axes3D(fig)
    cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])

    c = np.array(final_values_mustars)
    p = ax.scatter(final_xs, final_ys, final_zs, s=5, c = c)
    fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Mu star along a cone, final values only')

    fig = plt.figure()
    ax = Axes3D(fig)
    cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])



    c = np.array(final_values_integs)
    p = ax.scatter(final_xs, final_ys, final_zs, s=5, c = c)
    fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Theta integral term along a cone, final values only')


    nodes_final_dotq = []
    for i in range(0, len(final_xs_q)):
        nodes_final_dotq.append([final_xs_q[i], final_ys_q[i], final_zs_q[i]])

    WRITTEN_Q = WriteDataNodes(filename_dotq_out, final_values_dotqws, nodes_final_dotq)

    nodes_final_theta = []
    for i in range(0, len(final_xs)):
        nodes_final_theta.append([final_xs[i], final_ys[i], final_zs[i]])

    WRITTEN_THETA = WriteDataNodes(filename_theta_out, final_values_thetas, nodes_final_theta)


    nodes_final_skincf = []
    for i in range(0, len(final_xs)):
        nodes_final_skincf.append([final_xs[i], final_ys[i], final_zs[i]])

    WRITTEN_SKINCF = WriteDataNodes(filename_skincf_out, final_values_skincfs, nodes_final_skincf)



