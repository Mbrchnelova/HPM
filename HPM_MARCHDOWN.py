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
from HPM_marchdown_integrate import *

from HPM_SETUP import *
from HPM_NAMECONFIG import *



if VISCOUS:
















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

