from HPM_SETUP import *
from HPM_write import *

#########OUTPUT FILE MANAGEMENT

filename_infofile = str(path + 'INFO.txt')

filename_normals = str(path + 'normals.txt')
filename_centroids = str(path + 'centroids.txt')
filename_connectivity = str(path + 'connectivity.txt')
filename_centroids = str(path + 'centroids.txt')

filename_p_out = str(path + 'p_nodes.txt')
filename_u_out = str(path + 'u_nodes.txt')
filename_v_out = str(path + 'v_nodes.txt')
filename_w_out = str(path + 'w_nodes.txt')
filename_T_out = str(path + 'T_nodes.txt')
filename_rho_out = str(path + 'rho_nodes.txt')
filename_M_out = str(path + 'M_nodes.txt')
filename_mu_out = str(path + 'mu_nodes.txt')

filename_u_ele_out = str(path + 'u_elements.txt')
filename_v_ele_out = str(path + 'v_elements.txt')
filename_w_ele_out = str(path + 'w_elements.txt')
filename_T_ele_out = str(path + 'T_elements.txt')
filename_rho_ele_out = str(path + 'rho_elements.txt')
filename_mu_ele_out = str(path + 'mu_elements.txt')
filename_M_ele_out = str(path + 'M_elements.txt')

filename_stag_point = str(path + 'stg_point.txt')

filename_intersections = str(path + 'BTR_intersections.txt')
filename_nodepaths = str(path + 'BTR_paths.txt')
filename_nodecoords = str(path + 'BTR_nodes_coord.txt')
filename_nodesresolved = str(path + 'BTR_nodes_resolved.txt')
filename_epsilonnodes = str(path + 'BTR_epsilon_nodes.txt')

filename_dotq_out = str(path + 'dotq_nodes_out.txt')
filename_theta_out = str(path + 'theta_nodes_out.txt')
filename_skincf_out = str(path + 'skincf_nodes_out.txt')

data_to_write_info = [pinf, Tinf, Mach, rhoinf, muinf, R, cp, gamma, name,
                      USE_IDS, IDWP, COMPUTE_STAGNATION_POINT_FROM_GEOMETRY,
                      epsilon, T_w, PERZHIKAR, INTERPOLATE_N, beta_tol,
                      OWN_TRANSITION, max_x_user, dist_yz_min]
WRITTEN = WriteRunInfo(filename_infofile, data_to_write_info)

