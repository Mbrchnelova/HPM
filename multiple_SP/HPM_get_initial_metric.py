# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_get_initial_metric.py: this module contains all the function necessary to compute the initial metric coefficients for the Hamilton method.
# Right now, the Hamilton method is not preferred since it takes longer and oftentimes has convergence problems compared to Parzhikar. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *
from HPM_numerics import *
from HPM_mathematics import *



def FindInitialMetricCoefficientsCircle(betas, epsilon, epsilon_nodes, tol, VERIFY):
    hs = []
    ss = []
    sequence = []
    regions = []
    #This function looks over all the betas of the backtraced streamlines and clusters them within a certain beta tolerance
    for i in range(0, len(betas)):
        NEW_REGION = True
        for j in range(0, len(regions)):
            if abs(betas[i] - regions[j]) < tol:                                #If the stream has beta within tolerance of another existing region, do not create a new region
                NEW_REGION = False
   
        if NEW_REGION:                                                          #If there is no region which could take in the given backtraced streamline, create a new region
            regions.append(betas[i])
   
    if VERIFY:
        '''
        xs = []
        ys = []
        for i in range(0, len(regions)):
            x, y = 0.159 * math.cos(regions[i]), 0.159 * math.sin(regions[i])
            xs.append(x)
            ys.append(y)
        plt.scatter(xs, ys, c = 'k')
        plt.show()
        '''

    temp_regions = []
   
    #Create a temporary array to hold region beta information
    for i in range(0, len(regions)):
        temp_regions.append(regions[i])
   
   
    #This loop determines the order of the regions starting from the bottom of the epsilon line
    for i in range(0, len(regions)):
        min_region = 1e9
        for j in range(0, len(temp_regions)):
            region = temp_regions[j]
            if region < min_region:                                             #Finds the smallest beta existing among the remainder of the regions
                min_region = region
                min_regions_index = j

        #temp_regions[j] = 1e10
        temp_regions.remove(min_region)                                         #Removes the corresponding smallest beta region from the temporary region field
        sequence.append(regions.index(min_region))                              #Sequence array adds the index of the smallest best region

    #Initializes the array of metric coefficients
    for i in range(0, len(regions)):
        hs.append(epsilon)
        ss.append(epsilon)


    #This loop calculates the matric coefficient for each region according to the sequence previously determined
    for k in range(1, len(regions)):
        region = regions[sequence[k]]                                           #Take a region according to the found sequence
        region_old = regions[sequence[k-1]]                                     #Take the beta of the previous region considered (first, most bottom region has h of 0)
        delta_Beta = region - region_old                                        #Calculate the difference in beta between the two regions
        delta_Beta = region - regions[sequence[0]]                              #Calculate the difference in beta between the two regions
        delta_S = delta_Beta * epsilon                                          #The distance travelled along the line is equal to the arc length between the two regions
        index = sequence[k]                                                     #Absolute (original) index of the current region is found from the sequence array
        hs[index] = delta_S / delta_Beta                                        #The metric coefficient at this index is then computed
        ss[index] = delta_S
        #print regions[k], hs[index]


    hs_nodes = []
    ss_nodes = []
    #The array above only contains region information, not node information --> this has to be transformed to the node format
    for n in range(0, len(epsilon_nodes)):
        beta_node = betas[n]                                                    #Evaluate the beta of the current node in the right node order
        ASSIGNED = False
        for r in range(0, len(regions)):                                        #Scan over all regions defined above
            if abs(beta_node - regions[r]) < tol:                               #Find the region to which the current node belongs based on the tolerance                
                hs_nodes.append(hs[r])                                          #Assign the correct metric coefficient from the found 
                ss_nodes.append(ss[r])
                ASSIGNED = True
                break
        if not ASSIGNED:
            print "Metric coefficient assignment failure for node: ", n
    return hs_nodes, ss_nodes



###############################################################################



def FindInitial_wdydbetaz(hs, GradFs, Fxs, us, vs, ws, epsilon_nodes, list_of_crossed_elements):
    wdydbetazs = []
    for e in range(0, len(epsilon_nodes)):
        el = list_of_crossed_elements[e][-1]
        #print el
        V = math.sqrt(us[e]**2 + vs[e]**2 + ws[e]**2)
        #print hs[e], V, Fxs[el], GradFs[el]
        wdydbetaz = hs[e] * V * Fxs[el] / abs(GradFs[el])
        wdydbetazs.append(wdydbetaz)

    return wdydbetazs



###############################################################################



def FindInitial_wdydbetazStable(epsilon_nodes, nodes, T_w, T_nodes, mu_nodes, rho_nodes, M_nodes, centroids, u_elements, v_elements, w_elements, hs, ss, cp, k, R, gamma, stg_idx):
    wdydbetazs = []
    for e in range(0, len(epsilon_nodes)):
        node_index = stg_idx
        coordinates = epsilon_nodes[e]
        velocity_e, element = FindVelVectorNotOrdered(coordinates, u_elements, v_elements, w_elements, centroids)
        u_e = Mag(velocity_e)
        T_e = T_nodes[node_index]
        rho_e = rho_nodes[node_index]
        mu_e = mu_nodes[node_index]
        M_e = M_nodes[node_index]
        Pr = mu_e * cp / k
        V_e = math.sqrt(T_e * R * gamma) * M_e
        Re_e = mu_e * V_e / mu_e

        Re_star, rho_star, mu_star, T_star = ReturnEckert(False, T_w, T_e, M_e, gamma, Pr, mu_e, rho_e, Re_e)

        wdydbetaz = rho_star * mu_star * u_e * hs[e]**2. * ss[e]/4.
        wdydbetazs.append(wdydbetaz)

    return wdydbetazs



###############################################################################




def EpsilonCircTransformation(nodes_resolved):
    r_ns = []
    betas = []
    for n in range(0, len(nodes_resolved)):
        node = nodes_resolved[n]
        r_n = math.sqrt(node[1]**2 + node[2]**2)**0.5
        r_ns.append(r_n)
        beta = math.atan(-node[1]/node[2])
        betas.append(beta)

    return r_ns, betas


