






def LocalizeSPandLEP(centroids, normals, crit_angle, connectivity, u_velocities,
                     v_velocities, w_velocities, nodes, T_nodes, p_nodes, mu_nodes,
                     rho_nodes, u_nodes, v_nodes, w_nodes, Mach_nodes, nosetype,
                     Minf, beta, R):
    SPs = []
    xLEPs = []
    yLEPs = []
    zLEPs = []

    SPthermo = []
    xLEPthermo = []
    yLEPthermo = []
    zLEPthermo = []

    for c in range(0, len(centroids)):
        if c%100 == 0:
            print "Localizing at c: ", c

        c_c = centroids[c]
        norm_c = normals[c]

        x_crit = -GetShockZ(Minf, R, beta, c_c[1], nosetype)
        if abs(c_c[0]) > x_crit:
            u_c = u_velocities[c]
            v_c = v_velocities[c]
            w_c = w_velocities[c]

            v1 = connectivity[c][0]
            v2 = connectivity[c][1]
            v3 = connectivity[c][2]
            all_v_c = [v1, v2, v3] # this will look like [454, 453, 487] ... 

            c1 = nodes[v1]
            c2 = nodes[v2]
            c3 = nodes[v3]

            c_coords = [c1, c2, c3]

            for d in range(c+1, len(centroids)):
                XLEP_D = False
                YLEP_D = False
                ZLEP_D = False
                SP_D = False

                SIDECONNECTED = False

                v1b = connectivity[d][0]
                v2b = connectivity[d][1]
                v3b = connectivity[d][2]

                d1 = nodes[v1b]
                d2 = nodes[v2b]
                d3 = nodes[v3b]

                d_coords = [d1, d2, d3]

                match = 0

                matching_i = []
                for i in range(0, 3):
                    for j in range(0, 3):
                        if c_coords[i] == d_coords[j]:
                            match += 1
                            matching_i.append(all_v_c[i])
                            #matching_i will look like: [0, 1] or [0, 2] ...
                if match == 2:
                    SIDECONNECTED = True

                if SIDECONNECTED:
                    norm_d = normals[d]
                    d_c = centroids[d]
                    u_d = u_velocities[d]
                    v_d = v_velocities[d]
                    w_d = w_velocities[d]

                    #print "Found neighbouring cells. Velocities: ", round(u_d, 2),  round(u_c, 2),  '\t',  round(v_d, 2),  round(v_c, 2), '\t',  round(w_d, 2),  round(w_c, 2)


                    if u_d * u_c < 0.:
                        angle = math.acos((norm_c[0] * norm_d[0] + norm_c[1] * norm_d[1]  + norm_c[2] * norm_d[2]) / Mag(norm_c) / Mag(norm_d))
                        if angle > crit_angle:
                            XLEP_D = True
                            XLEP = [0.5 * (d_c[0] + c_c[0]), 0.5 * (d_c[1] + c_c[1]), 0.5 * (d_c[2] + c_c[2])]
                            xLEPs.append(XLEP)
                            print "Found XLEP at: ", XLEP

                            T1 = T_nodes[matching_i[0]]
                            T2 = T_nodes[matching_i[1]]
                            T_LEP  = 0.5 * (T1+T2)
                            rho1 = rho_nodes[matching_i[0]]
                            rho2 = rho_nodes[matching_i[1]]
                            rho_LEP  = 0.5 * (rho1+rho2)
                            mu1 = mu_nodes[matching_i[0]]
                            mu2 = mu_nodes[matching_i[1]]
                            mu_LEP  = 0.5 * (mu1+mu2)
                            p1 = p_nodes[matching_i[0]]
                            p2 = p_nodes[matching_i[1]]
                            p_LEP  = 0.5 * (p1+p2)
                            u1 = u_nodes[matching_i[0]]
                            u2 = u_nodes[matching_i[1]]
                            u_LEP  = 0.5 * (u1+u2)
                            v1 = v_nodes[matching_i[0]]
                            v2 = v_nodes[matching_i[1]]
                            v_LEP  = 0.5 * (v1+v2)
                            w1 = w_nodes[matching_i[0]]
                            w2 = w_nodes[matching_i[1]]
                            w_LEP  = 0.5 * (w1+w2)
                            M1 = Mach_nodes[matching_i[0]]
                            M2 = Mach_nodes[matching_i[1]]
                            M_LEP  = 0.5 * (M1+M2)
                            xLEPthermo.append([u_LEP, v_LEP, w_LEP, T_LEP, mu_LEP, rho_LEP, p_LEP, M_LEP])


                    if v_d * v_c < 0.:
                        angle = math.acos((norm_c[0] * norm_d[0] + norm_c[1] * norm_d[1]  + norm_c[2] * norm_d[2]) / Mag(norm_c) / Mag(norm_d))
                        if angle > crit_angle:
                            YLEP_D = True
                            YLEP = [0.5 * (d_c[0] + c_c[0]), 0.5 * (d_c[1] + c_c[1]), 0.5 * (d_c[2] + c_c[2])]
                            yLEPs.append(YLEP)
                            print "Found YLEP at: ", YLEP
                            T1 = T_nodes[matching_i[0]]
                            T2 = T_nodes[matching_i[1]]
                            T_LEP  = 0.5 * (T1+T2)
                            rho1 = rho_nodes[matching_i[0]]
                            rho2 = rho_nodes[matching_i[1]]
                            rho_LEP  = 0.5 * (rho1+rho2)
                            mu1 = mu_nodes[matching_i[0]]
                            mu2 = mu_nodes[matching_i[1]]
                            mu_LEP  = 0.5 * (mu1+mu2)
                            p1 = p_nodes[matching_i[0]]
                            p2 = p_nodes[matching_i[1]]
                            p_LEP  = 0.5 * (p1+p2)
                            u1 = u_nodes[matching_i[0]]
                            u2 = u_nodes[matching_i[1]]
                            u_LEP  = 0.5 * (u1+u2)
                            v1 = v_nodes[matching_i[0]]
                            v2 = v_nodes[matching_i[1]]
                            v_LEP  = 0.5 * (v1+v2)
                            w1 = w_nodes[matching_i[0]]
                            w2 = w_nodes[matching_i[1]]
                            w_LEP  = 0.5 * (w1+w2)
                            M1 = Mach_nodes[matching_i[0]]
                            M2 = Mach_nodes[matching_i[1]]
                            M_LEP  = 0.5 * (M1+M2)
                            yLEPthermo.append([u_LEP, v_LEP, w_LEP, T_LEP, mu_LEP, rho_LEP, p_LEP, M_LEP])



                    if w_d * w_c < 0.:
                        angle = math.acos((norm_c[0] * norm_d[0] + norm_c[1] * norm_d[1]  + norm_c[2] * norm_d[2]) / Mag(norm_c) / Mag(norm_d))
                        if angle > crit_angle:
                            ZLEP_D = True
                            ZLEP = [0.5 * (d_c[0] + c_c[0]), 0.5 * (d_c[1] + c_c[1]), 0.5 * (d_c[2] + c_c[2])]
                            zLEPs.append(ZLEP)
                            print "Found ZLEP at: ", ZLEP
                            T1 = T_nodes[matching_i[0]]
                            T2 = T_nodes[matching_i[1]]
                            T_LEP  = 0.5 * (T1+T2)
                            rho1 = rho_nodes[matching_i[0]]
                            rho2 = rho_nodes[matching_i[1]]
                            rho_LEP  = 0.5 * (rho1+rho2)
                            mu1 = mu_nodes[matching_i[0]]
                            mu2 = mu_nodes[matching_i[1]]
                            mu_LEP  = 0.5 * (mu1+mu2)
                            p1 = p_nodes[matching_i[0]]
                            p2 = p_nodes[matching_i[1]]
                            p_LEP  = 0.5 * (p1+p2)
                            u1 = u_nodes[matching_i[0]]
                            u2 = u_nodes[matching_i[1]]
                            u_LEP  = 0.5 * (u1+u2)
                            v1 = v_nodes[matching_i[0]]
                            v2 = v_nodes[matching_i[1]]
                            v_LEP  = 0.5 * (v1+v2)
                            w1 = w_nodes[matching_i[0]]
                            w2 = w_nodes[matching_i[1]]
                            w_LEP  = 0.5 * (w1+w2)
                            M1 = Mach_nodes[matching_i[0]]
                            M2 = Mach_nodes[matching_i[1]]
                            M_LEP  = 0.5 * (M1+M2)
                            zLEPthermo.append([u_LEP, v_LEP, w_LEP, T_LEP, mu_LEP, rho_LEP, p_LEP, M_LEP])

                    if YLEP_D and ZLEP_D:
                        SP_D = True
                        SP = [0.5 * (d_c[0] + c_c[0]), 0.5 * (d_c[1] + c_c[1]), 0.5 * (d_c[2] + c_c[2])]
                        SPs.append(SP)
                        SPthermo.append([[u_LEP, v_LEP, w_LEP, T_LEP, mu_LEP, rho_LEP, p_LEP, M_LEP]])

                    if SP_D:
                        print "Another stagnation point found at: ", SP, ". Check if it is indeed true."

    return SPs, xLEPs, yLEPs, zLEPs, SPthermo, xLEPthermo, yLEPthermo, zLEPthermo





def FilterOutSPs(SPs, stag_point, epsilon, thermo):
    thermo_new = []
    filtered_SPs = []

    for s in range(0, len(SPs)):
        dist_to_stag = ((SPs[s][0] - stag_point[0])**2 + (SPs[s][1] - stag_point[1])**2 + (SPs[s][2] - stag_point[2])**2 )**0.5
        if dist_to_stag > epsilon:
            filtered_SPs.append(SPs[s])
            thermo_new.append(thermo[s])

    return filtered_SPs, thermo_new




def GetShockZ(Minf, R, beta, y, nosetype):
    z_crit = 100.0
    #if nosetype == 'Sphere':
    #    Rc = R * 1.143 * math.exp(0.54/(Minf - 1.)**1.2)
    #    delta = Rc * 0.143 * math.exp(3.24/Minf**2)
    #    z_crit = R + delta - Rc * (1./math.tan(beta))**2 * ((1. + (y**2 * math.tan(beta)**2)/(Rc**2))**0.5 - 1.)
    #elif nosetype == 'Cylinder':   
    #    Rc = R * 1.386 * math.exp(1.8/(Minf - 1.)**0.75)
    #    delta = Rc * 0.386 * math.exp(4.67/Minf**2)
    #    z_crit = R + delta - Rc * (1./math.tan(beta))**2 * ((1. + (y**2 * math.tan(beta)**2)/(Rc**2))**0.5 - 1.)        
    #else:
    #    print "Nose geometry type not recognized." 

    return z_crit

