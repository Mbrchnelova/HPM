from HPM_import import *
from HPM_mathematics import *
from HPM_numerics import *
from HPM_viscous import *
from HPM_SETUP import *


def MarchDownstream(stg_coord, stg_cnds, wdydbetaz_ini, list_of_crossed_elements, list_of_coordinates, h, GradFs, Fxs, us_ele, vs_ele, ws_ele, T_ele, mu_ele, rho_ele, M_ele, w_fit, T_w, cp, k, wds, metric_coeffs, INTERPOLATE_N, OWN_TRANSITION, PERZHIKAR):

    d1s = wds[0]
    d2s = wds[1]
    d3s = wds[2]
    d4s = wds[3]
    d5s = wds[4]
    d6s = wds[5]
   
    wdydbetaz = wdydbetaz_ini
    coord_ini = list_of_coordinates[0]
   
    prev_Re_theta = 0.
    Pr = 1.0
    #stg_idx = T_ele.index(max(T_ele))
    #V_e_ini = (us_ele[stg_idx]**2 + vs_ele[stg_idx]**2  + ws_ele[stg_idx]**2)**0.5
    #T_e_ini = T_ele[stg_idx]
    #M_e_ini = M_ele[stg_idx]
    #rho_e_ini = rho_ele[stg_idx]
    #mu_e_ini = mu_ele[stg_idx]
   
    #stg_cnds = [Tstag_behindNSW, rhostag_behindNSW, pstag_behindNSW, mustag_behindNSW, Vstag_behindNSW]


    Tstag_behindNSW = stg_cnds[0]
    rhostag_behindNSW = stg_cnds[1]
    pstag_behindNSW = stg_cnds[2]
    mustag_behindNSW = stg_cnds[3]
    Vstag_behindNSW = stg_cnds[4]

    T_e_ini = Tstag_behindNSW
    rho_e_ini = rhostag_behindNSW
    p_e_ini = pstag_behindNSW
    V_e_ini = Vstag_behindNSW
    M_e_ini = Vstag_behindNSW / math.sqrt(1.4 * 287.05 * T_e_ini)
    T_e_ini = Tstag_behindNSW / (1. + (gamma - 1.0)/2. * M_e_ini**2 )
    S = 110.
    mu_e_ini = mustag_behindNSW / (1. + (gamma - 1.0)/2. * M_e_ini**2 )
    mu_e_ini = mustag_behindNSW * (Tstag_behindNSW/Tstag_behindNSW)**(3./2.) * (Tstag_behindNSW + S)/(Tstag_behindNSW + S)

    Re_e_ini = rho_e_ini * V_e_ini / mu_e_ini

    Re_star, rho_star, mu_star, T_star = ReturnEckert(False, T_w, T_e_ini, M_e_ini, gamma, Pr, mu_e_ini, rho_e_ini, Re_e_ini)
    rho_star_old = rho_star
    mu_star_old = mu_star
    rho_e_old = rho_e_ini
    V_e_old = V_e_ini
    u_e_old = 0.
    h_old = h
    coord_old = coord_ini


   
    all_points = []
    hs = []
    dotq_ws = []
    skin_cfs = []
    deltas = []
    thetas = []
    V_es = []
    rho_es = []
    mu_stars = []
    dss = []
    rho_stars = []
    integs = []
    transition_points = []
    transitioned = False
    NOTIFY_TRANSITION = True
    FIRST = True
    integ_tracing = 0.
    FIRST_turb = False
    #dwdy = 2d1y + d2 + d5z
    #dwdz = 2d3z + d4 + d5y

    for i in range(1, len(list_of_crossed_elements)):
        h_perzhikar = metric_coeffs[i]
        el = list_of_crossed_elements[i]
        el_old = list_of_crossed_elements[i-1]
        V = (us_ele[el]**2 + vs_ele[el]**2 + ws_ele[el]**2)**0.5
        V_old = (us_ele[el_old]**2 + vs_ele[el_old]**2 + ws_ele[el_old]**2)**0.5
        if V_old == 0.:
            V_old = 10.

        u_e = abs(us_ele[el])

        metric_coeff = metric_coeffs[i]

        coord = list_of_coordinates[i]
        x_coord = coord[0]
        ds = ((coord[0] - coord_old[0])**2 + (coord[1] - coord_old[1])**2 + (coord[2] - coord_old[2])**2)**0.5
        if V == 0.:
            V = 0.1




        if not PERZHIKAR:
            wdydbetaz_prev = wdydbetaz
            term_around_prev = 0.
            if (coord[1]**2 + coord[2]**2)**0.5 > 10.*epsilon:

                y_track += ds * vs_ele[el] / V
                z_track += ds * ws_ele[el] / V

                coord_old = list_of_coordinates[i-1]
                el_evenolder = list_of_crossed_elements[i-2]
                #dw_old =  ws_ele[el_old] - ws_ele[el_evenolder]
                #dy_old = list_of_coordinates[i-1][1] - list_of_coordinates[i-2][1]
                #dz_old = list_of_coordinates[i-1][2] - list_of_coordinates[i-2][2]

                #print d1s[el_old], coord[1]
                #dwdy = 2.*d1s[el_old]*coord[1] - d2s[el_old] - d5s[el_old]*coord[2]
                #dwdz = 2.*d3s[el_old]*coord[2] + d4s[el_old] - d5s[el_old]*coord[1]

                dwdy = 2.*d1s[el_old]*coord_old[1] + d2s[el_old] + d5s[el_old]*coord_old[2]
                dwdz = 2.*d3s[el_old]*coord_old[2] + d4s[el_old] + d5s[el_old]*coord_old[1]

                dwdy_new = 2.*d1s[el]*coord_old[1] + d2s[el] + d5s[el]*coord_old[2]
                dwdz_new = 2.*d3s[el]*coord_old[2] + d4s[el] + d5s[el]*coord_old[1]

                dw_old =  ws_ele[el_old] - ws_ele[el_evenolder]
                dy_old = list_of_coordinates[i-1][1] - list_of_coordinates[i-2][1]
                dz_old = list_of_coordinates[i-1][2] - list_of_coordinates[i-2][2]

                dw_new =  ws_ele[el] - ws_ele[el_old]
                dy_new = list_of_coordinates[i][1] - list_of_coordinates[i-1][1]
                dz_new = list_of_coordinates[i][2] - list_of_coordinates[i-1][2]

                dwdy = dw_old/dy_old
                dwdz = dw_old/dz_old

                dwdy = dw_new/dy_new
                dwdz = dw_new/dz_new

                dwdylh = dw_new / y_track
                dwdzlh = dw_new / z_track


                #term_around_prev =  1./V_old * (dw_old/dy_old + dw_old/dz_old) * ds 
                #term_around_prev =  1./V_old * (-dwdy + dwdz) * ds 


                #wdydbetaz = wdydbetaz_prev * (1. + (term_around_prev))
                '''
                
                c2 = 1./V * (-dwdy_new + dwdz_new)
                c1 = 1./V_old * (-dwdy + dwdz)
                
                #print (math.log(wdydbetaz_prev) + ((c2 - c1)*ds))
                avec = 0.5 * (c2+c1)
                wdydbetaz = math.exp(math.log(wdydbetaz_prev) + (avec*ds))
                '''


                wdydbetaz = wdydbetaz_prev + wdydbetaz_prev * 1./V_old * (dwdy + dwdz) * ds
                shittyterm1 = dwdy
                shittyterm2 = dwdz

            else:
                y_track += ds * vs_ele[el] / V
                z_track += ds * ws_ele[el] / V
                #if round(wdydbetaz_ini, 8) == 0.:
                #    wdydbetaz_ini =  (abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)) * V *  Fxs[el] / abs(GradFs[el])
                #wdydbetaz = wdydbetaz_ini
                wdydbetaz = (abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)) * V *  Fxs[el] / abs(GradFs[el])
                if i == 1:
                    wdydbetaz_ini = wdydbetaz
                    print "\n"
                    print "wdydbetaz_ini: ", wdydbetaz_ini
                term_around_prev = 0
                shittyterm1 = 0
                shittyterm2 = 0

            #wdydbetaz = max([wdydbetaz, wdydbetaz_ini])


            h1 = abs(GradFs[el]) / Fxs[el]/ V * (wdydbetaz)

            if abs(h1) > 3.:
                h1 = 3.
            #h = 0.1
            #h = abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)
            wdydbetaz_correct = h * V * Fxs[el] / abs(GradFs[el])
            print round(h, 4), round(h1, 4), '\t\t', round(wdydbetaz, 4), round(wdydbetaz_correct,4), shittyterm1, shittyterm2, ds, ws_ele[el]

        else:
            h = metric_coeff
            #h = abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)
        #h = abs(h1)
        #h = abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)
        if h <= 0.:
            break
        V_e = V
        T_e = T_ele[el]
        M_e = M_ele[el]
        rho_e = rho_ele[el]
        mu_e = mu_ele[el]
        Re_e = rho_e * V_e / mu_e
        Re_star, rho_star, mu_star, T_star = ReturnEckert(transitioned, T_w, T_e, M_e, gamma, Pr, mu_e, rho_e, Re_e)
        #if h < 0.:
        #    h = abs(h)

        dotq_w, delta, theta, Re_theta_e, integ_tracing, skin_cf = SolveThermal(transitioned,
                                            prev_Re_theta, cp, k, T_e, T_w, Pr, mu_e, mu_star,
                                            rho_star, V_e, u_e, h, ds, rho_e, mu_star_old,
                                            rho_star_old, V_e_old, u_e_old, h_old, rho_e_old,
                                            integ_tracing, FIRST, FIRST_turb, INTERPOLATE_N, x_coord)


        transitioned_prev = transitioned
        transitioned = Transition(M_e, Re_theta_e, theta, V_e, rho_e, mu_e, abs(x_coord), OWN_TRANSITION)
        #if h < 0.1:
        #transitioned = False
        if transitioned != transitioned_prev:
            #integ_tracing = 0.0
            FIRST_turb = True
        else:
            FIRST_turb = False
        if transitioned and NOTIFY_TRANSITION:
            print "Transition occurred at eq. body radius of:", h, " and xcoord of ", x_coord,  " with Re_T of: ", Re_theta_e
            transition_points.append(h)
            NOTIFY_TRANSITION = False

        #if not transitioned:
        prev_Re_theta = Re_theta_e
        mu_star_old = mu_star
        V_e_old = V_e
        u_e_old = u_e
        rho_e_old = rho_e
        rho_star_old = rho_star

        coord_old = coord

        h_old = h
        hs.append(h)
        V_es.append(V_e)
        rho_es.append(rho_e)
        mu_stars.append(mu_star)
        dss.append(ds)
        rho_stars.append(rho_star)
        integs.append(integ_tracing)
        dotq_ws.append(dotq_w)
        skin_cfs.append(skin_cf)
        deltas.append(delta)
        thetas.append(theta)
        all_points.append(list_of_coordinates[i])
        FIRST = False
    return hs, all_points, thetas, dotq_ws, deltas, transition_points, V_es, rho_es, mu_stars, rho_stars, integs, dss, skin_cfs


