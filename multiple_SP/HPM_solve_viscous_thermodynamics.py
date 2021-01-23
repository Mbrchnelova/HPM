# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_solve_viscous_thermodynamics.py: here, all functions required during marching down to compute the local thermodynamics are defined.
# These include the boundary layer momentum thickness and skin friction coefficient, depending on the flow state and Eckert's method.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *
from HPM_numerics import *
from HPM_mathematics import *




def SolveThermal(transitioned, prev_Re_theta, cp, k, T_e, T_w, Pr, mu_e, mu_star, rho_star, V_e, u_e, h, ds, rho_e, mu_star_old, rho_star_old, V_e_old, u_e_old, h_old, rho_e_old, integ_tracing, FIRST, INTERPOLATE_N, x):
    u_e = V_e
    u_e_old = V_e_old
    if V_e_old != 0.:
        theta_old = prev_Re_theta * mu_e / V_e_old / rho_e_old
    else:
        theta_old = 0.
    if not transitioned:
        mu_ref = mu_e
        T_ref = T_e
        S = 110.
        mu_w = mu_ref * (T_w/T_ref)**(3./2.) * (T_ref + S)/(T_w + S)
        Pr_w = mu_w * cp / k


        #This is only for fixing, remove this later
        #if rho_star > rho_star_old:
        #    rho_star = rho_star_old
        #if mu_star > mu_star_old:
        #    mu_star = mu_star_old

        if FIRST:
            integ = mu_star_old * rho_star_old * u_e * h**2 *ds/4.
            integ_tracing = integ
            tot_integ = integ
        else:
            f_old = mu_star_old * rho_star_old * u_e_old * h_old**2
            f_new = mu_star * rho_star * u_e * h**2

            integ = ReturnTrapezoid(f_old, f_new, ds)
            tot_integ = integ + integ_tracing
            integ_tracing += integ

        theta_L = 0.664 * tot_integ**0.5 / (rho_e * u_e * h)
        #theta_L = abs(0.002763/0.903293*x)
        H_e = cp * T_e
        H_w = cp * T_w
        R = Pr**0.5
        H_aw = H_e + R * V_e**2/2.
        H_w = H_aw/7.5
        Re_theta_e = theta_L * V_e * rho_e / mu_e

        dotq_w = 0.22 * 1./(Re_theta_e) * rho_star / rho_e * mu_star / mu_e * rho_e * u_e * (H_aw - H_w) * Pr_w**-0.6
        delta = 7.5 * theta_L
        theta = theta_L

        skin_cf = 2. * dotq_w * Pr**(2./3.) / rho_e / u_e / (H_aw - H_w)

    else:
        mu_ref = mu_e
        T_ref = T_e
        S = 110.
        mu_w = mu_ref * (T_w/T_ref)**(3./2.) * (T_ref + S)/(T_w + S)
        Pr_w = mu_w * cp / k

        if not INTERPOLATE_N:
            N = 12.67 - 6.5 * math.log(prev_Re_theta, 10) + 1.21 * math.log(prev_Re_theta**2, 10)
            if N < 0.:
                N = ReturnTurbN(Re_theta_e)
        else:
            N = ReturnTurbN(prev_Re_theta)

        m = 2./(N+1.)
        c5 = 2.2433 + 0.93 * N
        #print c5, N
        c1 = (1./c5)**((2. * N)/(N+1.))*(N/(N+1.)/(N+2.))**m
        c2 = (1. + m)*c1
        c3 = 1. + m
        c4 = 1./c3

        f_old = mu_star_old**m * rho_star_old * u_e_old * h_old**c3
        f_new = mu_star**m * rho_star * u_e * h**c3
        integ = ReturnTrapezoid(f_old, f_new, ds)
        integ_tot = integ + integ_tracing
        integ_tracing += integ

        theta_T = (c2* integ_tot)**c4 / (rho_e * u_e * h)
        #theta_T = abs(0.002763/0.903293*x)
        Re_theta_e = theta_T * V_e * rho_e / mu_e

        H_e = cp * T_e
        R = Pr**(1./3.)
        H_aw = H_e + R * V_e**2/2.
        H_w = cp * T_w
        H_w = H_aw/7.5

        dotq_w = c1 * (Re_theta_e)**(-m) * rho_star / rho_e * (mu_star / mu_e)**m * rho_e * u_e * (H_aw - H_w) * Pr_w**-0.4
        delta = theta_T * (N + 1. + (( (N+2.) /N ) * H_w/H_aw  + 1.) * (1. + 1.29 * Pr_w**0.333 * u_e**2/(2.*H_w)  )  )
        theta = theta_T

        if Re_theta_e == 0.:
            Re_theta_e += 0.00001

        skin_cf = 2. * c1 / Re_theta_e**m

    #print "Theta: ", theta, "Dotq_w: ", dotq_w, "delta: ", delta, "Pr_w: ", Pr_w
    return dotq_w, delta, theta, Re_theta_e, integ_tracing, skin_cf




###############################################################################



def ReturnTurbN(Re_theta_e):
    Ns = [4.00, 4.00, 4.01, 4.02, 4.05, 4.09, 4.11, 4.15, 4.19, 4.21, 4.26, 4.35,
          4.43, 4.53, 4.65, 4.75, 4.85, 4.93, 5.03, 5.11, 5.22, 5.33, 5.44, 5.58,
          5.68, 5.81, 5.90, 6.01, 6.08, 6.12, 6.18, 6.20, 6.21, 6.22, 6.25, 6.29,
          6.30, 6.31, 6.32, 6.34, 6.36, 6.39, 6.41, 6.42, 6.45, 6.47, 6.50, 6.52,
          6.55, 6.60, 6.64, 6.68, 6.73, 6.80, 6.83, 6.90, 6.98, 7.03, 7.12, 7.22,
          7.31, 7.40, 9.20, 11.10]
    Res = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000,
           3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000,
           9500, 10000, 11000, 18000, 19000, 20000, 21000, 22000, 23000, 24000,
           25000, 26000, 27000, 28000, 29000, 30000, 32000, 34000, 36000, 38000,
           40000, 44000, 47000, 50000, 55000, 60000, 65000, 70000, 80000, 90000,
           100000, 1000000, 10000000]

    CROSSED = False

    if Re_theta_e <= Res[0]:
        N = Ns[0]

    elif Re_theta_e >= Res[-1]:
        N = Ns[-1]

    else:
        for i in range(0, len(Res)):
            if Re_theta_e < Res[i]:
                if not CROSSED:
                    lower_i = i-1
                    upper_i = i
                    CROSSED = True
        N = Ns[lower_i] + (Re_theta_e - Res[lower_i]) * (Ns[upper_i] - Ns[lower_i])/(Res[upper_i] - Res[lower_i])

    return N


