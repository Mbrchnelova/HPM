# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_calculate_stagnation_heating.py: this function computes heating in the stagnation point. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be


from HPM_import import *
from HPM_solve_eckert import *
from HPM_SETUP import *



def GetStagnationHeating(stg_cnds, T_w, Pr, cp, pinf, epsilon):
    M_e = 0.
    T_e = stg_cnds[0]
    rho_e = stg_cnds[1]
    p_e = stg_cnds[2]
    mu_e = stg_cnds[3]
    V_e = stg_cnds[4]

    Re_e = V_e * rho_e/ mu_e
    Re_wall, rho_wall, mu_wall, T_wall = ReturnEckert(False, T_w, T_e, M_e, gamma, Pr, mu_e, rho_e, Re_e)

    h_w = cp * T_w
    h_s = cp * T_e

    du_dS = 1. / epsilon * math.sqrt(2. * (p_e - pinf)/rho_wall)
    q_stag = 0.767* (rho_e * mu_e)**0.43 * (rho_wall * mu_wall)**0.07 * du_dS**0.5 * (h_s - h_w) * Pr**-0.6

    return q_stag



