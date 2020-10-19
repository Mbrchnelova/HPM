from HPM_import import *


def ReturnEckert(transitioned, T_w, T_e, M_e, gamma, Pr, mu_e, rho_e, Re_e):
    if not transitioned:
        r = math.sqrt(Pr)                                                       #Laminar recovery
    else:
        r = Pr**(1./3.)                                                         #Turbulent recovery
    S = 110.                                                                    #Sutherland's constant for air
    T_star = T_e * ( 0.5 + 0.5 * T_w/T_e + 0.22*r*(gamma-1.)/2. * M_e**2)     #Eckert
    T_ref = T_e
    mu_ref = mu_e
    mu_star = mu_ref * (T_star/T_ref)**(3./2.) * (T_ref + S)/(T_star + S)       #Sutherland
    rho_star = rho_e * (T_e / T_star )**(gamma-1.0)                                        #Perfect gas
    Re_star = Re_e * (rho_star/rho_e) * (mu_e/mu_star)                          #Re number definition

    return Re_star, rho_star, mu_star, T_star



