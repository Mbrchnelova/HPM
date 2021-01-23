# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_get_transition.py: this function determines whether the flow has transitioned or not.
# Three conditions are currently implemented - Re(x) based (user), x-based (user) and ReT based.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *


def Transition(M_e, Re, theta_e, u_e, rho_e, mu_e, x_coord, OWN):
    ReT = 10.**(6.421 * math.exp(1.209e-4 * M_e**2.641))
    Rex = rho_e * abs(u_e) * x_coord / mu_e
    if OWN[0] != 'nan':
        if Rex >= OWN[0]:
            return True
    elif OWN[1] != 'nan':
        if x_coord >= OWN[1]:
            return True
    elif Rex >= ReT:                                                             #First transition condition
        return True
    else:
        Re_theta = rho_e * abs(u_e) * theta_e / mu_e
        Re_thetaT = 100. * M_e
        if Re_theta >= Re_thetaT:                                               #Second transition condition
            return False                                                         #SET TO FALSE JUST FOR TESTING, PUT IT BACK ONCE YOU ARE DONE
        else:
            return False


