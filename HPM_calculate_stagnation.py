from HPM_import import *
from HPM_SETUP import *

def GetStagPointConditions(Minf, pinf, rhoinf, Tinf):
    pstag = pinf * (1. + (gamma - 1.)/2. * Minf**2)**(gamma/(gamma - 1.))   
    pstag_behindNSW = pstag * (((gamma + 1.)*Minf**2)/ (2. + (gamma - 1.)*Minf**2))**(gamma/(gamma - 1.)) * ((gamma + 1.)/(2. * gamma * Minf**2 - (gamma - 1.)))**(1./(gamma - 1.))
    pinf_behindNSW = pinf * (2.0 * gamma * Minf**2 - gamma + 1.)/(gamma + 1.)
    Tstag = Tinf * (1. + (gamma - 1.)/2. * Minf**2)
    Tstag_behindNSW = Tstag
    Tinf_behindNSW = Tinf * ( 2. * gamma * Minf**2 - (gamma - 1.))*((gamma - 1.0)*Minf**2 + 2.0) / (gamma+1.0)**2 / Mach**2
    rhostag_behindNSW = pinf_behindNSW / (R * Tinf_behindNSW)
    
    if rhostag_behindNSW > rhoinf * 6.:
        rhostag_behindNSW = rhoinf * 6.
    
    rhostag_behindNSW = rhoinf * (gamma + 1.0)* Mach**2 / ((gamma - 1.0)* Mach**2 + 2.0)
    S = 110.
    mustag_behindNSW = muinf* (Tstag_behindNSW/Tinf)**(3./2.) * (Tinf + S)/(Tstag_behindNSW + S)
    Mach_behind_NSW = math.sqrt(((gamma - 1.)* Minf**2 +2.)/(2. * gamma * Minf**2 - gamma + 1.) )
    Vstag_behindNSW = Mach_behind_NSW * math.sqrt(gamma * R * Tinf_behindNSW)
    
    return Tstag_behindNSW, rhostag_behindNSW, pstag_behindNSW, mustag_behindNSW, Vstag_behindNSW

