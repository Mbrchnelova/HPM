# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_vector_mathematics.py: this module contains functions for basic operations with vectors.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be




def Cross(a, b):
    x1 = a[2] * b[3] - a[3] * b[2]
    x2 = a[3] * b[1] - a[1] * b[3]
    x3 = a[1] * b[2] - a[2] * b[1]
    return [x1, x2, x3]



############################################################################################



def Length(a):
    return (a[0]**2 + a[1]**2 + a[2]**2)**2



############################################################################################



def Dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]



############################################################################################



def Mag(v):
    return (v[0]**2 + v[1]**2 + v[2]**2)**0.5



