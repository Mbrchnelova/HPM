# This module contains operations to handle vectors



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


