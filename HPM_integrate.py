from HPM_import import *





def ReturnTrapezoid(fold, fnew, ds):
    return ds * 0.5 * (fnew + fold)





def RK4_3D(a, b, c, fa, fb, fc, hs):                                            #Fourth order RK in 3D from https://www.codeproject.com/Tips/792927/Fourth-Order-Runge-Kutta-Method-in-Python
    a1 = fa(a, b, c)*hs
    b1 = fb(a, b, c)*hs
    c1 = fc(a, b, c)*hs
    ak = a + a1*0.5
    bk = b + b1*0.5
    ck = c + c1*0.5
    a2 = fa(ak, bk, ck)*hs
    b2 = fb(ak, bk, ck)*hs
    c2 = fc(ak, bk, ck)*hs
    ak = a + a2*0.5
    bk = b + b2*0.5
    ck = c + c2*0.5
    a3 = fa(ak, bk, ck)*hs
    b3 = fb(ak, bk, ck)*hs
    c3 = fc(ak, bk, ck)*hs
    ak = a + a3
    bk = b + b3
    ck = c + c3
    a4 = fa(ak, bk, ck)*hs
    b4 = fb(ak, bk, ck)*hs
    c4 = fc(ak, bk, ck)*hs
    a = a + (a1 + 2.*(a2 + a3) + a4)/6.
    b = b + (b1 + 2.*(b2 + b3) + b4)/6.
    c = c + (c1 + 2.*(c2 + c3) + c4)/6.
    return a, b, c

