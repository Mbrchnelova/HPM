# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_matrix_mathematics.py: in this module, all functions which deal with basic matrix operations are defined.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be


from numpy.linalg import inv
import numpy as np


def idMatx(size):
    id = []
    for i in range(size):
        id.append([0]*size)
    for i in range(size):
        id[i][i] = 1
    return(id)



############################################################################################



def tranMtx(inMtx):
    tMtx = []
    for row in range(0, len(inMtx[0])):
        tRow = []
        for col in range(0, len(inMtx)):
            ele = inMtx[col][row]
            tRow.append(ele)
        tMtx.append(tRow)
    return(tMtx)



############################################################################################



def matxRound(matx, decPts=4):
    for col in range(len(matx)):
        for row in range(len(matx[0])):
            matx[col][row] = round(matx[col][row], decPts)



############################################################################################




def GaussJordanB(A, b=False, decPts=4):
    if not b == False:
        if (len(A) != len(b)):
            print 'A and b are not conformable'
            return
        Ab = A[:]
        Ab.append(b)
        m = tranMtx(Ab)
    else:
        ii = idMatx(len(A))
        Aa = A[:]
        for col in range(len(ii)):
            Aa.append(ii[col])
        tAa = tranMtx(Aa)
        m = tAa[:]

    (eqns, colrange, augCol) = (len(A), len(A), len(m[0]))

    for col in range(0, colrange):
        bigrow = col
        for row in range(col+1, colrange):
            if abs(m[row][col]) > abs(m[bigrow][col]):
                bigrow = row
                (m[col], m[bigrow]) = (m[bigrow], m[col])
    print 'm is ', m

    for rrcol in range(0, colrange):
        for rr in range(rrcol+1, eqns):
            cc = -(float(m[rr][rrcol])/float(m[rrcol][rrcol]))
            for j in range(augCol):
                m[rr][j] = m[rr][j] + cc*m[rrcol][j]

    for rb in reversed(range(eqns)):
        if ( m[rb][rb] == 0):
            if m[rb][augCol-1] == 0:
                continue
            else:
                print 'Error: system is inconsistent!'
                return
        else:
            for backCol in reversed(range(rb, augCol)):
                m[rb][backCol] = float(m[rb][backCol]) / float(m[rb][rb])
            if not (rb == 0):
                for kup in reversed(range(rb)):
                    for kleft in reversed(range(rb, augCol)):
                        kk = -float(m[kup][rb]) / float(m[rb][rb])
                        m[kup][kleft] += kk*float(m[rb][kleft])
    matxRound(m, decPts)
    if not b == False:
        return m
    else:
        mOut = []
        for row in range(len(m)):
            rOut = []
            for col in range(augCol/2, augCol):
                rOut.append(m[row][col])
            mOut.append(rOut)
    return mOut


############################################################################################



def GaussJordan(m, eps = 1.0/(10**10)):                                         #Gaussian Jordan elimination from https://elonen.iki.fi/code/misc-notes/python-gaussj/
     (h, w) = (len(m), len(m[0]))
     for y in range(0,h):
         maxrow = y
         for y2 in range(y+1, h):    # Find max pivot
             if abs(m[y2][y]) > abs(m[maxrow][y]):
                 maxrow = y2
         (m[y], m[maxrow]) = (m[maxrow], m[y])
         if abs(m[y][y]) <= eps:     # Singular?
             return False
         for y2 in range(y+1, h):    # Eliminate column y
             c = m[y2][y] / m[y][y]
             for x in range(y, w):
                 m[y2][x] -= m[y][x] * c
     for y in range(h-1, 0-1, -1): # Backsubstitute
         c  = m[y][y]
         for y2 in range(0,y):
             for x in range(w-1, y-1, -1):
                 m[y2][x] -=  m[y][x] * m[y2][y] / c
         m[y][y] /= c
         for x in range(h, w):       # Normalize row y
             m[y][x] /= c
     return True



############################################################################################



def GetJacobian(point_1, point_2, point_3, vel_1, vel_2, vel_3):

    x_dot_vec = [vel_1[0], vel_2[0], vel_1[0]]
    y_dot_vec = [vel_1[1], vel_2[1], vel_1[1]]

    matrix =  np.array([1., point_1[0], point_1[1]], [1., point_2[0], point_2[1]], [1., point_3[0], point_3[1]])

    inverse = inv(matrix)

    a_1 = inverse[0, 0] * x_dot_vec[0] + inverse[0, 1] * x_dot_vec[1] + inverse[0, 2] * x_dot_vec[2]
    b_1 = inverse[1, 0] * x_dot_vec[0] + inverse[1, 1] * x_dot_vec[1] + inverse[1, 2] * x_dot_vec[2]
    c_1 = inverse[2, 0] * x_dot_vec[0] + inverse[2, 1] * x_dot_vec[1] + inverse[2, 2] * x_dot_vec[2]

    a_2 = inverse[0, 0] * y_dot_vec[0] + inverse[0, 1] * y_dot_vec[1] + inverse[0, 2] * y_dot_vec[2]
    b_2 = inverse[1, 0] * y_dot_vec[0] + inverse[1, 1] * y_dot_vec[1] + inverse[1, 2] * y_dot_vec[2]
    c_2 = inverse[2, 0] * y_dot_vec[0] + inverse[2, 1] * y_dot_vec[1] + inverse[2, 2] * y_dot_vec[2]
   
    jacobian = np.array([b_1, c_1],[b_2, c_2])
   
    tr = b_1 + c_1
    det = b_1 * c_2 - b_2 * c_1


