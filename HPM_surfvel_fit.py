from HPM_import import *
from HPM_mathematics import *



#Format:    u[u1, u2, u3], v[v1, v2, v3], w[w1, w2, w3], of the three vertices
#           coord_1, coord_2, coord_3 = [x, y, z]
def CreateSurfaceFit(u, v, w, coord_1, coord_2, coord_3):             #CHECKED AND WORKS
    u_vec = [v, w, u]   
    coords = [coord_2, coord_3, coord_1]
    print u_vec, coords,
    matrix = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
    for i in range(0, 6):
        if i%2 == 0:
            for j in range(0, 6):
                if j%2 == 0:
                    matrix[i][j] = coords[i/2][j/2]**2
                else:
                    matrix[i][j] = coords[i/2][j/2]
        else:
            for j in range(0, 6):
                if j%2 == 0:
                    matrix[i][j] = 2. * u_vec[i/2][j/2] * coords[i/2][j/2]
                else:
                    matrix[i][j] = u_vec[i/2][j/2]
    b = [0.,0.,0.,0.,0.,0.]
    a_vector = list(solve(matrix, b))                                           #REDO THIS INTO A PURE PYTHON STYLE
    print matrix
    print '\n'
    #a_vector = GaussJordanB(matrix, [0.,0.,0.,0.,0.,0.])
    return a_vector



###############################################################################



#Format:    V_1 - V_6 are velocity components to be fitted (so u, v or w) or also pressure and enthalpy
#           coord_1, coord_2, coord_3 .., = [x, y, z]
def CreateVelFit(V_1, V_2_, V_3, V_4, V_5, V_6, coord_1, coord_2, coord_3, coord_4, coord_5, coord_6):
    y_coords = [coord_1[1], coord_2[1], coord_3[1], coord_4[1], coord_5[1], coord_6[1]]
    z_coords = [coord_1[2], coord_2[2], coord_3[2], coord_4[2], coord_5[2], coord_6[2]]
    v_vels = [V_1, V_2_, V_3, V_4, V_5, V_6]
    matrix = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],]
    for i in range(0, 6):
        for j in range(0, 6):
            if j == 0:
                matrix[i][j] = y_coords[i]**2
            if j == 1:
                matrix[i][j] = y_coords[i]
            if j == 2:
                matrix[i][j] = z_coords[i]**2
            if j == 3:
                matrix[i][j] = z_coords[i]
            if j == 4:
                matrix[i][j] = y_coords[i] * z_coords[i]
            if j == 5:
                matrix[i][j] = 1.
    b_vector = GaussJordanB(matrix, v_vels)
    return b_vector
