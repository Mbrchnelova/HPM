#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 12:52:00 2019

@author: brch
"""

#import HPM_3_ply
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
from scipy.linalg import solve
import copy
import time
from numpy.linalg import inv

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

############################################################################### MATHEMATICS

def idMatx(size):
    id = []
    for i in range(size):
        id.append([0]*size)
    for i in range(size):
        id[i][i] = 1
    return(id)
 
 
#def MatrixMultiplyVector(matrix, vector):
#    for i in range 
    
   
def tranMtx(inMtx):
    tMtx = []
    for row in range(0, len(inMtx[0])):
        tRow = []
        for col in range(0, len(inMtx)):
            ele = inMtx[col][row]
            tRow.append(ele)
        tMtx.append(tRow)
    return(tMtx)
 
    
    
    
def matxRound(matx, decPts=4):
    for col in range(len(matx)):
        for row in range(len(matx[0])):
            matx[col][row] = round(matx[col][row], decPts)
 
 
    
    
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

def Cross(a, b):
    x1 = a[2] * b[3] - a[3] * b[2]
    x2 = a[3] * b[1] - a[1] * b[3]
    x3 = a[1] * b[2] - a[2] * b[1]
    return [x1, x2, x3]

def Length(a):
    return (a[0]**2 + a[1]**2 + a[2]**2)**2

def Dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def ProjectTriangleFro3Dto2D(coord_1, coord_2, coord_3):
    x_1 = coord_1[0]   
    x_2 = coord_2[0]   
    x_3 = coord_3[0]   
    y_1 = coord_1[1]   
    y_2 = coord_2[1]   
    y_3 = coord_3[1]  
    z_1 = coord_1[2]   
    z_2 = coord_2[2]   
    z_3 = coord_3[2]  
    
    vector_12 = [x_2 - x_1, y_2 - y_1, z_2 - z_1 ]
    mag_vector_12 = Mag(vector_12)
    if mag_vector_12 == 0.:
        print "Error: vector of zero magnitude! Adjusted to very small numbers."
        mag_vector_12 += 0.000001
        
    vector_12 = [vector_12[0]/ mag_vector_12, vector_12[1]/ mag_vector_12, vector_12[2]/ mag_vector_12]
    
    vector_13 = [x_3 - x_1, y_3 - y_1, z_3 - z_1 ]
    #mag_vector_13 = Mag(vector_13)
    #vector_13 = [vector_13[0]/ mag_vector_13, vector_13[1]/ mag_vector_13, vector_13[2]/ mag_vector_13]
    
    cross_A = Cross(vector_12, vector_13)
    cross_A_mag = Mag(cross_A)
    if cross_A_mag == 0.:
        print "Error: vector of zero magnitude! Adjusted to very small numbers."
        cross_A_mag += 0.000001
    cross_A = [cross_A[0]/cross_A_mag, cross_A[1]/cross_A_mag, cross_A[2]/cross_A_mag]
    cross_B = Cross(cross_A, vector_12)
    
    x1 = 0.0
    y1 = 0.0
    x2 = Length(vector_12)
    y2 = 0.0
    x3 = Dot(vector_13, vector_12)
    y3 = Dot(vector_13, cross_B)
    
    coord_1_2D = [x1, y1]
    coord_2_2D = [x2, y2]
    coord_3_2D = [x3, y3]
    
    return [coord_1_2D, coord_2_2D, coord_3_2D]
    
    
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
 
    
def Mag(v):
    return (v[0]**2 + v[1]**2 + v[2]**2)**0.5

def Dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] 

def process_mesh(filename, VERIFY):
    with open(filename, "r") as ins:
        array = []
        for line in ins:
            array.append(line)
    
    alllines = []
    COUNT = False
    for line in array:
        if line[-1] == '\n':
            line = line[0:-1]
            
        line = line.strip(' ')
        
        if COUNT:
            alllines.append(line.split(' '))
        
        if line == 'end_header':
            COUNT = True
        
    nodes = []
    normals_nodes = []
    connectivity = []
    centroids = []
    normals = []
    test_dot_products = []
    areas = []

    for line in alllines:
        if len(line) == 6:
            nodes.append([float(line[0]), float(line[1]), float(line[2])])
            normals_nodes.append([float(line[3]), float(line[4]), float(line[5])])
            
        if len(line) == 4:
            connectivity.append([int(line[1]), int(line[2]), int(line[3])])
                
    elements_num = len(connectivity)
    
    if len(connectivity[0]) > 3:
        print "Incorrect mesh formar! Triangular elements required, please re-format and try again."
        
    for i in range(0, elements_num):
        first_point = connectivity[i][0]
        second_point = connectivity[i][1]
        third_point = connectivity[i][2]
        
        ctr_x = (nodes[first_point][0] + nodes[second_point][0] + nodes[third_point][0])/3.
        ctr_y = (nodes[first_point][1] + nodes[second_point][1] + nodes[third_point][1])/3.
        ctr_z = (nodes[first_point][2] + nodes[second_point][2] + nodes[third_point][2])/3.
        
        centroids.append([ctr_x, ctr_y, ctr_z])
    
    for i in range(0, elements_num):
        first_point = nodes[connectivity[i][0]]
        second_point = nodes[connectivity[i][1]]
        third_point = nodes[connectivity[i][2]]
        
        vector1 = [second_point[0] - first_point[0], second_point[1] - first_point[1], second_point[2] - first_point[2]]
        vector2 = [third_point[0] - first_point[0], third_point[1] - first_point[1], third_point[2] - first_point[2]]
    
        n_x = vector1[1] * vector2[2] - vector1[2] * vector2[1]
        n_y = vector1[2] * vector2[0] - vector1[0] * vector2[2]
        n_z = vector1[0] * vector2[1] - vector1[1] * vector2[0]
        
        normal = [n_x, n_y, n_z]
        
        dot_product = n_x * centroids[i][0] +  n_y * centroids[i][1] +  n_z * centroids[i][2]
        
        if dot_product < 0.:
            normal = [-n_x, -n_y, -n_z]    
        
        test_normal_1 = normals_nodes[connectivity[i][0]]
        test_normal_2 = normals_nodes[connectivity[i][1]]
        test_normal_3 = normals_nodes[connectivity[i][2]]
        
        ave_x = (test_normal_1[0] + test_normal_2[0] + test_normal_3[0])/3.
        ave_y = (test_normal_1[1] + test_normal_2[1] + test_normal_3[1])/3.
        ave_z = (test_normal_1[2] + test_normal_2[2] + test_normal_3[2])/3.
        
        test_normal_ave = [ave_x, ave_y, ave_z]
        #print test_normal_1, test_normal_2, test_normal_3, test_normal_4
        test_dot_product = normal[0] * test_normal_ave[0] +  normal[1] * test_normal_ave[1] +  normal[2] * test_normal_ave[2]
        
        test_dot_products.append(test_dot_product)
        
        if test_dot_product > 0:
            normal = [-n_x, -n_y, -n_z]    
            #print 'problem', test_dot_product
            
        normals.append(normal)
    
    for i in range(0, elements_num):
        first_point = nodes[connectivity[i][0]]
        second_point = nodes[connectivity[i][1]]
        third_point = nodes[connectivity[i][2]]
        
        vector1 = [second_point[0] - first_point[0], second_point[1] - first_point[1], second_point[2] - first_point[2]]
        vector2 = [third_point[0] - first_point[0], third_point[1] - first_point[1], third_point[2] - first_point[2]]
         
        a_x = vector1[1] * vector2[2] - vector1[2] * vector2[1]
        a_y = vector1[2] * vector2[0] - vector1[0] * vector2[2]
        a_z = vector1[0] * vector2[1] - vector1[1] * vector2[0]
        
        area = 0.5 * (a_x **2 + a_y**2 + a_z**2)*0.5
        areas.append(area)
    
    #RETURN node coordinates, node connectivity, element normals, element centroids, element areas
    
    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        
        x_coords_ele = []
        y_coords_ele = []
        z_coords_ele = []
        #u_vel_mag = []
        for c in range(0, len(centroids)):
            x_coords_ele.append(centroids[c][0])
            y_coords_ele.append(centroids[c][1])
            z_coords_ele.append(centroids[c][2])
            #u_vel_mag.append(velocities[c][2])
        c = np.array(areas)
        ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Nearest points resolved for parameter fitting')
        
        t_check = [96, 45, 234, 97, 232, 800]
        x_coords_neigh_check = [nodes[t_check[0]][0], nodes[t_check[1]][0], nodes[t_check[2]][0], nodes[t_check[3]][0], nodes[t_check[4]][0], nodes[t_check[5]][0]]
        y_coords_neigh_check = [nodes[t_check[0]][1], nodes[t_check[1]][1], nodes[t_check[2]][1], nodes[t_check[3]][1], nodes[t_check[4]][1], nodes[t_check[5]][1]]
        z_coords_neigh_check = [nodes[t_check[0]][2], nodes[t_check[1]][2], nodes[t_check[2]][2], nodes[t_check[3]][2], nodes[t_check[4]][2], nodes[t_check[5]][2]]
        
        ax.scatter(x_coords_neigh_check, y_coords_neigh_check, z_coords_neigh_check, s=15, c = 'r')
        
    return[nodes, connectivity, normals, centroids, areas] 


def Newton(nodes, connectivity, normals, centroids, areas, gamma, V_vector, Mach, VERIFY):
    x_coords = []
    y_coords = []
    z_coords = []
    
    for node in nodes:
        x_coords.append(node[0])
        y_coords.append(node[1])
        z_coords.append(node[2])
          
    x_coords_ele = []
    y_coords_ele = []
    z_coords_ele = []
    
    x_coords_ele_cp = []
    y_coords_ele_cp = []
    z_coords_ele_cp = []
    
    nonzero_cp = []
    shadowed = []
    
    for centroid in centroids:
        x_coords_ele.append(centroid[0])
        y_coords_ele.append(centroid[1])
        z_coords_ele.append(centroid[2])
        
    '''
    fig = plt.figure()
    ax = Axes3D(fig)
    '''
    
    V_mag = (V_vector[0]**2 + V_vector[1]**2 + V_vector[2]**2)**0.5 
    Cps = []
    
    if gamma == 1.:
        print "Gamma cannot be 1! Adjusted to prevent singularities in calculations."
        gamma += 0.001
            
    for i in range(0, len(connectivity)):
        normal = normals[i]
        dot_product = -(normal[0] * V_vector[0] + normal[1] * V_vector[1] + normal[2] * V_vector[2])
        normal_mag = (normal[0]**2 + normal[1]**2 + normal[2]**2)**0.5
        
            
        #Cp_max = (1. - gamma + 2. * Mach**2 * gamma)/(gamma + 1.) * ((gamma + 1.)**2 * Mach**2/(4. * gamma * Mach**2 - 2.*( gamma - 1.)))**(-gamma / (gamma - 1.))
        Cp_max = ((gamma + 1.)**2 / (4. * gamma))**(gamma / (gamma - 1.))  * (4. / (gamma + 1.))
        if dot_product < 0.:
            theta = math.acos(dot_product/(V_mag * normal_mag))
            Cp = Cp_max * math.sin(math.pi/2. - theta)**2
            shadowed.append('False')
            #print theta
            #print dot_product, normal_mag
        else:
            Cp = 0.
            shadowed.append('True')
            
        if Cp == 0:
            centroid = centroids[i]
            x_coords_ele_cp.append(centroid[0])
            y_coords_ele_cp.append(centroid[1])
            z_coords_ele_cp.append(centroid[2])
            nonzero_cp.append(Cp)
            
        Cps.append(Cp)
    #print Cp_max
    
    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(Cps)
        
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Cp on a cone using modified Newton technique')
        ax.axes.set_aspect('equal')

    return Cps, shadowed



def ReturnTrapezoid(fold, fnew, ds):
    return ds * 0.5 * (fnew + fold)


def ReadNullSpace(filename):
    j = 0.
    null_spaces = []
    with open(filename, "r") as ins:
        nullspace_array = []
        for line in ins:
            line = line.strip()
            line = line.split()
            if len(line) == 1:
                nullspace_array.append(float(line[0]))
            if len(line) > 1:
                null_spaces.append(nullspace_array)
                nullspace_array = []
            
            j += 1
    null_spaces.append(nullspace_array)
    
    null_spaces_formatted = []
    for i in range(0, len(null_spaces)):
        null_spaces_formatted.append(null_spaces[i])
        
    return null_spaces_formatted[1:]
    
def ReadInviscidData(filename):
    data = []
    nodes = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            
            nodes.append([x, y, z])
            data.append(float(line[3]))

    return data, nodes


def ReadCentroids(filename):
    centroids = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            
            centroids.append([x, y, z])

    return centroids


def ReadNormals(filename):
    normals = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            
            normals.append([x, y, z])

    return normals

def ReadStagPointData(filename):
    stag_points = []
    stg_idxs = []
    epss = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            stag_point = [float(line[0]), float(line[1]), float(line[2])]
            stg_idx = int(line[3])
            eps = float(line[4])
            
            stag_points.append(stag_point)
            stg_idxs.append(stg_idx)
            epss.append(eps)
    
    return stag_points, stg_idxs, epss
            
    

def ReadConnectivity(filename):
    connectivity = []
    with open(filename, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = int(line[0])
            y = int(line[1])
            z = int(line[2])
            
            connectivity.append([x, y, z])

    return connectivity


def ReadSolvedFitting(filename):
    j = 0.
    vectors = []
    with open(filename, "r") as ins:
        vector_array = []
        for line in ins:
            line = line.strip()
            line = line.split()
            if len(line) == 1:
                vector_array.append(float(line[0]))
            if len(line) > 1:
                vectors.append(vector_array)
                vector_array = []
            
            j += 1
    vectors.append(vector_array)
    
    solutions_formatted = []
    for i in range(0, len(vectors)):
        solutions_formatted.append(vectors[i])

    return solutions_formatted[1:]


def ReadBacktracingData(filename_intersections, filename_nodepaths, filename_nodecoords, filename_nodesresolved, filename_epsilonnodes):
    intersections = []
    node_paths_elem = []
    node_paths_coord = []
    epsilon_nodes = []
    nodes_resolved_idxs = []

    with open(filename_epsilonnodes, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            
            epsilon_nodes.append([x, y, z])
    
    with open(filename_intersections, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            
            intersections.append([x, y, z])
            
    with open(filename_nodesresolved, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            node_idx = int(line[0])
            
            nodes_resolved_idxs.append(node_idx)
            
        
    with open(filename_nodepaths, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            node_path = []
            for i in range(0, len(line)):
                node_path.append(int(line[i]))
            
            node_paths_elem.append(node_path)      
            
            
    with open(filename_nodecoords, "r") as ins:
        for line in ins:
            line = line.strip()
            line = line.split()
            node_path_coord = []
            coord = []
            for i in range(0, len(line)):
                #print line
                if line[i][0] != ',':
                    if line[i] != '' and line[i] != ' ' and line[i] != '\n':
                        coord_cur = float(line[i])
                        coord.append(coord_cur)
                else:
                    if line[i] != '' and line[i] != ' '  and line[i] != '\n':
                        node_path_coord.append(coord)  
                        if len(line[i][1:]) > 1:
                            new_coord = line[i][1:]
                            coord = [float(new_coord)]
                            
            node_paths_coord.append(node_path_coord)    
            
            
    return intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs 


############################################################################### NUMERICS
    
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


def GetVertexValues(nodes, centroids, centroid_values, p, VERIFY):
    values_nodes = []
    for n in range(0, len(nodes)):
        value = 0.
        D = 0.
        for c in range(0, len(centroids)):
            distance = ((nodes[n][0] - centroids[c][0])**2 + (nodes[n][1] - centroids[c][1])**2 + (nodes[n][2] - centroids[c][2])**2 )**0.5
            if distance == 0.:
                print "Node overlaps with a mesh cell centroid. Try to improve the quality of the mesh."
                distance = 0.001
            value += centroid_values[c] * 1./(distance)**p
            D += 1./(distance)**p
        value /= D
        values_nodes.append(value)
        if n%100 == 0.:
            print 'Node number...:', n,
    
    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        
        x_nodes= []
        y_nodes= []
        z_nodes= []
        
        for n in range(0, len(nodes)):
            x_nodes.append(nodes[n][0])
            y_nodes.append(nodes[n][1])
            z_nodes.append(nodes[n][2])
            
        c = np.array(values_nodes)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_nodes, y_nodes, z_nodes, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Node resolved from the ISD centroid interpolation')
    return values_nodes


def GetVertexValuesFaster(nodes, centroids, connectivity, u_centroids, v_centroids, w_centroids, p_centroids):
    u_nodes = []
    v_nodes = []
    w_nodes = []
    p_nodes = []
    how_many = []
    
    for i in range(0, len(nodes)):
        u_nodes.append(0.)
        v_nodes.append(0.)
        w_nodes.append(0.)
        p_nodes.append(0.)
        how_many.append(0)
    
    for i in range(0, len(centroids)):
        u_ele = u_centroids[i]
        v_ele = v_centroids[i]
        w_ele = w_centroids[i]
        p_ele = p_centroids[i]
        
        vertices = connectivity[i]
        v1 = vertices[0]
        v2 = vertices[1]
        v3 = vertices[2]
        
        u_nodes[v1] += u_ele
        v_nodes[v1] += v_ele
        w_nodes[v1] += w_ele
        p_nodes[v1] += p_ele
        
        u_nodes[v2] += u_ele
        v_nodes[v2] += v_ele
        w_nodes[v2] += w_ele
        p_nodes[v2] += p_ele
        
        u_nodes[v3] += u_ele
        v_nodes[v3] += v_ele
        w_nodes[v3] += w_ele
        p_nodes[v3] += p_ele
        
        how_many[v1] += 1
        how_many[v2] += 1
        how_many[v3] += 1
    
    for i in range(0, len(nodes)):
        how_many_at_i = how_many[i]
        if how_many_at_i == 0:
            print "Not all nodes included in the connectivity matrix. Either fix your mesh or use IDS for value interpolations."
            
        u_nodes[i] = u_nodes[i]/how_many_at_i
        v_nodes[i] = v_nodes[i]/how_many_at_i
        w_nodes[i] = w_nodes[i]/how_many_at_i
        p_nodes[i] = p_nodes[i]/how_many_at_i
    
    return u_nodes, v_nodes, w_nodes, p_nodes


############################################################################### HYPERSONICS, INVISCID

def GetThermodynamics(ps, Minf, nodes, gamma, R, muinf, Tinf, rhoinf, VERIFY):
    if gamma == 1.:
        print "Gamma cannot be 1! Adjusted to prevent singularities."
        gamma += 0.001
        
    if gamma == 0.:
        print "Gamma cannot be 0! Adjusted to prevent singularities. Results will be unphysical until you get your shit together with the input data."
        gamma += 0.001
        
    pstag = pinf * (1. + (gamma - 1.)/2. * Minf**2)**(gamma/(gamma - 1.))   
    pstag_behindNSW = pstag * (((gamma + 1.)*Minf**2)/ (2. + (gamma - 1.)*Minf**2))**(gamma/(gamma - 1.)) * ((gamma + 1.)/(2. * gamma * Minf**2 - (gamma - 1.)))**(1./(gamma - 1.))
    Tstag = Tinf * (1. + (gamma - 1.)/2. * Minf**2)
    Tstag_behindNSW = Tstag
    rhostag_behindNSW = pstag_behindNSW / (R * Tstag_behindNSW)
    if rhostag_behindNSW > rhoinf * 6.:
        print "Accuracy compromised, rho behind NSW had to be adjusted to the maximum value (6rho_inf). The resulting data are unreliable."
        rhostag_behindNSW = rhoinf * 6.
    T_nodes = []
    mu_nodes = []
    rho_nodes = []
    M_nodes = []
    for n in range(0, len(nodes)):
        p = ps[n]
        if p == 0.:
            "Zero pressure detected on the windward side! Check the suitability of the mesh and input data. Adjusted, but unphysical."
            p += 0.001
            
        T = Tstag_behindNSW * (p/pstag)**((gamma - 1.)/gamma)
        T = Tstag_behindNSW * (p/pstag_behindNSW)**((gamma - 1.)/gamma)
        
        T_nodes.append(T)
        if p > pstag_behindNSW:
            p = pstag_behindNSW
        M = math.sqrt( 2. / (gamma - 1.) * ((p/pstag_behindNSW)**((gamma - 1.)/(-gamma)) - 1.))
        M_nodes.append(M)
        S = 110.
        mu = muinf* (T/Tinf)**(3./2.) * (Tinf + S)/(T + S)
        mu_nodes.append(mu)
        #rho = p/(R * T)
        rho = (p/pstag_behindNSW)**(1./gamma) * rhostag_behindNSW
        rho_nodes.append(rho)
    
    if VERIFY:
        x_coords_ele = []
        y_coords_ele = []
        z_coords_ele = []
        for n in range(0, len(nodes)):
            x_coords_ele.append(nodes[n][0])
            y_coords_ele.append(nodes[n][1])
            z_coords_ele.append(nodes[n][2])
        
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(T_nodes)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('T of nodes on a cone using modified Newton technique')
        
        
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(rho_nodes)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Rho of nodes on a cone using modified Newton technique')
        
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(mu_nodes)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('mu of nodes on a cone using modified Newton technique')
        
    return T_nodes, mu_nodes, rho_nodes, M_nodes
        

def GetVelocities(V_vector, normals, centroids, Cps, pinf, Minf, Tinf, gamma, R, VERIFY):
    if gamma == 1.:
        print "Gamma cannot be 1! Adjusted to prevent singularities."
        gamma += 0.001
        
    if gamma == 0.:
        print "Gamma cannot be 0! Adjusted to prevent singularities. Results will be unphysical until you get your shit together with the input data."
        gamma += 0.001
        
    velocities = []                                                             
    pressures = []
    pstag = pinf * (1. + (gamma - 1.)/2. * Minf**2)**(gamma/(gamma - 1.))   
    pstag_behindNSW = pstag * (((gamma + 1.)*Minf**2)/ (2. + (gamma - 1.)*Minf**2))**(gamma/(gamma - 1.)) * ((gamma + 1.)/(2. * gamma * Minf**2 - (gamma - 1.)))**(1./(gamma - 1.))
    Tstag = Tinf * (1. + (gamma - 1.)/2. * Minf**2)
    Tstag_behindNSW = Tstag
    rhostag_behindNSW = pstag_behindNSW / (R * Tstag_behindNSW)
    if rhostag_behindNSW > rhoinf * 6.:
        rhostag_behindNSW = rhoinf * 6.

    for c in range(0, len(centroids)):
        if Cps[c] == 0.:
            Cps[c] = - 1./Minf**2
        n = normals[c]
        p = pinf * (Cps[c] * gamma * Minf**2 / 2. + 1.)
        if p == 0.:
            "Zero pressure detected on the windward side! Check the suitability of the mesh and input data. Adjusted, but unphysical."
            p += 0.001
        T = Tstag_behindNSW * (p/pstag_behindNSW)**((gamma - 1.)/gamma)
        pressures.append(p)
        #print 2. / (gamma - 1.) * ((p/pstag_behindNSW)**((gamma - 1.)/(-gamma)) - 1.), p/pstag_behindNSW
        if p/pstag_behindNSW > 1.:
            p = pstag_behindNSW
        M = math.sqrt( 2. / (gamma - 1.) * ((p/pstag_behindNSW)**((gamma - 1.)/(-gamma)) - 1.))
        #M = math.sqrt(2./(gamma * Cps[c]) * ((p/pinf) - 1.))
        if (gamma * R * T ) < 0.:
            print "The product of gamma * R * T is negative. Calculation will be terminated."
        V_mag = math.sqrt(gamma * R * T ) * M
        #print M, (gamma * Cps[c]),  ((p/pinf) - 1.), p, V_mag
        if Mag(n) == 0.:
            print "Normal of a zero magnitude detected. Adjusted to prevent singularities. Please fix your mesh for physical results."
            n[0] += 0.001
        projection = Dot(V_vector, n) / Mag(n) * 1. / Mag(n)
        v_proj = [n[0] * projection, n[1] * projection, n[2] * projection]
        v_rej = [V_vector[0] - v_proj[0], V_vector[1] - v_proj[1], V_vector[2] - v_proj[2] ]
        v_rej_mag = (v_rej[0]**2 + v_rej[1]**2 + v_rej[2]**2)**0.5
        if v_rej_mag > 0.:
            v_rej = [v_rej[0]/ v_rej_mag, v_rej[1]/ v_rej_mag, v_rej[2]/ v_rej_mag]
            V = [V_mag * v_rej[0], V_mag * v_rej[1], V_mag * v_rej[2]]
        else:
            V = [0., 0., 0.]
        velocities.append(V)
    
    if VERIFY:
        x_coords_ele = []
        y_coords_ele = []
        z_coords_ele = []
        u_vel_mag = []
        v_vel_mag = []
        w_vel_mag = []
        for c in range(0, len(centroids)):
            if centroids[c][0] > -2.:
                x_coords_ele.append(centroids[c][0])
                y_coords_ele.append(centroids[c][1])
                z_coords_ele.append(centroids[c][2])
                u_vel_mag.append(velocities[c][0])
                v_vel_mag.append(velocities[c][1])
                w_vel_mag.append(velocities[c][2])
        
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(u_vel_mag)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('U velocity on a cone using modified Newton technique')
        
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(v_vel_mag)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('V velocity on a cone using modified Newton technique')
        
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(w_vel_mag)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_ele, y_coords_ele, z_coords_ele, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('W velocity on a cone using modified Newton technique')
        
        
    return velocities, pressures


def SeparateCentroidsByDx(centroids, xmin, xmax, N):
    centroids_ordered = []
    centroids_ordered_dxs = []
    dx = abs(-xmin + xmax)/N
    for i in range(0, N+1):
        centroids_ordered_dxs.append(xmax - (i*dx))
        centroids_ordered.append([])
        
    for i in range(0, len(centroids)):
        min_x_dist = 1e9
        dx_region = -1
        for j in range(0, N+1):
            x_dist = abs(centroids_ordered_dxs[j] - centroids[i][0])
            if x_dist < min_x_dist:
                min_x_dist = x_dist
                dx_region = j
                
        centroids_ordered[dx_region].append([i, centroids[i]])
    
    centroids_ordered_final = []
    centroids_ordered_dxs_final = []
    for i in range(0, len(centroids_ordered)):
        if len(centroids_ordered[i]) > 0.:
            centroids_ordered_final.append(centroids_ordered[i])
            centroids_ordered_dxs_final.append(centroids_ordered_dxs[i])
            
    return centroids_ordered_final, centroids_ordered_dxs_final
        
    
def FindVelVectorNotOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids):
    dmin = 1e10
    for e in range(0, len(centroids)):
        centroid = centroids[e]
        d = ((coordinates[0] - centroid[0])**2 + (coordinates[1] - centroid[1])**2 + (coordinates[2] - centroid[2])**2)**0.5
        if d < dmin:
            dmin = d
            c = e
            found_element = e
            
    local_velocity = [us_mesh[c], vs_mesh[c], ws_mesh[c]]
    return local_velocity, found_element


def FindVelVectorOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids_ordered_local, centroids):
    dmin = 1e10
    for e in range(0, len(centroids_ordered_local)):
        centroid = centroids_ordered_local[e][1]
        d = ((coordinates[0] - centroid[0])**2 + (coordinates[1] - centroid[1])**2 + (coordinates[2] - centroid[2])**2)**0.5
        if d < dmin:
            dmin = d
            index = centroids_ordered_local[e][0]
            found_element = index
            
    local_velocity = [us_mesh[index], vs_mesh[index], ws_mesh[index]]
    return local_velocity, found_element


def DetermineDt(current_coord, u_vector, element, connectivity):
    element_nodes = connectivity[element]
    max_x = -1e9
    for i in range(0, len(element_nodes)):
        node = nodes[element_nodes[i]]
        if node[0] > max_x:
            max_x = node[0]
    dist_to_max_x = abs(max_x - current_coord[0])
    u_vec = abs(u_vector[0])
    if round(u_vec, 2) != 0.:
        dt = dist_to_max_x / u_vec
    else:
        dt = 1e-6
    return dt


def TraceStreamlinesBack(max_x_user, nodes, nodes_relevance, us, vs, ws, dt_end, stag_node, epsilon, us_mesh, vs_mesh, ws_mesh, centroids, centroids_ordered, centroids_ordered_dxs, normals, connectivity, xmin, xmax, ymin, ymax, zmin, zmax, VERIFY, PRINT):
    node_paths_coord = []
    node_paths_elem = []
    nodes_resolved = []
    nodes_resolved_idxs = []
    intersections = []
    streambacked_nodes = []
    #iter_number = int(length / (min(us)*dt))
    threshold_vel = 0.
    #nodes = nodes[:33]
    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Backtracing streamlines to SP from nodes')
        ax.scatter(stag_node[0], stag_node[1], stag_node[2], s=10, c = 'g')
    
    #nodes = nodes[:17]
    for n in range(0, len(nodes)):
        s_to_travel = epsilon/2.
        #if VERIFY:                                                              
        #    ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=1, c = 'b')
        RESOLVED = False
        FAILED = False
        if PRINT:
            print 'Node: ', n
        elif n % 100 == 0:
            print 'Node: ', n
        local_u_vector = [us[n], vs[n], ws[n]]
        
        coordinates = copy.deepcopy(nodes[n])
        dist_to_stag_point = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
        nit_crit = int(dist_to_stag_point*20000) 
        nit_crit = dist_to_stag_point / s_to_travel * 2
        if PRINT:
            print dist_to_stag_point, coordinates, stag_node, local_u_vector
        
        
        if local_u_vector[0] < threshold_vel and local_u_vector[1] < threshold_vel and local_u_vector[2] < threshold_vel:
            if PRINT:
                print 'Node in shadow: ', n                                         #THESE ARE ALL FAILURE OF TRACING POSSIBILITIES
            #if VERIFY:
            #    ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=4, c = 'r')
        elif nodes_relevance[n] == 0:
            if PRINT:
                print 'Node irrelevant: ', n                                        #THESE ARE ALL FAILURE OF TRACING POSSIBILITIES
            #if VERIFY:
            #    ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=4, c = 'r')
                
            
        elif nodes[n][0] > max_x_user:
            if PRINT:
                print 'Node irrelevant: ', n                                        #THESE ARE ALL FAILURE OF TRACING POSSIBILITIES
            #if VERIFY:
            #    ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=4, c = 'r')
            
                
        else: 
            crossed = CrossedStagnationLine(coordinates, stag_node, epsilon)[0]
            if dist_to_stag_point < epsilon:
                path_coord = [[stag_node[0], stag_node[1], stag_node[2]]]
                local_u_vector, element = FindVelVectorNotOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids)
                path_element = [element]                                        #THESE ALL INDICATE THAT THE POINT ALREADY IS IN A STAGNATION REGION
                
            elif crossed == True:
                path_coord = [[coordinates[0], coordinates[1], coordinates[1]]]
                local_u_vector, element = FindVelVectorNotOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids)
                path_element = [element]                                        #THESE ALL INDICATE THAT THE POINT ALREADY IS IN A STAGNATION REGION                                          #THESE ALL INDICATE THAT THE POINT ALREADY IS IN A STAGNATION REGION
            else:
                local_u_vector, element = FindVelVectorNotOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids)
                path_coord = [[coordinates[0], coordinates[1], coordinates[2]]]
                path_element = [element]
                
                crossed = CrossedStagnationLine(coordinates, stag_node, epsilon)[0]
                it = 0.
                dt = 1.*dt_end 
                dist_to_stag_point_previous = dist_to_stag_point
                while crossed == False:
                    if dist_to_stag_point_previous > 2.*epsilon:
                        s_to_travel = epsilon
                    else:
                        s_to_travel = epsilon/2.
                        
                    '''
                    if FIRST_LOOP:
                        dist_to_stag_point_previous = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
                        dist_to_stag_point_previous_x = (coordinates[0] - stag_node[0])
                        dist_to_stag_point_previous_y = (coordinates[1] - stag_node[1])
                        dist_to_stag_point_previous_z = (coordinates[2] - stag_node[2])
                        #dt = dist_to_stag_point_previous / us[n] 
                        FIRST_LOOP = False
                    else:
                        element = path_element[-1]
                        dist_to_stag_point_previous = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
                        dist_to_stag_point_previous_x = (coordinates[0] - stag_node[0])
                        dist_to_stag_point_previous_y = (coordinates[1] - stag_node[1])
                        dist_to_stag_point_previous_z = (coordinates[2] - stag_node[2])
                        
                        #dt = DetermineDt(coordinates, local_u_vector, element, connectivity, centroids, normals)
                        #dt = DetermineDt(coordinates, local_u_vector, element, connectivity)
                        #print dt
                    '''
                    #print  (local_u_vector)
                    '''
                    if max(abs_vector) != 0.:
                        if max(abs_vector) == abs_vector[0]:
                            dt = abs(dist_to_stag_point_previous_x) / (abs_vector[0])*0.25
                        if max(abs_vector) == abs_vector[1]:
                            dt = abs(dist_to_stag_point_previous_y) / (abs_vector[1])*0.25
                        if max(abs_vector) == abs_vector[2]:
                            dt = abs(dist_to_stag_point_previous_z) / (abs_vector[2])*0.25
                        
                    else:
                        dt = 0.001* dt
                    ''' 
                    #dt = dt_end*20.
                    total_u = (local_u_vector[0]**2  + local_u_vector[1]**2  + local_u_vector[2]**2)**0.5
                    if total_u == 0.:
                        dt = 0.00001
                    else:
                        dt = s_to_travel / total_u
                        
                    '''
                    if dist_to_stag_point_previous < 2.0 * dt * max(local_u_vector):
                        dt = 1.*dt_end
                        #dt = dist_to_stag_point_previous / max(local_u_vector)
                        
                    if dist_to_stag_point_previous < 1.0 * dt * max(local_u_vector):
                        dt = 0.25*dt_end
                        #dt = dist_to_stag_point_previous / max(local_u_vector)
                        
                    if dist_to_stag_point_previous < 0.1 * dt * max(local_u_vector):
                        dt = 0.1*dt_end
                        #dt = dist_to_stag_point_previous / max(local_u_vector)
                        
                    if dist_to_stag_point_previous < 0.01 * dt * max(local_u_vector):
                        dt = 0.01*dt_end
                        #dt = dist_to_stag_point_previous / max(local_u_vector)
                
                    '''
                    coordinates[0] = coordinates[0] - local_u_vector[0]*dt
                    coordinates[1] = coordinates[1] - local_u_vector[1]*dt
                    coordinates[2] = coordinates[2] - local_u_vector[2]*dt
                    

                    LOC = 0
                    min_loc = 0
                    max_loc = 1
                    for i in range(0, len(centroids_ordered_dxs)):
                        if coordinates[0] > centroids_ordered_dxs[i] and not LOC:
                            LOC = 1
                            min_loc = i-1
                            max_loc = i
                                
                    centroids_ordered_local = centroids_ordered[min_loc] + centroids_ordered[max_loc]
                    
                    local_u_vector, element = FindVelVectorOrdered(coordinates, us_mesh, vs_mesh, ws_mesh, centroids_ordered_local, centroids)
                    
                    normal = normals[element]
                    centroid = centroids[element]
                                        
                    normal = [normal[0]/Mag(normal), normal[1]/Mag(normal), normal[2]/Mag(normal)]
                    diff = [coordinates[0] - centroid[0], coordinates[1] - centroid[1], coordinates[2] - centroid[2]]
                    dot_scalar = Dot(diff, normal)
                    element_proj = [coordinates[0] - dot_scalar * normal[0], coordinates[1] - dot_scalar * normal[1], coordinates[2] - dot_scalar * normal[2] ]     
                    coordinates = element_proj
                    
                    if VERIFY:
                        if nodes[n][0] < 1.0:                                  #Only these streamlines will be shown for clarity
                            ax.scatter(coordinates[0], coordinates[1], coordinates[2], s=2, c = 'k')
                    

                    crossed = CrossedStagnationLine(coordinates, stag_node, epsilon)[0]
                    if crossed:
                        dist_to_stag_point_current = CrossedStagnationLine(coordinates, stag_node, epsilon)[1]
                        #coordinates = path_coord[-1]
                        dist_to_stag_point_previous = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
                        dt = dt * dist_to_stag_point_previous / (dist_to_stag_point_previous + dist_to_stag_point_current)
                        
                        normal = normals[element]
                        centroid = centroids[element]
                                            
                        normal = [normal[0]/Mag(normal), normal[1]/Mag(normal), normal[2]/Mag(normal)]
                        diff = [coordinates[0] - centroid[0], coordinates[1] - centroid[1], coordinates[2] - centroid[2]]
                        dot_scalar = Dot(diff, normal)
                        element_proj = [coordinates[0] - dot_scalar * normal[0], coordinates[1] - dot_scalar * normal[1], coordinates[2] - dot_scalar * normal[2] ]     
                        coordinates = element_proj
                        
                        if VERIFY:
                            ax.scatter(coordinates[0], coordinates[1], coordinates[2], s=4, c = 'g')
                            
                        nodes_resolved.append(coordinates)
                        nodes_resolved_idxs.append(n)
                        path_element.append(element)
                        path_coord.append([coordinates[0], coordinates[1], coordinates[2]])
                        RESOLVED = True
                        break

                    path_element.append(element)
                    path_coord.append([coordinates[0], coordinates[1], coordinates[2]])
                    it += 1.
                    
                    if coordinates[0] < xmin or coordinates[0] > xmax or coordinates[1] < ymin or coordinates[1] > ymax or coordinates[2] < zmin or coordinates[2] > zmax:
                        FAILED = True
                        if PRINT:
                            print 'Out of bounds: ', coordinates, 
                        break   
                    if it > nit_crit:
                        FAILED = True
                        if PRINT:
                            print 'Too many iterations: ', 
                        break
                if FAILED:
                    if PRINT:
                        print 'Could not find tracing for node:', n
                    if VERIFY:
                        ax.scatter(nodes[n][0], nodes[n][1], nodes[n][2], s=4, c = 'r')
            if RESOLVED:
                node_paths_coord.append(path_coord)
                node_paths_elem.append(path_element)
                intersections.append([coordinates[0], coordinates[1], coordinates[2]])                    #THIS WILL BE WRONG, YOU HAVE TO REDO THIS AFTER YOU TEST THE REST
                streambacked_nodes.append(nodes[n])
        
    return intersections, node_paths_elem, node_paths_coord, nodes_resolved, nodes_resolved_idxs
        


def EpsilonCircTransformation(nodes_resolved):
    r_ns = []
    betas = []
    for n in range(0, len(nodes_resolved)):
        node = nodes_resolved[n]
        r_n = math.sqrt(node[1]**2 + node[2]**2)**0.5
        r_ns.append(r_n)
        beta = math.atan(-node[1]/node[2])
        betas.append(beta)
        
    return r_ns, betas
        

def CrossedStagnationLine(coordinates, stag_node, epsilon):
    dist_to_stag_point = ((coordinates[0] - stag_node[0])**2 + (coordinates[1] - stag_node[1])**2 + (coordinates[2] - stag_node[2])**2)**0.5
    if dist_to_stag_point <= epsilon:
        return [True, dist_to_stag_point]
    else:
        return [False]


def FindLength(nodes):
    xmax = -1e9
    for n in range(0, len(nodes)):
        if nodes[n][0] > xmax:
            xmax = nodes[n][0]             
    return xmax

def FindBounds(nodes):
    xmax = -1e9
    xmin = 1e9
    ymax = -1e9
    ymin = 1e9
    zmax = -1e9
    zmin = 1e9
    for n in range(0, len(nodes)):
        if nodes[n][0] > xmax:
            xmax = nodes[n][0] 
        if nodes[n][0] < xmin:
            xmin = nodes[n][0]
        
        if nodes[n][1] > ymax:
            ymax = nodes[n][1] 
        if nodes[n][1] < ymin:
            ymin = nodes[n][1]
            
        if nodes[n][2] > zmax:
            zmax = nodes[n][2] 
        if nodes[n][2] < zmin:
            zmin = nodes[n][2]
            
    return xmin, xmax, ymin, ymax, zmin, zmax


    
def FindStagnationPointFromPressure(nodes, pressures):
    pmax = 0.
    for n in range(0, len(nodes)):
        p = pressures[n]
        if p > pmax:
            stag_point = nodes[n]
            stag_p = p
            pmax = p
            nstag = n
            
    return [stag_point], [stag_p], [nstag]
        
def FindStagnationPointFromGeometry(nodes, pressures):
    xmax = 1e9
    stag_points = []
    stag_ps = []
    nstags = []
    for n in range(0, len(nodes)):
        x = nodes[n][0]
        if x < xmax:
            stag_point = nodes[n]
            xmax = x
            
    for n in range(0, len(nodes)):
        x = nodes[n][0]
        if x == xmax:
            stag_point = nodes[n]
            stag_points.append(stag_point)
            stag_ps.append(pressures[n])
            nstags.append(n)
    return stag_points, stag_ps, nstags


def FindStagnationPointGeneral(centroids, u_velocity, v_velocity, w_velocity):
    return
    
    
    

def AssgnStagnationPoint(node, nodes, Minf, pinf, rhoinf, Tinf, p_nodes, rho_nodes, T_nodes, mu_nodes, u_nodes, v_nodes, w_nodes, VERIFY):
    gamma = 1.4
    R = 287.05
    for n in range(0, len(nodes)):
        if nodes[n] == node:
            #pstag = pinf * (1. + (gamma - 1.)/2. * Minf**2)**(gamma/(gamma - 1.))   
            #pstag_behindNSW = pstag * (((gamma + 1.)*Minf**2)/ (2. + (gamma - 1.)*Minf**2))**(gamma/(gamma - 1.)) * ((gamma + 1.)/(2. * gamma * Minf**2 - (gamma - 1.)))**(1./(gamma - 1.))
            #Tstag = Tinf * (1. + (gamma - 1.)/2. * Minf**2)
            #Tstag_behindNSW = Tstag
            #rhostag_behindNSW = pstag_behindNSW / (R * Tstag_behindNSW)
            #if rhostag_behindNSW > rhoinf * 6.:
            #    rhostag_behindNSW = rhoinf * 6.
            
            #S = 110.
            #mustag_behindNSW = muinf* (Tstag_behindNSW/Tinf)**(3./2.) * (Tinf + S)/(Tstag_behindNSW + S)
            
            #p_nodes[n] = pstag_behindNSW
            #T_nodes[n] = Tstag_behindNSW
            #mu_nodes[n] = mustag_behindNSW
            #rho_nodes[n] = rhostag_behindNSW
            Mach_behindNSW = math.sqrt(((gamma - 1.0)*Minf**2 + 2.0)/(2. * gamma * Minf**2 - gamma + 1.))
            u_nodes[n] = Mach_behindNSW * math.sqrt(gamma * R * T_nodes[n])
            v_nodes[n] = 0.
            w_nodes[n] = 0.
    '''            
    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(node[0], node[1], node[2], s=10, c = 'r')
    '''
    return p_nodes, T_nodes, mu_nodes, rho_nodes, u_nodes, v_nodes, w_nodes
        

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
    

def DetectIrrelevantNodes(nodes, p_nodes, pinf, Minf, gamma, VERIFY):
    relevant_nodes = []
    for n in range(0, len(nodes)):
        p = p_nodes[n]
        Cp = (p - pinf) / (0.5 * gamma * pinf * Minf**2)
        if Cp > 0.:
            relevant_nodes.append(1)
            
        else:
            relevant_nodes.append(0)
            
    if VERIFY:
        fig = plt.figure()
        ax = Axes3D(fig)
        for n in range(0, len(nodes)):
            node = nodes[n]
            if relevant_nodes[n] == 0:                
                ax.scatter(node[0], node[1], node[2], s=10, c = 'r')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Irrelevant nodes detected')
    return relevant_nodes    

    


############################################################################### HYPERSONICS, VISCOUS

def FindInitialMetricCoefficientsCircle(betas, epsilon, epsilon_nodes, tol, VERIFY):        
    hs = []    
    ss = []    
    sequence = []    
    regions = []
    #This function looks over all the betas of the backtraced streamlines and clusters them within a certain beta tolerance
    for i in range(0, len(betas)):
        NEW_REGION = True
        for j in range(0, len(regions)):
            if abs(betas[i] - regions[j]) < tol:                                #If the stream has beta within tolerance of another existing region, do not create a new region
                NEW_REGION = False
                
        if NEW_REGION:                                                          #If there is no region which could take in the given backtraced streamline, create a new region
            regions.append(betas[i])
            
    if VERIFY:
        xs = []
        ys = []
        for i in range(0, len(regions)):
            x, y = 0.159 * math.cos(regions[i]), 0.159 * math.sin(regions[i])
            xs.append(x)
            ys.append(y)
        plt.scatter(xs, ys, c = 'k')
        plt.show()

    temp_regions = []
    
    #Create a temporary array to hold region beta information
    for i in range(0, len(regions)):
        temp_regions.append(regions[i])
    
    
    #This loop determines the order of the regions starting from the bottom of the epsilon line
    for i in range(0, len(regions)): 
        min_region = 1e9
        for j in range(0, len(temp_regions)): 
            region = temp_regions[j]
            if region < min_region:                                             #Finds the smallest beta existing among the remainder of the regions
                min_region = region
                min_regions_index = j
                
        #temp_regions[j] = 1e10
        temp_regions.remove(min_region)                                         #Removes the corresponding smallest beta region from the temporary region field
        sequence.append(regions.index(min_region))                              #Sequence array adds the index of the smallest best region
    
    #Initializes the array of metric coefficients
    for i in range(0, len(regions)):
        hs.append(epsilon)
        ss.append(epsilon)
        
    #This loop calculates the matric coefficient for each region according to the sequence previously determined
    for k in range(1, len(regions)):
        region = regions[sequence[k]]                                           #Take a region according to the found sequence
        region_old = regions[sequence[k-1]]                                     #Take the beta of the previous region considered (first, most bottom region has h of 0)
        delta_Beta = region - region_old                                        #Calculate the difference in beta between the two regions
        delta_Beta = region - regions[sequence[0]]                              #Calculate the difference in beta between the two regions
        delta_S = delta_Beta * epsilon                                          #The distance travelled along the line is equal to the arc length between the two regions
        index = sequence[k]                                                     #Absolute (original) index of the current region is found from the sequence array
        hs[index] = delta_S / delta_Beta                                        #The metric coefficient at this index is then computed
        ss[index] = delta_S
        #print regions[k], hs[index]
        
    
    hs_nodes = []
    ss_nodes = []
    #The array above only contains region information, not node information --> this has to be transformed to the node format
    for n in range(0, len(epsilon_nodes)):                              
        beta_node = betas[n]                                                    #Evaluate the beta of the current node in the right node order
        ASSIGNED = False
        for r in range(0, len(regions)):                                        #Scan over all regions defined above
            if abs(beta_node - regions[r]) < tol:                               #Find the region to which the current node belongs based on the tolerance                
                hs_nodes.append(hs[r])                                          #Assign the correct metric coefficient from the found 
                ss_nodes.append(ss[r])
                ASSIGNED = True
                break
        if not ASSIGNED:
            print "Metric coefficient assignment failure for node: ", n
    return hs_nodes, ss_nodes



def FindInitial_wdydbetaz(hs, GradFs, Fxs, us, vs, ws, epsilon_nodes, list_of_crossed_elements):
    wdydbetazs = []
    for e in range(0, len(epsilon_nodes)):
        el = list_of_crossed_elements[e][-1]
        #print el
        V = math.sqrt(us[e]**2 + vs[e]**2 + ws[e]**2)
        #print hs[e], V, Fxs[el], GradFs[el]
        wdydbetaz = hs[e] * V * Fxs[el] / abs(GradFs[el])
        wdydbetazs.append(wdydbetaz)
        
    return wdydbetazs


def FindInitial_wdydbetazStable(epsilon_nodes, nodes, T_w, T_nodes, mu_nodes, rho_nodes, M_nodes, centroids, u_elements, v_elements, w_elements, hs, ss, cp, k, R, gamma, stg_idx):
    wdydbetazs = []
    for e in range(0, len(epsilon_nodes)):
        node_index = stg_idx
        coordinates = epsilon_nodes[e]
        velocity_e, element = FindVelVectorNotOrdered(coordinates, u_elements, v_elements, w_elements, centroids)
        u_e = Mag(velocity_e)
        T_e = T_nodes[node_index]
        rho_e = rho_nodes[node_index]
        mu_e = mu_nodes[node_index]
        M_e = M_nodes[node_index]
        Pr = mu_e * cp / k
        V_e = math.sqrt(T_e * R * gamma) * M_e
        Re_e = mu_e * V_e / mu_e
        
        Re_star, rho_star, mu_star, T_star = ReturnEckert(False, T_w, T_e, M_e, gamma, Pr, mu_e, rho_e, Re_e)
        
        wdydbetaz = rho_star * mu_star * u_e * hs[e]**2. * ss[e]/4.
        wdydbetazs.append(wdydbetaz)
        
    return wdydbetazs

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


def MarchDownstream(stg_coord, stg_cnds, wdydbetaz_ini, list_of_crossed_elements, list_of_coordinates, h, GradFs, Fxs, us_ele, vs_ele, ws_ele, T_ele, mu_ele, rho_ele, M_ele, w_fit, T_w, cp, k, wds, metric_coeffs, INTERPOLATE_N, OWN_TRANSITION, PERZHIKAR):
        
    d1s = wds[0]
    d2s = wds[1]
    d3s = wds[2]
    d4s = wds[3]
    d5s = wds[4]
    d6s = wds[5]
    
    wdydbetaz = wdydbetaz_ini
    coord_ini = list_of_coordinates[0]
    
    prev_Re_theta = 0.
    Pr = 1.0
    #stg_idx = T_ele.index(max(T_ele))
    #V_e_ini = (us_ele[stg_idx]**2 + vs_ele[stg_idx]**2  + ws_ele[stg_idx]**2)**0.5
    #T_e_ini = T_ele[stg_idx]
    #M_e_ini = M_ele[stg_idx]
    #rho_e_ini = rho_ele[stg_idx]
    #mu_e_ini = mu_ele[stg_idx]
    
    T_e_ini = Tstag_behindNSW
    rho_e_ini = rhostag_behindNSW
    p_e_ini = pstag_behindNSW
    V_e_ini = Vstag_behindNSW 
    M_e_ini = Vstag_behindNSW / math.sqrt(1.4 * 287.05 * T_e_ini)
    T_e_ini = Tstag_behindNSW / (1. + (gamma - 1.0)/2. * M_e_ini**2 )
    S = 110.
    mu_e_ini = mustag_behindNSW / (1. + (gamma - 1.0)/2. * M_e_ini**2 )
    mu_e_ini = mustag_behindNSW * (Tstag_behindNSW/Tstag_behindNSW)**(3./2.) * (Tstag_behindNSW + S)/(Tstag_behindNSW + S)
    
    Re_e_ini = rho_e_ini * V_e_ini / mu_e_ini 
        
    Re_star, rho_star, mu_star, T_star = ReturnEckert(False, T_w, T_e_ini, M_e_ini, gamma, Pr, mu_e_ini, rho_e_ini, Re_e_ini)
    rho_star_old = rho_star
    mu_star_old = mu_star
    rho_e_old = rho_e_ini
    V_e_old = V_e_ini
    u_e_old = 0.
    h_old = h
    coord_old = coord_ini
    
    all_points = []
    hs = []
    dotq_ws = []
    skin_cfs = []
    deltas = []
    thetas = []
    V_es = []
    rho_es = []
    mu_stars = []
    dss = []
    rho_stars = []
    integs = []
    transition_points = []
    transitioned = False
    NOTIFY_TRANSITION = True
    FIRST = True
    integ_tracing = 0.
    FIRST_turb = False
    #dwdy = 2d1y + d2 + d5z
    #dwdz = 2d3z + d4 + d5y

    for i in range(1, len(list_of_crossed_elements)):
        h_perzhikar = metric_coeffs[i]
        el = list_of_crossed_elements[i]
        el_old = list_of_crossed_elements[i-1]
        V = (us_ele[el]**2 + vs_ele[el]**2 + ws_ele[el]**2)**0.5
        V_old = (us_ele[el_old]**2 + vs_ele[el_old]**2 + ws_ele[el_old]**2)**0.5
        if V_old == 0.:
            V_old = 10.
            
        u_e = abs(us_ele[el])
        
        metric_coeff = metric_coeffs[i]
        
        coord = list_of_coordinates[i]
        x_coord = coord[0]
        ds = ((coord[0] - coord_old[0])**2 + (coord[1] - coord_old[1])**2 + (coord[2] - coord_old[2])**2)**0.5
        if V == 0.:
            V = 0.1
            
        if not PERZHIKAR:
            wdydbetaz_prev = wdydbetaz
            term_around_prev = 0.
            if (coord[1]**2 + coord[2]**2)**0.5 > 10.*epsilon:
                
                y_track += ds * vs_ele[el] / V
                z_track += ds * ws_ele[el] / V
                
                coord_old = list_of_coordinates[i-1]
                el_evenolder = list_of_crossed_elements[i-2]
                #dw_old =  ws_ele[el_old] - ws_ele[el_evenolder]
                #dy_old = list_of_coordinates[i-1][1] - list_of_coordinates[i-2][1]
                #dz_old = list_of_coordinates[i-1][2] - list_of_coordinates[i-2][2]
                
                #print d1s[el_old], coord[1]
                #dwdy = 2.*d1s[el_old]*coord[1] - d2s[el_old] - d5s[el_old]*coord[2]
                #dwdz = 2.*d3s[el_old]*coord[2] + d4s[el_old] - d5s[el_old]*coord[1]
            
                dwdy = 2.*d1s[el_old]*coord_old[1] + d2s[el_old] + d5s[el_old]*coord_old[2]
                dwdz = 2.*d3s[el_old]*coord_old[2] + d4s[el_old] + d5s[el_old]*coord_old[1]
                
                dwdy_new = 2.*d1s[el]*coord_old[1] + d2s[el] + d5s[el]*coord_old[2]
                dwdz_new = 2.*d3s[el]*coord_old[2] + d4s[el] + d5s[el]*coord_old[1]
                
                dw_old =  ws_ele[el_old] - ws_ele[el_evenolder]
                dy_old = list_of_coordinates[i-1][1] - list_of_coordinates[i-2][1]
                dz_old = list_of_coordinates[i-1][2] - list_of_coordinates[i-2][2]
                
                dw_new =  ws_ele[el] - ws_ele[el_old]
                dy_new = list_of_coordinates[i][1] - list_of_coordinates[i-1][1]
                dz_new = list_of_coordinates[i][2] - list_of_coordinates[i-1][2]
                
                dwdy = dw_old/dy_old
                dwdz = dw_old/dz_old
                
                dwdy = dw_new/dy_new
                dwdz = dw_new/dz_new
                
                dwdylh = dw_new / y_track
                dwdzlh = dw_new / z_track
                
                
                #term_around_prev =  1./V_old * (dw_old/dy_old + dw_old/dz_old) * ds 
                #term_around_prev =  1./V_old * (-dwdy + dwdz) * ds 
                
                
                #wdydbetaz = wdydbetaz_prev * (1. + (term_around_prev))
                '''
                
                c2 = 1./V * (-dwdy_new + dwdz_new)
                c1 = 1./V_old * (-dwdy + dwdz)
                
                #print (math.log(wdydbetaz_prev) + ((c2 - c1)*ds))
                avec = 0.5 * (c2+c1)
                wdydbetaz = math.exp(math.log(wdydbetaz_prev) + (avec*ds))
                '''
                
                wdydbetaz = wdydbetaz_prev + wdydbetaz_prev * 1./V_old * (dwdy + dwdz) * ds 
                shittyterm1 = dwdy
                shittyterm2 = dwdz
                
            else:
                y_track += ds * vs_ele[el] / V
                z_track += ds * ws_ele[el] / V
                #if round(wdydbetaz_ini, 8) == 0.:
                #    wdydbetaz_ini =  (abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)) * V *  Fxs[el] / abs(GradFs[el])
                #wdydbetaz = wdydbetaz_ini
                wdydbetaz = (abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)) * V *  Fxs[el] / abs(GradFs[el])
                if i == 1:
                    wdydbetaz_ini = wdydbetaz
                    print "\n"
                    print "wdydbetaz_ini: ", wdydbetaz_ini
                term_around_prev = 0
                shittyterm1 = 0
                shittyterm2 = 0
            
            #wdydbetaz = max([wdydbetaz, wdydbetaz_ini])
    
            
            h1 = abs(GradFs[el]) / Fxs[el]/ V * (wdydbetaz)
            
            if abs(h1) > 3.:
                h1 = 3.
            #h = 0.1
            #h = abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)
            wdydbetaz_correct = h * V * Fxs[el] / abs(GradFs[el])
            print round(h, 4), round(h1, 4), '\t\t', round(wdydbetaz, 4), round(wdydbetaz_correct,4), shittyterm1, shittyterm2, ds, ws_ele[el]
            
        else:
            h = metric_coeff
            #h = abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)
        #h = abs(h1)
        #h = abs(((coord[1]-stg_coord[1])**2 + (coord[2]-stg_coord[2])**2)**0.5)
        if h <= 0.:
            break
        V_e = V
        T_e = T_ele[el]
        M_e = M_ele[el]
        rho_e = rho_ele[el]
        mu_e = mu_ele[el]
        Re_e = rho_e * V_e / mu_e 
        Re_star, rho_star, mu_star, T_star = ReturnEckert(transitioned, T_w, T_e, M_e, gamma, Pr, mu_e, rho_e, Re_e)
        #if h < 0.:
        #    h = abs(h)

        dotq_w, delta, theta, Re_theta_e, integ_tracing, skin_cf = SolveThermal(transitioned, 
                                            prev_Re_theta, cp, k, T_e, T_w, Pr, mu_e, mu_star, 
                                            rho_star, V_e, u_e, h, ds, rho_e, mu_star_old, 
                                            rho_star_old, V_e_old, u_e_old, h_old, rho_e_old, 
                                            integ_tracing, FIRST, FIRST_turb, INTERPOLATE_N, x_coord)
        
        transitioned_prev = transitioned
        transitioned = Transition(M_e, Re_theta_e, theta, V_e, rho_e, mu_e, abs(x_coord), OWN_TRANSITION)
        #if h < 0.1:
        #transitioned = False
        if transitioned != transitioned_prev:
            #integ_tracing = 0.0
            FIRST_turb = True
        else:
            FIRST_turb = False
        if transitioned and NOTIFY_TRANSITION:
            print "Transition occurred at eq. body radius of:", h, " and xcoord of ", x_coord,  " with Re_T of: ", Re_theta_e
            transition_points.append(h)
            NOTIFY_TRANSITION = False
        
        #if not transitioned:
        prev_Re_theta = Re_theta_e
        mu_star_old = mu_star
        V_e_old = V_e
        u_e_old = u_e
        rho_e_old = rho_e
        rho_star_old = rho_star

        coord_old = coord 
        
        h_old = h
        hs.append(h)
        V_es.append(V_e)
        rho_es.append(rho_e)
        mu_stars.append(mu_star)
        dss.append(ds)
        rho_stars.append(rho_star)
        integs.append(integ_tracing)
        dotq_ws.append(dotq_w)
        skin_cfs.append(skin_cf)
        deltas.append(delta)
        thetas.append(theta)
        all_points.append(list_of_coordinates[i])
        FIRST = False
    return hs, all_points, thetas, dotq_ws, deltas, transition_points, V_es, rho_es, mu_stars, rho_stars, integs, dss, skin_cfs




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
            



def SolveThermal(transitioned, prev_Re_theta, cp, k, T_e, T_w, Pr, mu_e, mu_star, 
                 rho_star, V_e, u_e, h, ds, rho_e, mu_star_old, rho_star_old, V_e_old, 
                 u_e_old, h_old, rho_e_old, integ_tracing, FIRST, FIRST_turb, 
                 INTERPOLATE_N, x):        
    u_e = V_e
    u_e_old = V_e_old
    #FIRST_turb = True
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
        if rho_star > rho_star_old:
            rho_star = rho_star_old
        if mu_star > mu_star_old:
            mu_star = mu_star_old
        
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
            if prev_Re_theta <= 0.:
                print "prev_Re_theta: ", prev_Re_theta
            N = 12.67 - 6.5 * math.log(prev_Re_theta, 10) + 1.21 * (math.log(prev_Re_theta, 10))**2 
            if N < 0.:
                N = ReturnTurbN(prev_Re_theta)
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
        integ_tot = integ #+ integ_tracing
        integ_tracing += integ

        theta_T = (c2* integ_tot)**c4 / (rho_e * u_e * h)
        #theta_T = abs(0.002763/0.903293*x)
        if FIRST_turb:
            Re_theta_e = theta_T * V_e * rho_e / mu_e
            #Re_theta_e_first = Re_theta_e
            FIRST_turb = False
        #else:
        #    Re_theta_e = Re_theta_e_first
        
        H_e = cp * T_e
        R = Pr**(1./3.)
        H_aw = H_e + R * V_e**2/2.
        H_w = cp * T_w
        #H_w = H_aw/7.5
        Re_theta_e = theta_T * V_e * rho_e / mu_e
        
        dotq_w = c1 * (Re_theta_e)**(-m) * rho_star / rho_e * (mu_star / mu_e)**m * rho_e * u_e * (H_aw - H_w) * Pr_w**-0.4
        delta = theta_T * (N + 1. + (( (N+2.) /N ) * H_w/H_aw  + 1.) * (1. + 1.29 * Pr_w**0.333 * u_e**2/(2.*H_w)  )  )
        theta = theta_T
        
        if Re_theta_e == 0.:
            Re_theta_e += 0.00001
            
        skin_cf = 2. * c1 / Re_theta_e**m
        
    #print "Theta: ", theta, "Dotq_w: ", dotq_w, "delta: ", delta, "Pr_w: ", Pr_w
    return dotq_w, delta, theta, Re_theta_e, integ_tracing, skin_cf




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



def ComputeMetricPerzhikar(nodes, node_resolved_idxs, node_path_coordinates_orig, stag_point, dist_yz_min):
    node_path_coordinates = copy.deepcopy(node_path_coordinates_orig)
    #for n in range(0, len(node_path_coordinates)):
    #    first_point = nodes[node_resolved_idxs[n]]
    #    node_path_coordinates[n] = [first_point] + node_path_coordinates[n][:]
        
    total_ss_all = []
    #For each node that is traced
    for n in range(0, len(node_resolved_idxs)):
        if n%10 == 0:
            print "Node... ", n
        total_ss = []
        node_path_coordinates_curr = node_path_coordinates[n]
        
        #For each step made during tracing
        for m in range(1, len(node_path_coordinates_curr)):
            total_s = 0.
            
            #For each coordinate it has passed
            for k in range(m, len(node_path_coordinates_curr)):
                x_i = node_path_coordinates_curr[k][0]
                y_i = node_path_coordinates_curr[k][1]
                z_i = node_path_coordinates_curr[k][2]
                x_imin1 = node_path_coordinates_curr[k-1][0]
                y_imin1 = node_path_coordinates_curr[k-1][1]
                z_imin1 = node_path_coordinates_curr[k-1][2]
                total_s += ((x_i - x_imin1)**2 + (y_i - y_imin1)**2 + (z_i - z_imin1)**2)**0.5
                
            total_ss.append(total_s)
        #print "For node ", nodes[node_resolved_idxs[n]], " with ", len(node_path_coordinates_curr), " path elements, the ss are: ", total_ss
        #print '\n'
        if len(total_ss) == 0.:
            coords_zero_path = nodes[node_resolved_idxs[n]]
            radius = ((coords_zero_path[1]-stag_point[1])**2 + (coords_zero_path[2]-stag_point[2])**2)**0.5
            new_ss_extrapolated = radius
        
        elif len(total_ss) == 1:
             new_ss_extrapolated = total_ss[0]
             
        else:
            new_ss_extrapolated = (total_ss[0] - total_ss[1] ) + total_ss[0]
            
        total_ss = [new_ss_extrapolated] + total_ss[:]
        total_ss_all.append(total_ss)
        
    hs_perzhikar_all = []
    for n in range(0, len(node_resolved_idxs)):
        #if n%100 == 0:
        #    print "Node... cycle b ", n
        node_path_coordinates_curr = node_path_coordinates[n]
        hs_perzhikar = []
        streamline_dist_min = 1e9
        nearest_streamline_order = 0
        nearest_streamline_point_idx = 0
        
        #if n == 6522:
        #    print "We have a streamline starting from node: ", node_path_coordinates_curr[0]
        
        for l in range(0, len(node_resolved_idxs)):
            node_coord = nodes[node_resolved_idxs[l]]
            node_coord_curr = node_path_coordinates_curr[0]
            dist_yz = (0.*(node_coord[0] - node_coord_curr[0])**2 + (node_coord[1] - node_coord_curr[1])**2 + (node_coord[2] - node_coord_curr[2])**2)**0.5
            
            if dist_yz < streamline_dist_min and dist_yz > dist_yz_min and node_coord[0] >= node_coord_curr[0]:
                streamline_dist_min = dist_yz
                nearest_streamline_point_idx = node_resolved_idxs[l]
                nearest_streamline_order = l
           
        #At this point, the node index from which the nearest y-z streamline originates is known
        
        node_path_coordinates_nearest = node_path_coordinates[nearest_streamline_order]
        total_s_nearest = total_ss_all[nearest_streamline_order]
        total_s_current = total_ss_all[n]
        
        #if n == 6522:
        #    print "We found the nearest streamline to originate from: ",  node_path_coordinates_nearest[0], " with a y-z distance of ", streamline_dist_min
        
        #March through the points on the streamline
        for m in range(1, len(node_path_coordinates_curr)):
            current_s = total_s_current[m]
            s_dist_min = 1e9
            nearest_point = 0
            
            for k in range(0, len(node_path_coordinates_nearest)):
                s_nearest_at_k = total_s_nearest[k]
                #print current_s, s_nearest_at_k
                s_dist  = abs(current_s - s_nearest_at_k)
                if s_dist < s_dist_min:
                    s_dist_min = s_dist
                    nearest_point = k
            
            second_nearest_point = 0
            second_s_dist_min = 1e9
            for l in range(0, len(node_path_coordinates_nearest)):
                s_nearest_at_l = total_s_nearest[l]
                second_s_dist  = abs(current_s - s_nearest_at_l)
                if second_s_dist < second_s_dist_min and second_s_dist != s_dist_min:
                    second_s_dist_min = second_s_dist
                    second_nearest_point = l
            
            #Figure out linearly interpolated point between these two
            p_n_nearest = node_path_coordinates_nearest[nearest_point]
            p_n_secondnearest = node_path_coordinates_nearest[second_nearest_point]
            
            x_p_n1 = p_n_nearest[0]
            y_p_n1 = p_n_nearest[1]
            z_p_n1 = p_n_nearest[2]
            x_p_n2 = p_n_secondnearest[0]
            y_p_n2 = p_n_secondnearest[1]
            z_p_n2 = p_n_secondnearest[2]
            
            if second_s_dist_min == 0.:
                second_s_dist_min += 1e-6
            if s_dist_min ==0.:
                s_dist_min += 1e-6
                
            x_p_n = ((1./s_dist_min) * x_p_n1 + (1./ second_s_dist_min) * x_p_n2) / ((1./s_dist_min) + (1./ second_s_dist_min)) 
            y_p_n = ((1./s_dist_min) * y_p_n1 + (1./ second_s_dist_min) * y_p_n2) / ((1./s_dist_min) + (1./ second_s_dist_min)) 
            z_p_n = ((1./s_dist_min) * z_p_n1 + (1./ second_s_dist_min) * z_p_n2) / ((1./s_dist_min) + (1./ second_s_dist_min)) 
            
            #At this point, we know which of the points, k, on the nearest streamline is the closest one to the point of interest, m
            
            
            #p_n = node_path_coordinates_nearest[nearest_point]
            p_n = [x_p_n, y_p_n, z_p_n]
            p_i = node_path_coordinates_curr[m]
            p_imin1 = node_path_coordinates_curr[m-1]
            
            #if n == 6522:
            #    print "Staring with the second step ", p_i, " after ", p_imin1, " and taking ", p_n, "from the second streamline, with an error of ", s_dist_min
             
            hbetadbeta = ((p_n[0] - p_i[0])**2 + (p_n[1] - p_i[1])**2 + (p_n[2] - p_i[2])**2)**0.5
            xpipimin1 = ((p_imin1[0] - p_i[0])**2 + (p_imin1[1] - p_i[1])**2 + (p_imin1[2] - p_i[2])**2)**0.5
            xpimin1pn = ((p_imin1[0] - p_n[0])**2 + (p_imin1[1] - p_n[1])**2 + (p_imin1[2] - p_n[2])**2)**0.5
            
            #print xpipimin1, hbetadbeta
            
            term = -(xpimin1pn**2) + xpipimin1**2 + hbetadbeta**2 
            if hbetadbeta == 0. or xpipimin1 == 0.:
                print "Error for a node ", node_path_coordinates_curr[0], " using streamline from the node ", nodes[nearest_streamline_point_idx], " at coordinates pi, pimin1 and pn ", p_i, '\t', p_imin1, '\t', p_n
                print node_path_coordinates_curr
                print node_path_coordinates_nearest
                
                
            gamma = math.acos(0.5*term / hbetadbeta / xpipimin1 )
            
            
            hdbeta = hbetadbeta * math.sin(gamma)
            if m == 1:
                factor = hbetadbeta / ((p_i[1]-stag_point[1])**2 + (p_i[2]-stag_point[2])**2)**0.5
             
            if n == 6522:
                print "Comparison: ", hbetadbeta / factor, '\t' ,((p_i[1]-stag_point[1])**2 + (p_i[2]-stag_point[2])**2)**0.5
            
            #print "The found angle is ", gamma *180./3.14 , " and the metric ", round(hdbeta, 5), " while the actual metric is ", ((p_i[1]-stag_point[1])**2 + (p_i[2]-stag_point[2])**2)**0.5, "\n \n"
            
            hs_perzhikar.append(hdbeta)
        hs_perzhikar.append(hdbeta)    
        hs_perzhikar_all.append(hs_perzhikar)
    print "Pre-computing of metric coefficients has finished."
    return hs_perzhikar_all
        
            
            
            
            

            
def WriteRunInfo(filename, data):
    
    pinf = data[0]
    Tinf = data[1]
    Mach = data[2]
    rhoinf = data[3]
    muinf = data[4]
    R = data[5]
    cp = data[6]
    gamma = data[7]
    name = data[8]
    USE_IDS = data[9]
    IDWP = data[10]
    COMPUTE_STAGNATION_POINT_FROM_GEOMETRY = data[11]
    epsilon = data[12]
    T_w = data[13]
    PERZHIKAR = data[14]
    INTERPOLATE_N = data[15]
    beta_tol = data[16]
    OWN_TRANSITION = data[17]
    max_x_user  = data[18]
    dist_yz_min = data[19]

    f = open(filename, "w")
    
    f.write("Simulation at p(inf), T(inf), M(inf), rho(inf), mu(inf), R, cp and gamma: \n" + str(pinf) + ', ' + str(Tinf) + ', ' + str(Mach) + ', ' + str(rhoinf) + ', ' + str(muinf) + ', ' + str(R) + ', ' + str(cp) + ', ' + str(gamma) + '\n')
    f.write("Mesh used: " + str(name) + '\n \n')
    
    f.write("####### Inviscid computation data " + '\n')
    f.write("IDW used: " + str(USE_IDS) + " with IDWP: " + str(IDWP) + '\n' )
    f.write("SP computed from geometry: " + str(COMPUTE_STAGNATION_POINT_FROM_GEOMETRY) + '\n\n' )
            
    f.write("####### Viscous computation data " + '\n')
    f.write("Epsilon: " + str(epsilon) + '\n')
    f.write("T wall: " + str(T_w) + '\n')
    f.write("Perzhikar used: " + str(PERZHIKAR) + " " + '\n')
    f.write("Turbulent N was interpolated rather than computed: " + str(INTERPOLATE_N) + '\n')
    f.write("Beta tolerance while projecting (Hamilton): " + str(beta_tol) + '\n')
    f.write("Own transition vector: " + str(OWN_TRANSITION) + '\n')
    f.write("Data ignored beyond: " + str(max_x_user) + '\n')
    f.write("Min. distance between streamlines (Perzhikar): " + str(dist_yz_min) + '\n')
    
    f.close()

    return True



def WriteDataForSurfaceFit(filenamevel, filenamecoord, connectivity, centroids, u_nodes, v_nodes, w_nodes, nodes):
    f = open(filenamevel, "w")
    
    for c in range(0, len(centroids)):
        conn = connectivity[c]
        vertex_1 = conn[0]
        vertex_2 = conn[1]
        vertex_3 = conn[2]
        
        u = [u_nodes[vertex_1], u_nodes[vertex_2], u_nodes[vertex_3]]
        v = [v_nodes[vertex_1], v_nodes[vertex_2], v_nodes[vertex_3]]
        w = [w_nodes[vertex_1], w_nodes[vertex_2], w_nodes[vertex_3]]
        
        line = str(u[0]) + ' ' +  str(u[1]) + ' ' +  str(u[2]) + ' ' +  str(v[0]) + ' ' +  str(v[1]) + ' ' +  str(v[2]) + ' ' +  str(w[0]) + ' ' +  str(w[1]) + ' ' +  str(w[2]) + ' ' + '\n'
        f.write(line)
    f.close()
    
    f = open(filenamecoord, "w")
    
    for c in range(0, len(centroids)):
        conn = connectivity[c]
        vertex_1 = conn[0]
        vertex_2 = conn[1]
        vertex_3 = conn[2]
        
        coord_1 = nodes[vertex_1]
        coord_2 = nodes[vertex_2]
        coord_3 = nodes[vertex_3]
        
        line = str(coord_1[0]) + ' ' + str(coord_1[1]) + ' ' + str(coord_1[2]) + ' ' + str(coord_2[0]) + ' '+ str(coord_2[1]) + ' '+ str(coord_2[2]) + ' ' + str(coord_3[0]) + ' '+ str(coord_3[1]) + ' ' + str(coord_3[2]) + '\n'
        f.write(line)
    f.close()
    return True



def WriteDataNodes(filename, data, nodes):
    f = open(filename, "w")
    
    for n in range(0, len(nodes)):
        node = nodes[n]
        data_n = data[n]    
        line = str(node[0]) + ' ' +  str(node[1]) + ' ' +  str(node[2]) + ' ' + str(data_n) + '\n'
        f.write(line)
    f.close()

    return True


def WriteDataCentroids(filename, filenamecon, centroids, connectivity):
    f = open(filename, "w")
    
    for n in range(0, len(centroids)):
        centroid = centroids[n]
        line = str(centroid[0]) + ' ' +  str(centroid[1]) + ' ' +  str(centroid[2]) + '\n'
        f.write(line)
    f.close()
    
    f = open(filenamecon, "w")
    
    for n in range(0, len(connectivity)):
        conn = connectivity[n]    
        line = str(conn[0]) + ' ' +  str(conn[1]) + ' ' +  str(conn[2]) + '\n'
        f.write(line)
    f.close()

    return True


def WriteStagPointData(filename, stag_points, stag_idxs, epsilon):
    f = open(filename, "w")
    
    i = 0
    for stag_point in stag_points:
        stag_idx = stag_idxs[i]
        line = str(stag_point[0]) + ' ' +  str(stag_point[1]) + ' ' +  str(stag_point[2]) + ' ' + str(stag_idx) + ' ' + str(epsilon) + '\n'
        f.write(line)
        i += 1
    f.close()
    
    return True
    
def WriteNormals(filename, normals):
    f = open(filename, "w")
    
    for n in range(0, len(normals)):
        normal = normals[n]    
        line = str(normal[0]) + ' ' +  str(normal[1]) + ' ' +  str(normal[2]) + '\n'
        f.write(line)
    f.close()
    
def FilterNodes(nodes, cps, centroids, connectivity, shadowed):
    windward_nodes = []
    for n in range(0, len(nodes)):                                              #Check every node
        USELESS = True
        for c in range(0, len(centroids)):                                      #Check every centroid
            conn = connectivity[c]                                              #Find all the nodes of the element
            for i in range(0, len(conn)):
                node_of_c = conn[i]
                if n == node_of_c:                                              #If the node is a part of the element 
                    if shadowed[c] == 'False':                                  #if the elemen is not shadowed
                        USELESS = False
        if not USELESS:
            windward_nodes.append('True')
        else:
            windward_nodes.append('False')

    return windward_nodes
            
def WriteBackTracing(filename_intersections, filename_nodepaths, filename_nodecoords, filename_nodesresolved, filename_epsilonnodes, intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs):
    f = open(filename_intersections, "w")
    for n in range(0, len(intersections)):
        intersection = intersections[n]
        line = str(intersection[0]) + ' ' +  str(intersection[1]) + ' ' +  str(intersection[2]) + '\n'
        f.write(line)
    f.close()
    
    f = open(filename_nodesresolved, "w")
    for n in range(0, len(nodes_resolved_idxs)):
        line = str(nodes_resolved_idxs[n]) + '\n'
        f.write(line)
    f.close()
    
    f = open(filename_epsilonnodes, "w")
    for n in range(0, len(epsilon_nodes)):
        epsilon_node = epsilon_nodes[n]
        line = str(epsilon_node[0]) + ' ' +  str(epsilon_node[1]) + ' ' +  str(epsilon_node[2]) + '\n'
        f.write(line)
    f.close()
    
    f = open(filename_nodepaths, "w")
    for n in range(0, len(node_paths_elem)):
        node_path = node_paths_elem[n]
        line = ""
        for m in range(0, len(node_path)):
            line += str(node_path[m]) + ' ' 
        
        line += '\n'
        f.write(line)
    f.close()
    
    f = open(filename_nodecoords, "w")
    for n in range(0, len(node_paths_elem)):
        node_path_coord = node_paths_coord[n]
        line = ""
        for m in range(0, len(node_path_coord)):
            coord = node_path_coord[m]
            for l in range(0, len(coord)):
                line += str(coord[l]) + ' '
                
            line += ","
        
        line += '\n'
        f.write(line)
    f.close()
    
    return True
    
    

















































############################################################################### CODE BODY


#########INPUT FILE MANAGEMENT
name = "/Users/brch/HPM/cone_tri_longer_coarse.ply"
name = "/Users/brch/HPM/sphere_coarse_shifted.ply"
name = "/Users/brch/HPM/rounded_cylinder_coarser.ply"
name = "/Users/brch/HPM/projectile_fine.ply"

name = "/Users/brch/HPM/projectile_smoother_fine.ply"
name = "/Users/brch/HPM/blunted_cone.ply"
#name = "/Users/brch/HPM/blunted_cone_20deg_stg_2_inverted.ply"
#name = "/Users/brch/HPM/cone_tri_detal.ply"
#name = "/Users/brch/HPM/projectile_finer_translatedz.ply"
name = "/Users/brch/HPM/blunted_cone_translatedz.ply"
name = "/Users/brch/HPM/blunted_cone_inverted.ply"

#name = "/Users/brch/HPM/blunted_cone_20deg_stg_2_inverted.ply"
name = "/Users/brch/HPM/blunted_cone_20deg_correctly_scaled_tiny.ply"
#name = "/Users/brch/HPM/sphere_basic_0_deg.ply"
#name = "/Users/brch/HPM/sthlol.ply"
name = "/Users/brch/HPM/testcase.ply"


path = '/Users/brch/Downloads/Blunted_cone_M6/'
path = '/Users/brch/Downloads/correctly_scaled_20_deg_tiny/'
path = '/Users/brch/Downloads/HPMtestcase/'

#path = '/Users/brch/Downloads/sphere_basic_0_deg/'
#path = '/Users/brch/Downloads/sthlol/'





#########INVISCID PART
INVSCID                                         = 0
CHECK_MESH                                      = 1
VERIFY_NEWTON                                   = 1

USE_IDS                                         = 0
VERIFY_ISD_NODES                                = 0

VERIFY_THERMODYNAMICS                           = 1
VERIFY_STAG_POINT                               = 0
WRITE_FIELDS                                    = 0


#########VISCOUS PART
VISCOUS                                         = 1

VERIFY_NODE_RELEVANCE                           = 0
VERIFY_BACKTRACING                              = 1
BACKTRACE                                       = 1
WRITE_BACKTRACING                               = 1

MARCH                                           = 0
VERIFY_HS_CLUTTERING                            = 0
VERIFY_A_COMPUTATION                            = 0

##########GENERAL SETTING
COMPUTE_STAGNATION_POINT_FROM_GEOMETRY          = 1
PERZHIKAR                                       = 1



#########USER INPUT MANAGEMENT
dt = 0.00001
epsilon = 0.0125
IDWP = 4.
gamma = 1.4
R = 287.05
#pinf = 101000.
pinf = 580.29
#Tinf = 288.
Tinf = 62.870

#muinf = 1.716e-5
muinf = 4.196e-6
rhoinf = pinf / (Tinf * R)
V_vector = [1., 0., 0.]
V_mag = (V_vector[0]**2 + V_vector[1]**2 + V_vector[2]**2)**0.5
N = 100.
T_w = 300.
cp = 1006.
k = 0.02435
Mach = 6.
gamma = 1.4
V_vector = [1., 0., 0.]                                                         #While the AOA can be changed in this direction, it is preferred (for stability) to do it with the mesh
INTERPOLATE_N = False                                                            #For turbulent heating, either compute N (more accurate, less stable) = False or interpolate N = True
OWN_TRANSITION_Re = 'nan'                                                       #User defined transition Re number (x-based)
OWN_TRANSITION_x = 0.12                                                        #User defined transition coordinate (x-based)
OWN_TRANSITION = [OWN_TRANSITION_Re, OWN_TRANSITION_x]
max_x_user = 4.0                                                                #Distance beyond which the results are not desired (i.e too downstream to be relevant) for speedup
beta_tol = 0.01                                                                 #Resolution for point cluttering on the epsilon line, for Hamilton
dist_yz_min = 0.07                                                              #Min distance between two streamlines, for Perzhikar




#########OUTPUT FILE MANAGEMENT

filename_infofile = str(path + 'INFO.txt')

filename_normals = str(path + 'normals.txt')
filename_centroids = str(path + 'centroids.txt')
filename_connectivity = str(path + 'connectivity.txt')
filename_centroids = str(path + 'centroids.txt')

filename_p_out = str(path + 'p_nodes.txt')
filename_u_out = str(path + 'u_nodes.txt')
filename_v_out = str(path + 'v_nodes.txt')
filename_w_out = str(path + 'w_nodes.txt')
filename_T_out = str(path + 'T_nodes.txt')
filename_rho_out = str(path + 'rho_nodes.txt')
filename_M_out = str(path + 'M_nodes.txt')
filename_mu_out = str(path + 'mu_nodes.txt')

filename_u_ele_out = str(path + 'u_elements.txt')
filename_v_ele_out = str(path + 'v_elements.txt')
filename_w_ele_out = str(path + 'w_elements.txt')
filename_T_ele_out = str(path + 'T_elements.txt')
filename_rho_ele_out = str(path + 'rho_elements.txt')
filename_mu_ele_out = str(path + 'mu_elements.txt')
filename_M_ele_out = str(path + 'M_elements.txt')

filename_stag_point = str(path + 'stg_point.txt')
    
filename_intersections = str(path + 'BTR_intersections.txt')
filename_nodepaths = str(path + 'BTR_paths.txt')
filename_nodecoords = str(path + 'BTR_nodes_coord.txt') 
filename_nodesresolved = str(path + 'BTR_nodes_resolved.txt') 
filename_epsilonnodes = str(path + 'BTR_epsilon_nodes.txt')

filename_dotq_out = str(path + 'dotq_nodes_out.txt')
filename_theta_out = str(path + 'theta_nodes_out.txt')
filename_skincf_out = str(path + 'skincf_nodes_out.txt')

data_to_write_info = [pinf, Tinf, Mach, rhoinf, muinf, R, cp, gamma, name, 
                      USE_IDS, IDWP, COMPUTE_STAGNATION_POINT_FROM_GEOMETRY, 
                      epsilon, T_w, PERZHIKAR, INTERPOLATE_N, beta_tol, 
                      OWN_TRANSITION, max_x_user, dist_yz_min]
WRITTEN = WriteRunInfo(filename_infofile, data_to_write_info)



if INVSCID:
    print 'Processing mesh...'

    print 'Mesh Processor start time:', time.time()
    nodes, connectivity, normals, centroids, areas = process_mesh(name, CHECK_MESH)
    print 'Mesh Processor stop time:', time.time()
    
    print 'Inviscid Solver start time:', time.time()

    x_coords = []
    y_coords = []
    z_coords = []
    
    for node in nodes:
        x_coords.append(node[0])
        y_coords.append(node[1])
        z_coords.append(node[2])
        
        
    x_coords_ele = []
    y_coords_ele = []
    z_coords_ele = []
    
    for centroid in centroids:
        x_coords_ele.append(centroid[0])
        y_coords_ele.append(centroid[1])
        z_coords_ele.append(centroid[2])


    print 'Computing pressure coefficients... '
    Cps, shadowed = Newton(nodes, connectivity, normals, centroids, areas, gamma, V_vector, Mach, VERIFY_NEWTON)
    
    print 'Computing velocities and pressures at the mesh elements...'
    V, p = GetVelocities(V_vector, normals, centroids, Cps, pinf, Mach, Tinf, gamma, R, VERIFY_NEWTON)
    
    
    u_velocities = []
    v_velocities = []
    w_velocities = []
    
    for v in V:
        u_velocities.append(v[0])
        v_velocities.append(v[1])
        w_velocities.append(v[2])
    
    
    if USE_IDS:
        print 'Computing U velocities at the vertices by IDW interpolation with p of' , IDWP, '...'
        u_nodes = GetVertexValues(nodes, centroids, u_velocities, IDWP, VERIFY_ISD_NODES)
        print 'Computing V velocities at the vertices...'
        v_nodes = GetVertexValues(nodes, centroids, v_velocities, IDWP, VERIFY_ISD_NODES)
        print 'Computing W velocities at the vertices...'
        w_nodes = GetVertexValues(nodes, centroids, w_velocities, IDWP, VERIFY_ISD_NODES)
        print 'Computing pressures at the vertices...'
        p_nodes = GetVertexValues(nodes, centroids, p, IDWP, VERIFY_ISD_NODES)
    
    else:
        print 'Computing node velocities and pressures by cell averaging... '
        u_nodes, v_nodes, w_nodes, p_nodes = GetVertexValuesFaster(nodes, centroids, connectivity, u_velocities, v_velocities, w_velocities, p)
    
    
    print 'Computing therodynamics at nodes...'
    T_nodes, mu_nodes, rho_nodes, M_nodes = GetThermodynamics(p_nodes, Mach, nodes, gamma, R, muinf, Tinf, rhoinf, VERIFY_THERMODYNAMICS)
    
    print 'Computing therodynamics at elements...'
    T_elements, mu_elements, rho_elements, M_elements = GetThermodynamics(p, Mach, centroids, gamma, R, muinf, Tinf, rhoinf, VERIFY_THERMODYNAMICS)
    
    
    p_elements = p
    u_elements = u_velocities
    v_elements = v_velocities
    w_elements = w_velocities
    
    print 'Searching for stagnation point...'
    if COMPUTE_STAGNATION_POINT_FROM_GEOMETRY:
        stag_points, stag_ps, stag_idxs = FindStagnationPointFromGeometry(nodes, p_nodes)  
    else:    
        stag_points, stag_ps, stag_idxs = FindStagnationPointFromPressure(nodes, p_nodes)                 #CURRENTLY ONLY WORKS FOR ONE STAGNATION POINT, MUST BE REDONE IN CASE OF A WINGED GEOMETRY
    
    if VERIFY_STAG_POINT:
        fig = plt.figure()
        ax = Axes3D(fig)
        
    for stag_point in stag_points:
        if VERIFY_STAG_POINT:
            ax.scatter(stag_point[0], stag_point[1], stag_point[2], s=10, c = 'r')
            
        p_nodes, T_nodes, mu_nodes, rho_nodes, u_nodes, v_nodes, w_nodes = AssgnStagnationPoint(stag_point, nodes, Mach, pinf, rhoinf, Tinf, p_nodes, rho_nodes, T_nodes, mu_nodes, u_nodes, v_nodes, w_nodes, VERIFY_STAG_POINT)
    
    
    print 'Stagnation point found at,', stag_points, ' with stagnation pressure of', stag_ps
    #print 'Backtracing the streamlines...'
    
    #intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs = TraceStreamlinesBack(nodes, u_nodes, v_nodes, w_nodes, dt, stag_point, stag_idx, epsilon, u_velocities, v_velocities, w_velocities, centroids, normals, VERIFY_BACKTRACING)
    
    print 'Inviscid Solver stop time:', time.time()

    print 'Creating surface fits... '

if WRITE_FIELDS:
    
    WRITTEN = WriteDataForSurfaceFit(str(path + 'velocities.txt'), str(path + 'coordinates.txt'), connectivity, centroids, u_nodes, v_nodes, w_nodes, nodes)
    
    WRITTEN = WriteDataNodes(filename_p_out, p_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_u_out, u_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_v_out, v_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_w_out, w_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_T_out, T_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_rho_out, rho_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_M_out, M_nodes, nodes)
    WRITTEN = WriteDataNodes(filename_mu_out, mu_nodes, nodes)
    
    WRITTEN = WriteDataNodes(filename_u_ele_out, u_elements, centroids)
    WRITTEN = WriteDataNodes(filename_v_ele_out, v_elements, centroids)
    WRITTEN = WriteDataNodes(filename_w_ele_out, w_elements, centroids)
    WRITTEN = WriteDataNodes(filename_T_ele_out, T_elements, centroids)
    WRITTEN = WriteDataNodes(filename_rho_ele_out, rho_elements, centroids)
    WRITTEN = WriteDataNodes(filename_mu_ele_out, mu_elements, centroids)
    WRITTEN = WriteDataNodes(filename_M_ele_out, M_elements, centroids)
    
    WRITTEN = WriteStagPointData(filename_stag_point, stag_points, stag_idxs, epsilon)
    WRITTEN = WriteNormals(filename_normals, normals)
    WRITTEN = WriteDataCentroids(filename_centroids, filename_connectivity, centroids, connectivity)



if VISCOUS:
    print 'Viscous Solver backtrace start time:', time.time()

    connectivity = ReadConnectivity(filename_connectivity)
    centroids = ReadCentroids(filename_centroids)
    normals =  ReadNormals(filename_normals)
    
    p_nodes, nodes =   ReadInviscidData(filename_p_out)
    u_nodes, nodes =   ReadInviscidData(filename_u_out)
    v_nodes, nodes =   ReadInviscidData(filename_v_out)
    w_nodes, nodes =   ReadInviscidData(filename_w_out)
    T_nodes, nodes =   ReadInviscidData(filename_T_out)
    rho_nodes, nodes = ReadInviscidData(filename_rho_out)
    mu_nodes, nodes =  ReadInviscidData(filename_mu_out)
    M_nodes, nodes =   ReadInviscidData(filename_M_out)

    u_velocities, centroids = ReadInviscidData(filename_u_ele_out)
    v_velocities, centroids = ReadInviscidData(filename_v_ele_out)
    w_velocities, centroids = ReadInviscidData(filename_w_ele_out)
    T_elements, centroids =   ReadInviscidData(filename_T_ele_out)
    mu_elements, centroids =  ReadInviscidData(filename_mu_ele_out) 
    rho_elements, centroids = ReadInviscidData(filename_rho_ele_out)
    M_elements, centroids =   ReadInviscidData(filename_M_ele_out)
    u_elements = u_velocities
    v_elements = v_velocities
    w_elements = w_velocities
    
    stag_points, stg_idxs, epss = ReadStagPointData(filename_stag_point)
    if len(stag_points) == 1:
        stag_point = stag_points[0]
        stag_idx = stg_idxs[0]
        epsilon = epss[0]
    if len(stag_points) > 1:
        stag_point_x = 0.
        stag_point_y = 0.
        stag_point_z = 0.
        for i in range(0, len(stag_points)):
            stag_point_x += stag_points[i][0]
            stag_point_y += stag_points[i][1]
            stag_point_z += stag_points[i][2]
        stag_point = [stag_point_x/len(stag_points), stag_point_y/len(stag_points), stag_point_z/len(stag_points)] 
        stag_idx = stg_idxs[0]
        epsilon = epss[0]
        
    
    if BACKTRACE:

        nodes_relevance = DetectIrrelevantNodes(nodes, p_nodes, pinf, Mach, gamma, VERIFY_NODE_RELEVANCE)
        print "Finding intersections with epsilon line..."
        #epsilon = 2. * epsilon
        xmin, xmax, ymin, ymax, zmin, zmax = FindBounds(nodes)
        
        centroids_ordered, centroids_ordered_dxs = SeparateCentroidsByDx(centroids, xmin, xmax, 2)
        PRINT_BTR = True
        intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs = TraceStreamlinesBack(max_x_user, nodes, nodes_relevance, u_nodes, v_nodes, w_nodes, dt, stag_point, epsilon, u_velocities, v_velocities, w_velocities, centroids, centroids_ordered, centroids_ordered_dxs, normals, connectivity, xmin, xmax, ymin, ymax, zmin, zmax, VERIFY_BACKTRACING, PRINT_BTR)
        print "Writing intersection data..."
        
        print 'Viscous Solver backtrace stop time:', time.time()

    if WRITE_BACKTRACING:
        WRITTEN = WriteBackTracing(filename_intersections, filename_nodepaths, filename_nodecoords, filename_nodesresolved, filename_epsilonnodes, intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs)

    else:
        print "Reading intersection data..."
        intersections, node_paths_elem, node_paths_coord, epsilon_nodes, nodes_resolved_idxs = ReadBacktracingData(filename_intersections, filename_nodepaths, filename_nodecoords, filename_nodesresolved, filename_epsilonnodes)
    
    
    if not PERZHIKAR:
        filename_null = path + "nullspace.txt"
        #filename_p = path + "p_nodes_fitted.txt"
        filename_u = path + "u_nodes_fitted.txt"
        filename_v = path + "v_nodes_fitted.txt"
        filename_w = path + "w_nodes_fitted.txt"
        #filename_T = path + "T_nodes_fitted.txt"
        #filename_rho = path + "rho_nodes_fitted.txt"
        #filename_M = path + "M_nodes_fitted.txt"
        #filename_mu = path + "mu_nodes_fitted.txt"
        
        
        null_spaces = ReadNullSpace(filename_null)
        #p_fitted = ReadSolvedFitting(filename_p)
        u_fitted = ReadSolvedFitting(filename_u)
        v_fitted = ReadSolvedFitting(filename_v)
        w_fitted = ReadSolvedFitting(filename_w)
        #T_fitted = ReadSolvedFitting(filename_T)
        #rho_fitted = ReadSolvedFitting(filename_rho)
        #M_fitted = ReadSolvedFitting(filename_M)
        #mu_fitted = ReadSolvedFitting(filename_mu)
        
    else:
        u_fitted = [] 
        v_fitted = []
        w_fitted = []
        for i in range(0, len(nodes)):
            u_fitted.append([0.,0.,0.,0.,0.,0.])
            v_fitted.append([0.,0.,0.,0.,0.,0.])
            w_fitted.append([0.,0.,0.,0.,0.,0.])
        
        null_spaces = []
        for i in range(0, len(centroids)):
            null_spaces.append([0.,0.,0.,0.,0.,0.])
        
    wd1s = []
    wd2s = []
    wd3s = []
    wd4s = []
    wd5s = []
    wd6s = []
    for i in range(0, len(w_fitted)):
        wd1s.append(w_fitted[i][0])
        wd2s.append(w_fitted[i][1])
        wd3s.append(w_fitted[i][2])
        wd4s.append(w_fitted[i][3])
        wd5s.append(w_fitted[i][4])
        wd6s.append(w_fitted[i][5])
        
    wds = [wd1s, wd2s, wd3s, wd4s, wd5s, wd6s]
        
    Fxs = []
    GradFs = []
    
    if not PERZHIKAR:
        for c in range(0, len(centroids)):
            
            a_vector = null_spaces[c]
            centroid = centroids[c]
            x = centroid[0]
            y = centroid[1]
            z = centroid[2]
            Fx = 2. * a_vector[4] * x + a_vector[5]
            Fy = 2. * a_vector[0] * y + a_vector[1]
            Fz = 2. * a_vector[2] * z + a_vector[3]
            
            GradF = (Fx**2 + Fy**2 + Fz**2)**0.5
            
            '''
            if GradF > 4.:
                GradF = 100000
            if Fx > 0.5:
                Fx = 1000.
            if Fx < -0.5:
                Fx = -1000.
            '''
            Fxs.append(Fx)
        
            GradFs.append(GradF)
        
        if VERIFY_A_COMPUTATION:
            errors = []
            xs_err = []
            ys_err = []
            zs_err = []
            
            for c in range(0, len(centroids)):
                a_vector = null_spaces[c]
                x = centroids[c][0]
                y = centroids[c][1]
                z = centroids[c][2]
                
                u = u_elements[c] 
                v = v_elements[c] 
                w = w_elements[c] 
                centroid = centroids[c]
                Fx = 2. * a_vector[4] * centroid[0] + a_vector[5]
                Fy = 2. * a_vector[0] * centroid[1] + a_vector[1]
                Fz = 2. * a_vector[2] * centroid[2] + a_vector[3]
                if abs(Fx) > 0.:
                    u_th = (v * Fy + w * Fz )/  (-Fx)
                else:
                    u_th = 0.
                #print "Error in a vector: ", (u - u_th)/u_elements[c]]
                if u_elements[c] > 0.:
                    err = (u - u_th)/u_elements[c]
                else:
                    err = 0.
                if -1. < err < 1.:
                    xs_err.append(x)
                    ys_err.append(y)
                    zs_err.append(z)
                    if u_elements[c] > 0.:
                        errors.append((u - u_th)/u_elements[c])
                    else:
                        errors.append(0.)
            fig = plt.figure()
            ax = Axes3D(fig)
            c = np.array(errors)
            cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
            p = ax.scatter(xs_err, ys_err, zs_err, s=15, c = c)
            fig.colorbar(p, cax = cax, orientation = 'vertical')
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            ax.set_title('Errors in the matrix inversion')
    
            
        x_coords_elefourth = []
        y_coords_elefourth = []
        z_coords_elefourth = []
        GradFsfourth = []
        Fxsfourth = []
        for i in range(0, len(centroids)):
            x = centroids[i][0]
            y = centroids[i][1]
            z = centroids[i][2]
            
            #if x > -0.5:
            #    if y > 0.0:
            #        if z > 0.0:
            x_coords_elefourth.append(x)
            y_coords_elefourth.append(y)
            z_coords_elefourth.append(z)
            
            if GradFs[i] > 30:
                GradFs[i] = 30.
            if GradFs[i] < -30:
                GradFs[i] = -30.
            if Fxs[i] > 1:
                Fxs[i] = 1.
            if Fxs[i] < -1:
                Fxs[i] = -1.
            
            GradFsfourth.append(GradFs[i])
            Fxsfourth.append(Fxs[i])
        
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(Fxsfourth)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_elefourth, y_coords_elefourth, z_coords_elefourth, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Fxs on a cone using modified Newton technique')
        
        fig = plt.figure()
        ax = Axes3D(fig)
        c = np.array(GradFsfourth)
        cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
        p = ax.scatter(x_coords_elefourth, y_coords_elefourth, z_coords_elefourth, s=5, c = c)
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('GradFs on a cone using modified Newton technique')
    
    
    else:
        for i in range(0, len(centroids)):
            GradFs.append(0.0)
            Fxs.append(0.0)
            
        
    
    '''
    y_coords_stg = []
    z_coords_stg = []
    x_coords_stg = []
    for i in range(0, len(rn_s)):
        x_coords_stg.append(math.cos(betas[i]) * rn_s[i])
        y_coords_stg.append(math.sin(betas[i]) * rn_s[i])
        z_coords_stg.append(stag_point[2] - epsilon)
        
    fig = plt.figure()
    ax = Axes3D(fig)
    c = np.array(GradFs)
    cax = fig.add_axes([0.875, 0.1, 0.05, 0.8])
    p = ax.scatter(x_coords_stg, y_coords_stg, z_coords_stg, s=5, c = 'k')
    #fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Stagnation points resolved by backtracing')
    '''
    
    if not PERZHIKAR:    
        print 'Performing epsilon projecton and calculating matric coefficients at the stagnation region...'
    
        rn_s, betas = EpsilonCircTransformation(epsilon_nodes)
        hs, ss = FindInitialMetricCoefficientsCircle(betas, epsilon, epsilon_nodes, beta_tol, VERIFY_HS_CLUTTERING)
        wdydbetazs = FindInitial_wdydbetazStable(epsilon_nodes, nodes, T_w, T_nodes, mu_nodes, rho_nodes, M_nodes, centroids, u_elements, v_elements, w_elements , hs, ss, cp, k, R, gamma, stag_idx)
    
    else:
        hs = [] 
        ss = []
        rn_s = []
        betas = []
        wdydbetazs = []
        for i in range(0, len(epsilon_nodes)):
            hs.append(epsilon)
            ss.append(0.0)
            wdydbetazs.append(0.0)
            rn_s.append(0.0)
            betas.append(0.0)

            
    '''
    for node in intersections:
        x_coords_node.append(node[0])
        y_coords_node.append(node[1])
        z_coords_node.append(node[2])
        
    scatter_1 = ax.scatter(x_coords_node, y_coords_node, z_coords_node, s=2, c = 'k')
    '''
    #USE ISENTROPIC RELATIONS POSSIBLY TO DETERMINE THE TEMPERATURES AND DENSITIES AT THE VERTICES BEFORE YOU CONTINUE
    '''
    for i in range(0, len(hs)):
        hs[i] = abs(hs[i])
    
    wdydbetazs = FindInitial_wdydbetaz(hs, GradFs, Fxs, u_elements, v_elements, w_elements, epsilon_nodes, node_paths_elem)
    '''
    
    final_xs = []
    final_ys = []
    final_zs = []
    final_values_thetas = []
    final_values_dotqws = []
    final_values_hs = []
    final_values_Ves = []
    final_values_rhoes = []
    final_values_rhostars = []
    final_values_mustars = []
    final_values_integs = []
    final_values_skincfs = []
    
    print "Computing SP heat flux..."
    Tstag_behindNSW, rhostag_behindNSW, pstag_behindNSW, mustag_behindNSW, Vstag_behindNSW = GetStagPointConditions(Mach, pinf, rhoinf, Tinf) 
    stg_cnds = [Tstag_behindNSW, rhostag_behindNSW, pstag_behindNSW, mustag_behindNSW, Vstag_behindNSW]
    qstag = GetStagnationHeating(stg_cnds, T_w, 1.0, cp, pinf, epsilon)
        
    
    if PERZHIKAR and MARCH:
        print 'Viscous Solver Parzhikar start time:', time.time()

        print "Pre-computing Perzhikar's metric coefficients..."
        metriccoeffs = ComputeMetricPerzhikar(nodes, nodes_resolved_idxs, node_paths_coord, stag_point, dist_yz_min)

        print 'Viscous Solver Parzhikar stop time:', time.time()

    else:
        metriccoeffs = []
        for i in range(0, len(epsilon_nodes)):
            metriccoeffs.append(0.0)
            
            
    if MARCH:
        print 'Marching forward and computing solution...'
        
        print 'Viscous Solver marchdown start time:', time.time()

        for s in range(0, len(epsilon_nodes)):
            metric_coeffs_current = copy.deepcopy(list(reversed(metriccoeffs[s])))
            h_ini = hs[s]
            wdydbetaz_ini = wdydbetazs[s]
            list_of_crossed_elements = copy.deepcopy(list(reversed(node_paths_elem[s])))
            list_of_coordinates = copy.deepcopy(list(reversed(node_paths_coord[s])))
            
            hs_res, hs_coords, thetas, dotq_ws, deltas, transition_points, V_es, rho_es, mu_stars, rho_stars, integs, dss, skin_cfs = MarchDownstream(stag_point, stg_cnds, wdydbetaz_ini, list_of_crossed_elements, list_of_coordinates, h_ini, GradFs, Fxs, u_elements, v_elements, w_elements, T_elements, mu_elements, rho_elements, M_elements, w_fitted, T_w, cp, k, wds, metric_coeffs_current, INTERPOLATE_N, OWN_TRANSITION, PERZHIKAR)
            
            if len(hs_coords) > 0.:
                y_hs = []
                x_hs = []
                z_hs = []
                
                if hs_coords[-1][0] > -1.:
                    for i in range(0, len(hs_coords)):
                        x_hs.append(hs_coords[i][0])
                        y_hs.append(hs_coords[i][1])
                        z_hs.append(hs_coords[i][2])
                    
                    #c = np.array(dss)  
                    
                    
                    final_xs.append(x_hs[-1])
                    final_ys.append(y_hs[-1])
                    final_zs.append(z_hs[-1])
                    final_values_thetas.append(thetas[-1])
                    final_values_dotqws.append(dotq_ws[-1])
                    final_values_hs.append(hs_res[-1])
                    final_values_Ves.append(V_es[-1])
                    final_values_rhoes.append(rho_es[-1])
                    final_values_rhostars.append(rho_stars[-1])
                    final_values_mustars.append(mu_stars[-1])
                    final_values_integs.append(integs[-1])
                    final_values_skincfs.append(skin_cfs[-1])
                    
                    final_xs.append(x_hs[0])
                    final_ys.append(y_hs[0])
                    final_zs.append(z_hs[0])
                    final_values_thetas.append(thetas[0])
                    final_values_dotqws.append(dotq_ws[0])
                    final_values_hs.append(hs_res[0])
                    final_values_Ves.append(V_es[0])
                    final_values_rhoes.append(rho_es[0])
                    final_values_rhostars.append(rho_stars[0])
                    final_values_mustars.append(mu_stars[0])
                    final_values_integs.append(integs[0])
                    final_values_skincfs.append(skin_cfs[0])
                    
                    '''
                    
                    final_xs.append(x_hs[-1])
                    final_ys.append(y_hs[-1])
                    final_zs.append(z_hs[-1])
                    final_values_thetas.append(thetas[-1])
                    final_values_dotqws.append(0.)
                    final_values_hs.append(hs_res[-1])
                    final_values_Ves.append(V_es[-1])
                    final_values_rhoes.append(rho_es[-1])
                    final_values_rhostars.append(rho_stars[-1])
                    final_values_mustars.append(mu_stars[-1])
                    final_values_integs.append(integs[-1])
                    
                    final_xs.append(x_hs[0])
                    final_ys.append(y_hs[0])
                    final_zs.append(z_hs[0])
                    final_values_thetas.append(thetas[0])
                    final_values_dotqws.append(0.)
                    final_values_hs.append(hs_res[0])
                    final_values_Ves.append(V_es[0])
                    final_values_rhoes.append(rho_es[0])
                    final_values_rhostars.append(rho_stars[0])
                    final_values_mustars.append(mu_stars[0])
                    final_values_integs.append(integs[-1])
                    '''
                    #p = ax.scatter(x_hs, y_hs, z_hs, s=5, c = c)
        '''
        fig.colorbar(p, cax = cax, orientation = 'vertical')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Ds along a cone')
        '''
        
    print 'Viscous Solver marchdown stop time:', time.time()

    final_values_dotqws.append(qstag/10000.)
    final_xs_q = copy.deepcopy(final_xs)
    final_ys_q = copy.deepcopy(final_ys)
    final_zs_q = copy.deepcopy(final_zs)
    final_xs_q.append(stag_point[0])
    final_ys_q.append(stag_point[1])
    final_zs_q.append(stag_point[2])
    
    fig = plt.figure()
    ax = Axes3D(fig)
    cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])
    
    #for theta in range(0, len(final_values_thetas)):
    #    if final_values_thetas[theta] > 0.005:
    #        final_values_thetas[theta] = 0.00
    c = np.array(final_values_thetas)  
    p = ax.scatter(final_xs, final_ys, final_zs, s=5, c = c)
    fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Thetas along a cone, final values only, m')
    
    
    fig = plt.figure()
    ax = Axes3D(fig)
    cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])
    
    for k in range(0, len(final_xs_q)-1):
        dist_to_stag_point = ((final_xs_q[k] - final_xs_q[-1])**2 + (final_ys_q[k] - final_ys_q[-1])**2 + (final_zs_q[k] - final_zs_q[-1])**2)**0.5
        if final_values_dotqws[k] > qstag:
            final_values_dotqws[k] = qstag
        if dist_to_stag_point < 1. * epsilon:
            final_values_dotqws[k] = qstag / 10000.
        else:
            final_values_dotqws[k] /= 10000.
    
    c = np.array(final_values_dotqws)  
    p = ax.scatter(final_xs_q, final_ys_q, final_zs_q, s=5, c = c)
    ax.scatter(final_xs_q[-1], final_ys_q[-1], final_zs_q[-1], s=100, c = 'r')
    fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Dotqws along a cone, final values only, W/cm2')
    
    fig = plt.figure()
    ax = Axes3D(fig)
    cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])
    
    c = np.array(final_values_mustars)  
    p = ax.scatter(final_xs, final_ys, final_zs, s=5, c = c)
    fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Mu star along a cone, final values only')
    
    fig = plt.figure()
    ax = Axes3D(fig)
    cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])
    
    c = np.array(final_values_integs)  
    p = ax.scatter(final_xs, final_ys, final_zs, s=5, c = c)
    fig.colorbar(p, cax = cax, orientation = 'vertical')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Theta integral term along a cone, final values only')
    
    
    nodes_final_dotq = []
    for i in range(0, len(final_xs_q)):
        nodes_final_dotq.append([final_xs_q[i], final_ys_q[i], final_zs_q[i]])
        
    WRITTEN_Q = WriteDataNodes(filename_dotq_out, final_values_dotqws, nodes_final_dotq)
    
    nodes_final_theta = []
    for i in range(0, len(final_xs)):
        nodes_final_theta.append([final_xs[i], final_ys[i], final_zs[i]])
        
    WRITTEN_THETA = WriteDataNodes(filename_theta_out, final_values_thetas, nodes_final_theta)
    
    
    nodes_final_skincf = []
    for i in range(0, len(final_xs)):
        nodes_final_skincf.append([final_xs[i], final_ys[i], final_zs[i]])
        
    WRITTEN_SKINCF = WriteDataNodes(filename_skincf_out, final_values_skincfs, nodes_final_skincf)
    
