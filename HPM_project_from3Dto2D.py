# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_project_from3Dto2D: this is a basic module used anytime a 3D mesh element is to be converted to a plane in 2D (e.g. velocity calculation).
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be


from HPM_vector_mathematics import *
from HPM_matrix_mathematics import *
 

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
