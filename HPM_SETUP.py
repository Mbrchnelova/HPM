import subprocess
import math


#########INPUT FILE MANAGEMENT
name = "/Users/brch/HPM_geometries/geom2.ply"					#The name of the mesh file to be read
path = '/Users/brch/HPM_results/geom_2/'					#The path to the results folder (does not need to exist a priori)





#########INVISCID PART
INVSCID                                         = 1				#If selected, the Newton calculation and mesh processing will take place
CHECK_MESH                                      = 1				#If selected, the loaded mesh will be displayed (without transformations)
VERIFY_NEWTON                                   = 1				#If selected, the Newton results will be displayer (after transformations)

USE_IDS                                         = 0				#Use inverse distance weighing. More precise but takes a very long time.
VERIFY_ISD_NODES                                = 0				#If using inverse distance weighing, the process will be displayed for verification

VERIFY_THERMODYNAMICS                           = 1				#Display the computed inviscid thermodynamics for verification
VERIFY_STAG_POINT                               = 1				#Display the position of the located stagnation point(s)
WRITE_FIELDS                                    = 1				#Save the information from the inviscid processing (Newton + mesh processing)


#########VISCOUS PART
VISCOUS                                         = 1				#If selected, viscous calculations will take place (inviscid data required!)

VERIFY_NODE_RELEVANCE                           = 0				#There if a function to recognise nodes in shadow to be ignored. This displays these nodes.
VERIFY_BACKTRACING                              = 1				#Display the results of backtracing - red nodes were unsuccessful and green successful
BACKTRACE                                       = 1				#If selected, perform backtracing. If not, backtracing can be read from the results folder.
WRITE_BACKTRACING                               = 1				#Save backtracing data

MARCH                                           = 1				#If selected, perform downstream marching (inviscid & backtracing data required!)
VERIFY_HS_CLUTTERING                            = 0				#This is for stagnation line resultion for the Hamilton method, to verify line projection onto 2D
VERIFY_A_COMPUTATION                            = 0				#This option currently does not give any useful info. Orinally set up to verify the Hamilton method.

##########GENERAL SETTING
COMPUTE_STAGNATION_POINT_FROM_GEOMETRY          = 1				#If selected, the front-most point is used as SP. SP not determined from pressures.
PERZHIKAR                                       = 1				#If selected, Parzhikar's method is used instead of Hamilton. This is much preferred.



#########USER INPUT MANAGEMENT

#Numerical
dt = 0.00001									#Right now, this does not matter - dt is calculate on the fly. This can be changed in the backtracing module.
epsilon = 0.05									#The radius of the epsilon line. Adjust based on the geometry.
IDWP = 4.									#If IDW is used to compute node data, IDWP indicates the exponent in IDW.

#Aerodynamical
gamma = 1.4									#The ratio of specific heats
R = 287.05									#The gas constant
pinf = 580.29									#The static freestream pressure
Tinf = 62.870									#The static freestream temperature
muinf = 4.196e-6								#The static freestream viscosity
rhoinf = pinf / (Tinf * R)							#The ideal gas law is used to compute the freestream density
Mach = 6.									#The Mach number

N = 100. 									#This should not be doing anything but I am too afraid to remove it
T_w = 300.									#The vehicle wall temperature (needed to compute heat transfer)
cp = 1006.									#The heat capacity of air
k = 0.02435									#The heat conductivity of air
V_vector = [1., 0., 0.]                                                    	#This was originaly intended to manipulate alpha and beta, but using the mesh transformations below is safer 


#Mesh transformation
alpha = 0.0 #-30. * math.pi/180.						#The angle of attack to transform the mesh. For positive AOA, this angle should be negative.
beta = 0. * math.pi/180								#The angle of sideslip to transform the mesh
scale = 0.5									#The scaling to transform the mesh

#Transition settings
INTERPOLATE_N = False                                                    	#For turbulent heating, either compute N (more accurate, less stable) = False or interpolate N = True
OWN_TRANSITION_Re = 'nan'                                               	#User defined transition Re number (x-based)
OWN_TRANSITION_x = 0.12                                               		#User defined transition coordinate (x-based)
OWN_TRANSITION = [OWN_TRANSITION_Re, OWN_TRANSITION_x]


#User/ pre-set numericallimits
max_x_user = 4.0                                                                #Distance beyond which the results are not desired (i.e too downstream to be relevant) for speedup
beta_tol = 0.01                                                                 #Resolution for point cluttering on the epsilon line, for Hamilton
dist_yz_min = 0.07                                                              #Min distance between two streamlines, for Perzhikar




#########THE FOLLOWING MIGHT NOT WORK ON WINDOWS --> REPLACE IF NOW UNIX/OS
bashCmd = "mkdir " + str(path)
print bashCmd
process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)

output, error = process.communicate()

print output





V_vector = [1., 0., 0.]
V_mag = (V_vector[0]**2 + V_vector[1]**2 + V_vector[2]**2)**0.5






