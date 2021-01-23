# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_main_analyse.py: this is the file to be executed. 
# By importing the two main modules, they are being executed; first the inviscid solver and then the viscous one.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be


import math
from HPM_read_gnu_data import *
from HPM_ISA import *
from HPM_FE_to_Fa import *

infiles = [] #"geom1.ply", "geom1.ply", "geom1.ply", ]
outflies = []
conditions = []


conditions = []


#                         (  0) time                  1
#                         (  1) height                2
#                         (  2) longitude             3
#                         (  3) latitude              4
#                         (  4) velocity              5
#                         (  5) fl.path angle         6
#                         (  6) heading               7 wrt north (kept at 90 deg, so due east)
#                         (  7) drag force           17
#                         (  8) Twall                39
#                         (  9) heat load            40 stagnation point
#                         ( 10) heat flux            41 for stagnation-point nose radius
#                         ( 11) gload                49
#                         ( 12) atm. density         93 US76 model
#                         ( 13) Mach number          96 US76 speed of sound
#                         ( 14) Re number            97 based on vehicle length
#                         ( 15) qdyn                 98 dynamic pressure


for i in range(1, 28):
	
	path_def = "./Output/" + str(i) + "/"
	for j in range(1, 9):
		path = ""
		path = path_def + "out0" + str(j) + ".gnu"
		conditions = read_gnu_data(path)

		heights = []
		flight_paths = []
		headings = []
		velocities = []
		number_of_entries = len(conditions)
		for k in range(0, number_of_entries):
			heights.append(conditions[k][1])
			flight_paths.append(conditions[k][5])
			headings.append(conditions[k][6])
			velocities.append(conditions[k][4])

		number_of_entries = len(heights)

		print "heights: ", heights
		for k in range(0, number_of_entries):

			# Below are the default conditions same for all simulations
			dt = 0.00001
			epsilon = 0.05 
			IDWP = 4.
			T_w = 300.
			cp = 1000.6
        	        kap = 0.02435
			V_vector = [1., 0., 0.] 
 
			gamma = 1.4
			R = 287.05

			scale = 1.0

			# Next, the freestream conditions are found using ISA
			height = heights[k] * 1000. # conversion to meters
			pinf, Tinf, rhoinf, muinf = ISA_from_height(height)

			ainf = math.sqrt(gamma * R * Tinf)
			Mach = velocities[k] / ainf

			# Finally, the heading and flight path angle is transformed into AOA and Beta
			V_E = [velocities[k], 0, 0]
			flight_path_angle = flight_paths[k] * math.pi/180. # conversion to rad
			heading_angle = headings[k] * math.pi/180. #conversion to rad
			bank_angle = 0.0
			V_a, alpha, beta = FE_to_Fa(V_E, heading_angle, flight_path_angle, bank_angle)			


			# Balistic flight
			alpha = alpha * 0.0
			beta = beta * 0.0

			print "Running simulation: "
			print "H: ", height/1000., "km, V: ", V_E, "m/s, Mach: ", Mach, " AOA: ", alpha * 180./math.pi, "deg, BETA: ", beta * 180./math.pi, "deg."
			infile = "/Users/brch/HPM_geometries/geom" + str(i) + ".ply"
			outfile = "/Users/brch/HPM_results/geom_" + str(i) + "_" + str(j) + "_" + str(k) +  "/"

			OVERWRITE_SETUP = True 
			from HPM_INVISCID import *
			from HPM_BACKTRACE_MARCHDOWN import *


#for i in range(1, 28):

#	infiles.append("/Users/brch/HPM_geometries/geom" + str(i) + ".ply")
#	outfile.append("/Users/brch/HPM_results/geom_" + str(i) +  "/" )

#aoas = []

#for i in range(0, 16):
#	aoas.append(math.pi/180.0 * float(i) * 5.0)
 

#	for i in range(1, 28):

#        	infiles.append("/Users/brch/HPM_geometries/geom" + str(i) + ".ply")
#       		outfile.append("/Users/brch/HPM_results/geom_" + str(i) +  "/" )


#from HPM_INVISCID import *
#from HPM_BACKTRACE_MARCHDOWN import *













