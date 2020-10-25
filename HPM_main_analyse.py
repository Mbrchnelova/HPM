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
	path = "./Output/" + str(i)
	for j in range(1, 9):
		path = path + "out0" + str(j) + ".gnu"
		conditions = HPM_read_gnu_data(path)
		heights = conditions[1]
		flight_paths = conditions[5]
		headings = conditions[6]
		velocities = conditions[4]

		number_of_entries = len(heights)

		for k in range(0, number_of_entries):

			

#for i in range(1, 28):

#	infiles.append("/Users/brch/HPM_geometries/geom" + str(i) + ".ply")
#	outfile.append("/Users/brch/HPM_results/geom_" + str(i) +  "/" )

#aoas = []

#for i in range(0, 16):
#	aoas.append(math.pi/180.0 * float(i) * 5.0)
 

#	for i in range(1, 28):

#        	infiles.append("/Users/brch/HPM_geometries/geom" + str(i) + ".ply")
#       		outfile.append("/Users/brch/HPM_results/geom_" + str(i) +  "/" )


from HPM_INVISCID import *
from HPM_BACKTRACE_MARCHDOWN import *













