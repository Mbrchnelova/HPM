import sys
import math


geom = sys.argv[1]

no_h = int(sys.argv[2])

nose_end = float(sys.argv[3])

cone_1_end = float(sys.argv[4])


max_nose_heat_flux = []
max_cone_1_heat_flux = []
max_cone_2_heat_flux = []
height = []

temp_heat_flux_nose = []
temp_heat_flux_cone_1 = []
temp_heat_flux_cone_2 = []

path_def = "/scratch/leuven/338/vsc33811/HPM/HPM_results/geom_" + geom + "_2_"
path_erwin = "/scratch/leuven/338/vsc33811/HPM/Output/" + geom + "/out02.gnu"
for i in range(0, no_h):
	path = path_def + str(i) + "/" + "dotq_nodes_out.txt"
	f = open(path, "r")

	lines = f.readlines()
	temp_heat_flux_cone_2 = []
	temp_heat_flux_cone_1 = []
	temp_heat_flux_nose = []

	for line in lines:
		line_formatted = line.split(" ")
		x = float(line_formatted[0])
		y = float(line_formatted[1])
		z = float(line_formatted[2])
		q = float(line_formatted[3])

		if x < nose_end:
			temp_heat_flux_nose.append(q)
		elif x < cone_1_end:
			temp_heat_flux_cone_1.append(q)
		else:
			temp_heat_flux_cone_2.append(q)


	max_nose = max(temp_heat_flux_nose)
	max_cone_1 = max(temp_heat_flux_cone_1)
	max_cone_2 = max(temp_heat_flux_cone_2)


	max_nose_heat_flux.append(max_nose)
	max_cone_1_heat_flux.append(max_cone_1)
	max_cone_2_heat_flux.append(max_cone_2)
	
	height.append(i)

f.close()

print "Heatflux in kW/m2"

g = open(path_erwin, "r")
lines = g.readlines() 


for i in range(0, no_h):
	erwin_line = lines[i]
	erwin_line = erwin_line.split(" ")
	erwin_line_formatted = []
	for j in range(0, len(erwin_line)):
		if len(erwin_line[j]) > 0:
			erwin_line_formatted.append(erwin_line[j])


	#print erwin_line
	erwin_qnose = float(erwin_line_formatted[10])
	erwin_qcone_1 = float(erwin_line_formatted[17])
	erwin_qcone_2 = float(erwin_line_formatted[18])
	print height[i], '\t', 10*max_nose_heat_flux[i], '\t', erwin_qnose, '\t', 10*max_cone_1_heat_flux[i], '\t', erwin_qcone_1, '\t', 10*max_cone_2_heat_flux[i], '\t', erwin_qcone_2







