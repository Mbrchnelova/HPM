# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_transform_mesh_beta.py: this module contains a function to transform the imported mesh using an angle of sideslip.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



import math




def TransformByBeta(centroids, nodes, normals, beta):
	n_c = len(centroids)
	n_n = len(nodes)
	n_r = len(normals)

	#The transformation will be in x and z only.

	for i in range(0, n_c):
		centroid = centroids[i]
		c_x = centroid[0]
                c_y = centroid[1]
                c_z = centroid[2]

		c_x_new = math.cos(beta) * c_x - math.sin(beta) * c_y
                c_y_new = math.sin(beta) * c_x + math.cos(beta) * c_y

		centroids[i] = [c_x_new, c_y_new, c_z]


        for i in range(0, n_n):
                node = nodes[i]
                n_x = node[0]
                n_y = node[1]
                n_z = node[2]

                n_x_new = math.cos(beta) * n_x - math.sin(beta) * n_y
                n_y_new = math.sin(beta) * n_x + math.cos(beta) * n_y

                nodes[i] = [n_x_new, n_y_new, n_z]


	for i in range(0, n_r):
		normal = normals[i]
                n_x = normal[0]
                n_y = normal[1]
                n_z = normal[2]

                n_x_new = math.cos(beta) * n_x - math.sin(beta) * n_y
                n_y_new = math.sin(beta) * n_x + math.cos(beta) * n_y

                normals[i] = [n_x_new, n_y_new, n_z]

	return centroids, nodes, normals
