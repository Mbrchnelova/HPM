# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_transform_mesh_alpha.py: this module contains a function to transform the imported mesh using AOA.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



import math




def TransformByAOA(centroids, nodes, normals, alpha):
	n_c = len(centroids)
	n_n = len(nodes)
	n_r = len(normals)

	#The transformation will be in x and z only.

	for i in range(0, n_c):
		centroid = centroids[i]
		c_x = centroid[0]
                c_y = centroid[1]
                c_z = centroid[2]

		c_x_new = math.cos(alpha) * c_x - math.sin(alpha) * c_z
                c_z_new = math.sin(alpha) * c_x + math.cos(alpha) * c_z

		centroids[i] = [c_x_new, c_y, c_z_new]


        for i in range(0, n_n):
                node = nodes[i]
                n_x = node[0]
                n_y = node[1]
                n_z = node[2]

                n_x_new = math.cos(alpha) * n_x - math.sin(alpha) * n_z
                n_z_new = math.sin(alpha) * n_x + math.cos(alpha) * n_z

                nodes[i] = [n_x_new, n_y, n_z_new]


	for i in range(0, n_r):
		normal = normals[i]
                n_x = normal[0]
                n_y = normal[1]
                n_z = normal[2]

                n_x_new = math.cos(alpha) * n_x - math.sin(alpha) * n_z
                n_z_new = math.sin(alpha) * n_x + math.cos(alpha) * n_z

                normals[i] = [n_x_new, n_y, n_z_new]

	return centroids, nodes, normals
