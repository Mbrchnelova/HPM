# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_transform_mesh_scale.py: this module contains a function to transform the imported mesh using a constant scaling factor.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



import math




def TransformByScaling(centroids, nodes, normals, scale):
	n_c = len(centroids)
	n_n = len(nodes)
	n_r = len(normals)

	#The transformation will be in x and z only.

	for i in range(0, n_c):
		centroid = centroids[i]
		c_x = centroid[0]
                c_y = centroid[1]
                c_z = centroid[2]

		centroids[i] = [scale * c_x, scale * c_y, scale * c_z]


        for i in range(0, n_n):
                node = nodes[i]
                n_x = node[0]
                n_y = node[1]
                n_z = node[2]

                nodes[i] = [scale * n_x, scale * n_y, scale * n_z]


	#There is no normal transformation necessary for pure scaling
	#for i in range(0, n_r):
	#	normal = normals[i]
        #        n_x = normal[0]
        #        n_y = normal[1]
        #        n_z = normal[2]

        #        n_x_new = math.cos(alpha) * n_x - math.sin(alpha) * n_z
        #        n_z_new = math.sin(alpha) * n_x + math.cos(alpha) * n_z

        #        normals[i] = [n_x_new, n_y, n_z_new]

	return centroids, nodes, normals
