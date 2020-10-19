# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_vertex_values.py: this module contains function to extrapolate element data onto the nodes.
# Two methods are present; the first one uses inverse-distance weighing and while being accurate, it is extremely slow.
# The other method only uses the neighbouring elements for extrapolation. 
# The use of the second method is recommended as it is order of magnitude faster (no nested for loops). 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *
from HPM_mathematics import *





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



###############################################################################




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


