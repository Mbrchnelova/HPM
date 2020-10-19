# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_reading_results.py: this is a standalone code, not required to run the code.
# It helps reading the output files and with further processing and visualisation. 
# Can and should be adapted by the user depending on the way analysis should be made.
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np

    
def ReadData(filename):
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


path = '/Users/brch/Downloads/correctly_scaled_20_deg_tiny/'
filename_dotq_out = str(path + 'dotq_nodes_out.txt')
filename_theta_out = str(path + 'theta_nodes_out.txt')


final_values_dotqws, nodes_dotqs = ReadData(filename_dotq_out)
final_values_thetas, nodes_thetas = ReadData(filename_theta_out)

final_xs_q = []
final_ys_q = []
final_zs_q = []

for i in range(0, len(nodes_dotqs)):
    final_xs_q.append(nodes_dotqs[i][0])
    final_ys_q.append(nodes_dotqs[i][1])
    final_zs_q.append(nodes_dotqs[i][2])
 

final_xs = []
final_ys = []
final_zs = []

for i in range(0, len(nodes_thetas)):
    final_xs.append(nodes_thetas[i][0])
    final_ys.append(nodes_thetas[i][1])
    final_zs.append(nodes_thetas[i][2])


fig = plt.figure()
ax = Axes3D(fig)
cax = fig.add_axes([0.825, 0.1, 0.05, 0.8])

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

c = np.array(final_values_dotqws)  
p = ax.scatter(final_xs_q, final_ys_q, final_zs_q, s=5, c = c)
ax.scatter(final_xs_q[-1], final_ys_q[-1], final_zs_q[-1], s=100, c = 'r')
fig.colorbar(p, cax = cax, orientation = 'vertical')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.set_title('Dotqws along a cone, final values only, W/cm2')
