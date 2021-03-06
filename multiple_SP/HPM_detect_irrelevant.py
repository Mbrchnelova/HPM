# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_detect_irrelevant.py: this function distinguishes which nodes should not be computed (here e.g. based on Cp).
# This can be adjusted based on user's wish (e.g. geometry based constraints).
# Correct detection allows acceleration of backtracing and marching down. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be



from HPM_import import *



def DetectIrrelevantNodes(nodes, p_nodes, pinf, Minf, gamma, VERIFY):
    relevant_nodes = []
    for n in range(0, len(nodes)):
        p = p_nodes[n]
        Cp = (p - pinf) / (0.5 * gamma * pinf * Minf**2)
        if Cp > 0.:
            relevant_nodes.append(1)
            
        else:
            relevant_nodes.append(0)
            
    if VERIFY:
        '''
        fig = plt.figure()
        ax = Axes3D(fig)
        for n in range(0, len(nodes)):
            node = nodes[n]
            if relevant_nodes[n] == 0:                
                ax.scatter(node[0], node[1], node[2], s=10, c = 'r')
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Irrelevant nodes detected')
        '''
    return relevant_nodes    


