# This code if a part of the HPM library for rapid hypersonic modelling.
# Heads up! This software most likely still contains errors.
# It is therefore distributed without warranty of merchantability.
#
#
# HPM_filter_nodes: this function determined whether an element is completely shadowed or not, and thus excluded from calculation. 
#
# Developed/ made available:   19/10/2020 by M. Brchnelova
# Questions?                   michaela.brchnelova@kuleuven.be


from HPM_import import *




   
def FilterNodes(nodes, cps, centroids, connectivity, shadowed):
    windward_nodes = []
    for n in range(0, len(nodes)):                                              #Check every node
        USELESS = True
        for c in range(0, len(centroids)):                                      #Check every centroid
            conn = connectivity[c]                                              #Find all the nodes of the element
            for i in range(0, len(conn)):
                node_of_c = conn[i]
                if n == node_of_c:                                              #If the node is a part of the element 
                    if shadowed[c] == 'False':                                  #if the elemen is not shadowed
                        USELESS = False
        if not USELESS:
            windward_nodes.append('True')
        else:
            windward_nodes.append('False')

    return windward_nodes
   
   
   
   

