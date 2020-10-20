import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import dynetx as dn

# a = np.random.randint(0, 2, size=(10, 10))
# print(a)
# D = nx.DiGraph(a)
# print(type(D))

# D.clear()

# edgelist for directed graph with arrow facing last node label, weight given as 1
edgelist = [('RAN', 'CLS'),('RAN', 'CLS'),('GRS', 'RAN'), ('AMP', 'KTC')]

# D.add_weighted_edges_from(edgelist)
# D.add_nodes_from(nodelist)
# nx.draw(D,with_labels=True, font_weight='bold')
# plt.show()
# next step is to create dynamical network with time label and add vertices to different times in a manual way

# directed graph
Dn = dn.DynDiGraph()
# undirected graph
U = dn.DynGraph()

# U.add_nodes_from(nodelist)

# could use tiny loop to create nodelist
# nodelist = [('RAP', {'lat': 180, 'long': 240}),
# ('RAN', {'lat': 18, 'long': 150}),('CLS', {'lat': 18, 'long': 120}),('GRS', {'lat': 170, 'long': 27}),
# ('KTC', {'lat': 189, 'long': 154}),('AMP', {'lat': 189, 'long': 10}),('PPS', {'lat': 13, 'long': 130}),
# ('RAP', {'lat': 134, 'long': 132})]

nodelist=['RAN','RAP','KTC','AMP','CLS','PPS','GRS']

# Need to give nodes attributes not edges in order to later filter by lattitude

Dn.add_nodes_from(nodelist)

# Dn.add_node('RAL', lat=180,long=270, t=2)

# now add egdes for different time stamps, direction goes from e.g 'RAN' to 'CLS' in 1st case
# nodes do not need to be added if there is an edgelist

Dn.add_interaction('RAN', 'CLS', t=0)
Dn.add_interaction('RAN', 'KTC', t=0)
Dn.add_interaction('RAN', 'KTC', t=2)

# small code snippet to show how to create array of networks for data managment

# Dn1 = dn.DynDiGraph()
# Dn1.add_interaction('RAN', 'AMP', t=2)

# dna = [Dn1,Dn]
# print(dna)
# dna[0].add_interaction('RAN', 'GRS', t=2)
# nx.draw(Dn1,with_labels=True, font_weight='bold')
# plt.show()

# adding edges from list
Dn.add_interactions_from(edgelist, t=3)

print(dn.nodes(Dn))


print([e for e in Dn.interactions_iter(['RAN','KTC'])])

Dn = Dn(['RAN','CLS'])

nx.draw(Dn,with_labels=True, font_weight='bold')
plt.show()

# for i in [0,1,2,3]:
	
# 	Dn.time_slice(i,i+1)

# 	nx.draw(Dn,with_labels=True, font_weight='bold')
# 	plt.show()

# a = np.array([1,2,3,4])
# b = np.array([2,3,4,5])
# a = np.stack((a,b))
# print(a)
# print(np.shape(a))
