import networkx as nx
import matplotlib.pyplot as plt
import dynetx as dn
import itertools
import numpy as np
from dynetx.readwrite import json_graph
import json


# using dictionary to map strings to variables

nl =["na0.json","na1.json","na2.json","na3.json"]

dnl = ["dna0.json","dna1.json","dna2.json","dna3.json"]

na = [[],[],[],[]]

dna = [[],[],[],[]]

for i in range(len(dnl)):

	# netdata1 = open(nl[i],"rb")

	# netdata2 = open(dnl[i],"rb")

	with open(nl[i]) as json_file:

		netdata1 = json.load(json_file)

	with open(dnl[i]) as json_file:
    	
		netdata2 = json.load(json_file)

	# netdata1 = json_graph.node_link_data(nl[i])

	# netdata2 = json_graph.node_link_data(dnl[i])


	na[i] = dn.readwrite.json_graph.node_link_graph(netdata1)

	dna[i] = dn.readwrite.json_graph.node_link_graph(netdata2, directed=True)

# print(list(dna[1].stream_interactions())[-1][3])

print(list(dna[3].stream_interactions()))

for i in [1,2]:

	print(i)

	# gives time for last interaction

	epochs = list(dna[i].stream_interactions())[-1][3]

	for j in range(epochs):

		# plotting each time slice

		fig, axs = plt.subplots(figsize=(8, 8), nrows=2)

		axs = axs.flatten()

		# print(j)

		dns = dna[i].time_slice(j)

		ns = na[i].time_slice(j)

		axs[0].set_axis_off()

		axs[1].set_axis_off()

		nx.draw(dns, ax=axs[0],with_labels=True, font_weight='bold',arrowsize=20, edgecolor='red',width=1.2)

		nx.draw(ns, ax=axs[1],with_labels=True , font_weight='bold',edgecolor='orange',width=1.2)

		axs[0].set_title(f'digraph with window epoch {j} out of {epochs}, Pc{i+2}')

		axs[1].set_title(f'graph with window epoch {j} out of {epochs}, Pc{i+2}')

		plt.show()



avg_deg_matrix_dn = [[],[],[],[]]

avg_deg_matrix_n = [[],[],[],[]]

a = np.load('step_size_arr.npy')

print()


for i in [0,1,2,3]:

	fig, axs = plt.subplots(figsize=(8, 8), nrows=2)

	# obtains the last time stamp in the network

	epochs = list(dna[i].stream_interactions())[-1][3]

	# loop for slicing consecutive time stamps and calculating degree

	print(f'Pc{i}')

	for j in range(epochs):

		dns = dna[i].time_slice(j)

		ns = na[i].time_slice(j)

		# number of edges divided by number of nodes for directed graph

		print(dns.size(),dns.order(), 'edges','nodes')
		
		# code segment for dirnetwork
		if dns.order()!= 0:

			k = dns.size()/dns.order()

			print(k, j, 'avg deg, count rc')

			avg_deg_matrix_dn[i].append(k)

		else:

			avg_deg_matrix_dn[i].append(np.nan)

		# code segment for undirnetwork

		if ns.order()!= 0:

			k = ns.size()/ns.order()

			print(k, j, 'avg deg, count rc')

			avg_deg_matrix_n[i].append(k)

		else:

			avg_deg_matrix_n[i].append(np.nan)


	countdn = np.arange(len(avg_deg_matrix_dn[i]))

	countn = np.arange(len(avg_deg_matrix_n[i]))

	print(countdn,countn)

	axs[0].set_title('dirnetwork <k>')

	axs[0].plot(countdn,avg_deg_matrix_dn[i])

	axs[1].set_title('undirnetwork <k>')

	axs[1].plot(countn, avg_deg_matrix_n[i])

	plt.show()

		#can multiply this with stepsize/2, mid part of window, where stepsize defined in program

		# what to do with time at zero?

		# indexing starts at zero for enumerate zip, so zero is 1 in this case

		# len includes 0 in sum so len(min0lags)*count pretty close to len(ts1)

		# need to add 1 to all values in count before multiplication by count/2
			

		# plotting each time slice


		# print(j)








# Return an iterator for (node, degree) at time t, 

# G = dn.DynGraph()
# G.add_path([0,1,2,3], t=0)

# print(list(G.degree_iter(t=0)))

# 'need to use DynGraph.size([t])	Return the number of edges at time DynGraph.order([t]) Return the number of nodes in the t snpashot of a dynamic graph.'



# dirnet_dict = {dna0: 0, dna1: 1, dna2: 2, dna3: 3}

# net_dict = {na0: 0, na1: 1, na2: 2, na3: 3}


# for var in dirnet_dict:

# 	var = dn.read_interactions(var, nodetype=int, timestamptype=int)

