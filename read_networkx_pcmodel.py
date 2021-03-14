import networkx as nx
import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
import time
import functools
from collections import Counter
from nxviz import MatrixPlot
import datetime

def load_arrays():
	'''function to load text arrays to load 4 directed and 4 undirected netowrks from label string arrays nlc and dlc (both arrays with four strings) 
	to create dynamical network objects which can be used for plotting where magp is the cluster to be used and comp the component n,e,z'''
	
	nlc = [[],[],[],[]]

	dlc = [[],[],[],[]]

	for i in [0,1,2,3]:

		# nlc, dcl array of names of files to read into text arrays to load files

		nlc[i] = f'networks_data/na{i}test.txt'
		
		dlc[i] = f'networks_data/dna{i}test.txt'

		# netdata1 = open(nlc[i],"rb")

		# netdata2 = open(dlc[i],"rb")

		# nlc and dlc text array overwritted to contain networks

		nlc[i] = nx.read_edgelist(nlc[i])

		dlc[i] = nx.read_edgelist(dlc[i],create_using=nx.DiGraph)

	return dlc ,nlc


dna, na = load_arrays()

print(list(dna[1].edges(data=True)))


# G = nlc[1]

# print(G.edges(['A08']))

# # use G.in_edges(node) G.out_edges(node) for directed graphs

# print(G['A08'])

# print(G.nodes())
# print(list(G.nodes)[0]) #for accesing elemnets, same with edges

# # number of edges
# print(len(G.edges()))

# print(G.nodes(data=True)) #also shows meta data for nodes

# # list comprehension, counter counts the number of entires in a list, benifts include compact code (less lines and no new empty arrays need to be defined)

# mf_counts = Counter([d['attr_dict']['t_window'] for n1, n2, d in G.edges(data=True)])

# # print(mf_counts)

# # in networkX graphs are represented as nested dictonaries, key in node id and values dictonary of node attrbutes
# #  edges are part of the attribute of G.edge which is also a nested dictonary accessed using G.edge[n1,n2]['attr_name']
# # node any non hashable object (immutable) i.e not lists

# print([d['attr_dict']['t_window'] for n1,n2,d in G.edges(data=True)])

# print(type(list(G.edges(data=True))))

# t_slice = [(n1,n2) for n1,n2,d in G.edges(data=True) if d['attr_dict']['t_window']==180]

mltrange = {'dawn':[3,9],'noon':[9,15],'dusk':[15,21],'midnight':[21,3]}

mlt = [(u,v) for u,v,d in dna[2].edges(data=True) if d['attr_dict']['t_window']==180 and 3 <= (d['attr_dict']['MLT1'] and d['attr_dict']['MLT2'])<=9]

# print(t_slice,'tslice')

print('mlt',mlt)

# degrees = sum([val for (node, val) in G.degree()]) / len(G.nodes())

# print(degrees, 'avg_degree', degrees/len(G.nodes()))

def save_k_vs_t(dn, nn):
	
	''' function to calculate number of connections (k) at each time window
	in for all network arrays in dnn and nn to and returns file of k vs t to later use
	and speed up plotting of related graphs'''

	avg_deg_matrix_dn = [[],[],[],[]]

	avg_deg_matrix_n = [[],[],[],[]]

	times = [[],[],[],[]]

	for i in [0,1,2,3]:

		# obtains the ordered time stamps in the network

		time_stamps_dn = [d['attr_dict']['t_window'] for n1,n2,d in dn[i].edges(data=True)]

		time_stamps_nn = [d['attr_dict']['t_window'] for n1,n2,d in nn[i].edges(data=True)]

			# key= lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))
		# need to a

		# remove duplicates

		# statment to pick longest time stamp range to set for both arrays to allow for over plotting


		if len(time_stamps_dn) > len(time_stamps_nn):
			
			time_stamps = sorted(list(set(time_stamps_dn)))
		
		else:

			time_stamps = sorted(list(set(time_stamps_nn)))



		print(time_stamps, 'time_stamps')

		# can also save relevent utc ARRAY

		# loop for slicing consecutive time stamps and calculating degree

		print(f'Pc{i+2}')

		for j in time_stamps:

			ts_edge_list1 = [ (n1,n2) for n1,n2,d in dn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

			ts_edge_list2 = [ (n1,n2) for n1,n2,d in nn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

			mlt = [(u,v,d['attr_dict']['t_window'],d['attr_dict']['MLT1']) for u,v,d in dn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and 10 <= (d['attr_dict']['MLT1'] and d['attr_dict']['MLT2'])<=20]

			print('mlt',mlt)

			utc_times = [(j, d['attr_dict']['UTC2']) for n1,n2,d in dn[i].edges(data=True) if d['attr_dict']['t_window'] == j]

			# utc_times = sorted(list(set(utc_times)), key=lambda tup: tup[0]) 

			times[i].append(list(set(utc_times)))

			# count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
			
			# take a unique list of nodes in given time window

			# like double for loop, one loop for acessing tuple value and the next for unpacking it

			ts_node_list1 = list(set([val for tup in ts_edge_list1 for val in tup]))

			ts_node_list2 = list(set([val for tup in ts_edge_list2 for val in tup]))


			# number of edges divided by number of nodes for directed graph
		
			if len(ts_node_list1) == 0:

				avg_deg_matrix_dn[i].append(np.nan)

			else:

				avg_deg = len(ts_edge_list1) / len(ts_node_list1)

				# print(avg_deg, j, 'avg deg, count rc')

				avg_deg_matrix_dn[i].append(avg_deg)

				

			# code segment for undirnetwork

			if len(ts_node_list2) == 0:

				avg_deg_matrix_n[i].append(np.nan)

			else:

				avg_deg = len(ts_edge_list2) / len(ts_node_list2)

				print(avg_deg, j, 'avg deg, count rc')

				avg_deg_matrix_n[i].append(avg_deg)
			
			# avg_deg_matrix_dn[i].append(time_stamps)

		times[i] = sorted(times[i]) #key= lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S')


		print(avg_deg_matrix_dn[i], 'HERE')

		
		# print(avg_deg_matrix_n, 'HERE here')
		print(len(avg_deg_matrix_dn[i]),times[i], time_stamps,'len deg vs len utc dn vs timestamps')

	return avg_deg_matrix_dn, avg_deg_matrix_n, times


		# np.save(f'networks_data/avg_deg_{label}Pc{i+2}_dn.npy',avg_deg_matrix_dn)
		# np.save(f'networks_data/avg_deg_{label}Pc{i+2}_n.npy',avg_deg_matrix_n)



def plot_k_vs_t(avg_k_dn, avg_k_n, time):

	'''using file produced from function save_degree_k plots average degree with label the name of the cluster used to
	find file in directory and label plots'''

	# loading values

	# step_size = np.load('step_size_arr.npy')

	fig, ax = plt.subplots(nrows=4, figsize=(8, 15), facecolor='w', edgecolor='k',sharex='row')
	fig.subplots_adjust(hspace=0.8)


	for i in [0,1,2,3]:

		# avg_deg_matrix_dn = np.load(f'networks_data/avg_deg_{label}Pc{i+2}_dn.npy',allow_pickle=True, fix_imports=True)
		# avg_deg_matrix_n = np.load(f'networks_data/avg_deg_{label}Pc{i+2}_n.npy',allow_pickle=True, fix_imports=True)

		# # total time in seconds 
		# t = 4*3600

		# x is giving the time axis 

		# print(time[i],'time[i]')

		x = list(zip(*time[i]))

		x = [ datetime.datetime.strptime(n2, '%Y-%m-%dT%H:%M:%S') for tup in time[i] for n1,n2 in tup ]

		print('here',time[i],x, type(x[0]))

		ax[i].plot(x, avg_k_dn[i], color='r', label=f'Pc{i+2} directed network')

		ax[i].plot(x, avg_k_n[i], label= f'Pc{i+2} instantaneously directed network')

		ax[i].set_xlabel('time (UTC)')

		ax[i].set_ylabel('average degree')

		# ax[i].set_xlim(0,4)

		ax[i].legend(bbox_to_anchor=(0,1.02,1,0.2),loc="lower left",
		            mode="expand", borderaxespad=0, ncol=2)

		ax[i].grid()

	plt.show()

avg_deg__dn, avg_deg__n, time = save_k_vs_t(dna, na) 

# plot_k_vs_t(avg_deg__dn, avg_deg__n, time)


# nx.draw(G)

# m = MatrixPlot(G)
# m.draw()
# plt.show()

# print(list(G.edges(data=True))[0]['attr_dict']['t_window'])

# time_ordered_edge_list = sorted(G.edges(), key = lambda x: x['attr_dict']['t_window'] )

# print(time_ordered_edge_list)

# can create a sub graph by inputing nbunch


# for n, nbrs in G.items():
#    for nbr, eattr in nbrs.items():
#        wt = eattr['weight']
#        if wt < 0.5: print(f"({n}, {nbr}, {wt:.3})")

# need to work with time windows to need to use that to make many a subgraph at each time window just like with timeslice in dynetx
# for each time window subgraph take the metric (such as avg degree) and save it to an array also take the id
# can create a graph with both metic and id lists
# can resuse parts of old code to do this!!!!!

# for clusters first creat a subgraph with the given mlt then repeat above


