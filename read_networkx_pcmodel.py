import networkx as nx
import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
import time
import functools
from collections import Counter
#from nxviz import MatrixPlot
import datetime
import matplotlib.dates as mdates
import pdb


def utc_sort_reduce(utc_times):

	# takes an unordered list of utc times and removes repeated values and then orders based on tme value

	times_dn = list(set(utc_times))

	times_dn = sorted(times_dn, key = lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))

	return times_dn

def inter_cluster_pairs(mlta):
	
	# returns all adjecent pairs from a list, used for interscluster code
	
	mltar = np.roll(mlta,-1)

	adjclusts = list(zip(mlta, mltar))

	return adjclusts

def load_arrays(Pcbands, comp):
	'''function to load text arrays to load 4 directed and 4 undirected networks from label string arrays nlc and dlc (both arrays with four strings) 
	to create dynamical network objects which can be used for plotting where magp is the cluster to be used and comp the component n,e,z'''
	
	nlc = [[],[],[],[]]

	dlc = [[],[],[],[]]

	for i in Pcbands:

		# nlc, dcl array of names of files to read into text arrays to load files

		nlc[i] = f'networks_data/dna{i}_{comp}_spd.txt'
		
		dlc[i] = f'networks_data/na{i}_{comp}_spd.txt'

		print(f'loading Pc{i+2} networks')

		# nlc and dlc text array overwritted to contain networks

		# if i == 3 :

		# 	nlc[i] = nx.Graph()

		# else:
		# 	nlc[i] = nx.read_edgelist(nlc[i])

		nlc[i] = nx.read_edgelist(nlc[i])

		dlc[i] = nx.read_edgelist(dlc[i],create_using=nx.DiGraph)

	return dlc ,nlc


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

# mltrange = {'dawn':[3,9],'noon':[9,15],'dusk':[15,21],'midnight':[21,3]}

# mlt = [(u,v) for u,v,d in dna[2].edges(data=True) if d['attr_dict']['t_window']==180 and 3 <= (d['attr_dict']['MLT1'] and d['attr_dict']['MLT2'])<=9]

# print(t_slice,'tslice')

# degrees = sum([val for (node, val) in G.degree()]) / len(G.nodes())

# print(degrees, 'avg_degree', degrees/len(G.nodes()))

def save_k_vs_t_global(dn, nn):
	
	''' function to calculate number of connections (k) globally at each time window
	in for all network arrays in dnn and nn to and returns file of k vs t to later use
	and speed up plotting of related graphs'''

	# arrays used for plotting gaphs see cluster function for details

	dict_dn = {}

	dict_n ={}

	for i in ['times', 'avg_k']:

		dict_dn[i] = {0:[],1:[],2:[],3:[]}

		dict_n[i] = {0:[],1:[],2:[],3:[]}

	for i in [0,1,2,3]:

		# obtains only UTC times which have edges

		utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in dn[i].edges(data=True) ]

		# removing duplicates in utc_times_dn

		times_dn = list(set(utc_times_dn))

		# ordering via time order

		times_dn = sorted(times_dn, key = lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))

		dict_dn['times'][i].append(times_dn)

		utc_times_n = [ d['attr_dict']['UTC2'] for n1,n2,d in nn[i].edges(data=True) ]

		times_n = list(set(utc_times_nn))

		times_n = sorted(times_n, key = lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))

		dict_n['times'][i].append(times_n)
	
		utc_times_nn = [ d['attr_dict']['UTC2'] for n1,n2,d in nn[i].edges(data=True)]

		times_n = list(set(utc_times_nn))

		times_n = sorted(times_n, key = lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))

		dict_n['times'][i].append(times_n)

		# obtains the ordered time stamps in the network

		time_stamps_dn = [d['attr_dict']['t_window'] for n1,n2,d in dn[i].edges(data=True)]

		time_stamps_nn = [d['attr_dict']['t_window'] for n1,n2,d in nn[i].edges(data=True)]

		# remove duplicates
			
		time_stamps_dn = sorted(list(set(time_stamps_dn)))
		
		time_stamps_nn = sorted(list(set(time_stamps_nn)))

		# can also save relevent utc ARRAY

		# loop for slicing consecutive time stamps and calculating degree for directed network

		for j in time_stamps_dn:

			ts_edge_list1 = [ (n1,n2) for n1,n2,d in dn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

			# count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
			
			# take a unique list of nodes in given time window

			# like double for loop, one loop for acessing tuple value and the next for unpacking it

			ts_node_list1 = list(set(np.reshape(ts_edge_list1, len(ts_edge_list1)*2)))

		
			avg_deg = len(ts_edge_list1) / len(ts_node_list1)

			# print(avg_deg, j, 'avg deg, count rc')

			dict_dn['avg_k'][i].append(avg_deg)


		# loop for slicing consecutive time stamps and calculating degree for undirected network
			
		for j in time_stamps_nn:

			ts_edge_list2 = [ (n1,n2) for n1,n2,d in nn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

			# count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
			
			# take a unique list of nodes in given time window

			# like double for loop, one loop for acessing tuple value and the next for unpacking it

			ts_node_list2 = list(set([val for tup in ts_edge_list2 for val in tup]))
				
			avg_deg = len(ts_edge_list2) / len(ts_node_list2)

			dict_n['avg_k'][i].append(avg_deg)


	np.save(f'networks_data/spd_analysis/avg_deg_global_dn.npy',dict_dn)
	np.save(f'networks_data/spd_analysis/avg_deg_global_n.npy',dict_n)

	# files.download(f'networks_data/spd_analysis/avg_deg_Pc_dn.npy',avg_deg_dict_dn)
	# files.download(f'networks_data/spd_analysis/avg_deg_Pc_n.npy',avg_deg_dict_n)					


def save_k_vs_t_cluster(dn, nn, cluster, save_label):
	
	''' function to calculate number of connections (k) for a given MLT range given by window
	 at each time window in for all network arrays in dnn and nn to and returns file of k vs t to later use
	and speed up plotting of related graphs'''

	# magnetic local time 0 is midnight 6 dawn, 12 noon, 18 dusk 

    # dawn: 03:00-09:00 (ie 06:00 +/-3hrs)

    # noon: 09:00:15:00 (12:00 +/- 3hrs)

    # dusk:15:00-21:00

    # midnight (21:00-03:00)

	# arrays used for plotting gaphs use timesvalues in dict_dn['times'][i].app, for <k> dict_dn['avg_k'][i]
	# for nodes .app,dict_dn['n_nodes'][i].app in total each dict has 12 arrays to store values inside with shape [[[pc2t], [pc3t], [pc4t], [pc5t]], [[pc2k], [], [], []], [[pc2n], [], [], []]]


	# dict within dict for assings Pc values for all bands used for relevent parameter to store values

	dict_dn = {}

	dict_n ={}

	for i in ['times', 'avg_k', 'n_nodes']:

		dict_dn[i] = {0:[],1:[],2:[],3:[]}

		dict_n[i] = {0:[],1:[],2:[],3:[]}

	print('dict', dict_dn)
	
	# dictonary and values used to obtain mlt ranges to filter for

	mltr = {'dawn':[3,9],'noon':[9,15],'dusk':[15,21],'midnight':[21,3]}

	mlt1 = mltr[cluster][0]

	mlt2 = mltr[cluster][1]

	for i in [0,1,2,3]:

		# time stamps for looping, obtains only those timestamps which have edges

		time_stamps_dn = [d['attr_dict']['t_window'] for n1,n2,d in dn[i].edges(data=True) 
		if mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

		# could append none to empty degree array to use for plotting

		if len(time_stamps_dn) == 0:

			print(f'no values for Pc{i+2} {cluster} directed network array ')

			continue

		# obtains only UTC times which have edges

		utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in dn[i].edges(data=True) 
		if mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

		# removing duplicates in utc_times_dn

		times_dn = list(set(utc_times_dn))

		# ordering via time order

		times_dn = sorted(times_dn, key = lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))

		dict_dn['times'][i].extend(times_dn)

		# remove duplicates
			
		time_stamps_dn = sorted(list(set(time_stamps_dn)))
		
		# loop for slicing consecutive time stamps and calculating degree for directed network
		
		print(f'pc{i+2}, dir {cluster}')

		for j in time_stamps_dn:

			# print(f'pc{i+2}, dir {cluster}, {j} out of {max(time_stamps_dn)}')

			mlt_edge_list1 = [ [ n1, n2 ] for n1,n2,d in dn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

			# utc_check = [(n1, n2 ,j, d['attr_dict']['UTC2'], d['attr_dict']['MLT1'], d['attr_dict']['MLT2']) for n1,n2,d in dn[i].edges(data=True) 
			# if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]
			
			# print(f'checking dir {cluster} clusters pc{i+2}', utc_check)
			
			# reduced dim of edgelist and remove duplicated

			mlt_node_list1 = list(set(np.reshape(mlt_edge_list1, len(mlt_edge_list1)*2)))

			num_nodes = len(mlt_node_list1)

			avg_deg = len(mlt_edge_list1) / num_nodes

			# print(avg_deg, j, 'avg deg, count rc')

			dict_dn['avg_k'][i].append(avg_deg)

			dict_dn['n_nodes'][i].append(num_nodes)

	# loop for slicing consecutive time stamps and calculating degree for undirected network

	for i in [0,1,2,3]:

		print(f'pc{i+2}, undir {cluster}')

		time_stamps_nn = [d['attr_dict']['t_window'] for n1,n2,d in nn[i].edges(data=True) 
		if mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

		if len(time_stamps_nn) == 0:

			print(f'no values for Pc{i+2} {cluster} directed network array ')

			continue
		
		utc_times_nn = [ d['attr_dict']['UTC2'] for n1,n2,d in nn[i].edges(data=True) 

		if mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

		times_n = list(set(utc_times_nn))

		times_n = sorted(times_n, key = lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))

		dict_n['times'][i].extend(times_n)

		time_stamps_nn = sorted(list(set(time_stamps_nn)))

			
		for j in time_stamps_nn:

			# print(f'pc{i+2}, undir ,{cluster}, {j} out of {max(time_stamps_nn)}')

			mlt_edge_list2 = [ [ n1,n2 ] for n1,n2,d in nn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

			mlt_node_list2 = list(set(np.reshape(mlt_edge_list2, len(mlt_edge_list2)*2)))
				
			num_nodes2 = len(mlt_node_list2)

			avg_deg = len(mlt_edge_list2) / num_nodes2

			dict_n['avg_k'][i].append(avg_deg)

			dict_n['n_nodes'][i].append(num_nodes2)


		# print('times_dn',times_dn, len(times_dn))
		# print('times_n', times_n, len(times_n))

		# print(f'dict dn times {i}',dict_dn['times'][i], len(dict_dn['times'][i]))
		# print(f'dict n times {i}',dict_n['times'][i],len(dict_n['times'][i]))

		# print(dict_dn['avg_k'][i], len(dict_dn['avg_k'][i]))
		# print(dict_n['avg_k'][i],len(dict_n['avg_k'][i]))

	# for saving files more completley in the future can add the network array with comp value to output save file

	np.save(f'{save_label}_dn.npy', dict_dn)
	np.save(f'{save_label}_n.npy', dict_n)

	return dict_dn, dict_n


def save_k_vs_t_intercluster(dn, nn, Pcbands, cluster1, cluster2, save_label):
	
	''' function to calculate number of connections (k) for a between MLT ranges cluster1 and cluster2
	given by window  at each time window in for all network arrays in dnn and nn to and returns file of k vs t to later use
	and speed up plotting of related graphs'''

	# magnetic local time 0 is midnight 6 dawn, 12 noon, 18 dusk 

    # dawn: 03:00-09:00 (ie 06:00 +/-3hrs)

    # noon: 09:00:15:00 (12:00 +/- 3hrs)

    # dusk:15:00-21:00

    # midnight (21:00-03:00)

	# dict within dict for assings Pc values for all bands used for relevent parameter to store values

	dict_dn = {}

	dict_n ={}

	for i in ['times', 'times_in', 'times_out', 'avg_k', 'avg_k_in', 'avg_k_out' ,'n_nodes','n_nodes_in','n_nodes_out']:

		dict_dn[i] = {0:[],1:[],2:[],3:[]}

	for i in ['times', 'avg_k', 'n_nodes']:

		dict_n[i] = {0:[],1:[],2:[],3:[]}
	
	# dictonary and values used to obtain mlt ranges to filter for

	mltr = {'dawn':[3,9],'noon':[9,15],'dusk':[15,21],'midnight':[21,3]}

	mlt11 = mltr[cluster1][0]

	mlt12 = mltr[cluster1][1]

	mlt21 = mltr[cluster2][0]
	
	mlt22 = mltr[cluster2][1]

	for i in Pcbands:

		# time stamps for looping, obtains only those timestamps which have edges and
		# are within cluster1 or cluster2 mlt ranges

		time_stamps_dn = [d['attr_dict']['t_window'] for n1,n2,d in dn[i].edges(data=True) 
		if mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22 
		or mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22 ]

		# could append none to empty degree array to use for plotting

		if len(time_stamps_dn) == 0:

			print(f'no values for Pc{i+2} {cluster1}-{cluster2} directed network array')

			continue

		# obtains only UTC times which have edges

		utc_times_dn = [d['attr_dict']['UTC2'] for n1,n2,d in dn[i].edges(data=True) 
		if mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22 
		or mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22 ]

		utc_times_in = [d['attr_dict']['UTC2'] for n1,n2,d in dn[i].edges(data=True) 
		if mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22 ]

		utc_times_out = [d['attr_dict']['UTC2'] for n1,n2,d in dn[i].edges(data=True) 
		if mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22]

		# removing duplicates in utc_times_dn

		times_dn = utc_sort_reduce(utc_times_dn)

		times_in = utc_sort_reduce(utc_times_in)

		times_out = utc_sort_reduce(utc_times_out)

		dict_dn['times'][i].extend(times_dn)

		dict_dn['times_in'][i].extend(times_in)

		dict_dn['times_out'][i].extend(times_out)

		# remove duplicates
			
		time_stamps_dn = sorted(list(set(time_stamps_dn)))
		
		# loop for slicing consecutive time stamps and calculating degree for directed network
		
		print(f'pc{i+2}, dir {cluster1} {cluster2} out of 4 pairs')

		for j in time_stamps_dn:

			# print(f'pc{i+2}, dir {cluster}, {j} out of {max(time_stamps_dn)}')

			# MLT2 will always be the MLT value of the second node and MLT1 of the first node

			# first node (MLT1) in the second cluster and the second node (MLT2) in the first cluster B -> A in

			mlt_edge_list_in = [ [ n1, n2 ] for n1,n2,d in dn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22]

			# first node in the first cluster and the second node in the second cluster A -> B in

			mlt_edge_list_out = [ [ n1, n2 ] for n1,n2,d in dn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22]

			# will need to check intercluster code by ensuring each node mlt out of n1 and n2 in different cluster

			# utc_check = [(n1, n2 ,j, d['attr_dict']['UTC2'], d['attr_dict']['MLT1'], d['attr_dict']['MLT2']) for n1,n2,d in dn[i].edges(data=True) 
			# if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]
			
			# print(f'checking dir {cluster} clusters pc{i+2}', utc_check)
			
			# reduced dim of edgelist and remove duplicated

			mlt_node_list_in = list(set(np.reshape(mlt_edge_list_in, len(mlt_edge_list_in)*2)))

			mlt_node_list_out = list(set(np.reshape(mlt_edge_list_out, len(mlt_edge_list_out)*2)))

			num_nodes_in = len(mlt_node_list_in)

			num_nodes_out = len(mlt_node_list_out)

			# prevent double counting of nodes from both lists

			# print(mlt_node_list_in)

			# print(mlt_node_list_out)

			# need to do addition before list-set

			all_nodes_list = mlt_node_list_in + mlt_node_list_out

			all_nodes = list(set(all_nodes_list))

			# print(j,'j stamps',all_nodes)

			num_edge_in = len(mlt_edge_list_in)

			num_edge_out = len(mlt_edge_list_out)

			# might be an error

			avg_deg = (num_edge_in + num_edge_out) / len(all_nodes)

			# print(avg_deg, j, 'avg deg, count rc')

			dict_dn['avg_k'][i].append(avg_deg)

			if num_nodes_in != 0:

				dict_dn['avg_k_in'][i].append(num_edge_in/num_nodes_in)

				dict_dn['n_nodes_in'][i].append(num_nodes_in)

			if num_nodes_out != 0:

				dict_dn['avg_k_out'][i].append(num_edge_out/num_nodes_out)

				dict_dn['n_nodes_out'][i].append(num_nodes_out)

			dict_dn['n_nodes'][i].append(len(all_nodes))





	# loop for slicing consecutive time stamps and calculating degree for undirected network

	# two for loops due to continue statmetns to check for time-stamps for relevent Pc bands

	for i in Pcbands:

		print(f'pc{i+2}, undir {cluster1} {cluster2} out of 4 pairs')

		time_stamps_nn = [d['attr_dict']['t_window'] for n1,n2,d in nn[i].edges(data=True) 
		if mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22 
		or mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22 ]

		if len(time_stamps_nn) == 0:

			print(f'no values for Pc{i+2} {cluster1}-{cluster2} directed network array ')

			continue
		
		utc_times_nn = [ d['attr_dict']['UTC2'] for n1,n2,d in nn[i].edges(data=True) 
		if mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22 
		or mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22 ]

		times_n = utc_sort_reduce(utc_times_nn)

		dict_n['times'][i].extend(times_n)

		time_stamps_nn = sorted(list(set(time_stamps_nn)))

			
		for j in time_stamps_nn:

			# print(f'pc{i+2}, undir ,{cluster}, {j} out of {max(time_stamps_nn)}')

			mlt_edge_list2 = [ [ n1,n2 ] for n1,n2,d in nn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22 
			or mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22 ]

			mlt_node_list2 = list(set(np.reshape(mlt_edge_list2, len(mlt_edge_list2)*2)))
				
			num_nodes2 = len(mlt_node_list2)

			avg_deg = len(mlt_edge_list2) / num_nodes2

			dict_n['avg_k'][i].append(avg_deg)

			dict_n['n_nodes'][i].append(num_nodes2)


		# print('times_dn',times_dn, len(times_dn))
		# print('times_n', times_n, len(times_n))

		# print(f'dict dn times {i}',dict_dn['times'][i], len(dict_dn['times'][i]))
		# print(f'dict n times {i}',dict_n['times'][i],len(dict_n['times'][i]))

		# print(dict_dn['avg_k'][i], len(dict_dn['avg_k'][i]))
		# print(dict_n['avg_k'][i],len(dict_n['avg_k'][i]))

	# for saving files more completley in the future can add the network array with comp value to output save file

	np.save(f'{save_label}_dn.npy', dict_dn)
	np.save(f'{save_label}_n.npy', dict_n)

	return dict_dn, dict_n
				


def plot_k_vs_t(dict_dn, dict_n, Pcbands, comp, dykey1, dykey2, dtimeskey1, dtimeskey2, uykey1, uykey2, utimeskey1, utimeskey2, plots, label, save_plot = False):

	'''code to plot the average degree of networkx directed and undirected multi-dim network parameters dict_dn, dict_n 
	with label the name of single cluster or two cluster (for connection between clusters)

	dictonary keys to plot relevent values {ykey1} as avg_deg_in andf {ykey2} as nodes_in with relevent time keys
	for both directed and undirected networks denoted with d... or u...
	
	plots can be either dir, undir or both to specify number of network plots on a sungle subplot

	'''

	# check which Pc bands are not empty

	non_empty_pcbands = []
	
	for i in Pcbands:

		if len(dict_dn[f'{dykey1}'][i]) != 0:

			print(f' {label} non-empty plot Pc',i+2)

			non_empty_pcbands.append(i)

	# if all no Pc band values

	if len(non_empty_pcbands)==0:

		if plots == 'dir':

			print(f'no network values for {label}, {dykey1}, {comp}')

		else:

			print(f'no network values for {label}, {uykey1}, {comp}')


	else:

		print('relevent pc bands',len(non_empty_pcbands), non_empty_pcbands)

		fig, ax = plt.subplots(nrows= 6 + len(non_empty_pcbands), figsize=(8, 15), facecolor='w', edgecolor='k')
		fig.subplots_adjust(hspace=0.8)

		md = pd.read_csv('SME_SMR_spd.csv')

		md2 = pd.read_csv('spd_Bzgsm_p.csv')

		# print(md.head())

		index_date = md['Date_UTC']

		index_date2 = md2['Date_UTC']

		sme = md['SME']

		smr = md['SMR']

		print('start_time, end_time',index_date.iloc[0],index_date.iloc[-1])

		axsmer = [ datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') for x in index_date ]

		axsmer2 = [ datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') for x in index_date2 ]

		ax[0].plot(axsmer,sme, color='black')

		ax[0].set_xlabel('time (hh:mm:ss)')

		ax[0].set_ylabel('SME (nT)')

		ax[0].grid()

		datemin = np.datetime64(index_date.iloc[0])
		datemax = np.datetime64(index_date.iloc[-1])

		ax[0].set_xlim(datemin,datemax)

		formatter = mdates.DateFormatter("%H:%M:%S")

		ax[0].xaxis.set_major_formatter(formatter)

		ax[1].plot(axsmer,smr, color='g')

		# ax[1].set_xticks(indx_mask)

		# ax[1].set_xticklabels(datel)

		ax[1].set_xlabel('time (hh:mm:ss)')

		ax[1].set_ylabel('SMR (nT)')

		ax[1].grid()

		datemin = np.datetime64(index_date.iloc[0])
		datemax = np.datetime64(index_date.iloc[-1])

		ax[1].set_xlim(datemin,datemax)

		# need an array of non empty pc values

		# also need some kind of counting function

		dyn_pressure = md2['PDYN']

		# GSM field values

		bz, bx, by = md2['GSM_Bz'],md2['GSM_Bx'],md2['GSM_By']


		for i in [0,1,2,3]:

			y_vals = [bx,by,bz, dyn_pressure]

			y_labs = ['GSM_Bz','GSM_Bx','GSM_By','Dynamic pressure']

			c = ['red','black','green','orange']

			ax[i+2].xaxis.set_major_formatter(formatter)

			ax[i+2].plot(axsmer2, y_vals[i], label=y_labs[i], color = c[i])

			ax[i+2].set_xlim(datemin,datemax)

			ax[i+2].set_xlim(datemin,datemax)

			ax[i+2].set_ylim(-30,30)

			ax[i+2].legend()

			ax[i+2].legend()

			ax[i+2].grid()

			if i ==3:

				ax[i+2].set_ylabel('nPa')

			else:

				ax[i+2].set_ylabel('nT')


		# set title for relevent label for all plots

		if label and f'{dykey1}'=='avg_k_in':

			print(label)

			c1 , c2, c3 = label.split('_')

			if c1 == 'dir':

				fig.suptitle(f'{dykey1}  from {c3} to {c2}, directed networks, B_{comp} component')

			else:

				fig.suptitle(f' avg_deg from {c2} {c3}, undirected networks, B_{comp} component')


		elif label  and f'{dykey1}'=='avg_k_out':

			c1 , c2, c3 = label.split('_')
			
			if c1 == 'dir':

				fig.suptitle(f'{dykey1}  from {c2} to {c3}, directed networks, B_{comp} component')

			else:

				fig.suptitle(f' avg_deg from {c2} {c3}, undirected networks, B_{comp} component')

		elif label:

			fig.suptitle(f' {dykey1} {label} networks, B_{comp} component')

		else:

			fig.suptitle(f'global networks, B_{comp} component')

		# can do plotting on laptop however first need updated results matrices including node values

		# run code to obtain network parametsremotley without plotting

		# need to figure out way to plot degree on an independantly scaled axis

		for i in range(len(non_empty_pcbands)):

			pcind = non_empty_pcbands[i]

			if plots == 'all':

				xd = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_dn[f'{dtimeskey1}'][pcind] ]

				xn = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_n[f'{utimeskey1}'][pcind] ]

				ax[i+6].scatter(xd, dict_dn[f'{dykey1}'][pcind], color='black', s=6)

				ax[i+6].scatter(xn, dict_n[f'{uykey1}'][pcind], color='black', s=6)

				ax[i+6].plot(xd, dict_dn[f'{dykey1}'][pcind], color='r', label=f'Pc{i+2} directed network', lw=2)

				ax[i+6].plot(xn, dict_n[f'{uykey1}'][pcind], label= f'Pc{i+2} instantaneously directed network', lw=2)

			elif plots == 'undir':

				xn = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_n[f'{utimeskey1}'][pcind] ]

				xn2 = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_n[f'{utimeskey2}'][i] ]

				ax[i+6].scatter(xn, dict_n[f'{uykey1}'][pcind], color='black', s=6)

				ax[i+6].plot(xn, dict_n[f'{uykey1}'][pcind], label= f'Pc{i+2} instantaneously directed network', lw=2)

				yax2 = ax[i+6].twinx()

				yax2.plot(xn2, dict_n[f'{uykey2}'][i], label='number of nodes', color='green', linestyle='--' )

			elif plots == 'dir':

				xd = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_dn[f'{dtimeskey1}'][pcind] ]

				xd2 = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_dn[f'{dtimeskey2}'][i] ]

				# print statmetns here

				# print('i',i,'xd',dict_dn[f'{dtimeskey1}'][pcind], len(dict_dn[f'{dtimeskey1}'][pcind]))

				# print(f'{dykey1}',dict_dn[f'{dykey1}'][pcind], len(dict_dn[f'{dykey1}'][pcind]))

				# print('i',i,'xd2',dict_dn[f'{dtimeskey2}'][i], len(dict_dn[f'{dtimeskey2}'][i]))

				# print(f'{dykey2}',dict_dn[f'{dykey2}'][i], len(dict_dn[f'{dykey2}'][i]))
			
				ax[i+6].scatter(xd, dict_dn[f'{dykey1}'][pcind], color='orange', s=6)

				ax[i+6].plot(xd, dict_dn[f'{dykey1}'][pcind], color='r', label=f'Pc{i+2} directed network', lw=2)

				yax2 = ax[i+6].twinx()

				yax2.plot(xd2, dict_dn[f'{dykey2}'][i], label='number of nodes', color='orange', linestyle='--' )

			
			# if statmetn to remove empty subplots
			if len(dict_dn[f'{dykey1}'][pcind]) == 0:

				print(f' {label} empty plot Pc',i+2)

				ax[i+6].set_visible(False)

				yax2.set_visible(False)

			yax2.set_ylabel('num nodes')

			ax[i+6].set_xlabel('time (UTC)')

			ax[i+6].set_ylabel('average degree')

			# print(index_date[0],index_date.iloc[-1],type(datel[0]))

			datemin = np.datetime64(index_date.iloc[0])
			datemax = np.datetime64(index_date.iloc[-1])

			ax[i+6].set_xlim(datemin,datemax)

			formatter = mdates.DateFormatter("%H:%M:%S")

			ax[i+6].xaxis.set_major_formatter(formatter)


			ax[i+6].legend(bbox_to_anchor=(0,1.02,1,0.2),loc="lower left",
			            mode="expand", borderaxespad=0, ncol=2)

			ax[i+6].grid()


			if save_plot:

				# saves the plot

				plt.savefig(f'plots/{dykey1[i]}_{label}_net_{comp}.png')


		plt.show()




# # code for saving cluster and intercluster analysis file

# include commpand for plot

# need to check saveplot labels and commands, also

def save_plot_all_clusters(save_k, comp ,plot, option, save_plot):

	# option to chose to save all (option variable) intercluster or cluster analysis files or plot all analysis files for relevent dictonary keys
	# option given within the function using save or plot to true or false
	# with component (comp) for magnetic field 'n', 'e', 'z'
	# j can be wither 'dir', 'undir' or 'all'	

	# variable for chooseing pc-bands to compute 

	Pcband = [0,1,2,3]

	mltz = ['dawn','noon','dusk','midnight']

	if save_k:

		# # load network arrays!

		dna, na = load_arrays(Pcband, comp)

	# loop over network arrays

	for i in [0,1,2,3]:

		# labels for incluster savefiles and plots

		if option == 'intercluster':

			# adjacent mlt labels for analysis and labeling

			adjclusts = inter_cluster_pairs(mltz)

			# file label with file path, file 
			file_label = f'networks_data/spd_analysis/avg_deg_{adjclusts[i][0]}_{adjclusts[i][1]}_{comp}'

		else:

			file_label = f'networks_data/spd_analysis/avg_deg_{mltz[i]}_{comp}'


		# run and save dictonaries for plotting by filtering dynamical network array file

		if save_k:

			if option == 'cluster':

				save_k_vs_t_cluster(dna, na, mltz[i], file_label)

			elif option == 'intercluster':

				save_k_vs_t_intercluster(dna, na, Pcband, adjclusts[i][0], adjclusts[i][1], file_label)

		# plot the analysis files if specified

		if plot:

			if option == 'cluster':

				dict_dn = np.load(f'{file_label}_dn.npy',allow_pickle=True)[()]

				dict_n = np.load(f'{file_label}_n.npy',allow_pickle=True)[()]

				for j in ['dir', 'undir']:

					plot_label = '_'.join( [j, mltz[i]] )

					plot_k_vs_t(dict_dn, dict_n, Pcband, 'e', 'avg_k', 'n_nodes','times', 'times', 'avg_k', 
					'n_nodes','times', 'times', label = plot_label, plots= j, save_plot = False)

			
			elif option == 'intercluster':

				dict_dn = np.load(f'{file_label}_dn.npy', allow_pickle=True)[()]

				dict_n = np.load(f'{file_label}_n.npy', allow_pickle=True)[()]

				# labels for dictionary keys

				# 'avg_k_in', 'avg_k_out', 'n_nodes','n_nodes_in', 'times', 'times_in', 'times_out', 'n_nodes_out'

				# indexer in plotting function chooses Pc band and label a joining of strings from both clusters

				# def plot_k_vs_t(dict_dn, dict_n, Pcbands, comp, dykey1, dykey2, dtimeskey1, dtimeskey2, 
				# uykey1, uykey2, utimeskey1, utimeskey2, plots, label, save_plot = False):

				for j in ['dir', 'undir']:

					plot_label = '_'.join( [ j, adjclusts[i][0], adjclusts[i][1] ] )

					plot_k_vs_t(dict_dn, dict_n, Pcband, 'e', 'avg_k_in', 'n_nodes_in','times_in', 'times_in', 'avg_k', 
					'n_nodes','times', 'times', label = plot_label, plots= j, save_plot = False)



# code for reading and plotting network analysis files 

# save_plot_all_clusters(save_k, comp, ,plot, option, j, save_plot = False):

# save_plot_all_clusters(save_k = False, comp = 'e', plot=True,  option = 'cluster', save_plot = False)

for i in ['n','z','e']:

	for j in ['cluster', 'intercluster']:

		save_plot_all_clusters(save_k = True, comp = i, plot=False,  option = j, save_plot = False)

def simpfunc(k):

	if k:

		print('heeey')

simpfunc(k=False)


# for j in ['dir', 'undir']:

# 	for i in [0,1,2,3]:

# 		'''code for cluster plotting or intercluster ----------------

# 		(dict_dn, dict_n, ykey1, ykey2, timeskey1, timeskey2, plots, label=False):

# 		['avg_k', 'avg_k_in', 'avg_k_out', 'n_nodes','n_nodes_in', 'times', 'times_in', 'times_out', 'n_nodes_out']:'''

# 		# brackets at the end to allow numpy to store saved dictionary

# 		# cluster plotting, first loading the save files



# 		# intercluster plotting

# 		adjclusts = inter_cluster_pairs(mltz)

# 		dict_dn = np.load(f'networks_data/spd_analysis/avg_deg_{adjclusts[i][0]}_{ adjclusts[i][1]}_dn.npy', allow_pickle=True)[()]
		
# 		dict_n = np.load(f'networks_data/spd_analysis/avg_deg_{adjclusts[i][0]}_{ adjclusts[i][1]}_n.npy', allow_pickle=True)[()]

# 		# intercluster label

# 		# labels for dictionary keys

# 		# 'avg_k_in', 'avg_k_out', 'n_nodes','n_nodes_in', 'times', 'times_in', 'times_out', 'n_nodes_out'

# 		# indexer in plotting function chooses Pc band and label a joining of strings from both clusters

# 		plot_k_vs_t(dict_dn, dict_n, Pcband, 'e', 'avg_k_in', 'n_nodes_in','times_in', 'times_in', 'avg_k', 'n_nodes','times', 'times',  
# 			label = '_'.join( [ j, adjclusts[i][0], adjclusts[i][1] ] ), plots=j)

# # code for plotting global networks

# avg_deg_dn, avg_deg_n, time_dn, time_nn = save_k_vs_t_global(dna, na) 

# global network plotting

# avg_deg_dn = np.load('networks_data/spd_analysis/avg_deg_global_dn.npy',allow_pickle=True, fix_imports=True)
# avg_deg_n = np.load('networks_data/spd_analysis/avg_deg_global_n.npy',allow_pickle=True, fix_imports=True)

# plot_k_vs_t(avg_deg_dn, avg_deg_n, time_dn, time_nn, plots='dir')


# nx.draw(G)

# m = MatrixPlot(G)
# m.draw()
# plt.show()
