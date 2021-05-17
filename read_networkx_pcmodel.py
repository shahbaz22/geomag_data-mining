import networkx as nx
import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
import time
import functools
from collections import Counter, OrderedDict
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

def ulf_power( networks_filename, ulf_power_filename, times, band ,sec_times):

	''' function to take the filename of the networks analysis file and time labels then get 
	top 3 ULF power and store 'MLT and MAGLAT' values for each time stamps
	compare with the results of the networks file to ensure the timestamps are the same

	possible to speed up the code if not the same timestamps needed'''

	# first step ensure only the nodes in the network are being used for comparison with ulf wave power
	# by finding all stations with complete data from timeseries data set used to create networks 
	
	
	md = pd.read_csv(networks_filename)

	# print(md.head())

	# print(print(md.columns.tolist()))

	s_labels = md.drop_duplicates( subset=['IAGA'] )['IAGA']

	# checking all stations for incomplete data 4*3600 in the lenght of the data in seconds, keeping only complete data sets

	net_stations = [(lab, len( md[md['IAGA'] == lab]['Date_UTC'] ) / (sec_times)) for lab in s_labels]

	net_stations = [n1 for n1,n2 in net_stations if n2==1]

	# print(net_stations,len(net_stations))

	md1 = pd.read_csv(ulf_power_filename)

	# print(md1.head())

	# print(print(md1.columns.tolist()))

	# returns only those rows with the same labels as those in the timeseries data set used to make the networks

	ulf_powers = md1[ md1['IAGA'].isin(net_stations) ]
	
	ulf_powers_labels = ulf_powers.drop_duplicates( subset=['IAGA'] )['IAGA']

	dict_ulf ={}

	# highest values ulf power arrays to be recorded

	values = ['1st', '2nd', '3rd']

	for i in values:

		dict_ulf[i] = { 'MLT': [], 'MLAT':[], 'pc_power' : [] }



	# code for clusters, set cluster1 and cluster2 as function parameters


	# mltr = {'dawn':[3,9],'noon':[9,15],'dusk':[15,21],'midnight':[21,3]}

	# mlt11 = mltr[cluster1][0]

	# mlt12 = mltr[cluster1][1]

	# mlt21 = mltr[cluster2][0]

	# mlt22 = mltr[cluster2][1]

	# # if cluster 1 is equal to cluster two then single mlt zone is obtained

	# ulf_filt = ulf_powers.loc[ (ulf_powers['MLT'] >= mlt11) & (ulf_powers['MLT'] < mlt22) ]

	# print(ulf_filt)

	# print(ulf_filt['MLT'].max(), ulf_filt['MLT'].min())

	# if stamtnet to chose times for directed or undirected network

	# where band is Pc band - 2

	ulf_filt = ulf_powers

	# print(times,len(times))


	for i in np.squeeze(times):

		# for given time i, return the ulf values, then pick top three ulf power values, for those top three pick the MLT and MLAT

		# print('iiii',i)

		i = i.replace('T', ' ')


		ts = ulf_filt[ulf_filt['Date_UTC'] == i][[ 'Date_UTC', 'IAGA', 'MAGLAT', 'MLT', 'PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']] 
		
		# values are 1st 2nd and 3rd

		for j, lab in enumerate(values):

			# power labels for relevent band
		
			powers = ['PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']

			# number of values in array

			ts = ts.groupby(powers[band]).head()[0:len(values)]

			dict_ulf[lab]['MLT'].append(ts['MLT'].iloc[j])

			dict_ulf[lab]['MLAT'].append(ts['MAGLAT'].iloc[j])

			dict_ulf[lab]['pc_power'].append(ts[powers[band]].iloc[j])


	return dict_ulf

def top3append(reshaped_edge_list, dictt, i):
	
	# fuction to extract top 3 nodes using couting function with the reshaped edgelist then append the MLT, MLAT
	# and edge numbers  for relevent nodes  to dictonary dictt where i is the Pc band array for the values to be appended

	# note reshaped edgelist is a list without deliminator and has not been reduced using set() function

	counts = Counter(list(reshaped_edge_list)).most_common()

	#to add more nodes need an if statment, also two indices used to access tuple value

	nodevalues = [ counts[0][1], counts[1][1] ]

	# node label to obtain station information

	nodelabel1 = counts[0][0]

	nodelabel2 = counts[1][0]

	station1, timestamp1, mlt1, mlat1 = nodelabel1.split('_')

	station2, timestamp2, mlt2, mlat2 = nodelabel2.split('_')
	
	mlta = [float(mlt1), float(mlt2)]

	mlata = [float(mlat1), float(mlat2)]

	for ind1, k in enumerate(['1st','2nd']):

		dictt['avg_k'][k]['MLT'][i].append(mlta[ind1])

		dictt['avg_k'][k]['MLAT'][i].append(mlata[ind1])

		dictt['avg_k'][k]['e_num'][i].append(nodevalues[ind1])

	# loop for code for special case, 3rd node

	if len(counts) >= 3:

		nodevalue3 =  counts[2][1]

		nodelabel3 = counts[2][0]

		station3, timestamp3, mlt3, mlat3 = nodelabel3.split('_')

		dictt['avg_k']['3rd']['MLT'][i].append(float(mlt3))

		dictt['avg_k']['3rd']['MLAT'][i].append(float(mlat3))

		dictt['avg_k']['3rd']['e_num'][i].append(nodevalue3)

	else:

		dictt['avg_k']['3rd']['MLT'][i].append(np.nan)

		dictt['avg_k']['3rd']['MLAT'][i].append(np.nan)

		dictt['avg_k']['3rd']['e_num'][i].append(0)




def load_arrays(Pcbands, comp):
	'''function to load text arrays to load 4 directed and 4 undirected networks from label string arrays nlc and dlc (both arrays with four strings) 
	to create dynamical network objects which can be used for plotting where magp is the cluster to be used and comp the component n,e,z'''
	
	nlc = [[],[],[],[]]

	dlc = [[],[],[],[]]

	for i in Pcbands:

		# nlc, dcl array of names of files to read into text arrays to load files

		nlc[i] = f'networks_data/dna{i}_{comp}_spd.txt'
		
		dlc[i] = f'networks_data/na{i}_{comp}_spd.txt'

		print(f'loading Pc{i+2} {comp} network')

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

# mlt = [(u,v) for u,v,d in dna[2+n].edges(data=True) if d['attr_dict']['t_window']==180 and 3 <= (d['attr_dict']['MLT1'] and d['attr_dict']['MLT2'])<=9]

# print(t_slice,'tslice')

# degrees = sum([val for (node, val) in G.degree()]) / len(G.nodes())

# print(degrees, 'avg_degree', degrees/len(G.nodes()))

def save_k_vs_t_global(dn, nn, save_label, comp, pcbands):
	
	''' function to calculate number of connections (k) globally at each time window
	in for all network arrays in dnn and nn to and returns file of k vs t to later use
	and speed up plotting of related graphs'''

	# arrays used for plotting gaphs see cluster function for details

	dict_dn = {}

	dict_n ={}

	for i in ['times', 'n_nodes']:

		dict_dn[i] = {0:[], 1:[], 2:[], 3:[]}

	for i in ['times', 'n_nodes']:

		dict_n[i] = {0:[], 1:[], 2:[], 3:[]}

	for i in ['avg_k']:

		dict_n[i] = {0:[], 1:[], 2:[], 3:[]}

		for j in ['1st', '2nd', '3rd']:

			dict_n[i].update( { j :{ 'MLT':{0:[], 1:[], 2:[], 3:[]}, 'MLAT':{0:[], 1:[], 2:[], 3:[]}, 'e_num' : {0:[], 1:[], 2:[], 3:[]}} } )

	for i in ['avg_k']:

		dict_dn[i] = {0:[], 1:[], 2:[], 3:[]}

		for j in ['1st', '2nd', '3rd']:

			dict_dn[i].update( { j :{ 'MLT':{0:[], 1:[], 2:[], 3:[]}, 'MLAT':{0:[], 1:[], 2:[], 3:[]}, 'e_num' : {0:[], 1:[], 2:[], 3:[]}} } )


	for i in pcbands:

		# obtains only UTC times which have edges

		utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in dn[i].edges(data=True) ]

		# removing duplicates and sorting in utc_times_dn
		
		times_dn = utc_sort_reduce(utc_times_dn)

		dict_dn['times'][i].extend(times_dn)

		utc_times_n = [ d['attr_dict']['UTC2'] for n1,n2,d in nn[i].edges(data=True) ]

		times_n = utc_sort_reduce(utc_times_n)

		dict_n['times'][i].extend(times_n)

		# obtains the ordered time stamps in the network

		time_stamps_dn = [d['attr_dict']['t_window'] for n1,n2,d in dn[i].edges(data=True)]

		time_stamps_nn = [d['attr_dict']['t_window'] for n1,n2,d in nn[i].edges(data=True)]

		# remove duplicates
			
		time_stamps_dn = sorted(list(set(time_stamps_dn)))
		
		time_stamps_nn = sorted(list(set(time_stamps_nn)))

		# loop for slicing consecutive time stamps and calculating degree for directed network

		print(f'pc{i+2}, dir global {comp}')


		for j in time_stamps_dn:

			ts_edge_list = [ [n1,n2] for n1,n2,d in dn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

			# count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
			
			# take a unique list of nodes in given time window

			# like double for loop, one loop for acessing tuple value and the next for unpacking it

			reshaped_edge_list = np.reshape(ts_edge_list, len(ts_edge_list)*2)			

			nodes = list(set(reshaped_edge_list))
		
			avg_deg = len(ts_edge_list) / len(nodes)

			dict_dn['avg_k'][i].append(avg_deg)

			dict_dn['n_nodes'][i].append(len(nodes))

			top3append(reshaped_edge_list, dict_dn, i)


		# loop for slicing consecutive time stamps and calculating degree for undirected network

		print(f'pc{i+2}, undir global, {comp}')
			
		for j in time_stamps_nn:

			ts_edge_list = [ [n1,n2] for n1,n2,d in nn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

			# count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
			
			# take a unique list of nodes in given time window

			# like double for loop, one loop for acessing tuple value and the next for unpacking it

			reshaped_edge_list = np.reshape(ts_edge_list, len(ts_edge_list)*2)			

			nodes = list(set(reshaped_edge_list))
				
			avg_deg = len(ts_edge_list) / len(nodes)

			dict_n['avg_k'][i].append(avg_deg)

			dict_n['n_nodes'][i].append(len(nodes))	

			top3append(reshaped_edge_list, dict_n, i)


	np.save(f'{save_label}_dn.npy', dict_dn)
	np.save(f'{save_label}_n.npy', dict_n)			


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
	# d is the dictonary template used to store all pcband values (could no for used within for loops)
	# 1st 2nd and 3rd will contain arrays for top 3 most connected nodes at each timestamp

	dict_dn = {}

	dict_n ={}

	# each avg k in or out or total will have it's own top three nodes with connections

	for i in ['times', 'times_in', 'times_out', 'n_nodes','n_nodes_in','n_nodes_out']:

		dict_dn[i] = {0:[], 1:[], 2:[], 3:[]}

	for i in ['times', 'n_nodes']:

		dict_n[i] = {0:[], 1:[], 2:[], 3:[]}

	for i in ['avg_k']:

		dict_n[i] = {0:[], 1:[], 2:[], 3:[]}

		for j in ['1st', '2nd', '3rd']:

			dict_n[i].update( { j :{ 'MLT':{0:[], 1:[], 2:[], 3:[]}, 'MLAT':{0:[], 1:[], 2:[], 3:[]}, 'e_num' : {0:[], 1:[], 2:[], 3:[]}} } )

	for i in ['avg_k', 'avg_k_in', 'avg_k_out']:

		dict_dn[i] = {0:[], 1:[], 2:[], 3:[]}

		for j in ['1st', '2nd', '3rd']:

			dict_dn[i].update( { j :{ 'MLT':{0:[], 1:[], 2:[], 3:[]}, 'MLAT':{0:[], 1:[], 2:[], 3:[]}, 'e_num' : {0:[], 1:[], 2:[], 3:[]}} } )

	# print(dict_dn)
	
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

		# times_in = utc_sort_reduce(utc_times_in)

		# times_out = utc_sort_reduce(utc_times_out)

		dict_dn['times'][i].extend(times_dn)

		# dict_dn['times_in'][i].extend(times_in)

		# dict_dn['times_out'][i].extend(times_out)

		# remove duplicates
			
		time_stamps_dn = sorted(list(set(time_stamps_dn)))
		
		# loop for slicing consecutive time stamps and calculating degree for directed network
		
		print(f'pc{i+2}, dir {cluster1} {cluster2} out of 4 pairs')

		for ind, j in enumerate(time_stamps_dn):

			# MLT2 will always be the MLT value of the second node and MLT1 of the first node

			# first node (MLT1) in the second cluster and the second node (MLT2) in the first cluster B -> A in

			mlt_edge_list_in = [ [ n1, n2 ] for n1,n2,d in dn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22]

			# first node in the first cluster and the second node in the second cluster A -> B in

			mlt_edge_list_out = [ [ n1, n2 ] for n1,n2,d in dn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22]

			# will need to check intercluster code by ensuring each node mlt out of n1 and n2 in different cluster

			# utc_check_in = [(n1, n2 ,j, d['attr_dict']['UTC1'], d['attr_dict']['MLT1'], d['attr_dict']['MLT2']) for n1,n2,d in dn[i].edges(data=True) 
			# if d['attr_dict']['t_window']==j and mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22]
			
			# print('utc_check_out',utc_check_out)

			# reduced dimensionality list to use with counter function most counted node will be the node with the most edges

			all_edge_list = mlt_edge_list_in + mlt_edge_list_out

			reshaped_edge_list_in = np.reshape(mlt_edge_list_in, len(mlt_edge_list_in)*2)

			reshaped_edge_list_out = np.reshape(mlt_edge_list_out, len(mlt_edge_list_out)*2)

			reshaped_edge_list = np.reshape(all_edge_list, len(all_edge_list)*2)			
			
			# reduced dim of edgelist and remove duplicated, set and list also add commas inbetween spaced elments

			mlt_node_list_in = list(set(reshaped_edge_list_in))

			mlt_node_list_out = list(set(reshaped_edge_list_out))

			# prevent double counting of nodes from both lists
			# need to do addition before list-set

			all_nodes_list = mlt_node_list_in + mlt_node_list_out

			all_nodes = list(set(all_nodes_list))

			avg_deg = ( len(all_edge_list) ) / len(all_nodes)

			# counts for edge and MLT values

			# print(set(reshaped_edge_list))

			dict_dn['avg_k'][i].append(avg_deg)

			dict_dn['n_nodes'][i].append(len(all_nodes))

			# obtaining top nodes

			# obtaining top nodes and appending values

			top3append(reshaped_edge_list, dict_dn, i)



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

		for ind, j in enumerate(time_stamps_nn):

			mlt_edge_list2 = [ [ n1,n2 ] for n1,n2,d in nn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt11 <= d['attr_dict']['MLT1'] < mlt12 and mlt21 <= d['attr_dict']['MLT2'] < mlt22 
			or d['attr_dict']['t_window']==j and mlt11 <= d['attr_dict']['MLT2'] < mlt12 and mlt21 <= d['attr_dict']['MLT1'] < mlt22 ]

			# list set addes the commas between elements

			reshaped_edge_list = np.reshape(mlt_edge_list2, len(mlt_edge_list2)*2)

			mlt_node_list2 = list(set(reshaped_edge_list))
				
			num_nodes2 = len(mlt_node_list2)

			avg_deg = len(mlt_edge_list2) / num_nodes2

			dict_n['avg_k'][i].append(avg_deg)

			dict_n['n_nodes'][i].append(num_nodes2)

			top3append(reshaped_edge_list, dict_n, i)


		# for k in ['1st','2nd','3rd']:

		# 	# print('e_num dir',k,len(dict_dn['avg_k'][k]['e_num'][i]), dict_dn['avg_k'][k]['e_num'][i])

		# 	print('e_num undir',k,len(dict_n['avg_k'][k]['e_num'][i]), dict_n['avg_k'][k]['e_num'][i])

		# # print('num time stamps undir',len(time_stamps_nn))

		# # print('num time stamps dir',len(time_stamps_dn))

		# print('time undir', dict_n['times'][i], len(dict_n['times'][i]))

		# # print('time dir', dict_dn['times'][i], len(dict_dn['times'][i]))

		# # print('avg k dir',len(dict_dn['avg_k'][i]), dict_dn['avg_k'][i])

		# print('avg k undir',len(dict_n['avg_k'][i]), dict_n['avg_k'][i])


	# check undir network


	# for saving files more completley in the future can add the network array with comp value to output save file

	np.save(f'{save_label}_dn.npy', dict_dn)
	np.save(f'{save_label}_n.npy', dict_n)

	return dict_dn, dict_n
				


def plot_k_vs_t(dict_dn, dict_n, Pcbands, comp, dykey1, dykey2, dtimeskey1, uykey1, uykey2, utimeskey1, plots, label, save_plot = False):

	'''code to plot the average degree of networkx directed and undirected multi-dim network parameters dict_dn, dict_n 
	with label the name of single cluster or two cluster (for connection between clusters)

	# dkeys to be used in dictonary dn and ukeys to be used in undir dict n

	dictonary keys to plot relevent values {ykey1} as avg_deg_in andf {ykey2} as nodes_in with relevent time keys
	for both directed and undirected networks denoted with d... or u...
	
	plots can be either dir, undir or both to specify number of network plots on a sungle subplot

	relevents plots: SME/SMR, Bz/Be/Bn, Dynpressure, Pc power top 3, networks plots, MLT, MLAT

	'''
	# check which Pc bands are not empty

	# test_ulf('networks_data/spd_analysis/avg_deg_noon_e', dirr=True)


	non_empty_pcbands = []
	
	for i in Pcbands:

		if len(dict_dn[f'{dykey1}'][i]) != 0:

			print(f' {label} non-empty plot Pc',i+2)

			non_empty_pcbands.append(i)

	# If stamtent to prevent plotting empty network array

	if len(non_empty_pcbands)==0:

		if plots == 'dir':

			print(f'no network values for {label}, {dykey1}, {comp}')

		else:

			print(f'no network values for {label}, {uykey1}, {comp}')


	else:

		print('relevent pc bands, Pc',  [x+2 for x in non_empty_pcbands] )

		ref_plots = 6

		num_plots = ref_plots + len(non_empty_pcbands)

		fig, ax = plt.subplots(nrows = num_plots, figsize=(9, 17), facecolor='w', edgecolor='k')
		fig.subplots_adjust(hspace=0.9)

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

		ax[0].plot(axsmer,sme, color='black', label = 'SME')

		# ax[0].set_xlabel('time (hh:mm:ss)')

		ax[0].set_ylabel('SME (nT)')

		ax[0].grid()

		datemin = np.datetime64(index_date.iloc[0])
		datemax = np.datetime64(index_date.iloc[-1])

		ax[0].set_xlim(datemin,datemax)

		formatter = mdates.DateFormatter("%H:%M:%S")

		ax[0].xaxis.set_major_formatter(formatter)

		ax2 = ax[0].twinx()

		ax2.plot(axsmer,smr, color='g', label= 'SMR')

		ax2.xaxis.set_major_formatter(formatter)

		ax[0].legend()

		ax2.legend()

		# ax2.set_xticks(indx_mask)

		# ax2.set_xticklabels(datel)

		ax2.set_ylabel('SMR (nT)')

		ax2.grid()

		datemin = np.datetime64(index_date.iloc[0])
		datemax = np.datetime64(index_date.iloc[-1])

		ax2.set_xlim(datemin,datemax)

		# 1st 2nd and pc_power keys

		# need an array of non empty pc values

		# also need some kind of counting function

		dyn_pressure = md2['PDYN']

		# GSM field values

		bz, bx, by = md2['GSM_Bz'], md2['GSM_Bx'], md2['GSM_By']

		for i in [0,1,2,3]:

			y_vals = [bx, by, bz, dyn_pressure]

			y_labs = ['GSM_Bz','GSM_Bx','GSM_By','Dynamic pressure']

			c = ['red','black','green','orange']

			if i == 3:

				ax[2].set_ylabel('nPa')

				ax[2].set_ylim(-10,25)

				# get rid of peaks

				y = np.where( y_vals[i] == np.max(y_vals[i]) , np.nan , y_vals[i])

				ax[2].plot(axsmer2, y, label= y_labs[i], color = c[i])

				n = 1

			else:

				ax[1].plot(axsmer2, y_vals[i], label=y_labs[i], color = c[i])

				n = 0

			ax[1+n].xaxis.set_major_formatter(formatter)

			ax[1+n].set_xlim(datemin,datemax)

			ax[1+n].legend()

			ax[1+n].grid()

			ax[1].set_ylabel('nT')

			ax[1].set_ylim(-30,30)


		# set title for relevent label for all plots ----------------------------------------------------

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


		# -------------------------------------------------------------

		if plots == 'dirr':

			print('dir graphs:', plots)

			dictt = dict_dn[f'{dykey1}']

			tk1 = dict_dn[f'{dtimeskey1}']

		else:

			print('dir graphs:', plots)

			dictt = dict_n[f'{uykey1}']

			tk1 = dict_n[f'{utimeskey1}']

		# -------------------------------------------------------------


		for ind, num in enumerate(non_empty_pcbands):

			n = num_plots - 4

			xd = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in np.squeeze( tk1[num] ) ]

			# print(plots, len(xd), xd, len(dictt[f'{uykey1}'][num]), dictt[f'{utimeskey1}'][num])
	
			ax[ind+n].scatter(xd, dictt[num], label=f'Pc{num+2}  average degree', color='black', s=7)
			
			axs2 = ax[ind+n].twinx()

			# second yaxis-----------

			print('times',len(xd),xd)

			for order in ['1st','2nd','3rd']:

				print('connections',len(dictt[order]['e_num'][num]),dictt[order]['e_num'][num])

				axs2.scatter(xd, dictt[order]['e_num'][num], s=7, label= f'Pc{num+2} {order} most # connections')

			# y2range = range(int(np.min(dictt['3rd']['e_num'][num])), int(np.max(dictt['1st']['e_num'][num])),2)

			# print(np.min(dict_dn[f'{dykey1}']['2nd']['e_num'][num]), np.max(dict_dn[f'{uykey1}']['1st']['e_num'][num]), 'ranges min, max', i)

			# print(dict_dn[f'{dykey1}']['1st']['e_num'][num], '1st')
			
			# print(dict_dn[f'{dykey1}']['2nd']['e_num'][num], '2nd')

			axs2.set_ylabel('connections')

			# yax2.legend(bbox_to_anchor=(0.8,1.02,0.2,0.2), loc="lower right", borderaxespad=0, ncol=1)

			axs2.legend(bbox_to_anchor=(0,1.02,1,0.2),loc="lower left",
			            mode="expand", borderaxespad=0, ncol=4)

			axs2.grid()

			# axs2.set_yticks(y2range)

			# ax[ind+n].set_xlabel('time (UTC)')
     
			ax[ind+n].set_ylabel('avg. degree')

			# print(index_date[0],index_date.iloc[-1],type(datel[0]))

			datemin = np.datetime64(index_date.iloc[0])
			datemax = np.datetime64(index_date.iloc[-1])

			ax[ind+n].set_xlim(datemin,datemax)

			formatter = mdates.DateFormatter("%H:%M:%S")

			ax[ind+n].xaxis.set_major_formatter(formatter)

			ax[ind+n].legend()

			ax[ind+n]

			ax[ind+n].grid()

		# Loop for ULF power and MLT, MLAT plots
		
		point_style = ['.','^','>','<']

		# direction, c1, c2 = label.split('_')

		for ind, num in enumerate([num_plots-2, num_plots-1, 3]):

			# axes index stats from 0 to n-1

			for indpc, k in enumerate(non_empty_pcbands):

				# if statments to assign correct varablie MLT or MLAT for plot number in loop

				if ind == 0:

					y1 = lambda order: dictt[order]['MLT'][k]

					ylab = 'MLT'

				elif ind ==1:

					y1 = lambda order: dictt[order]['MLAT'][k]

					ylab = 'MLAT'

					ax[num].set_xlabel('time (UTC)')

				else:

					ylab = 'log(nT^2)'


				print('calculating coresponding highest ulf powers for Pc band', {k+2}, ylab)

				ulf = ulf_power('20201111-18-30-supermag.csv', 'networks_data/spd_ulfpower.csv', times = tk1[k], 
				band = k, sec_times = 4*3600)

				print('ulf power completed')
					
				xn = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in np.squeeze(tk1[k]) ]

				# if statments to perform correct operation for each plot

				if ind <=1:

					for m in ['1st','2nd','3rd']:

						print(label, m, 'xions', ylab, min(y1(m)), max(y1(m)))

						print(label, m, 'ulf', ylab, min(ulf[m][ylab]), max(ulf[m][ylab]))

						ax[num].scatter(xn, y1(m), s=7, label= f'Pc{k+2}, {m}, x', marker = point_style[indpc])

						ax[num].scatter(xn, ulf[m][ylab] , s=7, label= f'Pc{k+2}, {m}, power', marker = point_style[indpc])



				else:
					
					for m in ['1st','2nd','3rd']:
						
						ax[num].scatter(xn, ulf[m]['pc_power'] , s=7, label= f'Pc{k+2}, {m}, power', marker = point_style[indpc])


			# print('minmax',maxmin_vals)

			# print(min(maxmin_vals),max(maxmin_vals))

			ax[num].set_ylabel(ylab)

			formatter = mdates.DateFormatter("%H:%M:%S")

			ax[num].xaxis.set_major_formatter(formatter)


			ax[num].legend(bbox_to_anchor=(0,1.02,1,0.2),loc="lower left",
			            mode="expand", borderaxespad=0, ncol=6)

			ax[num].grid()	

			datemin = np.datetime64(index_date.iloc[0])
			datemax = np.datetime64(index_date.iloc[-1])

			ax[num].set_xlim(datemin,datemax)
					

			if save_plot:

				# saves the plot

				plt.savefig(f'plots/{dykey1}_{label}_net_{comp}3.png')


		plt.show()




# # code for saving cluster and intercluster analysis file

# include commpand for plot

# need to check saveplot labels and commands, also

def save_plot_all_clusters(save_k, comp , Pcband, option, plot, save_plot, nettype):

	# option to chose to save all (option variable) intercluster or cluster analysis files or plot all analysis files for relevent dictonary keys
	# option given within the function using save or plot to true or false
	# with component (comp) for magnetic field 'n', 'e', 'z'
	# nettype can be either 'dir', 'undir' or 'all given as an array'

	# variable for chooseing pc-bands to compute 

	mltz = ['dawn','noon','dusk','midnight']


	if save_k:

		# # load network arrays!

		dna, na = load_arrays(Pcband, comp)

	# loop over network arrays

	if option == 'global':

		file_label = f'networks_data/spd_analysis/avg_deg_global_{comp}'
		
		if save_k:

			if option == 'global':

				save_k_vs_t_global(dna, na, file_label, comp, Pcband)

		elif plot:

			dict_dn = np.load(f'{file_label}_dn.npy',allow_pickle=True)[()]

			dict_n = np.load(f'{file_label}_n.npy',allow_pickle=True)[()]

			for j in nettype:

				plot_label = f'global_{comp}_{j}'

				plot_k_vs_t(dict_dn, dict_n, Pcband, comp, 'avg_k', 'n_nodes','times', 'avg_k', 
				'n_nodes','times', label = plot_label, plots= j, save_plot=save_plot)
	else:

		# loop for calculating all clusters i.e 4 MLT cluters or 4 MLT cluster pairs

		for i, val in enumerate(mltz):

			# labels for incluster savefiles and plots

			if option == 'global':

				file_label = f'networks_data/spd_analysis/avg_deg_global_{comp}'

			elif option == 'intercluster':

				# adjacent mlt labels for analysis and labeling

				adjclusts = inter_cluster_pairs(mltz)

				# file label with file path, file 
				file_label = f'networks_data/spd_analysis/avg_deg_{adjclusts[i][0]}_{adjclusts[i][1]}_{comp}'

			elif option == 'cluster':

				file_label = f'networks_data/spd_analysis/avg_deg_{val}_{comp}'


			# run and save dictonaries for plotting by filtering dynamical network array file

			if save_k:

				if option == 'cluster':

					save_k_vs_t_cluster(dna, na, val, file_label)

				elif option == 'intercluster':

					# pc band option allows us to select a specific band for analysis and development

					save_k_vs_t_intercluster(dna, na, Pcband, adjclusts[i][0], adjclusts[i][1], file_label)

			# plot the analysis files if specified

			if plot:

				if option == 'cluster':

					dict_dn = np.load(f'{file_label}_dn.npy',allow_pickle=True)[()]

					dict_n = np.load(f'{file_label}_n.npy',allow_pickle=True)[()]

					# add a function here to clculate the ulf powers for the bands

					# will need to integrate it with the plotting function but should be easy once I have the arrays

					for j in nettype:

						plot_label = '_'.join( [j, mltz[i]] )

						plot_k_vs_t(dict_dn, dict_n, Pcband, comp, 'avg_k', 'n_nodes', 'times', 'avg_k', 
						'n_nodes','times',  label = plot_label, plots= j, save_plot=save_plot)

				
				elif option == 'intercluster':

					dict_dn = np.load(f'{file_label}_dn.npy', allow_pickle=True)[()]

					dict_n = np.load(f'{file_label}_n.npy', allow_pickle=True)[()]

					# labels for dictionary keys

					# indexer in plotting function chooses Pc band and label a joining of strings from both clusters

					# def plot_k_vs_t(dict_dn, dict_n, Pcbands, comp, dykey1, dykey2, dtimeskey1, 
					# dtimeskey2, uykey1, uykey2, utimeskey1, utimeskey2, plots, label, save_plot = False):


					for j in nettype:

						plot_label = '_'.join( [ j, adjclusts[i][0], adjclusts[i][1] ] )

						plot_k_vs_t(dict_dn, dict_n, Pcband, comp, 'avg_k', 'n_nodes', 'times', 'avg_k', 
						'n_nodes','times', label = plot_label, plots= j, save_plot=save_plot)




# code for reading and plotting network analysis files 

# save_plot_all_clusters(save_k, comp, ,plot, option, j, save_plot = False):

# after plot all other options are related to plottin

# for i in ['e','z','n']:

#     save_plot_all_clusters(save_k = True, comp = i, Pcband = [0,1,2,3], option = 'global', plot=False,  save_plot = False, nettype = ['dir','undir'])

#     print('done', i)


# save_plot_all_clusters(save_k = True, comp = 'e', Pcband = [1], option = 'global', plot=False,  save_plot = False, nettype = ['dir','undir'])


save_plot_all_clusters(save_k = False, comp = 'e', Pcband = [0,1,2,3], option = 'global', plot=True,  save_plot = False	
	, nettype = ['undir','dir'])











# nx.draw(G)

# m = MatrixPlot(G)
# m.draw()
# plt.show()

# label for the cluster
# f'networks_data/spd_analysis/avg_deg_{mltz[i]}_{comp}'


# test_ulf('networks_data/spd_analysis/avg_deg_noon_e', dirr=True)






