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

def load_arrays():
	'''function to load text arrays to load 4 directed and 4 undirected netowrks from label string arrays nlc and dlc (both arrays with four strings) 
	to create dynamical network objects which can be used for plotting where magp is the cluster to be used and comp the component n,e,z'''
	
	nlc = [[],[],[],[]]

	dlc = [[],[],[],[]]

	for i in [0,1,2,3]:

		# nlc, dcl array of names of files to read into text arrays to load files

		nlc[i] = f'networks_data/na{i}spd.txt'
		
		dlc[i] = f'networks_data/dna{i}spd.txt'

		print(f'loading Pc{i+2} netowrks')

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

	dict_dn = [[[]]*4]*2

	dict_n = [[[]]*4]*2

	for i in [0,1,2,3]:

		# obtains only UTC times which have edges

		utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in dn[i].edges(data=True) ]

		# removing duplicates in utc_times_dn

		times_dn = list(set(utc_times_dn))

		# ordering via time order

		times_dn = sorted(times_dn, key = lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))

		dict_dn['times'][i].append(times_dn)
	
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

			# utc_times = sorted(list(set(utc_times)), key=lambda tup: tup[0]) 

			dict_dn['times'][i].append(list(set(utc_times)))

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

			dict_n['times'][i].append(list(set(utc_times)))

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


def save_k_vs_t_cluster(dn, nn, cluster):
	
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

	# np.save(f'networks_data/spd_analysis/avg_deg_results_q_{cluster}_dn.npy', dict_dn)
	# np.save(f'networks_data/spd_analysis/avg_deg_results_q_{cluster}_n.npy', dict_n)

	return dict_dn, dict_n
			


def plot_k_vs_t(dict_dn, dict_n, plots, label=False,):

	'''code to plot the average degree of networkx directed and undirected multi-dim network parameters dict_dn, dict_n 
	with label the name of the cluster and to sepcify plots of undirected
	or directed netwrok or both using plots variable using 'all','dir','undir to determin which results to return on the plots'''

	# loading values

	# step_size = np.load('step_size_arr.npy')

	# [val for tup in mlt_edge_list1 for val in tup]
	
	# print('timedn!',np.squeeze(timedn))
	
	# code to give the correct number of subplots or otherwise not plot
	
	md = pd.read_csv('SME_SMR_spd.csv')

	# print(md.head())

	index_date = md['Date_UTC']

	sme = md['SME']

	smr = md['SMR']

	# print('start_time, end_time',index_date.iloc[0],index_date.iloc[-1])

	# sub divisions of time axis

	# n = 7

	# indx_mask = np.linspace(0,len(index_date)-1,n)

	# datel = dateflocal(index_date,indx_mask)

	axsmer = [ datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') for x in index_date ]

	fig, ax = plt.subplots(nrows=6, figsize=(8, 15), facecolor='w', edgecolor='k')
	fig.subplots_adjust(hspace=0.8)

	ax[0].plot(axsmer,sme, color='black')

	# ax[0].set_xticks(indx_mask)

	# ax[0].set_xticklabels(datel)

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

	ax[1].xaxis.set_major_formatter(formatter)

	# ax[0].set_xlim(indx_mask[0],indx_mask[-1])

	# ax[1].set_xlim(indx_mask[0],indx_mask[-1])


	# set title label if using cluster code

	if label:

		fig.suptitle(f'{label} network')

	else:

		fig.suptitle('global network')

	# can do plotting on laptop however first need updated results matrices including node values

	# run code to obtain network parametsremotley without plotting

	# need to figure out way to plot degree on an independantly scaled axis

	for i in [0,1,2,3]:

		print('times dir',i ,dict_dn['times'][i], len(dict_dn['times'][i]))

		print('k dir',i ,dict_dn['avg_k'][i], len(dict_dn['avg_k'][i]))

		if plots == 'all':

			xd = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_dn['times'][i] ]

			xn = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_n['times'][i] ]

			ax[i+2].scatter(xd, dict_dn['avg_k'][i], color='black', s=6)

			ax[i+2].scatter(xn, dict_dn['avg_k'][i], color='black', s=6)

			ax[i+2].plot(xd, avg_k_dn[i], color='r', label=f'Pc{i+2} directed network', lw=2)

			ax[i+2].plot(xn, avg_k_n[i], label= f'Pc{i+2} instantaneously directed network', lw=2)

		elif plots == 'undir':

			xn = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_n['times'][i] ]

			ax[i+2].scatter(xn, dict_n['avg_k'][i], color='black', s=6)

			ax[i+2].plot(xn, dict_n['avg_k'][i], label= f'Pc{i+2} instantaneously directed network', lw=2)

			yax2 = ax[i+2].twinx()

			yax2.plot(xn, dict_n['n_nodes'][i], color='green', linestyle='--' )

		elif plots == 'dir':

			xd = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in dict_dn['times'][i] ]

			print('xd',i, len(dict_dn['times'][i]))

			print('avg_k',i, len(dict_dn['avg_k'][i]))

			ax[i+2].scatter(xd, dict_dn['avg_k'][i], color='orange', s=6)

			ax[i+2].plot(xd, dict_dn['avg_k'][i], color='r', label=f'Pc{i+2} directed network', lw=2)

			yax2 = ax[i+2].twinx()

			yax2.plot(xd, dict_dn['n_nodes'][i], color='orange', linestyle='--' )



		yax2.set_ylabel('num nodes')

		ax[i+2].set_xlabel('time (UTC)')

		ax[i+2].set_ylabel('average degree')

		# print(index_date[0],index_date.iloc[-1],type(datel[0]))

		datemin = np.datetime64(index_date.iloc[0])
		datemax = np.datetime64(index_date.iloc[-1])

		ax[i+2].set_xlim(datemin,datemax)

		formatter = mdates.DateFormatter("%H:%M:%S")

		ax[i+2].xaxis.set_major_formatter(formatter)


		ax[i+2].legend(bbox_to_anchor=(0,1.02,1,0.2),loc="lower left",
		            mode="expand", borderaxespad=0, ncol=2)

		ax[i+2].grid()

		# if statment to save plots correctly 

		if label == False:

			plt.savefig(f'plots/global_net_{plots}.png')

		else:

			plt.savefig(f'plots/{label}_net_{plots}.png')


	# plt.show()


# load arrays!
# dna, na = load_arrays()

for j in ['dir','undir']:

	# # code for global network plotting-----------------

	# avg_deg_dn, avg_deg_n, time_dn, time_nn = save_k_vs_t_global(dna, na) 

	# avg_deg_dn = np.load('networks_data/spd_analysis/avg_deg_global_dn.npy',allow_pickle=True, fix_imports=True)
	# avg_deg_n = np.load('networks_data/spd_analysis/avg_deg_global_n.npy',allow_pickle=True, fix_imports=True)

	# time_dn = np.load('networks_data/spd_analysis/avg_deg_time_global_dn.npy',allow_pickle=True, fix_imports=True)
	# time_nn = np.load('networks_data/spd_analysis/avg_deg_time_global_n.npy',allow_pickle=True, fix_imports=True)

	# plot_k_vs_t(avg_deg_dn, avg_deg_n, time_dn, time_nn, plots=j)

	for i in ['dawn','noon','dusk','midnight']:

		# code for cluster plotting ----------------

		# dict_dn, dict_n = save_k_vs_t_cluster(dna, na, i)

		# only need two files moving forward

		dict_dn = np.load(f'networks_data/spd_analysis/avg_deg_results_q_{i}_dn.npy',allow_pickle=True)[()]
		
		dict_n = np.load(f'networks_data/spd_analysis/avg_deg_results_q_{i}_n.npy',allow_pickle=True)[()]

		plot_k_vs_t(dict_dn, dict_n, label = i, plots=j)



# nx.draw(G)

# m = MatrixPlot(G)
# m.draw()
# plt.show()
