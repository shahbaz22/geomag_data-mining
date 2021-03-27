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

def dateflocal(date,ind): 
	'''fuction to convert datetime utc formal to numpy array to be used for plotting 
	where ind is the total number of time points needed from the date-time series'''

	mask = date.index.isin(ind)

	date = date[mask]

	date = pd.to_datetime(date, format='%Y/%m/%d %H:%M:%S')

	# date = pd.to_datetime(date, format='%Y/%m/%d %H:%M:%S') - timedelta(hours=6.1)

	date = date.astype(str)

	date = date.str.extract(f'(\d\d:\d\d:\d\d)', expand=True)

	date = np.squeeze(date.values)

	return date


def load_arrays():
	'''function to load text arrays to load 4 directed and 4 undirected netowrks from label string arrays nlc and dlc (both arrays with four strings) 
	to create dynamical network objects which can be used for plotting where magp is the cluster to be used and comp the component n,e,z'''
	
	nlc = [[],[],[],[]]

	dlc = [[],[],[],[]]

	for i in [0,1,2,3]:

		# nlc, dcl array of names of files to read into text arrays to load files

		nlc[i] = f'networks_data/na{i}spd.txt'
		
		dlc[i] = f'networks_data/dna{i}spd.txt'

		# netdata1 = open(nlc[i],"rb")

		# netdata2 = open(dlc[i],"rb")

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

	matrix_dn = [[[]]*4]*2

	matrix_n = [[[]]*4]*2

	for i in [0,1,2,3]:

		# obtains the ordered time stamps in the network

		time_stamps_dn = [d['attr_dict']['t_window'] for n1,n2,d in dn[i].edges(data=True)]

		time_stamps_nn = [d['attr_dict']['t_window'] for n1,n2,d in nn[i].edges(data=True)]

		# print('timestamps_dn',time_stamps_dn)

		# print('timestamps_nn',time_stamps_nn)
			# key= lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))
		# need to a

		# remove duplicates
			
		time_stamps_dn = sorted(list(set(time_stamps_dn)))
		
		time_stamps_nn = sorted(list(set(time_stamps_nn)))

		# can also save relevent utc ARRAY

		# loop for slicing consecutive time stamps and calculating degree for directed network

		for j in time_stamps_dn:

			ts_edge_list1 = [ (n1,n2) for n1,n2,d in dn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

			utc_times = [(j, d['attr_dict']['UTC2']) for n1,n2,d in dn[i].edges(data=True) if d['attr_dict']['t_window'] == j]

			# utc_times = sorted(list(set(utc_times)), key=lambda tup: tup[0]) 

			matrix_dn[0][i].append(list(set(utc_times)))

			# count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
			
			# take a unique list of nodes in given time window

			# like double for loop, one loop for acessing tuple value and the next for unpacking it

			ts_node_list1 = list(set([val for tup in ts_edge_list1 for val in tup]))

		
			if len(ts_node_list1) != 0:

				avg_deg = len(ts_edge_list1) / len(ts_node_list1)

				# print(avg_deg, j, 'avg deg, count rc')

				matrix_dn[1][i].append(avg_deg)


		# loop for slicing consecutive time stamps and calculating degree for undirected network
			
		for j in time_stamps_nn:

			ts_edge_list2 = [ (n1,n2) for n1,n2,d in nn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

			utc_times = [(j, d['attr_dict']['UTC2']) for n1,n2,d in nn[i].edges(data=True) if d['attr_dict']['t_window'] == j]

			matrix_n[0][i].append(list(set(utc_times)))

			# count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
			
			# take a unique list of nodes in given time window

			# like double for loop, one loop for acessing tuple value and the next for unpacking it

			ts_node_list2 = list(set([val for tup in ts_edge_list2 for val in tup]))
				
			# code segment for undirnetwork

			if len(ts_node_list2) != 0:

				avg_deg = len(ts_edge_list2) / len(ts_node_list2)

				matrix_n[1][i].append(avg_deg)


	np.save(f'networks_data/spd_analysis/avg_deg_global_dn.npy',matrix_dn)
	np.save(f'networks_data/spd_analysis/avg_deg_global_n.npy',matrix_n)

	# files.download(f'networks_data/spd_analysis/avg_deg_Pc_dn.npy',avg_deg_matrix_dn)
	# files.download(f'networks_data/spd_analysis/avg_deg_Pc_n.npy',avg_deg_matrix_n)

	# files.download(f'networks_data/spd_analysis/avg_deg_time_Pc_dn.npy',times_dn)
	# files.download(f'networks_data/spd_analysis/avg_deg_time_Pc_n.npy',times_nn)
					


def save_k_vs_t_cluster(dn, nn, cluster):
	
	''' function to calculate number of connections (k) for a given MLT range given by window
	 at each time window in for all network arrays in dnn and nn to and returns file of k vs t to later use
	and speed up plotting of related graphs'''

	# magnetic local time 0 is midnight 6 dawn, 12 noon, 18 dusk 

    # dawn: 03:00-09:00 (ie 06:00 +/-3hrs)

    # noon: 09:00:15:00 (12:00 +/- 3hrs)

    # dusk:15:00-21:00

    # midnight (21:00-03:00)

	# arrays used for plotting gaphs use timesvalues in matrix_dn[0][i].app, for <k> matrix_dn[1][i]
	# for nodes .app,matrix_dn[2][i].app in total each matrix has 12 arrays to store values inside with shape [[[], [], [], []], [[], [], [], []], [[], [], [], []]]

	matrix_dn = [[[]]*4]*3

	matrix_n = [[[]]*4]*3

	# dictonary and values used to obtain mlt ranges to filter for

	mltr = {'dawn':[3,9],'noon':[9,15],'dusk':[15,21],'midnight':[21,3]}

	mlt1 = mltr[cluster][0]

	mlt2 = mltr[cluster][1]

	for i in [0,1,2,3]:

		print(f'pc{i+2}, {cluster} cluster {i}')

		# obtains the ordered time stamps in UTC from MLT range

		time_stamps_dn = [d['attr_dict']['t_window'] for n1,n2,d in dn[i].edges(data=True) 
		if mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]


		time_stamps_nn = [d['attr_dict']['t_window'] for n1,n2,d in nn[i].edges(data=True) 
		if mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

		# could append none to empty degree array to use for plotting

		if len(time_stamps_dn) == 0:

			print(f'no values for Pc{i+2} {cluster} directed network array ')

			continue

		elif len(time_stamps_nn) == 0:

			print(f'no values for Pc{i+2} {cluster} directed network array ')

			continue

		# remove duplicates
			
		time_stamps_dn = sorted(list(set(time_stamps_dn)))
		
		time_stamps_nn = sorted(list(set(time_stamps_nn)))

		# print(time_stamps_dn, f'timestamaps dir {i}')

		# can also save relevent utc ARRAY

		# loop for slicing consecutive time stamps and calculating degree for directed network

		for j in time_stamps_dn:

			mlt_edge_list1 = [ (n1,n2) for n1,n2,d in dn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

			utc_times = [(j, d['attr_dict']['UTC2']) for n1,n2,d in dn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

			# utc_check = [(n1, n2 ,j, d['attr_dict']['UTC2'], d['attr_dict']['MLT1'], d['attr_dict']['MLT2']) for n1,n2,d in dn[i].edges(data=True) 
			# if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]
			
			# print(f'checking dir {cluster} clusters pc{i+2}', utc_check)

			# appending time values
			matrix_dn[0][i].append(list(set(utc_times)))

			# count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
			
			# take a unique list of nodes in given time window

			# like double for loop, one loop for acessing tuple value and the next for unpacking it

			mlt_node_list1 = list(set([val for tup in mlt_edge_list1 for val in tup]))
		
			num_nodes = len(mlt_node_list1)

			if num_nodes != 0:

				avg_deg = len(mlt_edge_list1) / num_nodes

				# print(avg_deg, j, 'avg deg, count rc')

				matrix_dn[1][i].append(avg_deg)

				matrix_dn[2][i].append(num_nodes)


		# loop for slicing consecutive time stamps and calculating degree for undirected network
			
		for j in time_stamps_nn:

			mlt_edge_list2 = [ (n1,n2) for n1,n2,d in nn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

			utc_times = [(j, d['attr_dict']['UTC2']) for n1,n2,d in nn[i].edges(data=True) 
			if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]

			# utc_check2 = [(n1, n2 ,j, d['attr_dict']['UTC2'], d['attr_dict']['MLT1'], d['attr_dict']['MLT2']) for n1,n2,d in nn[i].edges(data=True) 
			# if d['attr_dict']['t_window']==j and mlt1 <= d['attr_dict']['MLT1'] < mlt2 and mlt1 <= d['attr_dict']['MLT2'] < mlt2]
			
			# print(f'checking undir {cluster} clusters pc{i+2}', utc_check2)

			matrix_n[0][i].append(list(set(utc_times)))

			mlt_node_list2 = list(set([val for tup in mlt_edge_list2 for val in tup]))
				
			num_nodes2 = len(mlt_node_list2)

			if num_nodes2 != 0:

				avg_deg = len(mlt_edge_list2) / num_nodes2

				matrix_dn[1][i].append(avg_deg)

				matrix_dn[2][i].append(num_nodes2)

	# for saving files more completley in the future can add the network array with comp value to output save file

	np.save(f'networks_data/spd_analysis/avg_deg_results_{cluster}_dn.npy', matrix_dn)
	np.save(f'networks_data/spd_analysis/avg_deg_results_{cluster}_n.npy', matrix_n)
			





def plot_k_vs_t(matrix_dn, matrix_n, plots, label=False,):

	'''code to plot the average degree of networkx directed and undirected multi-dim network parameters matrix_dn, matrix_n 
	with label the name of the cluster and to sepcify plots of undirected
	or directed netwrok or both using plots variable using 'all','dir','undir to determin which results to return on the plots'''

	# loading values

	# step_size = np.load('step_size_arr.npy')

	# [val for tup in mlt_edge_list1 for val in tup]
	
	# print('timedn!',np.squeeze(timedn))
	
	# code to give the correct number of subplots or otherwise not plot
	
	md = pd.read_csv('SME_SMR_spd.csv')

	print(md.head())

	index_date = md['Date_UTC']

	sme = md['SME']

	smr = md['SMR']

	print('start_time, end_time',index_date.iloc[0],index_date.iloc[-1])

	# sub divisions of time axis

	n = 7

	indx_mask = np.linspace(0,len(index_date)-1,n)

	datel = dateflocal(index_date,indx_mask)

	fig, ax = plt.subplots(nrows=6, figsize=(8, 15), facecolor='w', edgecolor='k')
	fig.subplots_adjust(hspace=0.8)

	ax[0].plot(index_date,sme, color='black')

	ax[0].set_xticks(indx_mask)

	ax[0].set_xticklabels(datel)

	ax[0].set_xlabel('time (hh:mm:ss)')

	ax[0].set_ylabel('SME (nT)')

	ax[0].grid()

	ax[1].plot(index_date,smr, color='g')

	ax[1].set_xticks(indx_mask)

	ax[1].set_xticklabels(datel)

	ax[1].set_xlabel('time (hh:mm:ss)')

	ax[1].set_ylabel('SMR (nT)')

	ax[1].grid()

	ax[0].set_xlim(indx_mask[0],indx_mask[-1])

	ax[1].set_xlim(indx_mask[0],indx_mask[-1])


	# set title label if using cluster code

	if label:

		fig.suptitle(f'{label} network')

	else:

		fig.suptitle('global network')

	# can do plotting on laptop however first need updated results matrices including node values

	# run code to obtain network parametsremotley without plotting

	# need to figure out way to plot degree on an independantly scaled axis

	for i in [0,1,2,3]:

		if plots == 'all':

			xd = [ datetime.datetime.strptime(n2, '%Y-%m-%dT%H:%M:%S') for tup in matrix_dn[0][i] for n1,n2 in tup]

			xn = [ datetime.datetime.strptime(n2, '%Y-%m-%dT%H:%M:%S') for tup in matrix_n[0][i] for n1,n2 in tup]

			ax[i+2].scatter(xd, matrix_dn[1][i], color='black', s=6)

			ax[i+2].scatter(xn, matrix_dn[1][i], color='black', s=6)

			ax[i+2].plot(xd, avg_k_dn[i], color='r', label=f'Pc{i+2} directed network', lw=2)

			ax[i+2].plot(xn, avg_k_n[i], label= f'Pc{i+2} instantaneously directed network', lw=2)

		elif plots == 'undir':

			xn = [ datetime.datetime.strptime(n2, '%Y-%m-%dT%H:%M:%S') for tup in timenn[i] for n1,n2 in tup]

			ax[i+2].scatter(xn, avg_k_n[i], color='black', s=6)

			ax[i+2].plot(xn, avg_k_n[i], label= f'Pc{i+2} instantaneously directed network', lw=2)

		elif plots == 'dir':

			xd = [ datetime.datetime.strptime(n2, '%Y-%m-%dT%H:%M:%S') for tup in timedn[i] for n1,n2 in tup]

			ax[i+2].scatter(xd, avg_k_dn[i], color='black', s=6)

			ax[i+2].plot(xd, avg_k_dn[i], color='r', label=f'Pc{i+2} directed network', lw=2)


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


	plt.show()


# load arrays!
dna, na = load_arrays()

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

		save_k_vs_t_cluster(dna, na, i)

		# only need two files moving forward

		# avg_deg_dn = np.load(f'networks_data/spd_analysis/avg_deg_{i}_dn.npy',allow_pickle=True, fix_imports=True)
		# print(avg_deg_dn)
		
		# avg_deg_n = np.load(f'networks_data/spd_analysis/avg_deg_{i}_n.npy',allow_pickle=True, fix_imports=True)
		# print(avg_deg_n)

		# time_dn = np.load(f'networks_data/spd_analysis/avg_deg_time_{i}_dn.npy',allow_pickle=True, fix_imports=True)
		# print(time_dn)
		
		# time_nn = np.load(f'networks_data/spd_analysis/avg_deg_time_{i}_n.npy',allow_pickle=True, fix_imports=True)
		# print(time_nn)

#		plot_k_vs_t(avg_deg_dn, avg_deg_n, time_dn, time_nn, label = i, plots=j)



# nx.draw(G)

# m = MatrixPlot(G)
# m.draw()
# plt.show()
