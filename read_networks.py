import networkx as nx
import matplotlib.pyplot as plt
import dynetx as dn
import itertools
import numpy as np
from dynetx.readwrite import json_graph
import pandas as pd
import time
import re
import itertools
import functools


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


def netvis():

	# Using loaded netowrk arrays to plot netowrk visualisation for each window dir and undir

	for i in [1,2]:

		print(i)

		# gives time for last interaction

		epochs = list(dna[i].stream_interactions())[-1][3]

		for j in [9097//step_size[i]]:

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

def save_k_vs_t(dnn,nn,label):
	
	''' function to calculate number of connections (k) at each time window
	in for all network arrays in dnn and nn to and returns file of k vs t to later use
	and speed up plotting of related graphs'''

	avg_deg_matrix_dn = [[],[],[],[]]

	avg_deg_matrix_n = [[],[],[],[]]

	for i in [0,1,2,3]:

		# obtains the last time stamp in the network

		epochs = list(dnn[i].stream_interactions())[-1][3]

		epochs1st = list(dnn[i].stream_interactions())[0][3] 

		print(epochs,'epochs',epochs1st,'epochs1st')

		# loop for slicing consecutive time stamps and calculating degree

		print(f'Pc{i+2}')

		for j in range(epochs):

			dns = dnn[i].time_slice(j,j+1)

			ns = nn[i].time_slice(j,j+1)

			# number of edges divided by number of nodes for directed graph

			# print(j,'time marker', dns.size(),dns.order(), 'edges','nodes')
		
			# code segment for dirnetwork
			if dns.order()!= 0:

				k = dns.size()/dns.order()

				# print(k, j, 'avg deg, count rc')

				avg_deg_matrix_dn[i].append(k)

			else:

				avg_deg_matrix_dn[i].append(np.nan)

			# code segment for undirnetwork

			if ns.order()!= 0:

				k = ns.size()/ns.order()

				# print(k, j, 'avg deg, count rc')

				avg_deg_matrix_n[i].append(k)

			else:

				avg_deg_matrix_n[i].append(np.nan)


		np.save(f'networks_data/avg_deg_{label}Pc{i+2}_dn.npy',avg_deg_matrix_dn)
		np.save(f'networks_data/avg_deg_{label}Pc{i+2}_n.npy',avg_deg_matrix_n)



def hist_degree(netowrk_arr, dirr):

	'''hist – plotting code for network arrays
	to plot either directred OR undirected (using dirr bool) network histogram subplots 
	at given timesnapshot'''

	fig, ax = plt.subplots(4, figsize=(11, 8), facecolor='w', edgecolor='k')

	# for j in [1000,2786,9097,2.8*3600]:
	for i in [0,1,2]:

		# time j of interest in seconds

		j= 2786

		t = j//step_size[i]

		if netowrk_arr[i].order()!= 0:

			# need to replot for only pc2 and pc3

			histdn = dn.classes.function.degree_histogram(netowrk_arr[i],t)

			print(histdn, f'Pc{i+2}')

			relevant_vals = np.array(histdn)

			# print(relevant_vals[relevant_vals>0])

			# bins in this case are the degree values!!

			bins = np.arange(0,len(histdn),1)

			print(bins,len(bins), 'should be num of stations 54')

			avg = np.sum(bins*relevant_vals)/np.sum(relevant_vals)

			# labels depending on network type

			if dirr == True:
				plabel = f'Pc{i+2} zero lag directed network'
				c = 'blue'
			else:
				plabel = f'Pc{i+2} zero lag undirected network'
				c = 'orange'

			ax[i].bar(bins,relevant_vals,alpha=0.5,edgecolor='black', label=plabel, color=c)

			# ax[i].set_xlabel('degree value')

			ax[i].vlines(avg, 0, np.max(relevant_vals), color='red', ls='--', label=f' mean value {np.round(avg,1)}')

			# ax[2].set_xticks(bins)

			# ax[i].set_ylabel('frequency')
			
			ax[i].legend()


	# plt.savefig(f'kvsf{j}undir.png')
	plt.show()

def plot_k_vs_t_SME_SMR():

	'''using file produced from function save_degree_k plots average degree as compared with SME, SMR and wave power 
	Pc indices'''

	# loading values

	step_size = np.load('step_size_arr.npy')

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

	print(datel)

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

	for i in [0,1,2,3]:

		avg_deg_matrix_dn = np.load(f'networks_data/avg_deg_{label}Pc{i+2}_dn.npy',allow_pickle=True, fix_imports=True)
		avg_deg_matrix_n = np.load(f'networks_data/avg_deg_{label}Pc{i+2}_n.npy',allow_pickle=True, fix_imports=True)


		# # total time in seconds 
		# t = 4*3600

		# x is giving the time axis 

		# FIX! time taken from leading edge of the window (mid is better) also need more data for higher bands Pc4-5

		x = (np.arange(len(avg_deg_matrix_dn[i]))+1)*step_size[i]/3600

		ax[i+2].plot(x, avg_deg_matrix_dn[i], color='r', label=f'Pc{i+2} directed network')

		# ax[i+2].plot(x, avg_deg_matrix_n[i], label= f'Pc{i+2} instantaneous directed network')

		ax[i+2].set_xlabel(f'epoch number (step size = {step_size[i]})')
		ax[i+2].set_xlabel('time (hrs)')

		ax[i+2].set_ylabel('average degree')

		ax[i+2].set_xlim(0,4)

		ax[i+2].legend(bbox_to_anchor=(0,1.02,1,0.2),loc="lower left",
		            mode="expand", borderaxespad=0, ncol=2)

		ax[i+2].grid()

	# plt.savefig('avg_kvst_undir.png')
	plt.show()

def plot_k_vs_t(label):

	'''using file produced from function save_degree_k plots average degree with label the name of the cluster used to
	find file in directory and label plots'''

	# loading values

	step_size = np.load('step_size_arr.npy')

	fig, ax = plt.subplots(nrows=4, figsize=(8, 15), facecolor='w', edgecolor='k')
	fig.subplots_adjust(hspace=0.8)


	for i in [0,1,2,3]:

		avg_deg_matrix_dn = np.load(f'networks_data/avg_deg_{label}Pc{i+2}_dn.npy',allow_pickle=True, fix_imports=True)
		avg_deg_matrix_n = np.load(f'networks_data/avg_deg_{label}Pc{i+2}_n.npy',allow_pickle=True, fix_imports=True)

		# # total time in seconds 
		# t = 4*3600

		# x is giving the time axis 

		# FIX! time taken from leading edge of the window (mid is better) also need more data for higher bands Pc4-5

		x = (np.arange(len(avg_deg_matrix_dn[i]))+1)*step_size[i]/3600

		ax[i].plot(x, avg_deg_matrix_dn[i], color='r', label=f'Pc{i+2} directed network, {label}')

		ax[i].plot(x, avg_deg_matrix_n[i], label= f'Pc{i+2} instantaneously directed network, {label}')

		ax[i].set_xlabel(f'epoch number (step size = {step_size[i]})')
		ax[i].set_xlabel('time (hrs)')

		ax[i].set_ylabel('average degree')

		ax[i].set_xlim(0,4)

		ax[i].legend(bbox_to_anchor=(0,1.02,1,0.2),loc="lower left",
		            mode="expand", borderaxespad=0, ncol=2)

		ax[i].grid()

	# plt.show()

	plt.savefig(f'plots/avg_kvst_e_{label}.png')


def cluster_net():
	'''funtion to take netowrk .txt arrays and longitude to calculate clusters for fixed geographical time (co-moving 
	with earth's roation'''

	# between longitude l1 and l2 create subgraph

	station_data = pd.read_csv('supermag-stations.csv')

	# create a dataset with MLT and station name

	station_data = pd.read_csv('20201111-18-30-supermag.csv')

	# magnetic local time 0 is midnight 6 dawn, 12 noon, 18 dusk 

	index = ['pc2','pc3','pc4','pc5']

	for i in [0,1,2,3]:

		#print(list(dna[i].stream_interactions()))
		nodes = list(dna[i].nodes())

		print(nodes)

		# just filter by longitude? need to addd lists of longitudes

		# here there should be code to filter relevent nodes

		'''when .combination used in previous code each calculation performed once to give direction 
		results in a .txt file line with nodes in either direction, hence .permutation used here'''

		labels = list(itertools.permutations(nodes,2)) 
		# strr = functools.reduce(operator.add, (labels)) 

		# convert tuple list to list 
		labels = [' '.join(x) for x in labels]


		print('labels',labels,type(labels), labels[0],type(labels[0]))

		# if within certain range pick station with this name | is the OR operator

		pattern = re.compile('|'.join(labels))

		# directory before name

		filename= f"networks_data/dirnet_clust_pc{i}.txt"

		# open file for writing ('w+') and create if it doesnt exist

		f = open(filename,'w+')

		# loop for openng global network file then filtering and saving filtered cluster network .txt file

		# if its too slow could use pandas.DataFrame.filter with regex and np.savetxt 

		for j, line in enumerate(open(dnl[i])):
			for match in re.finditer(pattern, line):
				print('found line',line)
				f.write(line)


# nlc =["networks_data/na0spd.txt","networks_data/na1spd.txt","networks_data/na2spd.txt","networks_data/na3spd.txt"]

# dlc = ["networks_data/dna0spd.txt","networks_data/dna1spd.txt","networks_data/dna2spd.txt","networks_data/dna3spd.txt"]

# load step size array of windows used for all networks




def load_arrays(magp, comp, savekvtfile=False):
	'''function to load text arrays to load 4 directed and 4 undirected netowrks from label string arrays nlc and dlc (both arrays with four strings) 
	to create dynamical network objects which can be used for plotting where magp is the cluster to be used and comp the component n,e,z'''
	
	nlc = [[],[],[],[]]

	dlc = [[],[],[],[]]

	for i in [0,1,2,3]:

		# nlc, dcl array of names of files to read into text arrays to load files

		nlc[i] = f'networks_data/na_spd_{magp}_Pc{i}{comp}.txt'
		
		dlc[i] = f'networks_data/dna_spd_{magp}_Pc{i}{comp}.txt'

		netdata1 = open(nlc[i],"rb")

		netdata2 = open(dlc[i],"rb")

		# nlc and dlc text array overwritted to contain networks

		nlc[i] = dn.readwrite.edgelist.read_interactions(netdata1, nodetype=str, timestamptype=int)

		dlc[i] = dn.readwrite.edgelist.read_interactions(netdata2, nodetype=str, timestamptype=int, directed=True)

	if savekvtfile == True:
		save_k_vs_t(dlc, nlc, magp)

	# 	# returns temporal ordered stream of interactions in form (node, node, op, timestamp)
	# 	# print(list(dna[i].stream_interactions()))
	# 	# print(list(na[i].stream_interactions()))

magp = ['dawn','noon','dusk','midnight']

component = 'e'

for i in magp:

	load_arrays(i,'e', savekvtfile = True)

	# # only needs to be run once
	# save_k_vs_t(dlc, nlc, td)

	plot_k_vs_t(i)


# net_read(dnl,nl)











