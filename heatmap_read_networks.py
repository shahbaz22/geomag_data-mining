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
import matplotlib.ticker as mticker
import pickle
import multiprocessing as mp
from multiprocessing.sharedctypes import Value, Array
from multiprocessing import Manager
import sys



def utc_sort_reduce(utc_times, deliminator='T'):

	# takes an unordered list of utc times and removes repeated values and then orders based on tme value

	times_dn = list(set(utc_times))

	times_dn = sorted(times_dn, key = lambda x: datetime.datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S'))

	times_dn = [ datetime.datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S') for x in np.squeeze(times_dn) ]

	return times_dn

def inter_cluster_pairs(mlta):
	
	# returns all adjecent pairs from a list, used for interscluster code
	
	mltar = np.roll(mlta,-1)

	adjclusts = list(zip(mlta, mltar))

	return adjclusts

def ulf_power(ulf_power_filename, times, band, clusters=False):

	''' function to take the ulf powder dataset and time labels from specific network band and calculate corresponding 
	top 3 ULF power and store 'MLT and MAGLAT' values for each time stamps and then append that to analysis dictonary'''

	# cluster is an an array containing two mlt cluster zones

	ulf_powers = pd.read_csv(ulf_power_filename)

	ulf_powers['Date_UTC'] = pd.to_datetime(ulf_powers['Date_UTC'])


	# print(md1.head())

	# print(print(md1.columns.tolist()))

	# returns only those rows with the same labels as those in the timeseries data set used to make the networks
	
	dict_ulf ={}

	# highest values ulf power arrays to be recorded

	values = ['1st', '2nd', '3rd']

	for i in values:

		dict_ulf[i] = { 'MLT': [], 'MLAT':[], 'pc_power' : [] }

	if clusters:

		mltr = {'dawn':[3,9],'noon':[9,15],'dusk':[15,21],'midnight':[21,3]}

		mlt11 = mltr[clusters[0]][0]

		mlt12 = mltr[clusters[0]][1]

		mlt21 = mltr[clusters[1]][0]

		mlt22 = mltr[clusters[1]][1]

		# if cluster 1 is equal to cluster two then single mlt zone is obtained

		ulf_filt = ulf_powers.loc[ (ulf_powers['MLT'] >= mlt11) & (ulf_powers['MLT'] < mlt22) ]

	else:

		ulf_filt = ulf_powers

	# print(times,len(times))


	for i in np.squeeze(times):

		# for given time i, return the ulf values, then pick top three ulf power values, for those top three pick the MLT and MLAT

		# .replace used to unify string format for times

		# i = i.replace('T', ' ')

		ts = ulf_filt[ulf_filt['Date_UTC'] == i][[ 'Date_UTC', 'IAGA', 'MAGLAT', 'MLT', 'PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']] 
		
		# values are 1st 2nd and 3rd

		powers = ['PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']

		
		for j, lab in enumerate(values):

			# power labels for relevent band		

			# oganise by values largest to smallest

			# if rows less than 3 then groupby doesn't work properly

			ts = ts.groupby(powers[band]).head(n=4)

			dict_ulf[lab]['MLT'].append(ts['MLT'].iloc[j])

			dict_ulf[lab]['MLAT'].append(ts['MAGLAT'].iloc[j])

			dict_ulf[lab]['pc_power'].append(ts[powers[band]].iloc[j])


	return dict_ulf

def ulf_power(ulf_power_filename, times, band, clusters=False):

	''' function to take the ulf powder dataset and time labels from specific network band and calculate corresponding 
	top 3 ULF power and store 'MLT and MAGLAT' values for each time stamps and then append that to analysis dictonary'''

	# cluster is an an array containing two mlt cluster zones

	ulf_powers = pd.read_csv(ulf_power_filename)

	ulf_powers['Date_UTC'] = pd.to_datetime(ulf_powers['Date_UTC'])


	# print(md1.head())

	# print(print(md1.columns.tolist()))

	# returns only those rows with the same labels as those in the timeseries data set used to make the networks
	
	dict_ulf ={}

	# highest values ulf power arrays to be recorded

	values = ['1st', '2nd', '3rd']

	for i in values:

		dict_ulf[i] = { 'MLT': [], 'MLAT':[], 'pc_power' : [] }

	if clusters:

		mltr = {'dawn':[3,9],'noon':[9,15],'dusk':[15,21],'midnight':[21,3]}

		mlt11 = mltr[clusters[0]][0]

		mlt12 = mltr[clusters[0]][1]

		mlt21 = mltr[clusters[1]][0]

		mlt22 = mltr[clusters[1]][1]

		# if cluster 1 is equal to cluster two then single mlt zone is obtained

		ulf_filt = ulf_powers.loc[ (ulf_powers['MLT'] >= mlt11) & (ulf_powers['MLT'] < mlt22) ]

	else:

		ulf_filt = ulf_powers

	# print(times,len(times))


	for i in np.squeeze(times):

		# print(i)

		# for given time i, return the ulf values, then pick top three ulf power values, for those top three pick the MLT and MLAT

		# .replace used to unify string format for times

		# i = i.replace('T', ' ')

		ts = ulf_filt[ulf_filt['Date_UTC'] == i][[ 'Date_UTC', 'IAGA', 'MAGLAT', 'MLT', 'PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']] 
		
		# values are 1st 2nd and 3rd

		powers = ['PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']

		
		for j, lab in enumerate(values):

			# power labels for relevent band		

			# oganise by values largest to smallest

			# if rows less than 3 then groupby doesn't work properly

			ts = ts.groupby(powers[band]).head(n=4)

			dict_ulf[lab]['MLT'].append(ts['MLT'].iloc[j])

			dict_ulf[lab]['MLAT'].append(ts['MAGLAT'].iloc[j])

			dict_ulf[lab]['pc_power'].append(ts[powers[band]].iloc[j])


	return dict_ulf


def load_arrays(Pcbands, comp, handel):
	'''function to load text arrays to load 4 directed and 4 undirected networks from label string arrays nlc and dlc (both arrays with four strings) 
	to create dynamical network objects which can be used for plotting where magp is the cluster to be used and comp the component n,e,z'''
	
	nlc = [[],[],[],[]]

	dlc = [[],[],[],[]]

	for i in Pcbands:

		nlc[i] = f'networks_data/na{i}_{handel}.txt'
		
		dlc[i] = f'networks_data/dna{i}_{handel}.txt'

		print(f'loading Pc{i+2} {comp} network dir and undir')

		nlc[i] = nx.read_edgelist(nlc[i])

		dlc[i] = nx.read_edgelist(dlc[i],create_using=nx.DiGraph)

	return dlc ,nlc

def append_degree_dist(net, dictt, time, delim):

	'function to use in heatmap_analysis_parallel_file function for parallel computing, dictt will be shared multiprocessing dict'
			
	edgelist = [ [n1,n2] for n1,n2,d in net[i].edges(data=True) if d['attr_dict']['UTC1'] == time ]

	reshaped_edgelist = np.reshape(edgelist, len(edgelist)*2)

	nodes = list(set(reshaped_edgelist))

	G = net[i].subgraph(nodes)
	
	# counts how many time a node name repeated in edgelist, i.e how many connections it has

	degree_list = [G.degree(n) for n in G.nodes()]

	# print(degrees)

	# counting frequency of degree values giving list of tuples with connections and frequency i.e [freq, degree]Ź
	freq = Counter(degree_list).most_common()

	degree, freq = zip(*freq)

	# print(counts)

	# print(freq)

	dtime = datetime.datetime.strptime(t, f'%Y-%m-%d{delim}%H:%M:%S')

	d = {dtime : {'degree': list(degree) ,'freq': list(freq)} } 

	dictt.update(d)




def heatmap_analysis_parallel_file(net, pcbands, comp, dirr, label, delim, ulf_filename):
	'''function to create degree distrobutions for a given number of plotted overlayed for 
	different times and bands, labeled and colour caded to their respective times'''

	# need to find way of picking relevent UTC times as close to each other as possible

	# take ordered UTC array data from of all bands and find closest value for both for given value of interest.

	# will get one heatmap for each pc band

	for indpc, i in enumerate(pcbands):

		utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in net[i].edges(data=True) ]

		# utc_times_dn =utc_times_dn

		# need to write below operation as a single function for a single time

		# initializing dict container dict

		lol
		manager = mp.Manager()

		shared_dict = manager.dict()

		pool = mp.Pool(mp.cpu_count())

		for ind, t in enumerate(utc_times_dn):

			pool.apply_async(append_degree_dist, args=(net, shared_dict, t, delim))

			sys.stderr.write('\rdone {0:%}'.format(ind/len(utc_times_dn)))



		pool.close()
		pool.join()


		df = pd.DataFrame.from_dict(shared_dict, orient='index')

		df = df.explode('degree')

		df = df.explode('freq')

		df.index.names = ['UTC_time']

		df = df.reset_index()

		df = df.groupby(['UTC_time','degree']).freq.first().unstack()

		pc_dict[f'Pc{indpc+2}_heatmap'] = df.T

		print(f'calculating pc {indpc + 2} power')

		pc_dict[f'pc{indpc+2}_power'] = ulf_power(ulf_filename, times = utc_times_dn, band = indpc)
		# add dir to filename
	
	a_file = open(f"heatmap_data_{dirr}_{label}.pkl", "wb")

	pickle.dump(pc_dict, a_file)

	# prevent file using up operating memeory
	a_file.close()

	# to read

	# a_file = open(f"heatmap_data_{dirr}.pkl", "rb")
	# output = pickle.load(a_file)
	# print(output)
	# a_file.close()


def heatmap_analysis_file(net, pcbands, comp, dirr, label, delim, ulf_filename):
	'''function to create degree distrobutions for a given number of plotted overlayed for 
	different times and bands, labeled and colour caded to their respective times'''

	# initializing dict

	pc_dict = {}


	# will get one heatmap for each pc band

	for indpc, i in enumerate(pcbands):

		container_dict = {}

		utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in net[i].edges(data=True) ]


		for ind, t in enumerate(utc_times_dn):

			print(ind, 'out of ', len(utc_times_dn))

			edgelist = [ [n1,n2] for n1,n2,d in net[i].edges(data=True) if d['attr_dict']['UTC1'] == t]

			reshaped_edgelist = np.reshape(edgelist, len(edgelist)*2)

			nodes = list(set(reshaped_edgelist))

			G = net[i].subgraph(nodes)
			# counts how many time a node name repeated in edgelist, i.e how many connections it has

			degree_list = [G.degree(n) for n in G.nodes()]

			# print(degrees)


			# counting frequency of degree values giving list of tuples with connections and frequency i.e [freq, degree]Ź
			freq = Counter(degree_list).most_common()

			degree, freq = zip(*freq)

			# print(counts)

			# print(freq)

			dtime = datetime.datetime.strptime(t, f'%Y-%m-%d{delim}%H:%M:%S')

			d = {dtime : {'degree': list(degree) ,'freq': list(freq)} } 

			container_dict.update(d)

		df = pd.DataFrame.from_dict(container_dict, orient='index')

		df = df.explode('degree')

		df = df.explode('freq')

		df.index.names = ['UTC_time']

		df = df.reset_index()

		df = df.groupby(['UTC_time','degree']).freq.first().unstack()

		pc_dict[f'Pc{indpc+2}_heatmap'] = df.T

		print(f'calculating pc {indpc + 2} power')

		pc_dict[f'pc{indpc+2}_power'] = ulf_power(ulf_filename, times = utc_times_dn, band = indpc)

		# add dir to filename
	
	a_file = open(f"heatmap_data_{dirr}_{label}.pkl", "wb")

	pickle.dump(pc_dict, a_file)

	# prevent file using up operating memeory
	a_file.close()

	# to read

	# a_file = open(f"heatmap_data_{dirr}.pkl", "rb")
	# output = pickle.load(a_file)
	# print(output)
	# a_file.close()


			

def plot_heatmap(dict_dn, dict_n, Pcbands, comp, dykey1, dykey2, dtimeskey1, uykey1, uykey2, utimeskey1, plots, label, save_plot = False):

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

		# non network plots include SME/SMR, B feilds, SP, all pc ULF power, all pc MLT, all pc MLAT (in this case more ulf band powers)

		ref_plots = 6

		num_plots = ref_plots + 2 #len(non_empty_pcbands)

		fig = plt.figure(figsize=(20, 17)) #__________________________________________________________________________________________________

		gs = fig.add_gridspec(num_plots, hspace=0)
		ax = gs.subplots(sharex=True, sharey=False)

		indices = pd.read_csv('networks_data/170313_indices.csv')

		print(indices.head())

		print(indices.columns.tolist())

		index_date = indices['Date_UTC']

		sme = indices['SME']

		axsmer = [ datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') for x in index_date ]

		ax[0].plot(axsmer,sme, color='black', label = 'SME')

		# ax[0].set_xlabel('time (hh:mm:ss)')

		ax[0].set_ylabel('SME (nT)')

		ax[0].grid()

		datemin = np.datetime64(index_date.iloc[0])
		datemax = np.datetime64(index_date.iloc[-1])

		# print('datemin',datemin)

		# print('datemax',datemax)

		ax[0].set_xlim(datemin,datemax)

		formatter = mdates.DateFormatter("%H:%M:%S")

		ax[0].xaxis.set_major_formatter(formatter)

		ax2 = ax[0].twinx()

		# code for plotting SMR sections for entire SMR yearly dataset

		indices['Date_UTC'] = pd.to_datetime(indices['Date_UTC']) 

		for i in ['SMR00','SMR06','SMR12','SMR18' ]:

			ax2.plot(axsmer,indices[(indices['Date_UTC'] >= datemin) & (indices['Date_UTC'] <= datemax)][i] , label= i)

		ax2.xaxis.set_major_formatter(formatter)

		ax[0].legend()

		ax2.legend()

		# ax2.set_xticks(indx_mask)

		# ax2.set_xticklabels(datel)

		ax2.set_ylabel('SMR (nT)')

		ax2.grid()

		datemin = np.datetime64(index_date.iloc[0])
		datemax = datetime.datetime.strptime('2013-3-17 12:00:00', '%Y-%m-%d %H:%M:%S')

		ax2.set_xlim(datemin,datemax)

		# 1st 2nd and pc_power keys

		# need an array of non empty pc values

		# also need some kind of counting function

		dyn_pressure = indices['PDYN']

		# GSM field values

		bz, bx, by = indices['GSM_Bz'], indices['GSM_Bx'], indices['GSM_By']

		for i in [0,1,2,3]:

			y_vals = [bx, by, bz, dyn_pressure]

			y_labs = ['GSM_Bz','GSM_Bx','GSM_By','Dynamic pressure']

			c = ['red','black','green','orange']

			if i == 3:

				ax[2].set_ylabel('nPa')

				ax[2].set_ylim(-10,25)

				# get rid of peaks

				y = np.where( y_vals[i] == np.max(y_vals[i]) , np.nan , y_vals[i])

				ax[2].plot(axsmer, y, label= y_labs[i], color = c[i])

				n = 1

			else:

				ax[1].plot(axsmer, y_vals[i], label=y_labs[i], color = c[i])

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

		if plots == 'dir':

			print('dir graphs:', plots)

			dictt = dict_dn[f'{dykey1}']

			tk1 = dict_dn[f'{dtimeskey1}']

			dict_all = dict_dn

		else:

			print('undir graphs:', plots)

			dictt = dict_n[f'{uykey1}']

			tk1 = dict_n[f'{utimeskey1}']

			dict_all = dict_n


		# -------------------------------------------------------------


		for ind, num in enumerate([0,1]):

			n = num_plots - 2-3  #len(non_empty_pcbands) -3

			# xd = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in np.squeeze( tk1[num] ) ]

			xd = tk1[num] 

			# print(plots, len(xd), xd, len(dictt[num]), dictt[num])

			if len(xd)!=0:
	
				markers, stemlines, baseline = ax[ind+n].stem(xd, dictt[num], label=f'Pc{num+2}  average degree')

				plt.setp(markers, color='purple', markersize=7, marker = '.')

				plt.setp(stemlines, color='grey', linestyle='--')

				nodes = np.array(dict_all['n_nodes'][num])

				# print('nodes',nodes)

				max_pos_conn =  nodes * (nodes -1) /2

				max_avg_deg = max_pos_conn / nodes

				ax[ind+n].scatter(xd, max_avg_deg, label=f'Pc{num+2} max average degree', marker='_', color='black')



			else:
				continue
			
			# axs2 = ax[ind+n].twinx()

			# second yaxis-----------

			# print('times',len(xd),xd)

			# for order in ['1st','2nd','3rd']:

			# 	# print('connections',len(dictt[order]['e_num'][num]),dictt[order]['e_num'][num])

			# 	axs2.scatter(xd, dictt[order]['e_num'][num], s=7, label= f'Pc{num+2} {order} \n # connections')

			# # y2range = range(int(np.min(dictt['3rd']['e_num'][num])), int(np.max(dictt['1st']['e_num'][num])),2)

			# # print(np.min(dict_dn[f'{dykey1}']['2nd']['e_num'][num]), np.max(dict_dn[f'{uykey1}']['1st']['e_num'][num]), 'ranges min, max', i)

			# # print(dict_dn[f'{dykey1}']['1st']['e_num'][num], '1st')
			
			# # print(dict_dn[f'{dykey1}']['2nd']['e_num'][num], '2nd')

			# axs2.set_ylabel('connections')

			# # yax2.legend(bbox_to_anchor=(0.8,1.02,0.2,0.2), loc="lower right", borderaxespad=0, ncol=1)

			# axs2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='x-small', ncol=1)

			# axs2.grid()

			# axs2.set_yticks(y2range)

			# ax[ind+n].set_xlabel('time (UTC)')
     
			ax[ind+n].set_ylabel('avg. degree')

			# print(index_date[0],index_date.iloc[-1],type(datel[0]))

			formatter = mdates.DateFormatter("%H:%M:%S")

			ax[ind+n].xaxis.set_major_formatter(formatter)

			ax[ind+n].legend()

			ax[ind+n]

			ax[ind+n].grid()

		# Loop for ULF power and MLT, MLAT plots
		
		point_style = ['.','^','>','<']

		# direction, c1, c2 = label.split('_')


		for ind, num in enumerate([num_plots-2, num_plots-1, num_plots-3]):

			# axes index stats from 0 to n-1

			for indpc, k in enumerate([0,1]):#enumerate(non_empty_pcbands):

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

		
				# xn = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in np.squeeze(tk1[k]) ]

				xn = tk1[k]

				labels = []

				# print(dict_all.keys())

				# if statments to perform correct operation for each plot

				if ind <=1:

					for m in ['1st','2nd','3rd']:


						ax[num].scatter(xn, y1(m), s=7, label= f'Pc{k+2}, {m}, x', marker = point_style[indpc])


						# print(ylab,m,len(dict_all[f'ulf_pc_{k+2}'][m][ylab]), dict_all[f'ulf_pc_{k+2}'][m][ylab])

						# print(len(xn),xn)


						ax[num].scatter(xn, dict_all[f'ulf_pc_{k+2}'][m][ylab] , s=7, label= f'Pc{k+2}, {m}, power', marker = point_style[indpc])

						labels.append(f'Pc{k+2}, {m}, x')

						labels.append(f'Pc{k+2}, {m}, power')

				# statment for speterate plot

				else:
					
					for m in ['1st','2nd','3rd']:
						
						ax[num].scatter(xn, dict_all[f'ulf_pc_{k+2}'][m]['pc_power'] , s=7, label= f'Pc{k+2}, {m}, power', marker = point_style[indpc])

						labels.append(f'Pc{k+2}, {m}, power')


			# print('minmax',maxmin_vals)

			# print(min(maxmin_vals),max(maxmin_vals))

			ax[num].set_ylabel(ylab)

			formatter = mdates.DateFormatter("%H:%M:%S")

			ax[num].xaxis.set_major_formatter(formatter)

			l1, l2 = np.array_split(labels,2)

			# colours automatically aligned between subplots can check

			if ind==1:

				ax[num].legend(loc='upper center', bbox_to_anchor=(0.5, -0.50), ncol=12, fontsize='small',title='MLAT and MLT')

			elif ind == 2:

				ax[num].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='x-small', ncol=1, title='Pc power')


			ax[num].grid()	

			ax[num].set_xlim(datemin,datemax)
					

			if save_plot:

				# saves the plot, label includes handel for the event

				plt.savefig(f'plots/{dykey1}_{label}_net_{comp}.png')


		plt.show()




# order of operations

# load networks (with proper window parameter), run analysis files, run plotting files



def save_plot( comp , pcbands,   ulf_filename, plot=True, save_plot=False, analysis_file=True):

	# option given within the function using save or plot to true or false
	# with component (comp) for magnetic field 'n', 'e', 'z'
	# nettype can be either 'dir', 'undir' or 'all given as an array'


	# handel for labeling anlysis files for each event, period multiple for window and 0.3 being peak height

	w_inc = '15'

	event_id = '170313'

	handel = f'test_{comp}_{event_id}_{w_inc}_0.3'


	if analysis_file:

		# # load network arrays for analysis

		dna, na = load_arrays(pcbands, comp, handel)

		# saving analysis files for both networks

		# heatmap_analysis_file(dna, pcbands, comp, 'dirr', f'window_inc_{w_inc}_peakh_0.3', ' ', ulf_filename)

		# heatmap_analysis_file(na, pcbands, comp, 'undir', f'window_inc_{w_inc}_peakh_0.3', ' ', ulf_filename)

		heatmap_analysis_parallel_file(dna, pcbands, comp, 'dirr', f'window_inc_{w_inc}_peakh_0.3', ' ', ulf_filename)
 
		heatmap_analysis_parallel_file(na, pcbands, comp, 'undir', f'window_inc_{w_inc}_peakh_0.3', ' ', ulf_filename)




save_plot('e', [0,1], 'networks_data/17march2013_4am.csv')





















