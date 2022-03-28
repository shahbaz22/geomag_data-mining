import networkx as nx
import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
import time
import functools
from collections import Counter
#from nxviz import MatrixPlot
from datetime import datetime 
import matplotlib.dates as mdates
import pdb
import matplotlib.ticker as mticker
import pickle
from NetDistr import NetDistr
import seaborn as sb
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec



# handel for labeling anlysis files for each event, period multiple for window and 0.3 being peak height

w_inc = '20'
peakh = 0.3
event_id = '170313'
ulf_filename = 'networks_data/17march2013_4am.csv'
path = 'networks_data/'

def utc_sort_reduce(utc_times, deliminator='T'):

	# takes an unordered list of utc times and removes repeated values and then orders based on tme value
	times_dn = list(set(utc_times))
	times_dn = sorted(times_dn, key = lambda x: datetime.datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S'))
	times_dn = [ datetime.datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S') for x in np.squeeze(times_dn) ]
	return times_dn


def ulf_power(ulf_power_filename, times, band):
	'''top two ulf powers for pc single band'''

	ulf_powers = pd.read_csv(ulf_power_filename)
	ulf_powers['Date_UTC'] = pd.to_datetime(ulf_powers['Date_UTC'])
	# print(md1.head())
	# print(print(md1.columns.tolist()))
	# returns only those rows with the same labels as those in the timeseries data set used to make the networks
	dict_ulf ={}
	# highest values ulf power arrays to be recorded
	values = ['1st', '2nd']
	for i in values:
		dict_ulf[i] = { 'pc_power' : [] }

	for i in np.squeeze(times):
		# for given time i, return the ulf values, then pick top three ulf power values, for those top three pick the MLT and MLAT
		# .replace used to unify string format for times
		# i = i.replace('T', ' ')

		ts = ulf_powers[ulf_powers['Date_UTC'] == i][[ 'Date_UTC', 'IAGA','PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']] 
		# values are 1st 2nd and 3rd
		powers = ['PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']

		for j, lab in enumerate(values):
			# power labels for relevent band		
			# oganise by values largest to smallest
			# if rows less than 3 then groupby doesn't work properly
			ts = ts.groupby(powers[band]).head(n=4)
			dict_ulf[lab]['pc_power'].append(ts[powers[band]].iloc[j])

	return dict_ulf


def plotting_files(filename_handel, ulf_filename, num_stations1):
    '''Func to create all plotting data including pivot talbe for undirected and directed networks
    also ulf power from edgelist network files'''
    global comp
    comp = 'z'
    mlt_lower_lim = 9
    pivot_table_dict = { 'dir_net':{0:[],1:[]}, 'undir_net_inphase':{0:[], 1:[]} }
    pivot_table_dict_mlt = { 'dir_net':{0:[],1:[]}, 'undir_net_inphase':{0:[], 1:[]} }
    pc_power_arr = [[[],[]],[[],[]]] # for pc2 and pc3
    time_arr = [[],[]]
    
    for label in ['dir_net','undir_net_inphase']:
        for num in [0,1]:
            print(label,num, comp)
            netdistr = NetDistr(f'{path}{label}{num}_{comp}_{filename_handel}.txt', num_stations1)
            # pt = netdistr.create_pivottable()
            # pt.index = pd.to_datetime(pt.index)
            # pivot_table_dict[label][num] = pt
            pt = netdistr.create_pivottable_by_mlt(mlt_lower_lim)
            pt.index = pd.to_datetime(pt.index)
            pivot_table_dict_mlt[label][num] = pt
            
            
#     for calculating ulf powers for corespending netowrk times
    for num in [0,1]:
        times = pivot_table_dict['dir_net'][num].index
        time_arr[num] = times
        ulf_dict = ulf_power(ulf_filename, times, num)
        for ind, m in enumerate(['1st', '2nd']):
            pc_power_arr[num][ind] = ulf_dict[m]['pc_power']

    return pivot_table_dict, pivot_table_dict_mlt, pc_power_arr, time_arr
        
pivot_table_dict, pivot_table_dict_mlt, pc_power_arr, times_arr = plotting_files('20-25_peakh_0.3_num_stations_122_170313', ulf_filename, 122)


def plot_heatmap(pivot_table_dict, pivot_table_dict_mlt, pc_power_arr, indices_filename, ulf_filename, label, save_plot = False):

    # test_ulf('networks_data/spd_analysis/avg_deg_noon_e', dirr=True)

    # ref plots include non network plots include SME/SMR, B feilds, SP, all pc ULF power, all pc MLT, all pc MLAT (in this case more ulf band powers)

    ref_plots = 4
    num_plots = ref_plots + 4 
    fig = plt.figure(figsize=(20, 17), constrained_layout=True)
    gs = fig.add_gridspec(num_plots, 3, width_ratios=[30,1,1])
    # fig.suptitle(f'Heatmap plots for {label} netowrks')
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 15})

    indices = pd.read_csv(indices_filename)
    print(indices.head())

    index_date = indices['Date_UTC']
    sme = indices['SME']
    axsmer = pd.to_datetime(index_date)

    print(type(axsmer[0]))
    
    f1 = fig.add_subplot(gs[0,0])
    f1.plot(axsmer,sme, color='black', label = 'SME')
    f1.set_ylabel('SME (nT)')
    f1.grid()
#     f1.set_xticklabels([])

    datemin = np.datetime64(index_date.iloc[0])
    datemax = np.datetime64(index_date.iloc[-1])

    f1.set_xlim(datemin,datemax)
    formatter = mdates.DateFormatter("%H:%M:%S")
    f1.xaxis.set_major_formatter(formatter)
    f1.legend()
    ax2 = f1.twinx()

    # code for plotting SMR sections for entire SMR yearly dataset
    indices['Date_UTC'] = pd.to_datetime(indices['Date_UTC']) 

    for i in ['SMR00','SMR06','SMR12','SMR18' ]:

        ax2.plot(axsmer,indices[(indices['Date_UTC'] >= datemin) & (indices['Date_UTC'] <= datemax)][i] , label= i)

    ax2.xaxis.set_major_formatter(formatter)
    ax2.legend()
    ax2.set_ylabel('SMR (nT)')
    ax2.grid()
    datemin = np.datetime64(index_date.iloc[0])
    datemax = datetime.strptime('2013-3-17 12:00:00', '%Y-%m-%d %H:%M:%S')
    ax2.set_xlim(datemin,datemax)
    ax2.set_xticklabels([])


    dyn_pressure = indices['PDYN']
    # GSM field values
    bz, bx, by = indices['GSM_Bz'], indices['GSM_Bx'], indices['GSM_By']
    f2 = fig.add_subplot(gs[1,0])
    f3 = fig.add_subplot(gs[2,0])
    for i in [0,1,2,3]:
        y_vals = [bx, by, bz, dyn_pressure]
        y_labs = ['GSM_Bz','GSM_Bx','GSM_By','Dynamic pressure']
        c = ['red','black','green','orange']

        if i == 3:
            f3.set_ylabel('nPa')
            f3.set_ylim(-10,25)
            # get rid of peaks
            y = np.where( y_vals[i] == np.max(y_vals[i]) , np.nan , y_vals[i])
            f3.plot(axsmer, y, label= y_labs[i], color = c[i])
            f = f3
        else:
            f2.plot(axsmer, y_vals[i], label=y_labs[i], color = c[i])
            f = f2

        f.xaxis.set_major_formatter(formatter)
        f.set_xlim(datemin,datemax)
        f.legend()
        f.grid()
        f.set_xticklabels([])
        f2.set_ylabel('nT')
        f2.set_ylim(-30,30)
        f2.set_xticklabels([])

    # plotting heatmaps 
    f4 = fig.add_subplot(gs[7,0])
    f4.grid()
    xfmt = mdates.DateFormatter('%H:%M:%S')
    for num in [0,1]:
        times = times_arr[num]
#         times =  [x.to_pydatetime() for x in pivot_arr[num].index ]
#         print(type(times),type(times[0]))
        for ind, m in enumerate(['1st', '2nd']):
            f4.scatter(times, pc_power_arr[num][ind] , s=7, label= f'Pc{ind+2}, {m}, power')
            f4.legend()
            f4.xaxis.set_major_formatter(formatter)
            f4.set_xlim(datemin,datemax)
            
    heat_ax1, colorbar_ax1, pt1 = fig.add_subplot(gs[3,0]), fig.add_subplot(gs[3,1]), pivot_table_dict['dir_net'][0]
    heat_ax2, colorbar_ax2, pt2 = fig.add_subplot(gs[4,0]), fig.add_subplot(gs[4,1]), pivot_table_dict['dir_net'][1]
    heat_ax3, colorbar_ax3, pt3 = fig.add_subplot(gs[5,0]), fig.add_subplot(gs[5,1]), pivot_table_dict['undir_net_inphase'][0]
    heat_ax4, colorbar_ax4, pt4 = fig.add_subplot(gs[6,0]), fig.add_subplot(gs[6,1]), pivot_table_dict['undir_net_inphase'][1]
    colour1 = 'Purples_r'
    colour2 = 'Oranges_r'
    colorbar_ax12, pt1_mlt =  fig.add_subplot(gs[3,2]), pivot_table_dict_mlt['dir_net'][0]
    colorbar_ax22, pt2_mlt =  fig.add_subplot(gs[4,2]), pivot_table_dict_mlt['dir_net'][1]
    colorbar_ax32, pt3_mlt =  fig.add_subplot(gs[5,2]), pivot_table_dict_mlt['undir_net_inphase'][0]
    colorbar_ax42, pt4_mlt =  fig.add_subplot(gs[6,2]), pivot_table_dict_mlt['undir_net_inphase'][1]
    
    
    ax1 = sb.heatmap(pt1.T, xticklabels = 3, yticklabels =10, ax = heat_ax1, cmap = colour1, 
               cbar_kws={'label': 'frequency'}, cbar_ax = colorbar_ax1)
    ax1.axes.get_xaxis().set_visible(False)
    heat_ax1.grid()
    heat_ax1.text(0.1, 1, f'directed_network_Pc2_{comp}', transform= heat_ax1.transAxes, bbox = dict(facecolor = 'white', alpha = 0.5))
    heat_ax1.invert_yaxis()
    heat_ax1.set_ylabel('degrees')
    heat_ax1.set_facecolor('grey')
    ax1 = sb.heatmap(pt1_mlt.T, xticklabels = 3, yticklabels = 10, ax = heat_ax1, cmap = colour2, 
               cbar_kws={'label': 'freq. MLT>9 hrs'}, cbar_ax = colorbar_ax12)
    ax1.axes.get_xaxis().set_visible(False)
    heat_ax1.grid()
    heat_ax1.invert_yaxis()
    heat_ax1.set_ylabel('degrees')
    heat_ax1.set_facecolor('grey')
    
    ax2 = sb.heatmap(pt2.T, xticklabels = 3, yticklabels = 10, ax = heat_ax2, cmap = colour1, 
               cbar_kws={'label': 'frequency'}, cbar_ax = colorbar_ax2)
    ax2.axes.get_xaxis().set_visible(False)
    heat_ax2.grid()
    heat_ax2.text(0.1, 1, f'directed_network_Pc3_{comp}', transform= heat_ax2.transAxes, bbox = dict(facecolor = 'white', alpha = 0.5))
    heat_ax2.invert_yaxis()
    heat_ax2.set_ylabel('degrees')
    heat_ax2.set_facecolor('grey')
    ax2 = sb.heatmap(pt2_mlt.T, xticklabels = 3, yticklabels = 10, ax = heat_ax2, cmap = colour2, 
               cbar_kws={'label': 'freq. MLT>9 hrs'}, cbar_ax = colorbar_ax22)
    ax2.axes.get_xaxis().set_visible(False)
    heat_ax2.grid()
    heat_ax2.invert_yaxis()
    heat_ax2.set_ylabel('degrees')
    heat_ax2.set_facecolor('grey')

    ax3 = sb.heatmap(pt3.T, xticklabels = 3, yticklabels = 10, ax = heat_ax3, cmap = colour1, 
               cbar_kws={'label': 'freq.'}, cbar_ax = colorbar_ax3)
    ax3.axes.get_xaxis().set_visible(False)
    heat_ax3.grid()
    heat_ax3.text(0.1, 1, f'undirected_network_Pc2_{comp}', transform= heat_ax3.transAxes, bbox = dict(facecolor = 'white', alpha = 0.5))
    heat_ax3.invert_yaxis()
    heat_ax3.set_ylabel('degrees')
    heat_ax3.set_facecolor('grey')
    ax3 = sb.heatmap(pt3_mlt.T, xticklabels = 3, yticklabels = 10, ax = heat_ax3, cmap = colour2, 
               cbar_kws={'label': 'freq. MLT>9 hrs'}, cbar_ax = colorbar_ax32)
    ax3.axes.get_xaxis().set_visible(False)
    heat_ax3.grid()
    heat_ax3.invert_yaxis()
    heat_ax3.set_ylabel('degrees')
    heat_ax3.set_facecolor('grey')

    ax4 = sb.heatmap(pt4.T, xticklabels = 3, yticklabels = 10, ax = heat_ax4, cmap = colour1, 
               cbar_kws={'label': 'frequency'}, cbar_ax = colorbar_ax4)
    ax4.axes.get_xaxis().set_visible(False)
    heat_ax4.grid()
    heat_ax4.text(0.1, 1, f'undirected_network_Pc3_{comp}', transform= heat_ax4.transAxes, bbox = dict(facecolor = 'white', alpha = 0.5))
    heat_ax4.invert_yaxis()
    heat_ax4.set_ylabel('degrees')
    heat_ax4.set_facecolor('grey')
    heat_ax4.set_ylim([0,20])
    ax4 = sb.heatmap(pt4_mlt.T, xticklabels = 3, yticklabels = 10, ax = heat_ax4, cmap = colour2, 
               cbar_kws={'label': 'freq. MLT>9 hrs'}, cbar_ax = colorbar_ax42)
    ax4.axes.get_xaxis().set_visible(False)
    heat_ax4.grid()
    heat_ax4.invert_yaxis()
    heat_ax4.set_ylabel('degrees')
    heat_ax4.set_facecolor('grey')
    heat_ax4.set_ylim([0,20])
        
    if save_plot:
        plt.savefig(f'plots/{label}_net.png', facecolor='w')

    plt.show()
    
label = f'event_{event_id}_{comp}_comp__windmult_{w_inc}_peakh_{0.3}_70_stations'

plot_heatmap(pivot_table_dict, pivot_table_dict_mlt, pc_power_arr, 'networks_data/170313_indices.csv', ulf_filename, label=f'B_{comp}')

