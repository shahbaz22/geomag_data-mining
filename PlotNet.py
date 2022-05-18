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
from matplotlib import rc
from tqdm import tqdm
import sys
import os

class PlotNet: 
    def __init__(self, fname_arr:list, ulf_fname:str, n:int, comp:str, pc:str, t_start:str, t_end:str):
        self.n = n
        self.fname_arr = fname_arr
        self.ulf_fname = ulf_fname
        self.pc = pc
        self.t_start = t_start
        self.t_end = t_end
        self.comp = comp
        # fname always ends in date and not .txt and includes path
        self.date = fname_arr[0].split('_')[-1]
        self.dict_labs = [label.split('/')[-1].split(f'_{self.comp}_')[0] for label in self.fname_arr]
        # create save path and save_path fname based on inputs

    def utc_sort_reduce(self, utc_times, deliminator='T'):
        # takes an unordered list of utc times and removes repeated values and then orders based on tme value
        times_dn = list(set(utc_times))
        times_dn = sorted(times_dn, key = lambda x: datetime.datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S'))
        times_dn = [ datetime.datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S') for x in np.squeeze(times_dn) ]
        return times_dn

    def ulf_power(self, ulf_power_filename, times, band):
        '''mean ulf powers for pc single band for corespondig times'''
        ulf_powers = pd.read_csv(ulf_power_filename)
        ulf_powers['Date_UTC'] = pd.to_datetime(ulf_powers['Date_UTC'])
        dict_ulf ={'mean_pc_power':[]}
        for i in np.squeeze(times):
            ts = ulf_powers[ulf_powers['Date_UTC'] == i][[ 'Date_UTC', 'IAGA','PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']] 
            powers = ['PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']
            dict_ulf['mean_pc_power'].append(ts[powers[band]].mean())
        return dict_ulf

    def plotting_files(self):
        '''Func to create all plotting data including pivot talbe for undirected and directed networks
        also ulf power from edgelist network files'''
        mlt_lower_lim = 9
        pivot_table_dict = {}
        pivot_table_dict_mlt = {}
        if self.pc == 'pc2':
            num = 0
        elif self.pc == 'pc3':
            num = 1
        for ind, val in enumerate(self.fname_arr):
            dict_lab = self.dict_labs[ind]
            pivot_table_dict[dict_lab] = []
            pivot_table_dict_mlt[dict_lab] = [] 
            netdistr = NetDistr(f'{val}.txt', self.n)
            print(dict_lab)
            print(val)
            if netdistr.num_edges<2:
                print(dict_lab, 'empty')
                continue
            pt = netdistr.create_pivottable(self.t_start, self.t_end, self.pc)
            pt.index = pd.to_datetime(pt.index)
            pivot_table_dict[dict_lab] = pt
            pt = netdistr.create_pivottable_by_mlt(mlt_lower_lim, self.t_start, self.t_end, self.pc)
            pt.index = pd.to_datetime(pt.index)
            pivot_table_dict_mlt[dict_lab] = pt
    #   for calculating ulf powers for corespending netowrk times
        times = pd.date_range(f'{netdistr.date} {self.t_start}', f'{netdistr.date} {self.t_end}', freq='30s')
        ulf_dict = self.ulf_power(ulf_fname, times, num)
        pc_power_arr = ulf_dict['mean_pc_power']

        return pivot_table_dict, pivot_table_dict_mlt, pc_power_arr, times


    def plot(self, pivot_table:dict, pivot_table_mlt:dict, pc_power:list, power_times:list, indices_fname:str, y_space:int, save_plot:bool):
        plt.rcParams.update({'font.size': 20})
        num_plots = len(self.fname_arr) + 4
        fig = plt.figure(figsize=(20, 17), constrained_layout=True)
        gs = fig.add_gridspec(num_plots, 3, width_ratios=[30,1,1])
        indices = pd.read_csv(indices_fname)
        index_date = indices['Date_UTC']
        sme = indices['SME']
        axsmer = pd.to_datetime(index_date)
        print('need to add date here to generalise for other events, along with other indices_fname')
        # datemin = datetime.strptime(f'{self.t_start}', '%H:%M:%S')
        # datemax = datetime.strptime(f'{self.t_end}', '%H:%M:%S')
        datemin = np.datetime64(index_date.iloc[0])
    #   can add t_Start and end
        datemax = datetime.strptime('2015-03-17T12:00:00', '%Y-%m-%dT%H:%M:%S')
        f1 = fig.add_subplot(gs[0,0])
        f1.plot(axsmer,sme, color='black', label = 'SME')
        f1.set_ylabel('SME (nT)')
        f1.grid()
    #   turn off labels for x axis while preserving ticks
        f1.set_xticklabels([])                           
        f1.set_xlim(datemin,datemax)
        formatter = mdates.DateFormatter("%H:%M:%S")
        f1.xaxis.set_major_formatter(formatter)
        f1.legend()
        f1.set_xlim(datemin,datemax)
        ax2 = f1.twinx()
        
        # code for plotting SMR sections for entire SMR yearly dataset
        indices['Date_UTC'] = pd.to_datetime(indices['Date_UTC']) 
        for i in ['SMR00','SMR06','SMR12','SMR18' ]:
            x = axsmer[(indices['Date_UTC'] >= datemin) & (indices['Date_UTC'] <= datemax)]
            y = indices[(indices['Date_UTC'] >= datemin) & (indices['Date_UTC'] <= datemax)][i]
            ax2.plot(x , y, label= i)
        ax2.xaxis.set_major_formatter(formatter)
        ax2.legend(loc=10)
        ax2.set_ylabel('SMR (nT)')
        ax2.grid()
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
                y = np.where( y_vals[i] == np.max(y_vals[i]) , np.nan , y_vals[i])
                f2.plot(axsmer, y, label=y_labs[i], color = c[i])
                f = f2

            f.xaxis.set_major_formatter(formatter)
            f.set_xlim(datemin,datemax)
            f.legend()
            f.grid()
            f.set_xticklabels([])
            f2.set_ylabel('nT')
            f2.set_ylim(-30,30)
            f2.minorticks_on()
            f2.set_xlim(datemin,datemax)
            #f2.set_xticklabels([])

        f4 = fig.add_subplot(gs[num_plots-1,0])
        f4.grid()
        f4.plot(power_times, pc_power, label= f'{pc} mean power')
        f4.legend(loc=4,markerscale=5)
        f4.xaxis.set_major_formatter(formatter)
        f4.set_xlim(datemin,datemax)
        f4.minorticks_on()
        f4.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.3)
        f4.set_ylabel('log(nT)^2')
        # plotting heatmaps 
        colour1 = 'Purples_r'
        colour2 = 'Oranges_r'
        y_tick_spacing = y_space 
        x_tick_spacing = 3
        for i, label in enumerate(self.dict_labs):
            i = i + 3 # fig position
            heat_ax1, colorbar_ax1, pt1 = fig.add_subplot(gs[i,0]), fig.add_subplot(gs[i,1]), pivot_table[label]    
            colorbar_ax12, pt1_mlt =  fig.add_subplot(gs[i,2]), pivot_table_mlt[label]
            ax1 = sb.heatmap(pt1.T, xticklabels = x_tick_spacing, yticklabels = y_tick_spacing, ax = heat_ax1, cmap = colour1, 
            cbar_kws={'label': 'frequency'}, cbar_ax = colorbar_ax1)
            ax1.axes.get_xaxis().set_visible(False)
            heat_ax1.grid()
            heat_ax1.text(0.1, 1, f'{label}_{pc}_{comp}', transform= heat_ax1.transAxes, bbox = dict(facecolor = 'white', alpha = 0.5))
            heat_ax1.invert_yaxis()
            heat_ax1.set_ylabel('degrees')
            heat_ax1.set_facecolor('grey')
            ax1 = sb.heatmap(pt1_mlt.T, xticklabels = x_tick_spacing, yticklabels = y_tick_spacing, ax = heat_ax1, cmap = colour2, 
            cbar_kws={'label': 'freq. MLT>9 hrs'}, cbar_ax = colorbar_ax12)
            ax1.axes.get_xaxis().set_visible(False)
            heat_ax1.grid()
            heat_ax1.invert_yaxis()
            heat_ax1.set_ylabel('degrees')
            heat_ax1.set_facecolor('grey')
            heat_ax1.set_ylim(0,80)
        
        if save_plot:
            s_path = f'/Users/sc/conda_envs/geomag/plots/{self.year}_nets/comp_{self.comp}/'
            path_exists = os.path.exists(s_path)
            if not path_exists:
                # Create a new directory because it does not exist 
                os.makedirs(s_path)
                print(f"The new directory {s_path} is created!")
                # sh = self.fname_arr[0]
            plt.savefig(f'{path}nets_try_{pc}_dist_ts.png', facecolor='w')

        plt.show()
        
# handel for labeling anlysis files for each event, period multiple for window and 0.3 being peak height
print('ideally, bellow code would be in a jupyter notebook')
comp ='e'
w_inc = 20
peakh = 0.3
event_id = '2015'
pc = 'pc2'
num_stations = 128
# dir_net_t10_e_128_0.3_0.25_2.5_2015
# dir_net_t10_e_128_0.3_0.25_2.5_2015
# dir_net_t11_e_128_0.3_0.25_2.5_2015
# dir_net_t11_e_128_0.3_0.25_2.5_2015
# fhandel = f'{comp}_{w_inc}_peakh_{peakh}_num_stations_{num_stations}_{event_id}_0.3'
# fhandel = f'surrogate_e_128_pc2_pc3_0.25_2.5_2015'
# fhandel = f'e_128_pc2_pc3_0.01_0.1_2015'
# fhandel = f'e_128_pc2_pc3_0.025_0.25_2015'
# fhandel = f'{comp}_{w_inc}_peakh_{peakh}_num_{num_stations}_{event_id}_0.25_2.5'
if pc == 'pc2':
    i=0
elif pc == 'pc3':
    i=1
fhandel = f'{i}_e_128_0.3_0.25_2.5_2015'
fhandel = f'{i}_surr_e_128_0.3_0.25_2.5_2015'
path = f'networks_data/{event_id}_nets/comp_{comp}/'
ulf_fname = 'networks_data/ulf_spd_2015.csv'
# fname_arr = ['dir_net']
fname_arr = ['dir_net_t1', 'dir_net_tn']

fname_arr = [path + x + fhandel for x in fname_arr]
indices_fname = 'networks_data/indices_spd_2015.csv'
t_start = '4:00:00'
t_end = '12:00:00'
# need to include date and 
pn = PlotNet( fname_arr, ulf_fname, num_stations, comp, pc, t_start, t_end)
pivot_table_dict, pivot_table_dict_mlt, pc_power_arr, p_times = pn.plotting_files()
pn.plot(pivot_table_dict, pivot_table_dict_mlt, pc_power_arr, p_times, indices_fname, y_space=10, save_plot=False)


