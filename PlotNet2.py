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
from matplotlib.colors import LogNorm
import os

class PlotNet: 
    def __init__(self, net_labs:list, ulf_fname:str, n:int, dt_start:str, dt_end:str, surrogate:bool):
        'comp needs to be input array and also bands'
        self.n = n
        self.net_labs= net_labs
        self.ulf_fname = ulf_fname
        self.dt_start = dt_start
        self.dt_end = dt_end
        # fname always ends in date and not .txt and includes path
        self.year = dt_start.split(' ')[0].split('-')[0]
        self.surrogate = surrogate
        # create save path and save_path fname based on inputs

    # can ajust file hadnels and paths from this function
    # passes file information to NetDistr class to then create dictionaries and dataframes for plotting
    def path_and_fhandel(self, comp, band):
        path = f'networks_data/{self.year}_nets/comp_{comp}/'
        if self.surrogate == True:
            fhandel = f'{band}_surr_{comp}_{self.n}_0.3_0.25_2.5_{self.year}.txt'
        else:
            fhandel = f'{band}_{comp}_{self.n}_0.3_0.25_2.5_{self.year}.txt'

        return path, fhandel

    def utc_sort_reduce(self, utc_times, deliminator='T'):
        # takes an unordered list of utc times and removes repeated values and then orders based on tme value
        times_dn = list(set(utc_times))
        times_dn = sorted(times_dn, key = lambda x: datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S'))
        times_dn = [ datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S') for x in np.squeeze(times_dn) ]
        return times_dn

    def ulf_power(self, band):
        '''mean ulf powers for pc single band for corespondig times'''
        ulf_powers = pd.read_csv(self.ulf_fname)
        ulf_powers['Date_UTC'] = pd.to_datetime(ulf_powers['Date_UTC'])
        dict_ulf ={'mean_pc_power':[]}
        power_times = pd.date_range(self.dt_start, self.dt_end, freq='30s')
        for i in np.squeeze(power_times):
            ts = ulf_powers[ulf_powers['Date_UTC'] == i][[ 'Date_UTC', 'IAGA','PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']] 
            powers = ['PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']
            dict_ulf['mean_pc_power'].append(ts[powers[band]].mean())
        return dict_ulf, power_times

    def plotting_files(self, net_label:str, band:int, comp:str, manual_label=False):
        '''Func to create all plotting data including pivot talbe for undirected and directed networks
        also ulf power from edgelist network files'''
        mlt_lower_lim = 9
        path, fhandel = self.path_and_fhandel(comp, band)

        if manual_label==False:
            fname = path + net_label + fhandel
        else: 
            fname = f'{path}{net_label}'
        
        try:
            netdistr = NetDistr(fname, self.n)
        except FileNotFoundError:
            print(f'File {fname} does not exist, change in self.path_and_fhandel')
        
        if netdistr.num_edges<2:
            print('pivot table empty')
            exit()
        pt = netdistr.create_pivottable(self.dt_start, self.dt_end, f'pc{band+2}')
        pt.index = pd.to_datetime(pt.index)
        pt_mlt = netdistr.create_pivottable_by_mlt(mlt_lower_lim, self.dt_start, self.dt_end, f'pc{band+2}')
        print(net_label, band, comp)
        return pt, pt_mlt

    def plotting_files_lags(self, net_label:str, band:int, comp:str, pm:int, mlt:bool, manual_label=False):
        '''Func to create all plotting data including pivot talbe for undirected and directed networks
        also ulf power from edgelist network files'''
        path, fhandel = self.path_and_fhandel(comp, band)
        if manual_label==False:
            fname = path + net_label + fhandel
        else: 
            fname = f'{path}{net_label}'
        
        try:
            netdistr = NetDistr(fname, self.n)
        except FileNotFoundError:
            print(f'File {fname} does not exist, change in self.path_and_fhandel')
            exit()

        if netdistr.num_edges<2:
            ('pivot table empty')
            exit()
        if mlt==True:
            pt = netdistr.create_pivottable_lags_mlt(9, self.dt_start, self.dt_end, f'pc{band+2}', period_mult=pm)
        else:
            pt = netdistr.create_pivottable_lags(self.dt_start, self.dt_end, f'pc{band+2}', period_mult=pm, maxlags=[20,20])
        pt.index = pd.to_datetime(pt.index)
        print(net_label, band, comp)
        return pt

    def combine_files(self, fnames:list, path:str, band:int):
        with open(f'{path}combined_net{band}.txt', 'w') as outfile:
            for fname in fnames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)            

    def combined_nets_conn(self, save_analysis_files:bool):
        'create file for plotting total number of connections for comined network for each component and pc band'
        num_conn = dict()
        for band in [0,1]:
            num_conn[band]={}
            for comp in ['n','e','z']:
                path, fhandel = self.path_and_fhandel(comp, band)
                if save_analysis_files==True:
                    fname_arr = [path + x + fhandel for x in self.net_labs]
                    self.combine_files(fname_arr, path, band)
                    print(f'combined networks file created! {band}, {comp}')
                netdistr = NetDistr(f'{path}combined_net{band}.txt', self.n)
                num_conn[band][comp] = netdistr.ts_all_conn()
                print(f'completed {band}, {comp}')
        return num_conn

    def all_nets_conn(self):
        num_conn = dict()
        for band in [0,1]:
            num_conn[band]={}
            for comp in ['n','e','z']:
                num_conn[band][comp] ={}
                for net_label in self.net_labs:
                    path, fhandel = self.path_and_fhandel(comp, band)
                    fname = path + net_label + fhandel
                    netdistr = NetDistr(fname, self.n)
                    num_conn[band][comp][net_label] = netdistr.ts_all_conn()
                    print(f'completed {band}, {comp}, {net_label}')
        return num_conn


    def all_nets_short_long_conn(self, band:int, comp:str):
        num_conn = dict()
        for net_label in self.net_labs:
            path, fhandel = self.path_and_fhandel(comp, band)
            fname = path + net_label + fhandel
            netdistr = NetDistr(fname, self.n)
            long_df, short_df = netdistr.ts_con_length_ratio_by_mlt(9)
            num_conn[net_label] = {'long': long_df, 'short': short_df} 
        return num_conn

    def single_net_all_comp_bands(self, net_type:str):
        num_conn = dict()
        for band in [0,1]:
            num_conn[band]={}
            for comp in ['n','e','z']:
                path, fhandel = self.path_and_fhandel(comp, band)
                fname = path + net_type + fhandel
                netdistr = NetDistr(fname, self.n)
                long_df, short_df = netdistr.ts_con_length_ratio_by_mlt(9)
                num_conn[band][comp] = {'long': long_df, 'short': short_df} 
        return num_conn

    def single_net_all_comp_and_bands_all_net_params(self, net_type:str, net_param:str):
        num_conn = dict()
        for comp in ['n','e','z']:
            path, fhandel = self.path_and_fhandel(comp, band = 0)
            fname = path + net_type + fhandel
            netdistr = NetDistr(fname, self.n)
            num_conn[comp] = {}
            
            if net_param == 'long-mlt_div_short-mlt':
                df = netdistr.ts_con_length_ratio_by_mlt(9)
                for ind, lab in enumerate(['ratio','top','bottom']):
                    num_conn[comp][lab] = df[ind]
            
            elif net_param == 'long_div_short':
                df = netdistr.ts_con_length_ratio_by_length()
                for ind, lab in enumerate(['ratio','top','bottom']):
                    num_conn[comp][lab] = df[ind]
            
            elif net_param == 'conj_div_north_south':                
                df = netdistr.ts_conj_div_north_south_conn() 
                for ind, lab in enumerate(['ratio','top','bottom']):
                    num_conn[comp][lab] = df[ind]            
            
            elif net_param == 'north_south_div_north':                
                df = netdistr.ts_north_south_div_north_conn() 
                for ind, lab in enumerate(['ratio','top','bottom']):
                    num_conn[comp][lab] = df[ind]            
            
            elif net_param == 'ew_div_we':                
                df = netdistr.ts_ew_div_we_conn() 
                for ind, lab in enumerate(['ratio','top','bottom']):
                    num_conn[comp][lab] = df[ind]            
            
            elif net_param == 'ns_div_sn':
                df = netdistr.ts_ns_div_sn_conn()
                for ind, lab in enumerate(['ratio','top','bottom']):
                    num_conn[comp][lab] = df[ind]            
            else:
                print('wrong net_param label used!') 
                print('choose from: conj_div_north_south, north_south_div_north, long_div_short, ew_div_we, ns_div_sn')
                break

        return num_conn
    
    def all_nets_conj_ns_conn(self, band:int, comp:str):
        num_conn = dict()
        for net_label in self.net_labs:
            path, fhandel = self.path_and_fhandel(comp, band)
            fname = path + net_label + fhandel
            netdistr = NetDistr(fname, self.n)
            df = netdistr.ts_conj_div_north_south_conn()
            num_conn[net_label] = df 
        return num_conn

    def single_net_conj_ns_all_comp_bands(self, net_type:str):
        num_conn = dict()
        for band in [0,1]:
            num_conn[band]={}
            for comp in ['n','e','z']:
                path, fhandel = self.path_and_fhandel(comp, band)
                fname = path + net_type + fhandel
                netdistr = NetDistr(fname, self.n)
                df = netdistr.ts_conj_div_north_south_conn()
                num_conn[band][comp] = df
        return num_conn
    
    def all_nets_n_ns_conn(self, band:int, comp:str):
        num_conn = dict()
        for net_label in self.net_labs:
            path, fhandel = self.path_and_fhandel(comp, band)
            fname = path + net_label + fhandel
            netdistr = NetDistr(fname, self.n)
            df = netdistr.ts_north_div_north_south_conn()
            num_conn[net_label] = df
        return num_conn

    def single_net_n_ns_all_comp_bands(self, net_type:str):
        num_conn = dict()
        for band in [0,1]:
            num_conn[band]={}
            for comp in ['n','e','z']:
                path, fhandel = self.path_and_fhandel(comp, band)
                fname = path + net_type + fhandel
                netdistr = NetDistr(fname, self.n)
                df = netdistr.ts_north_div_north_south_conn()
                num_conn[band][comp] = df
        return num_conn
    # can currenly run this code for all networks types and a single component and Pc band
    def plot_edge_ratio_ts(self, indices_fname:str):
#         plt.rcParams.update({'font.size': 20})
        num_plots = 5
        fig, ax = plt.subplots(num_plots, figsize=(12, 8))
        indices = pd.read_csv(indices_fname)
        index_date = indices['Date_UTC']
        sme = indices['SME']
        axsmer = pd.to_datetime(index_date)
        datemin = datetime.strptime(self.dt_start, '%Y-%m-%d %H:%M:%S')
        datemax = datetime.strptime(self.dt_end, '%Y-%m-%d %H:%M:%S')
        f1 = ax[0]
        f1.plot(axsmer,sme, color='black', label = 'SME')
        f1.set_ylabel('SME (nT)')
        f1.grid()
    #   turn off labels for x axis while preserving ticks
        f1.set_xticklabels([])                           
        f1.set_xlim(datemin,datemax)
        formatter = mdates.DateFormatter("%H:%M")
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
        # ax2.legend(loc=10)
        ax2.set_ylabel('SMR (nT)')
        ax2.grid()
        ax2.set_xlim(datemin,datemax)
        ax2.set_xticklabels([])

        dyn_pressure = indices['PDYN']
        # GSM field values
        bz, bx, by = indices['GSM_Bz'], indices['GSM_Bx'], indices['GSM_By']
        f2 = ax[1]
        f3 = ax[2]
        for j in [0,1,2,3]:
            y_vals = [bx, by, bz, dyn_pressure]
            y_labs = ['GSM_Bz','GSM_Bx','GSM_By','Dynamic pressure']
            c = ['red','black','green','orange']
            if j == 3:
                f3.set_ylabel('nPa')
                f3.set_ylim(-10,25)
                # get rid of peaks
                y = np.where( y_vals[j] == np.max(y_vals[j]) , np.nan , y_vals[j])
                f3.plot(axsmer, y, label= y_labs[j], color = c[j])
                f = f3
            else:
                y = np.where( y_vals[j] == np.max(y_vals[j]) , np.nan , y_vals[j])
                f2.plot(axsmer, y, label=y_labs[j], color = c[j])
                f = f2

            f.xaxis.set_major_formatter(formatter)
            f.set_xlim(datemin,datemax)
            f.legend()
            f.grid()
            ax3 = ax[2].twinx()
            for k in [0,1]:
                pc_power, power_times = self.ulf_power(k)
                print(pc_power)
                print(power_times)
                ax3.plot(power_times, pc_power, label= f'pc{k} mean p.')
                ax3.xaxis.set_major_formatter(formatter)
            ax3.legend()
            ax3.set_xlim(datemin,datemax)
            ax3.set_ylabel('log(nT)^2')
            f.set_xticklabels([])
            f2.set_ylabel('nT')
            f2.set_ylim(-30,30)
            f2.minorticks_on()
            f2.set_xlim(datemin,datemax)
            f2.set_xticklabels([])

        for band in [0,1]:
            ylims = dict()
            ylims[band]=[]
            for comp in ['n','e','z']:
                plot_num = 3+band
                vals = self.num_conn_dict[band][comp]
                ylims[band].append(max(vals['total']))
                t = pd.to_datetime(vals.index)
                ax[plot_num].plot(t ,vals['total'].values, label='Pc{band}, {comp}')
                ax[plot_num].legend(loc='lower left')
                ax[plot_num].xaxis.set_major_formatter(formatter)
                ax[plot_num].set_xlim(datemin,datemax)
            ax[plot_num].set_ylabel('all conn.')

        plt.show()   

    def plot_dist_ts(self, indices_fname:str, net_labs:list, comp:str, y_space:int, save_plot:bool, save_label:str):
        plt.rcParams.update({'font.size': 20})
        num_plots = len(self.net_labs) + 2
        fig = plt.figure(figsize=(20, 17), constrained_layout=True)
        gs = fig.add_gridspec(num_plots, 3, width_ratios=[30,1,1])
        indices = pd.read_csv(indices_fname)
        index_date = indices['Date_UTC']
        sme = indices['SME']
        axsmer = pd.to_datetime(index_date)
        datemin = datetime.strptime(self.dt_start, '%Y-%m-%d %H:%M:%S')
        datemax = datetime.strptime(self.dt_end, '%Y-%m-%d %H:%M:%S')
        f1 = fig.add_subplot(gs[num_plots-1,0])
        f1.plot(axsmer,sme, color='black', label = 'SME')
        f1.set_ylabel('SME (nT)')
        f1.grid()
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
        ax2.legend(loc='upper center',ncol=4, bbox_to_anchor=[0.5, 1.0])
        ax2.set_ylabel('SMR (nT)')
        ax2.grid()
        ax2.set_xlim(datemin,datemax)
        
        colour1 = 'Purples_r'
        colour2 = 'Oranges_r'
        y_tick_spacing = y_space 
        x_tick_spacing = 3
        step_arr = [0,2]
        for band in [0,1]:
            for i, label in enumerate(net_labs):
                i = step_arr[band] + i # fig position
                if band ==1:
                    i= len(self.net_labs)
                    label = 'combined_net'
                    pt1, pt1_mlt = self.plotting_files(f'combined_net{band}.txt', band, comp, manual_label=True)
                else:
                    pt1, pt1_mlt = self.plotting_files(label, band, comp)
                heat_ax1, colorbar_ax1, pt1 = fig.add_subplot(gs[i,0]), fig.add_subplot(gs[i,1]), pt1    
                colorbar_ax12, pt1_mlt =  fig.add_subplot(gs[i,2]), pt1_mlt
                ax1 = sb.heatmap(pt1.T, xticklabels = x_tick_spacing, yticklabels = y_tick_spacing, ax = heat_ax1, cmap = colour1, 
                cbar_kws={'label': 'frequency'}, cbar_ax = colorbar_ax1)
                ax1.axes.get_xaxis().set_visible(False)
                heat_ax1.grid()
                heat_ax1.text(0.1, 1, f'{label}_pc{band+2}_{comp}', transform= heat_ax1.transAxes, bbox = dict(facecolor = 'white', alpha = 0.5))
                heat_ax1.invert_yaxis()
                heat_ax1.set_ylabel('degrees')
                heat_ax1.set_facecolor('grey')
                ax1 = sb.heatmap(pt1_mlt.T, xticklabels = x_tick_spacing, yticklabels = y_tick_spacing, ax = heat_ax1, cmap = colour2, 
                cbar_kws={'label': 'freq. MLT>9 hrs'}, cbar_ax = colorbar_ax12)
                heat_ax1.grid()
                heat_ax1.invert_yaxis()
                heat_ax1.set_ylabel('degrees')
                heat_ax1.set_facecolor('grey')
                if band==1:
                    break

        
        if save_plot:
            year = self.dt_end.split('_')[0].split('-')[0]
            s_path = f'/Users/sc/conda_envs/geomag/plots/{year}_nets/comp_{comp}/'
            path_exists = os.path.exists(s_path)
            if not path_exists:
                # Create a new directory because it does not exist 
                os.makedirs(s_path)
                print(f"The new directory {s_path} is created!")
            plt.savefig(f'plots/{year}_nets/comp_{comp}/nets_{save_label}_dist_ts.png', facecolor='w')

        # plt.show()


    def plot_dist_ts_lags(self, indices_fname:str, comp:str, pm:int, save_label:str, save_plot:bool):
        plt.rcParams.update({'font.size': 18})
        num_plots = 3
        fig = plt.figure(figsize=(20, 17), constrained_layout=True)
        gs = fig.add_gridspec(num_plots, 3, width_ratios=[30,1,1])
        y_interval=[10,26]
        ylims = [[40,-40],[130,-130]]
        colour1 = 'seismic'
        for i in [0,1]:
            heat_ax1, colorbar_ax1 = fig.add_subplot(gs[i,0]), fig.add_subplot(gs[i,1])    
            label = 'combined_net'
            pt1 = self.plotting_files_lags(f'combined_net{i}.txt', i, comp, pm, mlt=False, manual_label=True)

            # ------ sorting by setting range of lags and limits
            # pt1 = pt1.reindex(sorted(pt1.columns), axis=1)
            try:
                pt1 = pt1[np.linspace(ylims[i][1],ylims[i][0],2*ylims[i][0]+1)]
            except KeyError:
                print(f'need to change column extent within {sorted(pt1.columns)}')
                raise KeyError()
            # # ---------
            # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            #             print(pt1[0])
            ax1 = sb.heatmap(pt1.T, xticklabels = 3, yticklabels = y_interval[i], ax = heat_ax1, cmap = colour1, 
            cbar_kws={'label': 'frequency'}, cbar_ax = colorbar_ax1)
            ax1.axes.get_xaxis().set_visible(False)
            heat_ax1.grid()
            heat_ax1.text(0.1, 1, f'{label}_pc{i+2}_{comp}', transform= heat_ax1.transAxes, bbox = dict(facecolor = 'white', alpha = 0.5))
            heat_ax1.invert_yaxis()
            heat_ax1.set_ylabel('lag')
            heat_ax1.set_facecolor('grey')
            heat_ax1.set_facecolor('xkcd:black')


        indices = pd.read_csv(indices_fname)
        index_date = indices['Date_UTC']
        sme = indices['SME']
        axsmer = pd.to_datetime(index_date)
        datemin = datetime.strptime(self.dt_start, '%Y-%m-%d %H:%M:%S')
        datemax = datetime.strptime(self.dt_end, '%Y-%m-%d %H:%M:%S')
        f1 = fig.add_subplot(gs[2,0])
        f1.plot(axsmer,sme, color='black', label = 'SME')
        f1.set_ylabel('SME (nT)')
        f1.grid()
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
        ax2.legend(loc='upper center',ncol=4, bbox_to_anchor=[0.5, 1.0])
        ax2.set_ylabel('SMR (nT)')
        ax2.grid()
        ax2.set_xlim(datemin,datemax)

        
        if save_plot:
            year = self.dt_end.split('_')[0].split('-')[0]
            s_path = f'/Users/sc/conda_envs/geomag/plots/{year}_nets/comp_{comp}/'
            path_exists = os.path.exists(s_path)
            if not path_exists:
                # Create a new directory because it does not exist 
                os.makedirs(s_path)
                print(f"The new directory {s_path} is created!")
            plt.savefig(f'plots/{year}_nets/comp_{comp}/nets_{save_label}_dist_lags_ts_0-3_lag.png', facecolor='w')
        plt.show()

# dn = NetDistr('networks_data/2015_nets/comp_e/combined_net0.txt',128)
# df = dn.ts_conj_conn()
# df.plot()
# plt.show()

# for i in [0]:
#     if i ==0:
#         year = '2015'
#         num_stations = 128
#         dt_start = '2015-03-17 4:00:00'
#         dt_end = '2015-03-17 12:00:00'
#         # dt_start ='2015-03-17 19:30:00'
#         # dt_end='2015-03-18 03:30:00'
#         ulf_fname = 'networks_data/ulf_spd_2015.csv'
#         indices_fname = 'networks_data/indices_spd_2015.csv'
#     if i ==1:
#         num_stations = 133
#         dt_start = '2012-01-21 20:00:00'
#         dt_end = '2012-01-22 16:00:00'
#         ulf_fname = 'networks_data/2012-01-21_22-ulf.csv'
#         indices_fname = 'networks_data/2012-01-21_22_indices.csv'
#     if i ==2:
#         num_stations = 122
#         dt_start = '2013-03-17 4:00:00'
#         dt_end = '2013-03-17 12:00:00'
#         ulf_fname = 'networks_data/17march2013_4am.csv'
#         indices_fname = 'networks_data/170313_indices.csv'

#     # hacky solution to include all these netwokr and combined Pc3 network
#     # net_labs = ['dir_net_t1', 'dir_net_tn', 'undir_net_inphase', 'undir_net_antiphase']
#     net_labs = ['dir_net_t1']
#     # pc2 t1, in_phase, anti_phase, pc3 combined

#     pn = PlotNet(net_labs, ulf_fname, num_stations, dt_start, dt_end)
#     for comp in ['n','e','z']:
#         lab = '_'.join(net_labs)
#         # pn.plot_dist_ts_lags(indices_fname, net_labs, comp, y_space=10, save_plot=False, save_label=f'{lab}_smr_min')
#         pn.plot_dist_ts_lags(indices_fname, comp , pm=2, save_plot=True, save_label=f'{lab}')     

