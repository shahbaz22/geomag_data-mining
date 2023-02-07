import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz, find_peaks, savgol_filter    
from PyAstronomy import pyaC # for zero crossing intrapolation
import networkx as nx
import matplotlib.pyplot as plt
import itertools
from datetime import datetime
from os import listdir
from os.path import isfile, join
from scipy import io
import time
from tqdm import tqdm
import os
import multiprocessing as mp
import matplotlib.dates as mdates
import matplotlib.ticker as mticker


def butter_bandpass(lowcut:float, highcut:float, fs:int, order:int)-> list:
    #function to plot butterworth bandpass filter coefficents
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    #band pass uses band freqs. as fraction of nyquist freq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data:list, lowcut:float, highcut:float, fs:int, order:int) -> list:
    #function to use bandpass coeffs to window then filter real data 
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    #creates the a digital frequency response and applies it to data
    #using the transfer function via a,b coefs in Z
    y = lfilter(b, a, signal.tukey(len(data))*data)
    return y


def consecutive_vals(y:list, num_nans_cutof, consec_val:float)->int:
    counts = 0
    consec_counts = 0
    for i in range(len(y)-1):
        if y[i] == consec_val:
            counts+=1
            if y[i]==y[i+1]:
                consec_counts+=1
        elif counts>num_nans_cutof:
            return counts, consec_counts
    return counts, consec_counts

def time_arr_from_unix(station_dict:dict ,comp:str, start_datetime:str, end_datetime:str)->list:
    '''retrun array of indices and times from unix start and end times from the begining of the year in seconds '''
    year = start_datetime.split('-')[0]
    unix_year = time.mktime(datetime.strptime(year, "%Y").timetuple())
    unix_start_time = unix_year + station_dict['begin_time'] 
    unix_end_time = unix_year + station_dict['begin_time'] + len(station_dict['n']) -1
    start_time_str = datetime.utcfromtimestamp(unix_start_time).strftime('%Y-%m-%d %H:%M:%S') 
    end_time_str = datetime.utcfromtimestamp(unix_end_time).strftime('%Y-%m-%d %H:%M:%S')
    all_times = pd.Series(pd.date_range(start_time_str, end_time_str, len(station_dict['n'])))
    relevent_times = all_times[all_times.between(start_datetime, end_datetime)]
    relevent_indices = relevent_times.index.tolist()
    return relevent_times, relevent_indices


def times(s1:np.array, times:list, window_mod:int) ->list:
    pc_period = [5,10,45,150,600]
    window_size_a = np.multiply(np.diff(pc_period), window_mod)
    step_size_a = np.zeros(len(window_size_a))
    step_size_a[0:2] = (window_size_a[0:2] * 1) // 2
    step_size_a = step_size_a.astype(np.int32)
    pc_ind = 0 
    step_size = step_size_a[pc_ind]
    times_arr=[]
    t_start = 0
    t_end = t_start + window_size_a[pc_ind]
    t_mid = t_end//2
    while t_end <= len(s1):
        t_start = t_start + step_size                
        t_end = t_end + step_size
        t_mid = t_mid + step_size
        time = times.iloc[t_mid]
        time = times.iloc[t_mid]   
        times_arr.append(time)
    return times_arr
    # return times_arr

def get_result(result):
    global results
    results.append(result)

def power_arr_append(s1:np.array, window_mod:list, gapfinder:float) ->None:
    # band-pass parameters 
    fs = 1
    order = 3
    pc_power_arr = []
    pc_period = [5,10,45,150,600]
    # custom window sizes for each band using the lower limit for each pc band period range
    window_size_a = np.multiply(np.diff(pc_period), window_mod)
    #custom step_size for different Pc bands multiplied by amount of window overlap // ensures whole number eg. *3//4 means 25% overlap
    step_size_a = np.zeros(len(window_size_a))
    step_size_a[0:2] = (window_size_a[0:2] * 1) // 2
    # ensure stepsize is an integer
    step_size_a = step_size_a.astype(np.int32)
    # loop to calculate values for each of four Pc bands
    pc_ind = 0 
    t_start = 0
    t_end = t_start + window_size_a[pc_ind]
    step_size = step_size_a[pc_ind]
    while t_end <= len(s1):
        y11 = s1[int(t_start):int(t_end)]
        counts11, consec_counts11 =  consecutive_vals(y11, pc_period[pc_ind]//4, consec_val = gapfinder)
        if consec_counts11>pc_period[pc_ind]//4: 
            t_start = t_start + step_size                
            t_end = t_end + step_size
            pc_power_arr.append(np.nan)
            continue 

        y11 = butter_bandpass_filter(y11, 1/(pc_period[pc_ind+1]), 1/pc_period[pc_ind], fs, order=order)
        pc_power11 = np.sum(y11**2)/len(y11)
        pc_power_arr.append(pc_power11)
        t_start = t_start + step_size                
        t_end = t_end + step_size

    return pc_power_arr
    # return times_arr

if __name__ == '__main__':
    comp='n'
    # save_path = f'/Users/sc/conda_envs/geomag/networks_data/{year}_nets/comp_{comp}'

    ind_f_names = [f for f in listdir('indices/')]
    print(ind_f_names)
    f = ind_f_names[5]
    # also need 10-2012 event
    fig, ax = plt.subplots(4)
    fig.set_figwidth(20)
    fig.set_figheight(8)
    indices = pd.read_csv(f'indices/{f}')
    index_date = indices['Date_UTC']
    print('filename', f)
    print('start_date from read-in index file', index_date.iloc[0])
    print('start_date from read-in index file', index_date.iloc[-1])
    date =  index_date[0].split(' ')[0]
    axsmer = pd.to_datetime(index_date)
    plt.suptitle(f)
    # copy and paste dates from print statments given in file.
    start_datetime = index_date.iloc[0]
    # end_datetime = '2012-10-09 23:11:00'

    end_datetime = index_date.iloc[-1]
    
    if date != start_datetime.split(' ')[0]:
        raise ValueError('dates not the same')

    splt_str= start_datetime.split('-')[0:2]
    year = splt_str[0]
    month_year = '-'.join([splt_str[1],splt_str[0]])
    print(month_year)

    datemin = datetime.strptime(start_datetime, '%Y-%m-%d %H:%M:%S')
    datemax = datetime.strptime(end_datetime, '%Y-%m-%d %H:%M:%S')

    formatter = mdates.DateFormatter("%d/%H:%M")

    for ind, lab in enumerate(['GSM_Bz','GSM_By','PDYN']):
        c = ['red','blue','black','green']
        if ind == 2:
            ax2 = ax[0].twinx()
            ax2.set_ylabel('nPa')
            # get rid of peaks
            y = np.where( indices[lab] == np.max(indices[lab]) , np.nan , indices[lab])
            ax2.plot(axsmer, y, label=lab, color = c[ind])
        else:
            y = np.where( indices[lab] == np.max(indices[lab]) , np.nan , indices[lab])
            ax[0].plot(axsmer, y, label=lab, color = c[ind])

    ax[0].xaxis.set_major_formatter(formatter)
    ax[0].set_ylabel('nT')
    # ax[0].set_xlim(datemin,datemax)
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax[0].get_legend_handles_labels()
    ax2.legend(h1+h2, l1+l2,bbox_to_anchor=[0.47, 1.15],ncol=4,loc='upper center')
    ax[0].grid()
    ax[0].set_xticklabels([])
    ax[0].set_xlim(datemin,datemax)

    #     code for plotting SMR sections for entire SMR yearly dataset
    sme = indices['SME']
    ax[1].plot(axsmer,sme, color='black', label = 'SME')
    ax[1].set_ylabel('SME (nT)')
    ax[1].grid()
    #   turn off labels for x axis while preserving ticks
    ax[1].set_xticklabels([]) 
    ax[1].set_xlim(datemin,datemax)

    # ax[1].xaxis.set_major_formatter(formatter)
    ax3 = ax[1].twinx()
    indices['Date_UTC'] = pd.to_datetime(indices['Date_UTC']) 
    for i in ['SMR00','SMR06','SMR12','SMR18' ]:
        x = axsmer #[(indices['Date_UTC'] >= datemin) & (indices['Date_UTC'] <= datemax)]
        y = indices[i] #[(indices['Date_UTC'] >= datemin) & (indices['Date_UTC'] <= datemax)][i]
        ax3.plot(x , y, label= i)
    ax3.xaxis.set_major_formatter(formatter)
    ax[1].legend(bbox_to_anchor=[0.1,1.1])
    ax3.legend(bbox_to_anchor=[0.7, 1.1],ncol=4,loc='upper center')
    ax3.set_ylabel('SMR (nT)')
    ax2.set_xlim(datemin,datemax)                          
    ax[1].set_xlim(datemin,datemax)
    ax3.set_xlim(datemin,datemax)
    ax3.set_xticklabels([])
    # plt.savefig(f'plots/{f}.png')

    # Pc2 power calculator------------------------
    results = []
    fpath = f'networks_data/new_events/{month_year}'
    f_xml_read = lambda x: scipy.io.readsav(x, idict=None, python_dict=True, uncompressed_file_name=None, verbose=False)
    station_f_names = [f for f in listdir(fpath)]
    num_stations = len(station_f_names)
    print(num_stations)
    print(station_f_names[0])
    # 6-21-12_2012_A08_1s_final
    time_reference_data = f_xml_read(f'{fpath}/{station_f_names[0]}')
    # print(time_reference_data)
    time_vals, time_indices = time_arr_from_unix(time_reference_data, comp, start_datetime, end_datetime)
    gapidentifier = 999999
    # need to indentify quality of the data
    bool_arr = []
    for ind, label in enumerate(station_f_names):
        # print(ind, label)
        dict_station_0 = f_xml_read(f'{fpath}/{label}')
        s0_data = dict_station_0[f'{comp}'][time_indices]
        data_gaps = (s0_data == gapidentifier)
        bool_arr.append(data_gaps)
    gaps_ts = np.count_nonzero(bool_arr, axis=0)/len(station_f_names)
    ax[2].plot(time_vals, gaps_ts, color = 'purple')
    ax[2].xaxis.set_major_formatter(formatter)
    ax[2].set_xlim(datemin,datemax)
    ax[2].set_ylabel('num. data gaps')


    pool = mp.Pool(mp.cpu_count())
    count =0
    for ind, label in enumerate(station_f_names):
        # print(ind, label)
        dict_station_0 = f_xml_read(f'{fpath}/{label}')
        s0_data = dict_station_0[f'{comp}'][time_indices]
        perc_gaps = np.count_nonzero(s0_data == gapidentifier)/len(s0_data)
        # print('perc of data gaps', perc_gaps)
        if perc_gaps>0.001:
            count=+1
            continue
        s0_data = dict_station_0[f'{comp}'][time_indices]
        pool.apply_async(power_arr_append, args=(s0_data, 10, gapidentifier), callback=get_result)
    pool.close()
    pool.join()
    average_power = np.log10(np.nanmean(results, axis=0))
    # smoothing function
    average_power = scipy.signal.savgol_filter(average_power,  window_length=70, polyorder=2)
    print('minimum average power', min(average_power))
    times = times(s0_data, time_vals, 10)
    print(len(station_f_names)-count)

    ax[3].plot(times, average_power, color = 'blue',label= 'pc2 mean power')
    ax[3].xaxis.set_major_formatter(formatter)
    ax[3].set_xlim(datemin,datemax)
    ax[3].set_ylabel('log(nT^2)')
    # ax[3].set_ylim(-1.2,1.5)
    ax[3].legend(bbox_to_anchor=[0.8, 0.45])
    plt.show()


# code for supermag ULF power
# datemin = '2012-09-30 00:00:00'
# datemax = '2012-09-30 23:59:59'
# def ulf_power(band, ulf_fname, dt_start, dt_end):
#     '''mean ulf powers for pc single band for corespondig times'''
#     ulf_powers = pd.read_csv(ulf_fname)
#     print(ulf_powers['Date_UTC'].iloc[0],ulf_powers['Date_UTC'].iloc[-1])
#     ulf_powers['Date_UTC'] = pd.to_datetime(ulf_powers['Date_UTC'])
#     dict_ulf ={'mean_pc_power':[]}
#     power_times = pd.date_range(dt_start, dt_end, freq='30s')
#     for i in np.squeeze(power_times):
#         ts = ulf_powers[ulf_powers['Date_UTC'] == i][[ 'Date_UTC', 'IAGA','PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']] 
#         powers = ['PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']
#         dict_ulf['mean_pc_power'].append(ts[powers[band]].mean())
#     return dict_ulf, power_times
# ulf_fname = f'networks_data/30_09_12supermag.csv'
# pc_power, power_times = ulf_power(0, ulf_fname, datemin, datemax)
# plt.plot(power_times, pc_power['mean_pc_power'], color = 'black',label= f'pc{2} mean power')
# plt.ylabel('log(nT^2)')
# plt.legend(bbox_to_anchor=[0.8, 0.45])
# plt.show()
