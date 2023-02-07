import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz, find_peaks
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

def plot_wf(s1:list, s2:list ,lags: list, xc: list, px:float, py:float, pcind:int, label:str)->None:
    # code to be used within if statments to plot wave forms if needed
    fig, ax = plt.subplots(nrows=2, facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=0.8)
    fig.suptitle(f'{label}')
    ax[0].plot(np.arange(len(s1)),s1, color='red', linestyle='--',label =f'time series 1, Pc{pcind+2}')
    ax[0].plot(np.arange(len(s2)),s2, label = f'time series 2, Pc{pcind+2}')
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('B (nT)')
    ax[1].plot(lags,xc)
    # ax[1].scatter(px, py, color = 'red')
    ax[1].set_xlabel('Time lags (s)')
    ax[1].set_ylabel('Normalised cross-correlation.')
    ax[1].axvline(x=0,c='red',ls='--')
    ax[1].grid()
    ax[0].grid()
    ax[0].legend()
    # plt.savefig(f'plots/waveform_{label}.png')
    plt.show()
    plt.cla()

def xcorr(x:list, y:list, lags:list, wparr:float, mode='full') -> np.array:
    '''function to window time series using a tukey window with window parameter alpha 
    then to perform a time laged cross correlation with normalistation norm1'''
    
    y = np.array(y-np.mean(y))
    x = np.array(x-np.mean(x))
    y = signal.tukey(len(y),alpha=wparr)*y
    x = signal.tukey(len(x),alpha=wparr)*x
    sd1 = np.std(x) 
    sd2 = np.std(y)
    # norm2 = (x.std() * y.std() * x.size)
    # if std is 0 signal has no power (average amplitude squared)
    if sd1 == 0 or sd2 == 0:
        return np.array([]) 
    elif sd1 == np.nan or sd2 == np.nan:
        return np.array([]) 
    elif len(y)!=len(x):
        return np.array([])
    else:
        corr = signal.correlate(x, y, mode=mode, method='fft')
        lnorm1 = len(x)-np.absolute(lags)
        return corr/(sd1 * sd2 * lnorm1)

def closest_to_0(x:list)->float:
    '''function to find the value closest to 0 in an array x'''
    clv = np.min(np.abs(x))
    # used for finding the sign of the peaks +/- 
    cli = np.argmin(np.abs(x))
    clv = np.sign(x[cli])*clv
    return clv

def c0peaks(y:list, x:list, cutofpeak:float)->list:
    '''function to find x value for (maximums) peaks in data closest to 0
    use y -> -y for minimums with hight above cutofpeak which we set to 0.2 from noise analysis'''
    # gives index position in in xcor array y of peakss satifying conditions 
    maxl = find_peaks(y, height=cutofpeak )[0]
    # if cutofpeak not satisfied and peakss not found then waveform must be noise
    # lablled with flag 'noise'
    if len(maxl)==0:
        return 'noise'
    else:
        clv = closest_to_0(x[maxl])
        return clv, x[maxl]

def fperiod(y: list, cutofpeak=0.25, viapeakss=False)-> np.array:
    '''finds period of discrete point signal y, using intrapolation or to find 0 crossings
     or peaks finding, then finds the average crossing distance and mult. by 2 to find period
     set to True to find period using peakss,however must be two or more values satifying find_peakss
     pyaC is the intrapol. function, first and last zero bypassed'''
    if viapeakss:
        ppts = find_peakss(y, height=cutofpeak, threshold = 0.08 )[0]
    else:
        ppts = pyaC.zerocross1d(np.arange(len(y)), y, getIndices=False)
    # np.diff gives interpoint spacing array
    s = np.diff(ppts)
    # print('fperiod',s)
    if len(s)<2:
        return 0
    else:
        # average of zero crossing seperation taken and doubled to obtain period
        return 2*np.mean(s)

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

def linear_interp(y:np.array, anom_val:float)->np.array:

    nans_bool_arr = y==anom_val
    f_nan_indcies = lambda z: z.nonzero()[0]
    nan_indcies = f_nan_indcies(nans_bool_arr)
    non_nan_indices = f_nan_indcies(~nans_bool_arr)
    y[nans_bool_arr] = np.interp(nan_indcies, non_nan_indices, y[~nans_bool_arr])
    return y

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

def edge_attr(ind:int, direction:str, times:str, lag:int)->dict:
    # function to create dictonary to label network with station names and loop index j
    d = { 't_window':ind, 'UTC': times, 'dir':direction, 'lag':lag}
    return d

# determines directions for directed edges
def dir_edge_assigner(net:nx.Graph, s1:str, s2:str,  j:int, delta:int, time:str, xcval:float, lag:int):
    if xcval>0 and lag > delta:
        net.add_edge(f'{s1}_{j}', f'{s2}_{j}', attr_dict = edge_attr(j,'A-B', time, lag))
    elif xcval>0 and lag < -delta:
        net.add_edge(f'{s2}_{j}', f'{s1}_{j}', attr_dict = edge_attr(j,'B-A', time, lag))
    elif xcval<0 and lag < -delta:
        net.add_edge(f'{s1}_{j}', f'{s2}_{j}', attr_dict = edge_attr(j,'A-B', time, lag))
    elif xcval<0 and lag > delta:
        net.add_edge(f'{s2}_{j}', f'{s1}_{j}', attr_dict = edge_attr(j,'B-A', time, lag))

def network_append(s1name:str, s1:np.array, s2name:str, s2:np.array, times:list, window_mod:list, cutofpeak:float, gapfinder:float, wave_amp_thresh:list) ->None:
    '''function to produce rolling window time lagged cross 
    correleation of two singals with peaks finding routine for signals with labels s1 and s2 
    returning array (of array) of xcor vals, lag at peaks closest to zero (irrespective of 
    aplitude as long as condtions in func c0peaks satisfied) and to append name strings values 
    to a network object for each window
    algorithm filters data to obtain wave-like signals with labels: 'noise', 'non-wave-like' and 'wave-like' for band pass filter 
    applied on each windowsample rate fs, order of butterworth filter order 
    componenet used can be n, e , or z '''
    # provides the whole data set of times given the name of the station, need to change
    # ts (and s1) are the time-series and utc and mlt are the UTC and MLT labels
    # band-pass parameters 
    fs = 1
    order = 3
    delta = 1
    #Pc wave period ranges and labels
    label = ['Pc2','Pc3','Pc4','Pc5']
    pc_period = [5,10,45,150,600]
    # custom window sizes for each band using the lower limit for each pc band period range
    window_size_a = np.multiply(np.diff(pc_period), window_mod)
    #custom step_size for different Pc bands multiplied by amount of window overlap // ensures whole number eg. *3//4 means 25% overlap
    step_size_a = np.zeros(len(window_size_a))
    step_size_a[0:2] = (window_size_a[0:2] * 1) // 2
    step_size_a[2:5] = (window_size_a[2:5] * 1) // 2
    # ensure stepsize is an integer
    step_size_a = step_size_a.astype(np.int32)
    # loop to calculate values for each of four Pc bands
    for i in [0]:
        # minimum window size 2xperiod for xcor
        # first step always whole window size
        t_start = 0
        t_end = t_start + window_size_a[i]
        t_mid = t_end//2
        step_size = step_size_a[i]
        lag_of_extr_close_0_arr = []
        xc_of_extr_close_0_arr = []
        pc_power1 = []
        pc_power2 = []
        # indices for time values
        t_mid_arr = []
        t_mid_arr.append(t_mid)
        while t_end <= len(s1):
            y11 = s1[int(t_start):int(t_end)]
            y21 = s2[int(t_start):int(t_end)]

            counts11, consec_counts11 =  consecutive_vals(y11, pc_period[i]//4, consec_val = gapfinder)
            counts21, consec_counts21 =  consecutive_vals(y21, pc_period[i]//4, consec_val = gapfinder)
            if consec_counts11>pc_period[i]//4 or consec_counts21>pc_period[i]//4:
                t_start = t_start + step_size                
                t_end = t_end + step_size
                t_mid = t_mid + step_size
                t_mid_arr.append(t_mid)
                # needed to keep timing of j indices
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                pc_power1.append(0)
                pc_power2.append(0)
                continue 
            if counts11>0:
                 y11 = linear_interp(y11, anom_val=gapfinder)
            if counts21>0:
                y21 = linear_interp(y21, anom_val=gapfinder)

            y11 = butter_bandpass_filter(y11, 1/(pc_period[i+1]), 1/pc_period[i], fs, order=order)
            y21 = butter_bandpass_filter(y21, 1/(pc_period[i+1]), 1/pc_period[i], fs, order=order)

            # add pcpower
            pc_power11 = np.sum(y11**2)/len(y11)
            pc_power22 = np.sum(y21**2)/len(y21)
            pc_power1.append(pc_power11)
            pc_power2.append(pc_power22)

            lags0 = np.arange(-(len(y21) - 1), len(y21))
            tlxcor = xcorr(y11,y21,lags0, wparr=0.6)
            peaks = c0peaks(tlxcor,lags0, cutofpeak)
            troughs = c0peaks(-tlxcor,lags0, cutofpeak)

            if (troughs == 'noise' and peaks =='noise') or np.max(abs(y11))< wave_amp_thresh[i] or np.max(abs(y21))< wave_amp_thresh[i]:
                
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                # updating time values
                t_start = t_start + step_size                
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                t_mid_arr.append(t_mid)
                continue 

            elif troughs == 'noise' and peaks!='noise':

                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                t_start = t_start + step_size
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                t_mid_arr.append(t_mid)
                continue

            elif troughs != 'noise' and peaks=='noise':
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                t_start = t_start + step_size
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                t_mid_arr.append(t_mid)
                continue

            # statments to check if data is strictly wave-like 
            # so ps approx. pxc cannot be smaller, nature of xcor)
            # xcor cannot have small period then average period of both signals
            # if the period of either signal is zero then cannot find a phase difference
            if fperiod(y11)==0 or fperiod(y21)==0:
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)     
                t_start = t_start + step_size             
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                t_mid_arr.append(t_mid)
                continue

            ps = 0.5 * (fperiod(y11) + fperiod(y21))
            # at this point in the algorithm period of xc pxc will always have a value
            pxc = fperiod(tlxcor)
            if pxc > ps*1.56:

                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                t_start = t_start + step_size
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                t_mid_arr.append(t_mid)
                continue

            lag_peak_closest_0 = peaks[0]
            lag_trough_closest_0 = troughs[0]
            lag_of_extr_close_0 = closest_to_0([lag_peak_closest_0,lag_trough_closest_0])
            xc_extr_close_0 = tlxcor[np.where(lags0==lag_of_extr_close_0)]
            lag_of_extr_close_0_arr.append(lag_of_extr_close_0)
            xc_of_extr_close_0_arr.append(xc_extr_close_0)
            t_start = t_start + step_size
            t_end = t_end + step_size
            t_mid = t_mid +step_size
            t_mid_arr.append(t_mid)

        for j, (xcval, lag) in enumerate(zip(xc_of_extr_close_0_arr, lag_of_extr_close_0_arr)):
            
            t_mid_window =  int(t_mid_arr[j])
            time = str(times.iloc[t_mid_window])   
            # np.nan values should be rejected in below if statments
            # if lags outside of delta range, make directed network
            if lag > delta or lag < -delta:
                # detmerining period of directed edge, within one period wither side of 0 lag
                if lag < pc_period[i+1] and lag > - pc_period[i+1]:
                    dir_edge_assigner(dir_net_t1[i], s1name, s2name, j, delta, time, xcval, lag)
                else:
                    dir_edge_assigner(dir_net_tn[i], s1name, s2name, j,  delta, time, xcval, lag)
            elif xcval>0:
                undir_net_inphase[i].add_edge(f'{s1name}_{j}', f'{s2name}_{j}', attr_dict = edge_attr(j,'NA', time, lag))
            elif xcval<0:
                undir_net_antiphase[i].add_edge(f'{s1name}_{j}', f'{s2name}_{j}', attr_dict = edge_attr(j,'NA', time, lag))

def network_global(save_path:str, comp:str, start_datetime:str, end_datetime:str, wave_amp_cutof:list, tlcc_thresh:float)-> list:
    '''function to create both a directed and undirected global netowrk from data for componenet comp n,e or z,
    Geomagnetic coordinates> (magnetic north (n), magnetic east (e) and vertical down (z)) 
    appling the xcor-peaks-finding and network appending function, network_append different pairs of stations combinatorially
    for all pairs and returning text files for all directed and undirected netowrks'''
    # MLT gives the magnetic longitude (E AND W, relative to IGRF model), MGALAT is the magnetic latitude is the (N and S)
    # initialising directed and undirected network arrays to store values for all Pc band
    print(start_datetime, end_datetime)
    # to locate analysis files
    file_path = f'networks_data/{year}'
    global undir_net_inphase
    global undir_net_antiphase
    global dir_net_t1
    global dir_net_tn

    dir_net_t1 = [nx.DiGraph(),nx.DiGraph()]
    dir_net_tn = [nx.DiGraph(),nx.DiGraph()]
    undir_net_inphase = [nx.Graph(),nx.Graph()]
    undir_net_antiphase =[nx.Graph(),nx.Graph()]

    splt_str= start_datetime.split('-')[0:2]
    year = splt_str[0]
    month_year = '-'.join([splt_str[1],splt_str[0]])

    fpath = f'networks_data/new_events/{month_year}'
    f_xml_read = lambda x: scipy.io.readsav(x, idict=None, python_dict=True, uncompressed_file_name=None, verbose=False)
    station_f_names = [f for f in listdir(fpath)]
    nc2_station_f_names = list(itertools.combinations(station_f_names,2))
    num_stations = len(station_f_names)
    print(num_stations)
    print(station_f_names[0])
    # 6-21-12_2012_A08_1s_final
    time_reference_data = f_xml_read(f'{fpath}/{station_f_names[0]}')
    # print(time_reference_data)
    time_vals, time_indices = time_arr_from_unix(time_reference_data, comp, start_datetime, end_datetime)
    gapidentifier = 999999
    # need to indentify quality of the data
    count =0
    for ind, label in tqdm(enumerate(nc2_station_f_names)):
        # print(ind, label)
        dict_station_0 = f_xml_read(f'{fpath}/{label[0]}')
        dict_station_1 = f_xml_read(f'{fpath}/{label[1]}')
        s0_data = dict_station_0[f'{comp}'][time_indices]
        s1_data = dict_station_1[f'{comp}'][time_indices]
        perc_gaps0 = np.count_nonzero(s0_data == gapidentifier)/len(s0_data)
        perc_gaps1 = np.count_nonzero(s1_data == gapidentifier)/len(s0_data)
        # print('perc of data gaps', perc_gaps)
        if perc_gaps0>0.04 or perc_gaps1>0.04:
            count=+1
            continue
        network_append(label[0].split('_')[2], s0_data, label[1].split('_')[2], s1_data, time_vals, [20], 
            0.3, gapidentifier, wave_amp_cutof)
    print('count',count)
    s_path = f'/home/space/phrmfl/Shared/disk/_Users_sc_conda_envs_geomag/networks_data/{month_year}_nets/comp_{comp}'
    path_exists = os.path.exists(s_path)
    if not path_exists:
      os.makedirs(s_path)
    print(f"The new directory {s_path} is created!")
    lab = f'{num_stations}_{month_year}'
    # need to add new labels for surrogate fuction
    for i in [0]:
        power_lab = f'{wave_amp_cutof[i]}'
        nx.write_edgelist(dir_net_t1[i], f'{s_path}/dir_net_t1{i}_{comp}_{power_lab}_{lab}.txt')
        nx.write_edgelist(dir_net_tn[i], f'{s_path}/dir_net_tn{i}_{comp}_{power_lab}_{lab}.txt')            
        nx.write_edgelist(undir_net_inphase[i], f'{s_path}/undir_net_inphase{i}_{comp}_{power_lab}_{lab}.txt')
        nx.write_edgelist(undir_net_antiphase[i], f'{s_path}/undir_net_antiphase{i}_{comp}_{power_lab}_{lab}txt')


# print('Enter event year, 2012, 2013, 2015:')
# year = input()
# print('Enter event component:')
comp = input()
# comp='n'
# print('Enter event network wave amp cutoff Pc2 (0.2,0.25):')
wave_amp = [0.28]
# s_path = f'/Users/sc/conda_envs/geomag/networks_data/{year}_nets/comp_{comp}'

start_dt = '2012-09-30 08:00:00'
end_dt = '2012-10-01 08:00:00'

network_global(save_path = s_path ,comp= comp, 
    start_datetime=start_dt,end_datetime=end_dt,
    wave_amp_cutof=wave_amp, tlcc_thresh = 0.3)

