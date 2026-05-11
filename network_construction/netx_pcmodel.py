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
from collections import Counter
import datetime
import time


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

def plot_wf(s1:str, s2:str ,lags: list, xc: list, px:float, py:float, pcind:int, label:str)->None:
    # code to be used within if statments to plot wave forms if needed
    fig, ax = plt.subplots(nrows=2, facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=0.8)
    fig.suptitle(f'{label}')
    ax[0].plot(np.arange(len(s1)),s1, color='red', linestyle='--',label =f'time series 1, Pc{pcind+2}')
    ax[0].plot(np.arange(len(s2)),s2, label = f'time series 2, Pc{pcind+2}')
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('B (nT)')
    ax[1].plot(lags,xc)
    ax[1].scatter(px, py, color = 'red')
    ax[1].set_xlabel('Time lags (s)')
    ax[1].set_ylabel('Normalised cross-correlation.')
    ax[1].axvline(x=0,c='red',ls='--')
    ax[1].grid()
    ax[0].grid()
    ax[0].legend()
    plt.savefig(f'plots/waveform_{label}.png')
    # plt.show()
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

def c0peaks(y:list, x:list, cutoffh:float)->list:
    '''function to find x value for (maximums) peaks in data closest to 0
    use y -> -y for minimums with hight above cutoffh which we set to 0.2 from noise analysis'''
    # gives index position in in xcor array y of peakss satifying conditions 
    maxl = find_peaks(y, height=cutoffh )[0]
    # if cutoffh not satisfied and peakss not found then waveform must be noise
    # lablled with flag 'A'
    if len(maxl)==0:
        return 'A'
    else:
        clv = closest_to_0(x[maxl])
        return clv, x[maxl]

def fperiod(y: list, cutoffh=0.25, viapeakss=False)-> np.array:
    '''finds period of discrete point signal y, using intrapolation or to find 0 crossings
     or peaks finding, then finds the average crossing distance and mult. by 2 to find period
     set to True to find period using peakss,however must be two or more values satifying find_peakss
     pyaC is the intrapol. function, first and last zero bypassed'''
    if viapeakss:
        ppts = find_peakss(y, height=cutoffh, threshold = 0.08 )[0]
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

def network_append(s1name:str, s2name:str, md:pd.DataFrame, comp:str, window_mod:int, cutoffh:float) ->None:
    '''function to produce rolling window time lagged cross 
    correleation of two singals with peaks finding routine for signals with labels s1 and s2 
    
    returning array (of array) of xcor vals, lag at peaks closest to zero (irrespective of 
    aplitude as long as condtions in func c0peaks satisfied) and to append name strings values 
    to a network object for each window
    algorithm filters data to obtain wave-like signals with labels: 'A', 'B' and 'C'
    for noise, non-wave-like correlated and wavelike for band pass filter 
    applied on each windowsample rate fs, order of butterworth filter order 
    componenet used can be n, e , or z '''
    # provides the whole data set of times given the name of the station, need to change
    # ts (and s1) are the time-series and utc and mlt are the UTC and MLT labels

    ts1 = md[md['IAGA'] == s1name][f'db{comp}_nez']
    ts2 = md[md['IAGA'] == s2name][f'db{comp}_nez']

    utc1 =  md[md['IAGA'] == s1name]['Date_UTC']
    utc2 =  md[md['IAGA'] == s2name]['Date_UTC']

    mlt1 =  md[md['IAGA'] == s1name]['MLT']
    mlt2 =  md[md['IAGA'] == s2name]['MLT']

    mlat1 =  md[md['IAGA'] == s1name]['MAGLAT']
    mlat2 =  md[md['IAGA'] == s2name]['MAGLAT']

    pc2_power1 =  md[md['IAGA'] == s1name]['PC2_IPOW']
    pc2_power2 =  md[md['IAGA'] == s2name]['PC2_IPOW']

    pc3_power1 =  md[md['IAGA'] == s1name]['PC3_IPOW']
    pc3_power2 =  md[md['IAGA'] == s2name]['PC3_IPOW']

    s1 = np.array(ts1)
    s2 = np.array(ts2)
    # band-pass parameters 
    fs = 1
    order = 3
    #Pc wave period ranges and labels
    label = ['Pc2','Pc3','Pc4','Pc5']
    tf = [5,10,45,150,600]
    ws = np.diff(tf)
    # custom window sizes for each band
    window_size_a = np.multiply(ws, window_mod)
    #custom step_size for different Pc bands multiplied by amount of window overlap // ensures whole number eg. *3//4 means 25% overlap
    step_size_a = np.zeros(len(window_size_a))
    step_size_a[0:2] = window_size_a[0:2] * 1 // 2
    step_size_a[2:5] = window_size_a[2:5] * 1 // 2
    # ensure stepsize is an integer
    step_size_a = step_size_a.astype(np.int32)
    # loop to calculate values for each of four Pc bands
    for i in [0,1]:
        # minimum window size 2xperiod for xcor
        #  first step always whole window size
        t_start = 0
        t_end = t_start + window_size_a[i]
        t_mid = t_end//2
        step_size = step_size_a[i]
        counter = 0
        lag_of_extr_close_0_arr = []
        xc_of_extr_close_0_arr = []
        flags = []
        # indices for time values
        t_mid_arr = []
        t_mid_arr.append(t_mid)
        while t_end <= len(s1):
            y11 = butter_bandpass_filter(s1[int(t_start):int(t_end)], 1/(tf[i+1]), 1/tf[i], fs, order=order)
            y21 = butter_bandpass_filter(s2[int(t_start):int(t_end)], 1/(tf[i+1]), 1/tf[i], fs, order=order)
            lags0 = np.arange(-(len(y21) - 1), len(y21))
            tlxcor = xcorr(y11,y21,lags0,wparr=0.6)
            peaks = c0peaks(tlxcor,lags0, cutoffh)
            troughs = c0peaks(-tlxcor,lags0, cutoffh)

            if (troughs == 'A' and peaks =='A') or (np.max(abs(y11))<0.1 or np.max(abs(y21))<0.1) :
                
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                flags.append('A')
                # updating time values
                t_start = t_start + step_size                
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                counter = +1
                t_mid_arr.append(t_mid)
                continue 

            elif troughs == 'A' and peaks!='A':
                # uncomment if 'B' non-wave like correlated needed in network
                # lag_of_extr_close_0_arr.append(peaks[0])
                # xc_of_extr_close_0_arr.append(rs[np.where(lags0==peaks[0])])
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                flags.append('B')
                t_start = t_start + step_size
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                counter = +1
                t_mid_arr.append(t_mid)
                continue

            elif troughs != 'A' and peaks=='A':
                # uncomment if 'B' non-wave like correlated needed in network
                # lag_of_extr_close_0_arr.append(troughs[0])
                # xc_of_extr_close_0_arr.append(rs[np.where(lags0==troughs[0])])
                flags.append('B')
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                t_start = t_start + step_size
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                counter = +1
                t_mid_arr.append(t_mid)
                continue

            # statments to check if data is strictly wave-like 
            # so ps approx. pxc cannot be smaller, nature of xcor)
            # xcor cannot have small period then average period of both signals
            # so only > bound used for range 
            # if the period of either signal is zero then cannot find a phase difference
            if fperiod(y11)==0 or fperiod(y21)==0:
                # uncomment if 'B' non-wave like correlated needed in network
                # lag_of_extr_close_0_arr.append(extr_0l)
                # xc_of_extr_close_0_arr.append(xc0l)
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                flags.append('B')
                t_start = t_start + step_size             
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                t_mid_arr.append(t_mid)
                counter = +1
                continue

            # ps is average singal period and pxc xcor period for use in flagging
            # if fp(y11) 
            ps = 0.5 * (fperiod(y11) + fperiod(y21))
            # at this point in the algorithm period of xc pxc will always have a value
            pxc = fperiod(tlxcor)
            if pxc > ps*1.56:
                # uncomment if 'B' non-wave like correlated needed in network
                # lag_of_extr_close_0_arr.append(extr_0l)
                # xc_of_extr_close_0_arr.append(xc0l)
                lag_of_extr_close_0_arr.append(np.nan)
                xc_of_extr_close_0_arr.append(np.nan)
                flags.append('B')
                t_start = t_start + step_size
                t_end = t_end + step_size
                t_mid = t_mid +step_size
                t_mid_arr.append(t_mid)
                counter = +1
                continue

            # now we can find case 'C' wave like which is what is needed
            # gives lag value at peaks closest to 0
            lag_peak_closest_0 = peaks[0]
            lag_trough_closest_0 = troughs[0]
            lag_of_extr_close_0 = closest_to_0([lag_peak_closest_0,lag_trough_closest_0])
            xc_extr_close_0 = tlxcor[np.where(lags0==lag_of_extr_close_0)]
            lag_of_extr_close_0_arr.append(lag_of_extr_close_0)
            xc_of_extr_close_0_arr.append(xc_extr_close_0)
            flags.append('C')
            t_start = t_start + step_size
            t_end = t_end + step_size
            t_mid = t_mid +step_size
            t_mid_arr.append(t_mid)
            # for loop for parallel interation to create networks
            # np.nan conditionals always return false
        for j, (xcval, lag) in enumerate(zip(xc_of_extr_close_0_arr, lag_of_extr_close_0_arr)):
            # print(t_mid_arr[j], j, 'tind,j')
            # print(utc1.iloc[t_mid_arr[j]],utc2.iloc[t_mid_arr[j]])
            # print(mlt1.iloc[t_mid_arr[j]],mlt2.iloc[t_mid_arr[j]])
            def dict_a2b(ind):
                # function to create dictonary to label network with station names and loop index j
                tind =  int(t_mid_arr[ind])
                d = { 't_window':ind, f'UTC1': utc1.iloc[tind],
                f'UTC2': utc2.iloc[tind], f'MLT1': mlt1.iloc[tind], f'MLAT1': mlat1.iloc[tind],
                f'MLT2': mlt2.iloc[tind], f'MLAT2': mlat2.iloc[tind] }
                return d
            def dict_b2a(ind):
                # same as dict_a2b however with reversed time and mlt labels to match directionality of network
                # e.g the first node now had time values from UTC2 and MLT2
                tind =  int(t_mid_arr[ind])
                d = { 't_window':ind, f'UTC1': utc2.iloc[tind],
                f'UTC2': utc1.iloc[tind], f'MLT1': mlt2.iloc[tind], f'MLAT1': mlat2.iloc[tind],
                f'MLT2': mlt1.iloc[tind], f'MLAT2': mlat1.iloc[tind] }
                return d
            def lab(station , ind, value1, value2):
                'label for edges in network used as timestamp j, mlt1, mlat1, eg. station_j_mlt1_mlat1'
                tind =  int(t_mid_arr[ind])
                return f'{station}_{ind}_{value1.iloc[tind]}_{value2.iloc[tind]}'
            
            if xcval>0 and lag >0:
                dir_net[i].add_edge(lab(s1name,j,pc2_power1, pc3_power1), lab(s2name, j,pc2_power2, pc3_power2) , attr_dict = dict_a2b(j) )
            elif xcval>0 and lag <0:
                dir_net[i].add_edge(lab(s2name,j,pc2_power2, pc3_power2), lab(s1name, j,pc2_power1, pc3_power1) , attr_dict = dict_b2a(j) )
            elif xcval<0 and lag <0:
                dir_net[i].add_edge( lab(s1name, j,pc2_power1, pc3_power1), lab(s2name, j,pc2_power2, pc3_power2) , attr_dict = dict_a2b(j) )
            elif xcval<0 and lag>0:
                dir_net[i].add_edge(lab(s2name,j,pc2_power2, pc3_power2), lab(s1name,j,pc2_power1, pc3_power1) , attr_dict = dict_b2a(j) )
            elif xcval>0 and lag==0:
                undir_net_inphase[i].add_edge(lab(s1name, j,pc2_power1,pc3_power1), lab(s2name, j,pc2_power2,pc3_power2) , attr_dict = dict_a2b(j))
            elif xcval<0 and lag==0:
                undir_net_antiphase[i].add_edge(lab(s1name, j,pc2_power1,pc3_power1), lab(s2name, j,pc2_power2,pc3_power2) , attr_dict = dict_a2b(j))
    
def network_global(data:pd.DataFrame, comp:str, times_sec:int)-> list:
    '''function to create both a directed and undirected global netowrk from data for componenet comp n,e or z,
    Geomagnetic coordinates> (magnetic north (n), magnetic east (e) and vertical down (z)) 
    appling the xcor-peaks-finding and network appending function, network_append different pairs of stations combinatorially
    for all pairs and returning text files for all directed and undirected netowrks'''
    # MLT gives the magnetic longitude (E AND W, relative to IGRF model), MGALAT is the magnetic latitude is the (N and S)
    # initialising directed and undirected network arrays to store values for all Pc bands
    global dir_net
    global undir_net_inphase
    global undir_net_antiphase
    dir_net = [nx.DiGraph(),nx.DiGraph(),nx.DiGraph(),nx.DiGraph()]
    undir_net_inphase = [nx.Graph(),nx.Graph(),nx.Graph(),nx.Graph()]
    undir_net_antiphase =[nx.Graph(),nx.Graph(),nx.Graph(),nx.Graph()]

    md = pd.read_csv(data)
    print(md.head())
    print(md.columns.tolist())
    print(md.tail())
    # data set with only station names to be used as a label
    s_labels = md.drop_duplicates(subset=['IAGA'])['IAGA']
    # checking all stations for incomplete data 4*3600 in the lenght of the data in seconds, keeping only complete data sets
    filt_stations = [(lab, len( md[md['IAGA'] == lab]['Date_UTC'] ) / (times_sec)) for lab in s_labels]
    filt_stations = [n1 for n1,n2 in filt_stations if n2==1]
    print(filt_stations, len(filt_stations))
    # scb lists station pair permutations as Nc2 
    scb = list(itertools.combinations(filt_stations,2))
    print('scb',scb,len(scb))
    num_stations = len(filt_stations)
    # for loop over permutation pairs to calculate network connections between station pairs and append to network containers
    s = len(scb)
    for w_inc in [20]:
        for i in range(len(scb)):
            # k is the output from the 2 station peaks-finding and netowrk appending algo
            print(scb[i][0],scb[i][1],'station pair out of', num_stations,'stations with', i,'out of',len(scb),'operations', comp)
            network_append(scb[i][0],scb[i][1],md, comp, [w_inc, w_inc,1.5,1.5], 0.3)
    for i in [0,1,2,3]:
        nx.write_edgelist(dir_net[i], f'networks_data/dir_net{i}_{comp}_{w_inc}_peakh_0.3_num_stations_{num_stations}_170313.txt')
        nx.write_edgelist(undir_net_inphase[i], f'networks_data/undir_net_inphase{i}_{comp}_{w_inc}_peakh_0.3_num_stations_{num_stations}_170313.txt')
        nx.write_edgelist(undir_net_antiphase[i], f'networks_data/undir_net_antiphase{i}_{comp}_{w_inc}_peakh_0.3_num_stations_{num_stations}_170313.txt')


# time series data from 24/03/2012 starting at 3am and lasting for 4 hours with second resolution containing 7 stations
for i in ['e','n','z']:
    network_global('networks_data/17march2013_4am.csv', i, times_sec = 8*3600)
    print(f'done {i} networks')
# run analysis file
# import read_networkx_pcmodel
# code to save all cluster networks for given component
# network_global('20201111-18-30-supermag.csv', 'e',times_sec = 4*3600)
