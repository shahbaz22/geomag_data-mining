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
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import seaborn as sns




# from dynetx.readwrite import json_graph
# import json

def plot_tslice_dist(net, pcbands, relevent_times, comp, dirr, handel, day_date, delim, title):
    '''function to create degree distrobutions for a given number of plotted overlayed for 
    different times and bands, labeled and colour caded to their respective times'''

    # DONT NEED ANALYSIS CODE FOR SPECIFIC TIME DEGREE DISTROBUTIONS    

    # need to find way of picking relevent UTC times as close to each other as possible

    # take ordered UTC array data from of all bands and find closest value for both for given value of interest.

    # also make sure all values are integers

    # double check dates!!!

    # handel ='spd'

    # day_date = '2015-03-17'

    # delim = 'T'

    # day date as string in format YYYY-MM-DD

    # reformating times into datetime for comutiation, make sure correct delimiter used for times

    # spd day date

    relevent_times = [f'{day_date}{delim}{x}' for x in relevent_times]

    relevent_times = [ datetime.datetime.strptime(x, f'%Y-%m-%d{delim}%H:%M:%S') for x in np.squeeze(relevent_times) ]

    # relevent_times = f'{day_date}{delim}{relevent_times}'

    # relevent_times = [datetime.datetime.strptime(relevent_times, f'%Y-%m-%d{delim}%H:%M:%S')]

    print('relevent times',relevent_times)

    fig, ax = plt.subplots(ncols=len(pcbands), figsize=(14,10)) #----------------------------------------------------

    if dirr == 'dir':
        fig.suptitle(f'directed network {title}')
    
    else:
        fig.suptitle(f'undirected network {title}')

    c=['blue','orange','grey','yellow','brown','red','black']

    for indpc, i in enumerate(pcbands):

        utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in net[i].edges(data=True) ]

        # removing duplicates and sorting in utc_times_dn
        
        times_dn = utc_sort_reduce(utc_times_dn, deliminator =  delim)

        print('times_dn[1]',times_dn[1])

        for ind, target_t in enumerate(relevent_times):

            # time which minimises the timedelta between the chosen time and times presnet in our analysis array

            close_t = min(times_dn, key=lambda d: abs(d - target_t))

            difference_seconds = abs(close_t - target_t).total_seconds()

            # if difference_seconds > 2000:
            #     raise Exception(f'date/time format needs to be checked, day date {day_date}, or target time {target_t} and closest time {close_t} \n with seconds difference {difference_seconds} and time format of net times{times_dn[0]}')


            # close_t back to string
            close_t = close_t.strftime(f'%Y-%m-%d{delim}%H:%M:%S')

            print(f'Pc{i+2} time and target time, dir ', target_t, type(target_t), close_t)

            edgelist = [ [n1,n2] for n1,n2,d in net[i].edges(data=True) if d['attr_dict']['UTC1'] == close_t ]

            reshaped_edgelist = np.reshape(edgelist, len(edgelist)*2)

            nodes = list(set(reshaped_edgelist))

            # counts how many time a node name repeated in edgelist, i.e how many connections it has

            counts = Counter(list(reshaped_edgelist)).most_common()

            # taking degree values for each node and storing then

            degrees_list = [v2 for v1,v2 in [v for v in counts]]

            # counting frequency of degree values giving list of tuples with connections and frequency

            freq = Counter(degrees_list).most_common()

            print('counts',counts)

            print('frequency',freq)

            degree = [v2 for v1,v2 in [v for v in freq]]
            
            freq = [v1 for v1,v2 in [v for v in freq]]

            ax[indpc].bar(degree, freq, label = f'Pc{i+2} {dirr} network, {close_t} UTC', color = c[ind])

            avg_deg = np.average(degree, weights=freq)

            ax[indpc].vlines(np.round(avg_deg,2), ymin = 0, ymax=max(freq), linestyle = '--',color=c[4+ind],label = f'{close_t} UTC, weighted average degree {avg_deg}, {close_t}')

            # set the spacing of points on the x and y axis

            ax[indpc].yaxis.set_major_locator(mticker.MultipleLocator(1))

            ax[indpc].xaxis.set_major_locator(mticker.MultipleLocator(2))

            ax[indpc].grid()
            print(f'Pc{i+2}')

            ax[indpc].set_ylabel('frequency')


        ax[indpc].legend()

        ax[indpc].set_xlabel('degree')

    plt.savefig(f'plots/degree_dist_net_{handel}_{comp}_{title}_pre_post_onset.png')

    # plt.show()

    plt.cla()

def plot_t_heat_dist(net, pcbands, comp, dirr, title):
    '''function to create degree distrobutions for a given number of plotted overlayed for 
    different times and bands, labeled and colour caded to their respective times'''

    # DONT NEED ANALYSIS CODE FOR SPECIFIC TIME DEGREE DISTROBUTIONS    

    # need to find way of picking relevent UTC times as close to each other as possible

    # take ordered UTC array data from of all bands and find closest value for both for given value of interest.




    fig, ax = plt.subplots(ncols=1, figsize=(14,10)) #----------------------------------------------------

    if dirr == 'dir':
        fig.suptitle(f'directed network {title}')
    
    else:
        fig.suptitle(f'undirected network {title}')

    # will get one heatmap for each pc band

    for indpc, i in enumerate(pcbands):

        dictt = {}

        utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in net[i].edges(data=True) ]
    
        for ind, t in enumerate(utc_times_dn):


            edgelist = [ [n1,n2] for n1,n2,d in net[i].edges(data=True) if d['attr_dict']['UTC1'] == t]

            reshaped_edgelist = np.reshape(edgelist, len(edgelist)*2)

            nodes = list(set(reshaped_edgelist))

            # counts how many time a node name repeated in edgelist, i.e how many connections it has

            counts = Counter(list(reshaped_edgelist)).most_common()

            # taking degree values for each node and storing then

            degrees_list = [v2 for v1,v2 in [v for v in counts]]

            # counting frequency of degree values giving list of tuples with connections and frequency i.e [freq, degree]

            freq = Counter(degrees_list).most_common()

            # print(counts)

            # print(freq)

            freqq, degree = zip(*freq)

            # print('counts',counts)

            # print('frequency',freq)

            # degree = [v2 for v1,v2 in [v for v in freq]]
            
            # freq = [v1 for v1,v2 in [v for v in freq]]

            # print(list(degree))

            # print(list(freqq))

            # degree = [v2 for v1,v2 in [v for v in freq]]
            
            # freq = [v1 for v1,v2 in [v for v in freq]]

            # print(degree)

            # print(freq)

            print(t)

            date, time = t.split(' ')

            print(date)

            print(time)

            time = datetime.datetime.strptime(time, '%H:%M:%S').time()

            print(time,type(time))

            d = {time: {'degree': list(degree) ,'freq': list(freqq)} } 

            dictt.update(d)

        df = pd.DataFrame.from_dict(dictt,orient='index' ,columns=['degree','freq'])

        df = df.explode('degree')

        df = df.explode('freq')

        df.index.names = ['time']

        df = df.reset_index()

        df = df.groupby(['time','degree']).freq.first().unstack()

        print(df)


    sns.heatmap(df.T, annot=True, fmt="g", cmap='YlOrBr',  xticklabels = 6, cbar_kws={'label': 'frequency'}, linewidths = 1, linecolor = 'white')


    plt.xticks(rotation=90) 


    plt.show()

    # plt.show()


    print(df)


    # print('00',df['degree'])


    # print('111',df['freq'])

    # print('111',df['date'])

    # print(df.index)
    

    print(ll)
   
    # sb.heatmap(df)

    # plt.show()



def utc_sort_reduce(utc_times, deliminator='T'):

    # takes an unordered list of utc times and removes repeated values and then orders based on tme value

    times_dn = list(set(utc_times))

    times_dn = sorted(times_dn, key = lambda x: datetime.datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S'))

    times_dn = [ datetime.datetime.strptime(x, f'%Y-%m-%d{deliminator}%H:%M:%S') for x in np.squeeze(times_dn) ]

    return times_dn

def save_k_vs_t_global(dn, nn, comp, pcbands):
    
    ''' function to calculate number of connections (k) globally at each time window
    in for all network arrays in dnn and nn to and returns file of k vs t to later use
    and speed up plotting of related graphs'''

    # arrays used for plotting gaphs see cluster function for details

    # magd anf ulf files are magnetometer and 

    delim = ' '

    dict_dn = {}

    dict_n ={}

    for i in ['times', 'n_nodes']:

        dict_dn[i] = {0:[], 1:[], 2:[], 3:[]}

    for i in ['times', 'n_nodes']:

        dict_n[i] = {0:[], 1:[], 2:[], 3:[]}

    for i in ['avg_k']:

        dict_n[i] = {0:[], 1:[], 2:[], 3:[]}

        for j in ['1st', '2nd', '3rd']:

            dict_n[i].update( { j :{ 'MLT':{0:[], 1:[], 2:[], 3:[]}, 'MLAT':{0:[], 1:[], 2:[], 3:[]}, 'e_num' : {0:[], 1:[], 2:[], 3:[]}} } )

    for i in ['avg_k']:

        dict_dn[i] = {0:[], 1:[], 2:[], 3:[]}

        for j in ['1st', '2nd', '3rd']:

            dict_dn[i].update( { j :{ 'MLT':{0:[], 1:[], 2:[], 3:[]}, 'MLAT':{0:[], 1:[], 2:[], 3:[]}, 'e_num' : {0:[], 1:[], 2:[], 3:[]}} } )


    for i in pcbands:

        # obtains only UTC times which have edges

        utc_times_dn = [ d['attr_dict']['UTC2'] for n1,n2,d in dn[i].edges(data=True) ]

        # removing duplicates and sorting in utc_times_dn
        
        times_dn = utc_sort_reduce(utc_times_dn, deliminator=delim)

        dict_dn['times'][i].extend(times_dn)

        utc_times_n = [ d['attr_dict']['UTC2'] for n1,n2,d in nn[i].edges(data=True) ]

        times_n = utc_sort_reduce(utc_times_n, deliminator=delim)

        dict_n['times'][i].extend(times_n)

        # obtains the ordered time stamps in the network

        time_stamps_dn = [d['attr_dict']['t_window'] for n1,n2,d in dn[i].edges(data=True)]

        time_stamps_nn = [d['attr_dict']['t_window'] for n1,n2,d in nn[i].edges(data=True)]

        # remove duplicates
            
        time_stamps_dn = sorted(list(set(time_stamps_dn)))
        
        time_stamps_nn = sorted(list(set(time_stamps_nn)))

        # loop for slicing consecutive time stamps and calculating degree for directed network

        # need to create a function that finds network paramets for given timestamp and computes in order

        # Or stores the value by using enumerate

        for ind, j in enumerate(time_stamps_dn):

            print(f'pc{i+2}, dir global {comp}, {ind} out of {len(time_stamps_dn)}')

            ts_edge_list = [ [n1,n2] for n1,n2,d in dn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

            # print(times_dn[ind], ts_edge_list)

            # count all edges in network at specific time (t_slice_node_list) divided by all unique nodes present, dont use g.degree
            
            # take a unique list of nodes in given time window

            # like double for loop, one loop for acessing tuple value and the next for unpacking it

            reshaped_edge_list = np.reshape(ts_edge_list, len(ts_edge_list)*2)          

            nodes = list(set(reshaped_edge_list))
        
            avg_deg = len(ts_edge_list) / len(nodes)

            dict_dn['avg_k'][i].append(avg_deg)

            dict_dn['n_nodes'][i].append(len(nodes))

        # loop for slicing consecutive time stamps and calculating degree for undirected network

        print(f'pc{i+2}, undir global, {comp}')
            
        for ind, j in enumerate(time_stamps_nn):

            print(f'pc{i+2}, undir global {comp}, {ind} out of {len(time_stamps_nn)}')

            ts_edge_list = [ [n1,n2] for n1,n2,d in nn[i].edges(data=True) if d['attr_dict']['t_window'] == j ]

            # print(times_n[ind], ts_edge_list)

            reshaped_edge_list = np.reshape(ts_edge_list, len(ts_edge_list)*2)          

            nodes = list(set(reshaped_edge_list))
                
            avg_deg = len(ts_edge_list) / len(nodes)

            dict_n['avg_k'][i].append(avg_deg)

            dict_n['n_nodes'][i].append(len(nodes)) 

    return dict_dn , dict_n

       


def plot_k_vs_t(dict_dn, dict_n, Pcbands, num_s, comp, dykey1,  dtimeskey1, uykey1, utimeskey1, plots, label):

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


    # non network plots include SME/SMR, B feilds, SP, all pc ULF power, all pc MLT, all pc MLAT (in this case more ulf band powers)


    num_plots = 2


    fig, ax = plt.subplots(nrows=2, figsize=(14,8))


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


    fig.suptitle(f'{label}')


    for ind, num in enumerate([0,1]):

        n = 0

        # xd = [ datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S') for x in np.squeeze( tk1[num] ) ]

        xd = tk1[num] 

        # print(plots, len(xd), xd, len(dictt[num]), dictt[num])

        if len(xd)!=0:

            markers, stemlines, baseline = ax[ind+n].stem(xd, dictt[num], label=f'Pc{num+2}  average degree')

            plt.setp(markers, color='purple', markersize=7, marker = '.')

            plt.setp(stemlines, color='grey', linestyle='--')

            nodes = np.array(dict_all['n_nodes'][num])

            # print('nodes',nodes)

            # max avg degree for available node

            max_pos_conn =  nodes * (nodes -1) /2

            max_avg_deg = max_pos_conn / nodes

            ax[ind+n].scatter(xd, max_avg_deg, label=f'Pc{num+2} max average degree', marker='_', color='black')

            # theoretical max avg degree

            max_deg = (num_s*(num_s-1))/(2*num_s)

            ax[ind+n].axhline(y=max_deg, xmin=0, xmax=1, linestyle='--',color='red', label = 'max possible avg. degree')        

 
        ax[ind+n].set_ylabel('avg. degree')

        # print(index_date[0],index_date.iloc[-1],type(datel[0]))

        formatter = mdates.DateFormatter("%H:%M:%S")

        ax[ind+n].xaxis.set_major_formatter(formatter)

        ax[ind+n].legend()

        ax[ind+n]

        ax[ind+n].grid()
        

    plt.savefig(f'avg_k_{label}.png')

    # plt.show()

    plt.cla()




def butter_bandpass(lowcut, highcut, fs, order):

    #function to plot butterworth bandpass filter coefficents
    nyq = 0.5 * fs

    low = lowcut / nyq

    high = highcut / nyq
    #band pass uses band freqs. as fraction of nyquist freq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order):

    #function to use bandpass coeffs to window then filter real data 
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    #creates the a digital frequency response and applies it to data
    #using the transfer function via a,b coefs in Z
    y = lfilter(b, a, signal.tukey(len(data))*data)
    return y

def plot_wf(s1,s2,lags, xc, px, py, pcind, label_title, label_save):

     # code to be used within if statments to plot wave forms if needed

     fig, ax = plt.subplots(nrows=2, facecolor='w', edgecolor='k')
     fig.subplots_adjust(hspace=0.8)

     fig.suptitle(label_title)

     ax[0].plot(np.arange(len(s1)),s1, color='red', linestyle='--',label = f'time series 1, Pc{pcind+2}')

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

     # want to have window size, wave period, and peak finding paramenter i.e peak hight

     # do this for a window of increasing size and at each window parameter increase the hight req.

     # also claculate entire network for these changes taking a representative sample of say 20 stations

     # record results in overleaf

     plt.savefig(f'plots/waveform_{label_save}.png')

     # plt.show()

     plt.cla()

def xcorr(x, y, lags, wparr, mode='full'):

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


def closest_to_0(x):

    '''function to find the value closest to 0 in an array x'''

    clv = np.min(np.abs(x))

    # used for finding the sign of the peak +/- 
    cli = np.argmin(np.abs(x))

    clv = np.sign(x[cli])*clv

    return clv


def c0peak(y, x, cutoffh):

    '''function to find x value for (maximums) peak in data closest to 0
    use y -> -y for minimums with hight above cutoffh which we set to 0.2 from noise analysis'''

    # gives index position in in xcor array y of peaks satifying conditions 
    maxl = find_peaks(y, height=cutoffh )[0]
    
    # if cutoffh not satisfied and peaks not found then waveform must be noise
    # lablled with flag 'A' 

    if len(maxl)==0:
        
        return 'A'
    else:

        clv = closest_to_0(x[maxl])

        return clv, x[maxl]

def fperiod(y, cutoffh=0.25, viapeaks=False):

    '''finds period of discrete point signal y, using intrapolation or to find 0 crossings
     or peak finding, then finds the average crossing distance and mult. by 2 to find period
     set to True to find period using peaks,however must be two or more values satifying find_peaks
     pyaC is the intrapol. function, first and last zero bypassed'''


    if viapeaks:
        ppts = find_peaks(y, height=cutoffh, threshold = 0.08 )[0]

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




def network_append( s1name, s2name, md, comp, window_mod, cutoffh ):

    '''function to produce rolling window time lagged cross 
    correleation of two singals with peak finding routine for signals with labels s1 and s2 
    
    returning array (of array) of xcor vals, lag at peak closest to zero (irrespective of 
    aplitude as long as condtions in func c0peak satisfied) and to append name strings values 
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

        # intial index to get mid slice values for time

        t_mid = t_end//2

        # step sizes

        step_size = step_size_a[i]

        counter = 0

        # keep empty arrays otside of loop(s) if values needed to be returned outside of function

        min0lags = []

        min0vals = []

        flags = []

        # indices for time values

        t_inda = []

        t_inda.append(t_mid)

        # t_ind array above while loop so can have one additional value before while loop terminates

        while t_end <= len(s1):

            # counter = t_end//window_size_a[i]

            # print(counter)

            # y11 , y21 filtering time series to use with TLXC

            y11 = butter_bandpass_filter(s1[int(t_start):int(t_end)], 1/(tf[i+1]), 1/tf[i], fs, order=order)
            y21 = butter_bandpass_filter(s2[int(t_start):int(t_end)], 1/(tf[i+1]), 1/tf[i], fs, order=order)

            # need to add code here to take utc time and 

            lags0 = np.arange(-(len(y21) - 1), len(y21))

            rs = xcorr(y11,y21,lags0,wparr=0.6)

            # gives values of lags at peaks and troughs satifsfying condition in C0peaks, and peak closest to 0
            # returns 'A' if no peaks found

            pp = c0peak(rs,lags0, cutoffh)

            pn = c0peak(-rs,lags0, cutoffh)

            # if statment bellow to catogrise noise into flag array for all cases and append nan into other arrays to keep time dimensionality
            # continue statment restarts the while loop rather than exiting out of it with break
            # first if statment, after or, to prevent xcor with a flat line at zero, otherwise period finding function needs modifcation
            # preventing xcor with noise

            # could speed up code by checking second condition of first if then pp and pn 

            if (pn == 'A' and pp =='A') or (np.max(abs(y11))<0.1 or np.max(abs(y21))<0.1) :
                
                min0lags.append(np.nan)

                min0vals.append(np.nan)

                flags.append('A')

                # updating time values

                t_start = t_start + step_size
                
                t_end = t_end + step_size

                t_mid = t_mid +step_size

                counter +=1

                t_inda.append(t_mid)

                continue           

            elif pn == 'A' and pp !='A':

                # uncomment if 'B' non-wave like correlated needed in network

                # min0lags.append(pp[0])

                # min0vals.append(rs[np.where(lags0==pp[0])])

                min0lags.append(np.nan)

                min0vals.append(np.nan)

                flags.append('B')

                t_start = t_start + step_size
                
                t_end = t_end + step_size

                t_mid = t_mid +step_size

                counter +=1

                t_inda.append(t_mid)

                continue

            elif pn != 'A' and pp =='A':

                # uncomment if 'B' non-wave like correlated needed in network

                # min0lags.append(pn[0])

                # min0vals.append(rs[np.where(lags0==pn[0])])

                flags.append('B')

                min0lags.append(np.nan)

                min0vals.append(np.nan)

                t_start = t_start + step_size
                
                t_end = t_end + step_size

                t_mid = t_mid +step_size

                counter +=1

                t_inda.append(t_mid)

                continue

            # statments to check if data is strictly wave-like 
            # so ps approx. pxc cannot be smaller, nature of xcor)
            # xcor cannot have small period then average period of both signals
            # so only > bound used for range 

            # if the period of either signal is zero then cannot find a phase difference

            if fperiod(y11)==0 or fperiod(y21)==0:

                # uncomment if 'B' non-wave like correlated needed in network

                # min0lags.append(extr_0l)

                # min0vals.append(xc0l)

                min0lags.append(np.nan)

                min0vals.append(np.nan)

                flags.append('B')

                t_start = t_start + step_size
                
                t_end = t_end + step_size

                t_mid = t_mid +step_size

                t_inda.append(t_mid)

                counter += 1

                continue

            # ps is average singal period and pxc xcor period for use in flagging

            # if fp(y11) 

            ps = 0.5 * (fperiod(y11) + fperiod(y21))

            # at this point in the algorithm period of xc pxc will always have a value

            pxc = fperiod(rs)

            if pxc > ps*1.56:

                # uncomment if 'B' non-wave like correlated needed in network

                # min0lags.append(extr_0l)

                # min0vals.append(xc0l)

                min0lags.append(np.nan)

                min0vals.append(np.nan)

                flags.append('B')

                t_start = t_start + step_size
                
                t_end = t_end + step_size

                t_mid = t_mid +step_size

                t_inda.append(t_mid)

                counter += 1

                continue

            # now we can find case 'C' wave like which is what is needed

            # lag values at peaks

            # make cuttoff an hour before



            mlagsp = pp[1]

            # gives lag value at peak closest to 0

            mlagp = pp[0]

            # now finding lag values at troughs

            # lag values at troughs

            mlagsn = pn[1]

            # print(mlagsn)

            # gives lag value of trough closest to 0

            mlagn = pn[0]

            # print(mlagn)

            # extremum (peak or trough) lag closest to zero

            extr_0l = closest_to_0([mlagp,mlagn])

            # print(extr_0l)

            # value of xcor at extr_0

            xc0l = rs[np.where(lags0==extr_0l)]

            # print('peak value 1111',xc0l, 'lag closest to 0,', extr_0l)

            # function to plot waveforms -------------

            counter +=1

            # waveform finding code --------------------------------------------

            # change times as required

            # if t_mid < len(utc2):

            #     t1 = datetime.datetime.strptime(utc2.iloc[int(t_mid)],'%Y-%m-%d %H:%M:%S')

            #     t2 = datetime.datetime(2013,3,17,6,30,0)

            # else:

            #     break

            # if t1 > t2:


            #     print(t1, t2, 'times')

            #     lab1 = f'Pc{i+2}_peakh_{cutoffh}_window_period_mult_{window_mod[i]}_post_onset_{utc2.iloc[int(t_mid)]}'

            #     lab2 = f'Pc{i+2}_peakh_{cutoffh}_window_period_mult_{window_mod[i]}_post_onset'


            #     plot_wf(y21, y11, lags0, rs, extr_0l, xc0l, i, label_title = lab1, label_save = lab2 )

            #     t_mid = t_mid +step_size

            #     break


            # t_mid = t_mid +step_size

            # continue





                # if xc0l < 0 and extr_0l == 0:

                #     print('peak value 2222',xc0l, 'lag closest to 0,', extr_0l)
                #     plot_wf(y21, y11, lags0, rs, extr_0l, xc0l, i, label = f'Pc{i+2}_peakh_{cutoffh}_window_period_multiple_{window_mod[i]}_neg_zero')

                #     print(counter, 'counts', f'Pc{i+2}_peakh_{cutoffh}_window_period_multiple_{window_mod[i]}')

                #     break

                # ----------------------------------------
                # print(counter,'peak value 11',xc0l, 'lag closest to 0,', extr_0l)


            min0lags.append(extr_0l)

            min0vals.append(xc0l)

            flags.append('C')

            t_start = t_start + step_size
            
            t_end = t_end + step_size

            t_mid = t_mid +step_size

            t_inda.append(t_mid)

            

            # for loop for parallel interation to create networks

            # np.nan conditionals always return false

        for j, (xcval, lag) in enumerate(zip(min0vals,min0lags)):

            # print(t_inda[j], j, 'tind,j')

            # print(utc1.iloc[t_inda[j]],utc2.iloc[t_inda[j]])

            # print(mlt1.iloc[t_inda[j]],mlt2.iloc[t_inda[j]])

            # both below dictonaries used to assign correct time values to node1 (MLT1) or node two (MLT2)
            # depending on directionality of the network, s.t MLT1 will always have the time for the node1
            # MLT2 will always have the correct value for the node2

            # still useful to have dictonaries attached as it makes filtering easier

            def dict_a2b(ind):
                # function to create dictonary to label network with station names and loop index j

                # utc is the same for both stations however need to check for clusters '2015-03-17T05:18:43' '2015-03-17T05:24:44'

                # strings ending in 1 represent first node and ending in 2 the second node of edge

                tind =  int(t_inda[ind])

                d = { 't_window':ind, f'UTC1': utc1.iloc[tind],

                f'UTC2': utc2.iloc[tind], f'MLT1': mlt1.iloc[tind], f'MLAT1': mlat1.iloc[tind],

                f'MLT2': mlt2.iloc[tind], f'MLAT2': mlat2.iloc[tind] }

                return d

            def dict_b2a(ind):
                # same as dict_a2b however with reversed time and mlt labels to match directionality of network

                # e.g the first node now had time values from UTC2 and MLT2

                tind =  int(t_inda[ind])

                d = { 't_window':ind, f'UTC1': utc2.iloc[tind],

                f'UTC2': utc1.iloc[tind], f'MLT1': mlt2.iloc[tind], f'MLAT1': mlat2.iloc[tind],

                f'MLT2': mlt1.iloc[tind], f'MLAT2': mlat1.iloc[tind] }

                return d

            def lab(station , ind, coordinate1, coordinate2):

                'label for edges in network used as timestamp j, mlt1, mlat1, eg. station_j_mlt1_mlat1'

                tind =  int(t_inda[ind])

                return f'{station}_{ind}_{coordinate1.iloc[tind]}_{coordinate2.iloc[tind]}'


            if xcval>0 and lag >0:
                dna[i].add_edge( lab(s1name,j,mlt1,mlat1) ,lab(s2name, j,mlt2,mlat2) , attr_dict = dict_a2b(j) )

            elif xcval>0 and lag <0:
                dna[i].add_edge( lab(s2name,j,mlt2,mlat2) ,lab(s1name, j,mlt1,mlat1) , attr_dict = dict_b2a(j) )

            elif xcval<0 and lag <0:
                dna[i].add_edge( lab(s1name, j,mlt1,mlat1) ,lab(s2name, j,mlt2,mlat2) , attr_dict = dict_a2b(j) )

            elif xcval<0 and lag>0:
                dna[i].add_edge( lab(s2name,j,mlt2,mlat2) ,lab(s1name,j,mlt1,mlat1) , attr_dict = dict_b2a(j) )

            # maybe add else: instead of elif lag==0 to speed up code
            elif lag==0:
                na[i].add_edge( lab(s1name, j,mlt1,mlat1) ,lab(s2name, j,mlt2,mlat2) , attr_dict = dict_a2b(j) )



    

def network_global(data, comp, times_sec):
    '''function to create both a directed and undirected global netowrk from data for componenet comp n,e or z,
    Geomagnetic coordinates> (magnetic north (n), magnetic east (e) and vertical down (z)) 
    appling the xcor-peak-finding and network appending function, network_append different pairs of stations combinatorially
    for all pairs and returning text files for all directed and undirected netowrks'''

    # MLT gives the magnetic longitude (E AND W, relative to IGRF model), MGALAT is the magnetic latitude is the (N and S)


    md = pd.read_csv(data)

    print(md.head())

    print(md.columns.tolist())

    print(md.tail())

    # data set with only station names to be used as a label
    s_labels = md.drop_duplicates(subset=['IAGA'])['IAGA'][0:4]

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

    # for i in range(len(scb)):

    #     # k is the output from the 2 station peak-finding and netowrk appending algo

    #     print(scb[i][0],scb[i][1],'station pair out of', num_stations,'stations with', i,'out of',len(scb),'operations', comp)

    #     network_append(scb[i][0],scb[i][1],md, comp, [10,10,1.5,1.5], 0.4)

    for w_ind, w_inc in enumerate([10]):

        for peak_height in [0.3]:

            global dna

            global na

            dna = [nx.DiGraph(),nx.DiGraph(),nx.DiGraph(),nx.DiGraph()]

            na = [nx.Graph(),nx.Graph(),nx.Graph(),nx.Graph()]

            for i in range(len(scb)):

                # k is the output from the 2 station peak-finding and netowrk appending algo

                print(scb[i][0],scb[i][1],'station pair out of', num_stations,'stations with', i,'out of',len(scb),'operations', comp)

                network_append(scb[i][0],scb[i][1],md, comp, [w_inc, w_inc, 1.5, 1.5], peak_height)

            # for i in [0,1,2,3]:

            #     nx.write_edgelist(dna[i], f'networks_data/dna{i}_{comp}_170313_params_{w_inc}_{peak_height}.txt')

            #     nx.write_edgelist(na[i], f'networks_data/na{i}_{comp}_170313_params_{w_inc}_{peak_height}.txt')

            # analysis of networks

            # d_an , un_a = save_k_vs_t_global(dna, na, comp, [0,1])

            for ind, dirr in enumerate(['undir', 'dir']):

                # histogram code

                netarr = [dna, na]


                # plot_tslice_dist(['04:30:00','07:30:00'], [0,1], time_arr[ind], comp, 
                #   dirr=dirr, handel = '170313', day_date = '2013-03-17' , delim = ' ', title = f'wind_mult_{w_inc}_peakh_{peak_height}_{dirr}_t_{ind}')
                
                plot_t_heat_dist(netarr[ind], [0,1], comp, dirr=dirr, title = f'wind_mult_{w_inc}_peakh_{peak_height}_{dirr}_t_{ind}')

                # average degree
                # plot_k_vs_t(d_an , un_a, [0,1], num_stations, comp, 'avg_k',  'times', 'avg_k',  'times', dirr, f'window_mult_{w_inc}_peakh_{peak_height}_{dirr}_num_s_{num_stations}')






    # for i in [0,1,2,3]:

    #     nx.write_edgelist(dna[i], f'networks_data/dna{i}_{comp}_170313.txt')

    #     nx.write_edgelist(na[i], f'networks_data/na{i}_{comp}_170313.txt')


# time series data from 24/03/2012 starting at 3am and lasting for 4 hours with second resolution containing 7 stations

for i in ['e']:

    network_global('networks_data/17march2013_4am.csv', i,times_sec = 8*3600)

    # code to reinitilize network arrays

    print(f'done {i} networks')




# run analysis file

# import read_networkx_pcmodel

# code to save all cluster networks for given component

# network_global('20201111-18-30-supermag.csv', 'e',times_sec = 4*3600