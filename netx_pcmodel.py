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

# from dynetx.readwrite import json_graph
# import json

# initialising directed and undirected network arrays to store values for all Pc bands

dna = [nx.DiGraph(),nx.DiGraph(),nx.DiGraph(),nx.DiGraph()]

na = [nx.Graph(),nx.Graph(),nx.Graph(),nx.Graph()]

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


def xcorr(x, y, lags, wparr, mode='full'):

    '''function to window time series using a tukey window with window parameter alpha 
    then to perform a time laged cross correlation with normalistation norm1'''
    

    y = signal.tukey(len(y),alpha=wparr)*y
    x = signal.tukey(len(x),alpha=wparr)*x
    y = np.array(y)
    x = np.array(x)

    corr = signal.correlate(x, y, mode=mode, method='fft')

    lnorm1 = len(x)-np.absolute(lags)

    sd1 = np.std(x) 

    sd2 = np.std(y)

    norm1 =  sd1 * sd2 * lnorm1

    # norm2 = (x.std() * y.std() * x.size)

    # print(lnorm1,'lnorm')

    if (sd1 or sd2) == 0:

        return np.array([]) 

    elif len(y)!=len(x):

        return np.array([])
    else:

        return corr/norm1


def closest_to_0(x):

    '''function to find the value closest to 0 in an array x'''

    clv = np.min(np.abs(x))

    # used for finding the sign of the peak +/- 
    cli = np.argmin(np.abs(x))

    clv = np.sign(x[cli])*clv

    return clv


def c0peak(y, x, cutoffh = 0.3):

    '''function to find x value for (maximums) peak in data closest to 0
    use y -> -y for minimums with hight above cutoffh which we set to 0.2 from noise analysis'''

    # gives index position in in xcor array y of peaks satifying conditions 
    maxl = find_peaks(y, height=cutoffh, prominence = 0.1 )[0]
    
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
        ppts = find_peaks(y, height=cutoffh, prominence = 0.1 )[0]

    else:
        ppts = pyaC.zerocross1d(np.arange(len(y)), y, getIndices=False)

    # np.diff gives interpoint spacing array
    s = np.diff(ppts)
    # if s empty return then the reulting xcor from 
    # print(s)

    # average of zero crossing seperation taken and doubled to obtain period
    return 2*np.mean(s)


def network_append( s1name, s2name, md, comp):

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


    s1 = np.array(ts1)
    s2 = np.array(ts2)

    # band-pass parameters 

    fs = 1

    order = 3

    #Pc wave period ranges and labels
    label = ['Pc2','Pc3','Pc4','Pc5']

    tf = [5,10,45,150,600]

    nm = 3

    num_windows = []

    window_size_a = []

    step_size_array = []

    # loop to calculate values for each of four Pc bands

    for i in [0,1,2,3]:

        # minimum window size 2xperiod for xcor

        ws = (tf[i+1] - tf[i])

        window_size = nm*ws

        # initialising time values

        t_start = 0
        t_end = t_start + window_size

        # intial index to get mid slice values for time

        t_mid = t_end//2

        # step_size for different Pc bands

        step_size = (window_size * 3)//4 # amount of window overlap

        # keep empty arrays otside of loop(s) if values needed to be returned outside of function

        min0lags = []

        min0vals = []

        flags = []

        # indices for time values

        t_inda = []

        t_inda.append(t_mid)

        # t_ind array above while loop so can have one additional value before while loop terminates

        step_size_array.append(step_size)

        while t_end <= len(s1):

            counter = t_end//window_size

            # print(counter)

            # y11 , y21 filtering time series to use with TLXC

            y11 = butter_bandpass_filter(s1[t_start:t_end], 1/(tf[i+1]), 1/tf[i], fs, order=order)
            y21 = butter_bandpass_filter(s2[t_start:t_end], 1/(tf[i+1]), 1/tf[i], fs, order=order)

            # need to add code here to take utc time and 

            lags0 = np.arange(-(len(y21) - 1), len(y21))

            rs = xcorr(y11,y21,lags0,wparr=1)

            # gives values of lags at peaks and troughs satifsfying condition in C0peaks, and peak closest to 0
            # returns 'A' if no peaks found

            pp = c0peak(rs,lags0)

            pn = c0peak(-rs,lags0)

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

                t_inda.append(t_mid)

                continue

            # now we can filter for case 'C' wave like which is what is needed

            # lag values at peaks

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

            # ps is average singal period and pxc xcor period for use in flagging

            ps = 0.5 * (fperiod(y11) + fperiod(y21))

            pxc = fperiod(rs)

            # print(ps,pxc)

            # statments to check if data is strictly wave-like 
            # so ps approx. pxc cannot be smaller, nature of xcor)
            # xcor cannot have small period then average period of both signals
            # so only > bound used for range 

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

                continue

            else:

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

            def dict_a2b(ind):
                # function to create dictonary to label network with station names and loop index j

                # utc is the same for both stations however need to check for clusters '2015-03-17T05:18:43' '2015-03-17T05:24:44'

                # strings ending in 1 represent first node and ending in 2 the second node of edge

                d = { 't_window':ind, f'UTC1': utc1.iloc[t_inda[ind]],

                f'UTC2': utc2.iloc[t_inda[ind]], f'MLT1': mlt1.iloc[t_inda[ind]],

                f'MLT2': mlt2.iloc[t_inda[ind]] }

                return d

            def dict_b2a(ind):
                # same as dict_a2b however with reversed time and mlt labels to match directionality of network

                # e.g the first node now had time values from UTC2 and MLT2

                d = { 't_window':ind, f'UTC1': utc2.iloc[t_inda[ind]],

                f'UTC2': utc1.iloc[t_inda[ind]], f'MLT1': mlt2.iloc[t_inda[ind]],

                f'MLT2': mlt1.iloc[t_inda[ind]] }

                return d
            
            if xcval>0 and lag >0:
                dna[i].add_edge(f'{s1name}_{j}', f'{s2name}_{j}', attr_dict = dict_a2b(j))

            elif xcval>0 and lag <0:
                dna[i].add_edge(f'{s2name}_{j}',f'{s1name}_{j}', attr_dict = dict_b2a(j))

            elif xcval<0 and lag <0:
                # print(i,'i','cond3','j',j)
                dna[i].add_edge(f'{s1name}_{j}', f'{s2name}_{j}', attr_dict = dict_a2b(j))

            elif xcval<0 and lag>0:
                dna[i].add_edge(f'{s2name}_{j}',f'{s1name}_{j}', attr_dict = dict_b2a(j))

            # maybe add else: instead of elif lag==0 to speed up code
            elif lag==0:
                na[i].add_edge(f'{s1name}_{j}', f'{s2name}_{j}', attr_dict = dict_a2b(j))


    # gives num of windows in network, last value gets updates over

    # in order to find Pc ranges require only diference values between indices hence one redundant value below arrays

    # t = np.arange(len(min0lags))*step_size

    # # print(window_size,'window_size')

    # print(t, len(t), 'time array')

    # s = np.roll(t_inda,-1) - t_inda

    # print(t_inda, len(t_inda),'mid time array')

    # print( s,'mid time spacing array')

    # print(len(min0lags),'0lags')

    num_windows.append(len(min0lags))

    window_size_a.append(window_size)

    # # can comment out when running large code, always be the same

    # np.save('step_size_arr.npy',step_size_array)


    return num_windows , window_size_a

    

def network_global(data, comp):
    '''function to create both a directed and undirected global netowrk from data for componenet comp n,e or z,
    Geomagnetic coordinates> (magnetic north (n), magnetic east (e) and vertical down (z)) 
    appling the xcor-peak-finding and network appending function, network_append different pairs of stations combinatorially
    for all pairs and returning text files for all directed and undirected netowrks'''

    md = pd.read_csv(data)

    print(md.head())

    print(md.tail())

    # data set with only station names to be used as a label
    s_labels = md.drop_duplicates(subset=['IAGA'])['IAGA']

    # checking all stations for incomplete data 4*3600 in the lenght of the data in seconds, keeping only complete data sets

    filt_stations = [(lab, len(md[md['IAGA'] == lab]['Date_UTC'])/(4*3600)) for lab in s_labels]

    filt_stations = [n1 for n1,n2 in filt_stations if n2==1]

    # scb lists station pair permutations as Nc2 

    scb = list(itertools.combinations(filt_stations,2))

    num_stations = len(filt_stations)

    # for loop over permutation pairs to calculate network connections between station pairs and append to network containers

    for i in range(len(scb)):

        # k is the output from the 2 station peak-finding and netowrk appending algo

        print(scb[i][0],scb[i][1],'station pair out of', num_stations,'stations with', i,'out of',len(scb),'operations', {comp})

        network_append(scb[i][0],scb[i][1],md, comp)

    # saving network below, plotting displays networks, i is the Pc index of interest 0,1,2,3
    # k[0] is the num of windows (which can overlap) (timestamps) for each network k[1] is the window size

    for i in [0,1,2,3]:

        nx.write_edgelist(dna[i], f'networks_data/dna{i}_{comp}_spd.txt')

        nx.write_edgelist(na[i], f'networks_data/na{i}_{comp}_spd.txt')


# time series data from 24/03/2012 starting at 3am and lasting for 4 hours with second resolution containing 7 stations

for i in ['z','n']:

    network_global('20201111-18-30-supermag.csv', i)


# code to save all cluster networks for given component


# test code to try labels

# md = pd.read_csv('20201111-18-30-supermag.csv')

# ts1 = md[md['IAGA'] == 'C10'][f'dbn_nez']
# ts2 = md[md['IAGA'] == 'C11'][f'dbn_nez']

# utc1 = md[md['IAGA'] == 'C10']['Date_UTC']
# utc2 = md[md['IAGA'] == 'C11']['Date_UTC']

# x = 4723

# # print('utc1',utc1.iloc[x-10:x+10])
# # print('utc1',utc2.iloc[x-10:x+10])

# # # utc1, utc2 = np.array(utc1), np.array(utc2)

# # print('utc1',len(utc1),utc1)
# # print('utc2',len(utc2),utc2)

# # 
# # dif =  sorted(np.setdiff1d(utc2,utc1), key= lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))

# # print('DIFFA', len(dif),dif[0],dif[-1],dif)

# # print('ticks',utc1[ticks1])
# # x=8
# # plt.scatter(utc1.iloc[::x],ts1.iloc[::x], s=2)
# # plt.scatter(utc2.iloc[::x],ts2.iloc[::x], s=2)

# s_labels = md.drop_duplicates(subset=['IAGA'])['IAGA']

# k = [(lab, len(md[md['IAGA'] == lab]['Date_UTC'])/(4*3600)) for lab in s_labels]

# print(k)

# remove arrays with missing data and carry on with analysis
# as algorithm cuts data by index value


# plt.show()

# network_cluster('20201111-18-30-supermag.csv','e',cluster='dusk')


# code to be used within if statments to plot wave forms if needed

# fig, ax = plt.subplots(nrows=2, facecolor='w', edgecolor='k')
# fig.subplots_adjust(hspace=0.8)

# ax[0].plot(np.arange(len(y21)),y21, color='red', linestyle='--',label ='time series 1')

# ax[0].plot(np.arange(len(y11)),y11, label ='time series 2')

# ax[0].set_xlabel('Time (s)')

# ax[0].set_ylabel('B (nT)')

# ax[1].plot(lags0,rs)

# ax[1].set_xlabel('Time lags (s)')

# ax[1].set_ylabel('Normalised cross-correlation.')

# ax[1].axvline(x=0,c='red',ls='--')

# ax[1].grid()

# ax[0].grid()

# ax[0].legend()

# plt.show()

