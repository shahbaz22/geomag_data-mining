import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz, find_peaks
# from geopacktest import apdvalue
from datetime import timedelta
import time
from PyAstronomy import pyaC # for zero crossing intrapolation
import networkx as nx
import matplotlib.pyplot as plt
import dynetx as dn
import itertools


# initialising directed and undirected network arrays for all Pc bands

dna = [dn.DynDiGraph(),dn.DynDiGraph(),dn.DynDiGraph(),dn.DynDiGraph()]

na = [dn.DynGraph(),dn.DynGraph(),dn.DynGraph(),dn.DynGraph()]

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

def apexd(lat, lon, date):
    '''function to return apex distance and convert dateimte to unix datetime UTC for datetime in supermag datetime 
    format, lon is the longitude and date value needed for 
    conversion to calculate apex distance for magnetic fields using geopack (function apdvalue)
    use unix date for apd values then compare to local time
    date foramt: YYYY/MM/DD hh/mm/ss
    input: unix utc date, longitude
    return: apex distance in Re multiple'''

    t0 = np.datetime64(str(date))

    t0 = t0.astype('datetime64[s]').astype('int')

    d0 = np.round(apdvalue(lat, lon, int(t0)),decimals=1)

    return d0

def dateftutc(ind, date):
    '''function to filter UTC dates to calculate apexd values for plotting'''

    mask = date.index.isin(ind)

    date = date[mask]

    return date

def dateflocal(lon, date, ind):
    ''' fuction to convert datetime utc formal to local time 
   selecting values from ind array to be used for plotting'''

    mask = date.index.isin(ind)

    date = date[mask]

    date = pd.to_datetime(date, format='%Y/%m/%d %H:%M:%S')

    if lon > 180:

        date = date - timedelta(hours=np.round((360-lon)/15))

    else:

        date = date + timedelta(hours=np.round(lon/15))

    # date = pd.to_datetime(date, format='%Y/%m/%d %H:%M:%S') - timedelta(hours=6.1)

    date = date.astype(str)

    date = date.str.extract(f'(\d\d:\d\d:\d\d)', expand=True)

    date = np.squeeze(date.values)

    return date

def xcorr(x, y, lags, wparr, mode='full'):
    '''function to window time series using a tukey window with window parameter alpha 
    then to perform a time laged cross correlation with normalistation norm1'''
    

    y = signal.tukey(len(y),alpha=wparr)*y
    x = signal.tukey(len(x),alpha=wparr)*x
    y = np.array(y)
    x = np.array(x)

    corr = signal.correlate(x, y, mode=mode, method='fft')

    lnorm1 = len(y)-abs(lags)

    norm1 = ( x.std() * y.std() * lnorm1)

    norm2 = (x.std() * y.std() * x.size)

    if norm2 == 0:
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

def maxlagdata( s1name, s2name, magdata, comp):
    '''function to produce rolling window time lagged cross correleation of two singals with peak finding routine for signals s1 and s2 
    returning array (of array) of xcor vals, lag at peak closest to zero and to append values to a network object, values obtained for each window
    plots a graph with time lags on the yaxis (2*Pc_wave_period) and windowed epochs on x axis (xmax=total_time/window size)
    sample rate fs, order of butterworth filter order '''

    '''fuction should reach into above function to grab names and time series in a combinatorial loop'''
    md = pd.read_csv(magdata)

    ts1 = md[md['IAGA'] == s1name][f'db{comp}_nez']
    ts2 = md[md['IAGA'] == s2name][f'db{comp}_nez']
    
    s1 = np.array(ts1)
    s2 = np.array(ts2)

    fs = 1

    order = 3

    #Pc wave period ranges and labels
    label = ['Pc2','Pc3','Pc4','Pc5']

    tf = [5,10,45,150,600]

    nm = 4

    #----setting up tlxc

    # initialising multi-dim arrays to store values
    # each returned element from function will be one of these arrays to later stack

    # mdxc = [[],[],[],[]]

    # mdl = [[],[],[],[]]

    # flags = [[],[],[],[]]

    num_windows = []

    window_size_a = []

    for i, num in enumerate(tf):
        if i<len(tf)-1:
        # if i==1:
        # could speed this loop ? 

            # minimum window size 2xperiod for xcor

            ws = (tf[i+1] - tf[i])
            print(ws)
            window_size = nm*ws

            t_start = 0
            t_end = t_start + window_size
            step_size = (window_size * 3)//4 # amount of window overlap

            # print('windsize',window_size,'stepsize', step_size)

            # keep empty arrays otside of loop(s) if values needed to be returned outside of function

            min0lags = []

            min0vals = []

            flags = []

            while t_end <= len(s1):

                counter = t_end//window_size

                # print(counter)

                y11 = butter_bandpass_filter(s1[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)
                y21 = butter_bandpass_filter(s2[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)

                lags0 = np.arange(-(len(y21) - 1), len(y21))

                rs = xcorr(y11,y21,lags0,wparr=1)

                # gives values of lags at peaks and troughs satifsfying condition in C0peaks, and peak closest to 0

                pp = c0peak(rs,lags0)

                pn = c0peak(-rs,lags0)

                # if statment bellow to catogrise noise into flag array for all cases and append nan into other arrays to keep time dimensionality
                # continue statment restarts the while loop rather than exiting out of it with break
                # first if statment, after or, to prevent xcor with a flat line at zero, otherwise period finding function needs modifcation
                # preventing xcor with noise

                if (pn == 'A' and pp =='A') or (np.max(abs(y11))<0.1 or np.max(abs(y21))<0.1) :
                    
                    min0lags.append(np.nan)

                    min0vals.append(np.nan)

                    flags.append('A')

                    t_start = t_start + step_size
                    
                    t_end = t_end + step_size

                    continue
                
                elif pn == 'A' and pp !='A':

                    min0lags.append(pp[0])

                    min0vals.append(rs[np.where(lags0==pp[0])])

                    flags.append('B')

                    plt.plot

                    t_start = t_start + step_size
                    
                    t_end = t_end + step_size

                    continue

                elif pn != 'A' and pp =='A':

                    min0lags.append(pn[0])

                    min0vals.append(rs[np.where(lags0==pn[0])])

                    flags.append('B')

                    t_start = t_start + step_size
                    
                    t_end = t_end + step_size

                    continue


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

                # statments to check if data is strictly wave-like so ps approx. pxc (cannot be smaller, nature of xcor)

                if pxc > ps*1.56:

                    min0lags.append(extr_0l)

                    min0vals.append(xc0l)

                    flags.append('B')

                    t_start = t_start + step_size
                    
                    t_end = t_end + step_size

                    continue

                else:

                    min0lags.append(extr_0l)

                    min0vals.append(xc0l)

                    flags.append('C')

                    t_start = t_start + step_size
                    
                    t_end = t_end + step_size

            # will loop over final arrays to creat the network
            # also need to stack tlxc, tau and flag arrays? could always filter the tlxc and tau arrays before hand
            # add arrays to check network i guess unless it's too much work
            # could return indcices for where a certain condtion is met i.e flag[i]='B' which gives time stamps

            # range of min0lags, min0vals and flags the same for each band
            # index i will be the window number

            # for loop for parallel interation to create networks

            for j, (xcval, lag) in enumerate(zip(min0vals,min0lags)):

                # print(i, xcval, lag)
                
                if xcval>0 and lag >0:
                    dna[i].add_interaction(s1name, s2name, t=j)

                if xcval>0 and lag <0:
                    dna[i].add_interaction(s2name,s1name, t=j)

                if xcval<0 and lag <0:
                    # print(i,'i','cond3','j',j)
                    dna[i].add_interaction(s1name, s2name, t=j)

                if xcval>0 and lag >0:
                    dna[i].add_interaction(s2name,s1name, t=j)

                if lag==0:
                    na[i].add_interaction(s1name, s2name, t=j)

        # gives num of windows in network, last value gets updates over

        # in order to find Pc ranges require only diference values between indices hence one redundant value below arrays
        print(i)

        num_windows.append(len(min0lags))

        window_size_a.append(window_size)

    return num_windows , window_size_a

    

def network(data, component):
    '''function to create both a directed and undirected netowrk from data for component n,e or z, using the xcor-peak-finding function 
    for two stations'''

    # will loop over final arrays to creat the network
    # also need to stack tlxc, tau and flag arrays? could always filter the tlxc and tau arrays before hand
    # add arrays to check network i guess unless it's too much work
    # could return indcices for where a certain condtion is met i.e flag[i]='B' which gives time stamps


    md = pd.read_csv(data)

    print(md.head())

    # data set with only one set of station names to be used as a label
    s_labels = md.drop_duplicates(subset=['IAGA'])['IAGA']

    print(s_labels[0],s_labels[1]) 

    # loop to create plot with all the time series if needed

    # fig, ax = plt.subplots(len(labels), figsize=(11, 8), facecolor='w', edgecolor='k')
    # fig.subplots_adjust(hspace = .5, wspace=.7)
    # ax = ax.ravel()

    # for i, txt in enumerate(labels):

    #     print(i,txt)
    #     ts1 = md[md['IAGA'] == txt]['dbn_nez']

    #     ts2 = md[md['IAGA'] == txt]['dbe_nez']

    #     ts3 = md[md['IAGA'] == txt]['dbz_nez']

    #     ax[i].plot(np.arange(len(ts1)),ts1, label='Bn')

    #     ax[i].plot(np.arange(len(ts2)),ts2, label='Be')

    #     ax[i].plot(np.arange(len(ts3)),ts3,label='Bz')

    #     ax[i].legend(loc=2)
    # plt.show()



    # y1 = md['N']
    # y2 = md['E']
    # y3 = md['Z']

    # y11 = mdxc['dbn_nez']
    # y22 = mdxc['dbe_nez']
    # y33 = mdxc['dbz_nez']


    # print(list(itertools.combinations(range(1,7),2)))
    scb = list(itertools.combinations(s_labels,2))
    for i in range(len(scb)):

        # print(scb[i][0],scb[i][1])
        k = maxlagdata(scb[i][0],scb[i][1],data, component)

    print(list(itertools.chain(scb)))


    # two station test for networks
    # k = maxlagdata(s_labels[0],s_labels[1],data, component)



    
    for i in [1,2,3]:

        for j in range(k[0][i]):

            fig, axs = plt.subplots(figsize=(8, 8), nrows=2)

            axs = axs.flatten()

            # print(j)

            dns = dna[i].time_slice(j,j+1)

            ns = na[i].time_slice(j,j+1)

            axs[0].set_axis_off()

            axs[1].set_axis_off()
           
            nx.draw(dns, ax=axs[0],with_labels=True, font_weight='bold',arrowsize=20, edgecolor='red',width=1.2)

            nx.draw(ns, ax=axs[1],with_labels=True , font_weight='bold',edgecolor='orange',width=1.2)

            axs[0].set_title(f'digraph with window epoch {j} out of {k[0][i]}, with window size {k[1][i]}, Pc{i+2}')

            axs[1].set_title(f'graph with window epoch {j} out of {k[0][i]}, with window size {k[1][i]}, Pc{i+2}')

            plt.show()

    # at this point network should have some interaction values, so print so check if everything is working
    # need to print interaction list of both directed and undirected graphs and compare with lags plots before moving forward
    # strange that plot for each time slice, easier to check for i>0.


# time series data from 24/02/2012 starting at 3am and lasting for 4 hours with second resolution containing 7 stations
network('20201020-13-39-supermag.csv','n')







# longg = float(ld['RAN'].iloc[0])
# latt = float(ld['RAN'].iloc[1])

# longg2 = float(ld['GIM'].iloc[0])
# latt2 = float(ld['GIM'].iloc[1])

# x = md.index

# y1 = md['N']
# y2 = md['E']
# y3 = md['Z']
# dateutc= md['Date_UTC']

# y11 = mdxc['dbn_nez']
# y22 = mdxc['dbe_nez']
# y33 = mdxc['dbz_nez']

# plt.plot(signal.tukey(len(y1),alpha=1.2))
# plt.plot(signal.slepian(len(y1),width=0.0001))
# plt.show()

# global factor to be multiplied by the period range
# must be aleast T_Pc for osclating signals


# maxlagdata(y2,y22)
