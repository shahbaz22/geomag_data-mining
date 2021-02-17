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
import dynetx as dn
import itertools
# from dynetx.readwrite import json_graph
# import json

# initialising directed and undirected network arrays to store values for all Pc bands

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

def dateflocal(date,ind):

    ''' fuction to convert datetime utc formal to numpy array to be used for plotting.
    where ind is the total number of time points needed from the date-time series'''

    mask = date.index.isin(ind)

    date = date[mask]

    date = pd.to_datetime(date, format='%Y/%m/%d %H:%M:%S')

    # date = pd.to_datetime(date, format='%Y/%m/%d %H:%M:%S') - timedelta(hours=6.1)

    date = date.astype(str)

    date = date.str.extract(f'(\d\d:\d\d:\d\d)', expand=True)

    date = np.squeeze(date.values)

    return date

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

    # md is type pandas dataframe


    ts1 = md[md['IAGA'] == s1name][f'db{comp}_nez']
    ts2 = md[md['IAGA'] == s2name][f'db{comp}_nez']

    s1 = np.array(ts1)
    s2 = np.array(ts2)

    fs = 1

    order = 3

    #Pc wave period ranges and labels
    label = ['Pc2','Pc3','Pc4','Pc5']

    tf = [5,10,45,150,600]

    nm = 3

    #----setting up tlxc

    # initialising multi-dim arrays to store values
    # each returned element from function will be one of these arrays to later stack

    # mdxc = [[],[],[],[]]

    # mdl = [[],[],[],[]]

    # flags = [[],[],[],[]]

    num_windows = []

    window_size_a = []

    step_size_array = []


    for i, num in enumerate(tf):
        if i<len(tf)-1:

            # minimum window size 2xperiod for xcor

            ws = (tf[i+1] - tf[i])

            window_size = nm*ws

            t_start = 0
            t_end = t_start + window_size
            step_size = (window_size * 3)//4 # amount of window overlap

            # print('windsize',window_size,'stepsize', step_size)
            print(step_size, 'step_size')
            # keep empty arrays otside of loop(s) if values needed to be returned outside of function

            min0lags = []

            min0vals = []

            flags = []

            step_size_array.append(step_size)

            while t_end <= len(s1):

                counter = t_end//window_size

                # print(counter)

                # y11 , y21 filtering time series to use with TLXC

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

                    fig, ax = plt.subplots(nrows=2, facecolor='w', edgecolor='k')
                    fig.subplots_adjust(hspace=0.8)

                    ax[0].plot(np.arange(len(y21)),y21, color='red', linestyle='--',label ='time series 1')

                    ax[0].plot(np.arange(len(y11)),y11, label ='time series 2')

                    ax[0].set_xlabel('Time (s)')

                    ax[0].set_ylabel('B (nT)')

                    ax[1].plot(lags0,rs)

                    ax[1].set_xlabel('Time lags (s)')

                    ax[1].set_ylabel('Normalised cross-correlation.')

                    ax[1].axvline(x=0,c='red',ls='--')

                    ax[1].grid()

                    ax[0].grid()

                    ax[0].legend()

                    plt.show()

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

                # statments to check if data is strictly wave-like 
                # so ps approx. pxc cannot be smaller, nature of xcor)
                # xcor cannot have small period then average period of both signals
                # so only > bound used for range 

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

                print(i, j, xcval, lag, 'i, j, xcval, lag')
                
                if xcval>0 and lag >0:
                    dna[i].add_interaction(s1name, s2name, t=j)

                elif xcval>0 and lag <0:
                    dna[i].add_interaction(s2name,s1name, t=j)

                elif xcval<0 and lag <0:
                    # print(i,'i','cond3','j',j)
                    dna[i].add_interaction(s1name, s2name, t=j)

                elif xcval<0 and lag>0:
                    dna[i].add_interaction(s2name,s1name, t=j)

                # maybe add else: instead of elif lag==0 to speed up code
                elif lag==0:
                    na[i].add_interaction(s1name, s2name, t=j)


        # gives num of windows in network, last value gets updates over

        # in order to find Pc ranges require only diference values between indices hence one redundant value below arrays

        t = len(min0lags)*step_size

        # print(window_size,'window_size')

        print(t,'time in seconds')

        num_windows.append(len(min0lags))

        window_size_a.append(window_size)

        # can comment out when running large code, always be the same

        np.save('step_size_arr.npy',step_size_array)


    return num_windows , window_size_a

    

def network_global(data, component):
    '''function to create both a directed and undirected global netowrk from data for component n,e or z, appling the xcor-peak-finding
    and network appending function, network_append different pairs of stations combinatorially
    for all pairs and returning text files for all directed and undirected netowrks'''

    md = pd.read_csv(data)

    print(md.head())

    # data set with only station names to be used as a label
    s_labels = md.drop_duplicates(subset=['IAGA'])['IAGA']

    num_stations = len(s_labels)

    s_labels = s_labels[0:3]

    # scb lists station pair permutations as Nc2 

    scb = list(itertools.combinations(s_labels,2))

    print(scb)

    # for loop over permutation pairs to calculate network connections between station pairs and append to network containers

    for i in range(len(scb)):

        # k is the output from the 2 station peak-finding and netowrk appending algo

        print(scb[i][0],scb[i][1],'station pair out of', num_stations,'stations with', i,'out of',len(scb),'operations')

        k = network_append(scb[i][0],scb[i][1],md, component)

    # print list for displaying perm pairs from scb 
    # print(list(itertools.chain(scb)))

    # print(list(dna[3].stream_interactions()))



    # saving network below, plotting displays networks, i is the Pc index of interest 0,1,2,3
    # k[0] is the num of windows (which can overlap) (timestamps) for each network k[1] is the window size

    for i in [0,1,2,3]:

        dn.readwrite.edgelist.write_interactions(dna[i], f'dna{i}.txt')

        dn.readwrite.edgelist.write_interactions(na[i], f'na{i}.txt')



# time series data from 24/03/2012 starting at 3am and lasting for 4 hours with second resolution containing 7 stations
# network_global('20201111-18-30-supermag.csv','e')

def network_cluster(data, component, magp):
    '''function to create both a directed and undirected (clust fixed in MLT) netowrk from data for component n,e or z, appling the xcor-peak-finding
    and network appending function, network_append different pairs of stations combinatorially
    for all pairs and returning text files for all directed and undirected netowrks'''

    # magnetic local time 0 is midnight 6 dawn, 12 noon, 18 dusk 

    # add parameters for dawn, dusk, midnight and noon

    # dawn: 03:00-09:00 (ie 06:00 +/-3hrs)

    # noon: 09:00:15:00 (12:00 +/- 3hrs)

    # dusk:15:00-21:00

    # midnight (21:00-03:00)

    md = pd.read_csv(data)

    # print(md.head())

    # data set with only station names to be used as a label

    s_labels = md.drop_duplicates(subset=['IAGA'])['IAGA']

    num_stations = len(s_labels)

    # s_labels = s_labels[0:3]

    # creating data set to use in analysis to filter MLT values

    if magp == 'dawn':
        
        mltpts = [3,9]

    elif magp == 'noon':
        
        mltpts = [9,15]

    elif magp == 'dusk':
       
        mltpts = [15,21]

    elif magp == 'midnight':
        
        mltpts = [21,3]       

    mltf = md[(md['MLT'] >= mltpts[0]) & (md['MLT'] <= mltpts[1])]

    # ts1 = mltf[mltf['IAGA'] == s_labels[0]]
    # ts2 = mltf[mltf['IAGA'] == 'PG2']

    # print(mltf)

    # print(ts2)

    # print(mltfiltered[(mltfiltered['IAGA'].isin(['C04','C12']))])

    # scb lists station pair permutations as Nc2 

    scb = list(itertools.combinations(s_labels,2))

    # print(scb)

    # for loop over permutation pairs to calculate network connections between station pairs and append to network containers

    for i in range(len(scb)):

        # k is the output from the 2 station peak-finding and netowrk appending algo

        print(scb[i][0],scb[i][1],'station pair out of', num_stations,'stations with', i,'out of',len(scb),'operations')

        # the next algortihm will then select entire station times based on the filtered mlt values between date ranges

        k = network_append(scb[i][0],scb[i][1],mltf, component)


    # saving network below, plotting displays networks, i is the Pc index of interest 0,1,2,3
    # k[0] is the num of windows (which can overlap) (timestamps) for each network k[1] is the window size

    # the time compenent for this dynamical network will measure how much the eath has moved 'time' 
    # underneath the floating section

    # for i in [0,1,2,3]:

    #     dn.readwrite.edgelist.write_interactions(dna[i], f'networks_data/dna_spd_{magp}_Pc{i+1}{component}.txt')

    #     dn.readwrite.edgelist.write_interactions(na[i], f'networks_data/na_spd_{magp}_Pc{i+1}{component}.txt')


# code to save all cluster networks for given component

# for i in ['dusk','midnight']:#'dawn','noon',
#     network_cluster('20201111-18-30-supermag.csv','e',magp=i)
network_cluster('20201111-18-30-supermag.csv','e',magp='dusk')


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

# plt.grid()

# plt.show()

