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

    norm1 = 1/( x.std() * y.std() * lnorm1)

    norm2 = 1/(x.std() * y.std() * x.size) 

    return corr*norm1

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

def fperiod(y, cutoffh=0.2, viapeaks=False):
    '''finds period of discrete point signal y, using intrapolation or to find 0 crossings
     or peak finding, then finds the average crossing distance and mult. by 2 to find period
     set to True to find period using peaks,however must be two or more values satifying find_peaks
     pyaC is the intrapol. function, first and last zero bypassed'''


    if viapeaks:
        ppts = find_peaks(y, height=cutoffh, prominence = 0.1 )[0]

    else:
        ppts = pyaC.zerocross1d(np.arange(len(y)), y, getIndices=False)


    s = []

    # loop to calculate serperation between zero crossings or peaks
    for i in range(len(ppts)-1):
        d = ppts[i+1] - ppts[i]
        s.append(d)

    # average of zero crossing seperation taken and doubled to obtain period
    return 2*np.mean(s)

def maxlagdata( s1, s2 ):
    # function to produce rolling window time lagged cross correleation of two singals, s1 and s2
    # plots a graph with time lags on the yaxis (2*Pc_wave_period) and windowed epochs on x axis (xmax=total_time/window size)
    # sample rate fs, order of butterworth filter order 
    
    fs = 1

    order = 3

    #Pc wave period ranges and labels
    label = ['Pc2','Pc3','Pc4','Pc5']

    tf = [5,10,45,150,600]
    
    # index values to select timestamps and filter time values
    # n = 5

    # indx = np.linspace(0,len(x)-1,n)
    # # local time filtered date
    # datel = dateflocal(longg, dateutc, indx)

    # dateu = dateftutc(indx, dateutc)

    # # loop to calculate apexd values sepatarley to speed up code
    
    # apexda = np.zeros(4)
    
    # for i in [1,2,3]:
    #     apexda[i] = apexd(latt,longg, dateu.iloc[i])
    # print(apexda)

    nm = 8

    #----setting up tlxc
    # 
    fig, ax = plt.subplots(2,2, figsize=(11, 8), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.7)
    ax = ax.ravel()

    for i, num in enumerate(tf):
        if i<len(tf)-1:
        # if i==1:        

            # minimum window size 2xperiod for xcor

            ws = (tf[i+1] - tf[i])
            window_size = nm*ws


            t_start = 0
            t_end = t_start + window_size
            step_size = (window_size * 2)//6 # amount of window overlap

            print('windsize',window_size,'stepsize', step_size)

            min0lags = []

            min0vals = []

            flags = []

            while t_end <= len(x):

                y11 = butter_bandpass_filter(s1[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)
                y21 = butter_bandpass_filter(s2[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)

                lags0 = np.arange(-(len(y21) - 1), len(y21))

                rs = xcorr(y11,y21,lags0,wparr=1)

                # gives values of lags at peaks and troughs satifsfying condition in C0peaks, and peak closest to 0

                pp = c0peak(rs,lags0)

                pn = c0peak(-rs,lags0)

                # if statment bellow to catogrise noise into flag array for all cases and append nan into other arrays to keep time dimensionality
                # continue statment restarts the while loop rather than exiting out of it with break
                # 

                if pn == 'A' and pp =='A':
                    
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

                print(pp,'pp')

                mlagsp = pp[1]

                print(mlagsp)

                # gives lag value at peak closest to 0

                mlagp = pp[0]

                # print(mlagp)

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

                # statments to check if data is strictly wave-like so ps approx. pxc

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

                # period of pxc cannot be less than that of signal ps

                # --- plotting to examine waveforms

                # if pxc > ps*1.56:
                # # if (xc0l <0 and len(mlagsn) <3):

                #     fig, ax = plt.subplots(2, figsize=(11, 8), facecolor='w', edgecolor='k')
                #     fig.subplots_adjust(hspace = .5, wspace=.7)
                #     ax = ax.ravel()
                    
                #     # ax[0].set_title(f'{label[i]} wave signals, T range {tf[i]} - {tf[i+1]} s, epoch # {t_end//window_size}')
                #     ax[0].set_title(f'period {ps}')

                #     ax[0].plot(np.arange(len(y11)),y11, linestyle='--', color='r', label='s1')
                #     ax[0].plot(np.arange(len(y11)),y21, label='s2')
                #     ax[0].grid()

                #     ax[1].plot(lags0,rs)

                #     # ax[1].set_title(f'Tlxcor plot with peak sep. {abs(nl)+pl}, xcor_T range {2*tf[i]} - {2*tf[i+1]} s ')
                #     ax[1].set_title(f'period {pxc}')

                #     ax[1].scatter(extr_0l,xc0l, color='red')

                #     ax[1].vlines(0,-1,1, linestyle='--',color='r')

                #     ax[1].grid()

                #     ax[0].legend()

                #     # plt.savefig(f'{label[i]}waveform_epoch{t_end//window_size}')

                # plt.show()

                # # ---------------


       


# time series data from 24/02/2012 starting at 3am and lasting for 4 hours with second resolution containing 7 stations
# header = 0 to remove title from dataset
# md = pd.read_csv('20201020-13-39-supermag.csv', header=0, delimiter=',')

md = pd.read_csv('20201020-13-39-supermag.csv')

# gmd = np.genfromtxt('20201020-16-03-supermag2.txt', dtype=None)

print(md.head())

# data set with only one set of station names to be used as a label
labels = md.drop_duplicates(subset=['IAGA'])['IAGA'] 

# loc selects rows and columns based on their name iloc used index position
# print(labels)
# print(labels.iloc[0], type(labels.iloc[0]))

# print(md.columns)

k = labels[0]

ts1 = md[md['IAGA'] == k][['IAGA','dbn_nez']]

# ts1 = md.loc[['RAN','dbn_nez']]
# ts1 = md.filter(like=k, axis=0)

print(ts1)

# loop to create plot with all the time series

print(len(labels))

fig, ax = plt.subplots(len(labels), figsize=(11, 8), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.7)
ax = ax.ravel()
for i, txt in enumerate(labels):
    print(i,txt)
    ts1 = md[md['IAGA'] == txt]['dbn_nez']
    ts2 = md[md['IAGA'] == txt]['dbe_nez']
    ts3 = md[md['IAGA'] == txt]['dbz_nez']
    ax[i].plot(np.arange(len(ts1)),ts1, label='Bn')
    ax[i].plot(np.arange(len(ts2)),ts2, label='Be')
    ax[i].plot(np.arange(len(ts3)),ts3,label='Bz')
    ax[i].legend(loc=2)
plt.show()












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
