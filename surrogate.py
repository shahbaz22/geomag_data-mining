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
import datetimegit
import time
import random
from random import shuffle
def _compute_AAFT_surrogates(a):
    i, data, angle = a
    # create Gaussian data
    gaussian = np.random.randn(len(data))
    # sorts the guassian data in size order
    gaussian.sort()
    # rescale data
    # data.argsort() indicies of sorted array
    ranks = data.argsort().argsort()
    rescaled_data = gaussian[ranks]
    # transform the time series to Fourier domain
    xf = np.fft.rfft(rescaled_data)
    # randomise the time series with random phases     
    cxf = xf * np.exp(1j * angle)
    # return randomised time series in time domain
    ft_surr_ts = np.fft.irfft(cxf, n = len(data) )
    # rescale back to amplitude distribution of original data
    sorted_original = data.copy()
    sorted_original.sort()
    ranks = ft_surr_ts.argsort().argsort()
    rescaled_data = sorted_original[ranks]
    
    return (i, rescaled_data)

def get_single_AAFT_surrogate(ts, seed = None):
    """
    Returns single amplitude adjusted FT surrogate.
    Seed / integer : when None, random seed, else fixed seed (e.g. for multivariate AAFT surrogates).
    """

    if seed is None:
        np.random.seed()
    else:
        np.random.seed(seed)

    xf = np.fft.rfft(ts, axis = 0)
    angle = np.random.uniform(0, np.pi, len(ts)//2+1)
    del xf

    return _compute_AAFT_surrogates([None, ts, angle])[-1]

x = np.linspace(1,20,20000)
# ts = [1 if i<5 else 0 for i in x]
fig, ax =  plt.subplots(2,2)
ts= np.sin(x)
# ts = signal.tukey(len(ts))*ts

# ts_fourier  = np.fft.rfft(ts)
# print(np.random.uniform(0,np.pi,len(ts)//2+1))
# random_phases = np.exp(np.random.uniform(0,np.pi,len(ts)//2+1)*1.0j)
# print(len(random_phases),random_phases)
# ts_fourier_new = ts_fourier*random_phases
# ts_new  = np.fft.irfft(ts_fourier_new)
# ts_new = get_single_AAFT_surrogate(ts)

# ax[0,0].plot(x,ts)
# ax[0,0].set_title('ts')
# ax[0,1].plot(x,ts_new)
# ax[0,1].set_title('surrogate ts random phase')
# print(len(ts),len(ts)//2+1)
# tsc = np.correlate(ts, ts, mode='full')
# new_tsc = np.correlate(ts_new, ts_new, mode='full')
# print(type(ts))
# lags = np.arange(-(len(ts) - 1), len(ts))
# ax[1,0].plot(lags,tsc)
# ax[1,1].plot(lags,new_tsc)
# ax[1,0].set_title('ts xcor')
# ax[1,1].set_title('surrogate xcor')
# plt.show()


# print(np.random.uniform(np.pi,len(ts)//2+1))

# new_ts = get_single_AAFT_surrogate(ts)


# ts_fourier  = np.fft.rfft(ts)
# print(np.average(np.random.uniform(0,np.pi,len(ts)//2+1)),np.average(ts))
# random_phases = np.exp(np.random.uniform(0,np.pi,len(ts)//2+1)*1.0j)
# ts_fourier_new = ts_fourier*random_phases
# new_ts = np.fft.irfft(ts_fourier_new)



def surrogates(x, ns, tol_pc=5., verbose=True, maxiter=10, sorttype="quicksort"):
    """
    Returns iAAFT surrogates of given time series.
    Parameter
    ---------
    x : numpy.ndarray, with shape (N,)
        Input time series for which IAAFT surrogates are to be estimated.
    ns : int
        Number of surrogates to be generated.
    tol_pc : float
        Tolerance (in percent) level which decides the extent to which the
        difference in the power spectrum of the surrogates to the original
        power spectrum is allowed (default = 5).
    verbose : bool
        Show progress bar (default = `True`).
    maxiter : int
        Maximum number of iterations before which the algorithm should
        converge. If the algorithm does not converge until this iteration
        number is reached, the while loop breaks.
    sorttype : string
        Type of sorting algorithm to be used when the amplitudes of the newly
        generated surrogate are to be adjusted to the original data. This
        argument is passed on to `numpy.argsort`. Options include: 'quicksort',
        'mergesort', 'heapsort', 'stable'. See `numpy.argsort` for further
        information. Note that although quick sort can be a bit faster than 
        merge sort or heap sort, it can, depending on the data, have worse case
        spends that are much slower.
    Returns
    -------
    xs : numpy.ndarray, with shape (ns, N)
        Array containing the IAAFT surrogates of `x` such that each row of `xs`
        is an individual surrogate time series.
    See Also
    --------
    numpy.argsort
    """
    # as per the steps given in Lancaster et al., Phys. Rep (2018)
    nx = x.shape[0]
    xs = np.zeros((ns, nx))
    ii = np.arange(nx)

    # get the fft of the original array
    x_amp = np.abs(np.fft.fft(x))
    x_srt = np.sort(x)
    r_orig = np.argsort(x)

    # loop over surrogate number
    for k in range(ns):

        # 1) Generate random shuffle of the data
        count = 0
        r_prev = np.random.permutation(ii)
        r_curr = r_orig
        z_n = x[r_prev]
        percent_unequal = 100.

        # core iterative loop
        while (percent_unequal > tol_pc) and (count < maxiter):
            r_prev = r_curr

            # 2) FFT current iteration yk, and then invert it but while
            # replacing the amplitudes with the original amplitudes but
            # keeping the angles from the FFT-ed version of the random
            y_prev = z_n
            fft_prev = np.fft.fft(y_prev)
            phi_prev = np.angle(fft_prev)
            e_i_phi = np.exp(phi_prev * 1j)
            z_n = np.fft.ifft(x_amp * e_i_phi)

            # 3) rescale zk to the original distribution of x
            r_curr = np.argsort(z_n, kind=sorttype)
            z_n[r_curr] = x_srt.copy()
            percent_unequal = ((r_curr != r_prev).sum() * 100.) / nx

            # 4) repeat until number of unequal entries between r_curr and 
            # r_prev is less than tol_pc percent
            count += 1

        if count >= (maxiter - 1):
            print("maximum number of iterations reached!")

        xs[k] = np.real(z_n)

    return xs
# x = np.linspace(1,20,2000)
# # ts = [1 if i<5 else 0 for i in x]
# fig, ax =  plt.subplots(2,2)
# ts= np.sin(x)
# ts = signal.tukey(len(ts))*ts

# ts_new = surrogates(ts,1)[0]

# ax[0,0].plot(x,ts)
# ax[0,0].set_title('ts')
# ax[0,1].plot(x,ts_new)
# ax[0,1].set_title('surrogate ts random phase')
# print(len(ts),len(ts)//2+1)
# tsc = np.correlate(ts, ts, mode='full')
# new_tsc = np.correlate(ts_new, ts_new, mode='full')
# print(type(ts))
# lags = np.arange(-(len(ts) - 1), len(ts))
# ax[1,0].plot(lags,tsc)
# ax[1,1].plot(lags,new_tsc)
# ax[1,0].set_title('ts xcor')
# ax[1,1].set_title('surrogate xcor')
# plt.show()


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
 


def two_station_surrogate_xcor(s1name:str, s2name:str, md:pd.DataFrame, pc_type:str, comp:str) ->None:
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
    fs = 1
    order = 3
    
    ts1 = md[md['IAGA'] == s1name][f'db{comp}_nez']
    ts2 = md[md['IAGA'] == s2name][f'db{comp}_nez']
    utc1 =  md[md['IAGA'] == s1name]['Date_UTC']
    utc2 =  md[md['IAGA'] == s2name]['Date_UTC']
    window_size = 5*20
    t_start = 0
    t_end = t_start + window_size
    count = 0

    pc_periods = [5,10,45,150,600]

    if pc_type == 'pc2':
        t1 = pc_periods[0]
        t2 = pc_periods[1]
    elif pc_type == 'pc3':
        t1 = pc_periods[1]
        t2 = pc_periods[2]

    while t_end <= len(ts1):
        
        y1 = ts1[int(t_start):int(t_end)]
        y2 = ts2[int(t_start):int(t_end)]
        s1 = np.array(y1)
        s2 = np.array(y2)
        num_nans_s1 = np.count_nonzero(np.isnan(s1))
        num_nans_s2 = np.count_nonzero(np.isnan(s2))

        if num_nans_s1>5 or num_nans_s2>5:
            t_start = t_start + window_size                
            t_end = t_end + window_size
            count += 1
            continue

        y11 = butter_bandpass_filter(s1, 1/t2, 1/t1, fs, order=order)
        y21 = butter_bandpass_filter(s2, 1/t2, 1/t1, fs, order=order)
        # print('y11',len(y11))
        # print('y21',len(y21))     
        if (np.max(abs(y11))<0.25 and np.max(abs(y21))<0.25):
            t_start = t_start + window_size                
            t_end = t_end + window_size
            count += 1
            continue

        lags = np.arange(-(len(y21) - 1), len(y21))
        xc = xcorr(y11,y21,lags,wparr=0.6)
        num_nans_xc = np.count_nonzero(np.isnan(xc))
        
        if num_nans_xc>10:
            t_start = t_start + window_size                
            t_end = t_end + window_size
            count += 1
            continue        

        # ss1= get_single_AAFT_surrogate(s1)
        # ss2 = get_single_AAFT_surrogate(s2)
        ss1 = s1
        ss2 = s2

        fig, ax = plt.subplots(nrows=3,ncols=3, facecolor='w', edgecolor='k',figsize=(13,8))
        fig.subplots_adjust(wspace=0.5, hspace=None)


        labs = [f'{s1name} Pc2 wf autocorr',f'{s2name} Pc2 wf auto corr', ]
        
        for ind, arr in enumerate([y11, y21]):
            # plotting waveform autocorrs
            auto_xc = xcorr(arr,arr,lags,wparr=0.6)
            ax[ind,2].plot(lags, auto_xc, label=f'{labs[ind]}')
            ax[ind,2].axvline(x=0,c='red',ls='--')
            ax[ind,2].grid()
            ax[ind,2].legend()

        ax[0,0].plot(np.arange(len(y1)),y1, color='red', linestyle='--',label =f'{s1name} time series 1')
        ax1 = ax[0,0].twinx()
        ax1.grid()
        ax1.set_ylabel(r'$B_{2} (nT)$')
        ax1.plot(np.arange(len(y2)),y2, label = f'{s2name}  time series 2')
        ax1.legend()
        
        ax[0,0].set_ylabel(r'$B_{1} (nT)$')
        ax[0,0].grid()
        ax[0,0].legend()
        ax[1,0].plot(np.arange(len(y11)),y11, color='red', linestyle='--',label =f'{s1name} Pc2 waveform')
        ax[1,0].plot(np.arange(len(y21)),y21, label = f'{s2name} Pc2 waveform')

        ax[1,0].set_ylabel('B (nT)')
        ax[1,0].grid()
        ax[1,0].legend()
        ax[2,0].plot(lags,xc)
        ax[2,0].set_xlabel('Time (s)')
        ax[2,0].set_ylabel('Normalised xcor.')
        ax[2,0].axvline(x=0,c='red',ls='--')
        ax[2,0].grid()
        ax[2,0].set_yticks(np.arange(-1,1,0.2))

        ax[0,1].plot(np.arange(len(ss1)),ss1, color='red', linestyle='--',label =f'{s1name} sur. time series 1')
        ax0 = ax[0,1].twinx()
        ax0.grid()
        ax0.set_ylabel(r'$B_{2} (nT)$')
        ax0.plot(np.arange(len(ss2)),ss2, label = f'{s2name} sur. time series 2')
        ax0.legend()

        y21_sur = np.random.normal(np.mean(y21), np.std(y21), size=len(y21))
        y11_sur = y11

        auto_xc = xcorr(y11_sur, y11_sur, lags,wparr=0.6)
        ax[2,2].plot(lags, auto_xc, label=f'{s2name} Pc2 wf shuffle autocorr')
        ax[2,2].axvline(x=0,c='red',ls='--')
        ax[2,2].grid()
        ax[2,2].legend()

        xc_sur = xcorr(y11_sur, y21_sur,lags,wparr=0.6)
        
        ax[0,1].set_ylabel(r'$B_{1} (nT)$')
        ax[0,1].grid()
        ax[0,1].legend()
        ax[1,1].plot(np.arange(len(y11_sur)),y11_sur, color='red', linestyle='--',label =f'{s1name} sur. Pc2 waveform')
        ax[1,1].plot(np.arange(len(y21_sur)),y21_sur, label = f'{s2name} sur. Pc2 waveform')
        
        ax[1,1].set_ylabel('B (nT)')
        ax[1,1].grid()
        ax[1,1].legend()
        ax[2,1].plot(lags,xc_sur, label='sur. xcor')
        ax[2,1].set_xlabel('Time (s)')
        # ax[2,1].set_ylabel('Normalised xcor.')
        ax[2,1].axvline(x=0,c='red',ls='--')
        ax[2,1].grid()
        ax[2,1].set_yticks(np.arange(-1,1,0.2))
        ax[2,1].legend()

        # plt.savefig(f'plots/waveform_{label}.png')
        plt.show()

        t_start = t_start + window_size                
        t_end = t_end + window_size
        count += 1
        print(count)

        # break
    


md = pd.read_csv('networks_data/17march2013_4am.csv')
print(md.head())
print(md.columns.tolist())
print(md.tail())
# data set with only station names to be used as a label
s_labels = md.drop_duplicates(subset=['IAGA'])['IAGA']
two_station_surrogate_xcor(s_labels[2], s_labels[10], md, pc_type = 'pc2', comp='e')

