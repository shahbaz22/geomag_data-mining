import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal


md=pd.read_csv('20200121-17-21-supermag.csv')
md.head()

md.dtypes
# index can be used as seconds count


x=md.index
y1=md['N']
y2=md['E']
y3=md['Z']
plt.plot(x,y1,label='N')
plt.plot(x,y2,label='E')
plt.plot(x,y3,label='Z')
plt.legend()
plt.ylabel('nT')
plt.xlabel('seconds')
plt.show()
# print(x)


#Filter implemntation 
from scipy.signal import butter, lfilter, freqz

#Pc wave period ranges and labels
label=['Pc2','Pc3','Pc4','Pc5']
tf=[5,10,45,150,600]

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

#function (tranforms) for T,f duel axis

def forward(x):
    return 1 /x


def inverse(x):
    return 1 /x


def filtdata(yi):
    #function to plot frequency responce of butterworth filter and plot bandpass spectrograms
    # sample rate fs, order of butterworth filter order 
    fs = 1
    order=3
    fig, ax = plt.subplots(constrained_layout=True)
    
    #for loop to plot frequency responce for all bandpass freqs
    for i, num in enumerate(tf):
        if i<len(tf)-1:
            b, a = butter_bandpass(1/(tf[i+1]),1/num, fs, order)
            #worN is the num of point in the bandpass visualisation
            #2d array output (x,y) from butter coefs
            w, h = freqz(b, a, worN=len(x))
            ax.plot((fs * 0.5 / np.pi) * w, abs(h), label=f'{label[i]}')
    
    ax.set_xscale('log')
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel('Gain')
    #ax.set_title(f' Bandpass filters for {label} waves')
    # error raised as 1/x uses 0 from bandpass function

    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel('period [s]')
    ax.legend(loc='best')
    plt.show()

    #creates multi-dim y array to store values for later use
    #y will have 0 to tf values down and y[i] values accross for each band of time series data 
    y=np.zeros((len(tf),len(yi)))
    fig, ax = plt.subplots(nrows=1, figsize=(7,7))
    #[x,y ,width, hight]
    #positioning of subplots and axis for subplots
    # ax1 = fig.add_axes([0.21,0.65,0.25,0.2])
    # ax2 = fig.add_axes([0.22,0.17,0.25,0.2])

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            #windowed time series data from butter_bandpass_filter function
            y[i] = butter_bandpass_filter(yi, 1/(tf[i+1]),1/num, fs, order=order)
            ax.plot(x, y[i], label=label[i])
            # ax1.plot(x, y[i], label=label[i])
            # ax2.plot(x, y[i], label=label[i])
    #plt.title(f'Band-pass Filtered signal for {label} waves')
    #ax1.set_xlim(3,tl[-1])
    ax.set_xlabel('time (seconds)')
    ax.set_ylabel('nT')
    # ax1.set_xlim(375,575)
    # ax1.set_ylim(-0.04,0.04)
    # ax2.set_xlim(475,575)
    # ax2.set_ylim(-0.007,0.005)
    ax.legend(loc='best')
    plt.show()

    # setting up empty subplots to add in below loop
    # note that the the number of subplots has to be accounted for in the loop size
    fig, axs = plt.subplots(4,1, figsize=(12, 7), facecolor='w', edgecolor='k')
    # fig.subplots_adjust(hspace = .5, wspace=.001)
    # ax.ravel similar to ax.flatten and combines arrays
    axs = axs.ravel()


    for i, num in enumerate(tf):
        if i<len(tf)-1:
            #spectrogram plots
            ts=((tf[i+1]) -num)
            print(ts)
            print(tf[i]+1)

            f, t, sxx = signal.spectrogram(y[i], 1, nperseg=round(tf[i+1]), window='hamming')

            # fa= np.append(fa, f)
            # ta= np.append(ta, t)
            # sxxa= np.append(sxxa, sxx)


            # f=f[1:-1]
            # t=t[1:-1]
            # sxx=sxx[1:,1:]
            # mask = np.logical_and(f >=float( 1/tf[i+1]) , f <= float(1/num))
            # print('mask',mask)
            # f=f[mask]
            # print(f)
            # t=t[mask]
            # sxx=sxx[[mask],[mask]]


            axs[i].pcolormesh(t, f, sxx)
            #axs[i].set_title(label[i])
            # print('ax ranges',1/tf[i+1],1/num)
            axs[i].set_ylim(1/tf[i+1],1/num)
            secax = axs[i].secondary_yaxis('right', functions=(forward, inverse))
            # secax.set_yticks([np.linspace(num,tf[i+1],2)])
            # secax.set_yticklabels([f"{num}",f"{tf[i+1]}"])
            secax.set_ylabel('period [s]')

            #axs[i].set_yscale('log')

            axs[i].set_ylabel('freq [Hz]')




    axs[3].set_xlabel('time [s]')
    plt.show()

filtdata(y1)

    