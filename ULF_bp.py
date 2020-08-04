import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal
# fuctions for plotting
from crosscor_bp2 import dateftutc, dateflocal, apexd

md=pd.read_csv('20200121-17-21-supermag.csv')
md.head()

md.dtypes

# label data
ld = pd.read_csv('supermag-stations.csv', delimiter=',')

longg = float(ld['RAN'].iloc[0])
latt = float(ld['RAN'].iloc[1])

dateutc= md['Date_UTC']


x=md.index
y1=md['N']
y2=md['E']
y3=md['Z']
fig, ax = plt.subplots(constrained_layout=True)

ax.plot(x,y1,label='N')
ax.plot(x,y2,label='E')
ax.plot(x,y3,label='Z')

ax.legend()
ax.set_ylabel('[nT]')
ax.set_xlabel('time [s]')
ax.set_xlim(0,len(x))

ax2 = ax.twiny() # Remark: twiny will create a new axes 

n = 5

indx = np.linspace(0,len(x)-1,n)
# local time filtered date
datel = dateflocal(longg, dateutc, indx)

ax2.set_xticks(indx) # copy over the locations of the x-ticks from the first axes
ax2.set_xticklabels(datel, rotation=0)
ax2.set_xlabel('timestamp [24 hour]')

plt.show()



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
    # fig, ax = plt.subplots(constrained_layout=True)
    
    #for loop to plot frequency responce for all bandpass freqs
    for i, num in enumerate(tf):
        if i<len(tf)-1:
            b, a = butter_bandpass(1/(tf[i+1]),1/num, fs, order)
            #worN is the num of point in the bandpass visualisation
            #2d array output (x,y) from butter coefs
            w, h = freqz(b, a, worN=len(x))
            # ax.plot((fs * 0.5 / np.pi) * w, abs(h), label=f'{label[i]}')
    
    # ax.set_xscale('log')
    # ax.set_xlabel('f [Hz]')
    # ax.set_ylabel('Gain')
    # #ax.set_title(f' Bandpass filters for {label} waves')
    # # error raised as 1/x uses 0 from bandpass function

    # secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    # secax.set_xlabel('period [s]')
    # ax.legend(loc='best')
    # plt.show()

    #creates multi-dim y array to store values for later use
    #y will have 0 to tf values down and y[i] values accross for each band of time series data 
    y0=np.zeros((len(tf),len(yi)))
    y1=np.zeros((len(tf),len(yi)))

    # fig, axs = plt.subplots(4,1, figsize=(12, 7), facecolor='w', edgecolor='k')

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            ts=((tf[i+1]) -num)

            # c=['b', 'g', 'r', 'm']

            #windowed time series data (yi input) from butter_bandpass_filter function to creat mutiple bands from
            y0[i] = butter_bandpass_filter(yi, 1/(tf[i+1]),1/num, fs, order=order)
            
            # y1[i] = butter_bandpass_filter(y2, 1/(tf[i+1]),1/num, fs, order=order)
            
    #         axs[i].plot(x, y0[i], label=label[i], color=c[i])
    #         axs[i].plot(x, y1[i], label=label[i])

    #         axs[i].set_ylabel('Magnetic field [nT]')
    #         axs[i].set_xlabel('time [s]')
    #         axs[i].set_title(label[i], color=c[i])
    #         # axs[i].set_xlim(len(x)/4, len(x)/4 + ts*20)
    #         ax2 = ax.twiny() # Remark: twiny will create a new axes 

    #         n = 5

    #         indx = np.linspace(0,len(x)-1,n)
    #         # local time filtered date
    #         datel = dateflocal(longg, dateutc, indx)

    #         ax2.set_xticks(indx) # copy over the locations of the x-ticks from the first axes
    #         ax2.set_xticklabels(datel, rotation=20)
    #         ax2.set_xlabel('timestamp [24 hour]')

    # plt.show()






    # axs[3].set_xlabel('time [s]')
    # axs[3].set_title(label[3])

    # plt.show()

    # setting up empty subplots to add in below loop
    # note that the the number of subplots has to be accounted for in the loop size
    fig, axs = plt.subplots(4,1, figsize=(12, 7), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = 1.5, wspace=.5)

    # fig.subplots_adjust(hspace = .5, wspace=.001)
    # ax.ravel similar to ax.flatten and combines arrays
    axs = axs.ravel()


    for i, num in enumerate(tf):
        if i<len(tf)-1:
            #spectrogram plots

            f, t, sxx = signal.spectrogram(y0[i], 1, nperseg=round(tf[i+1]), window='hamming')

            axs[i].pcolormesh(t, f, sxx)
            # print('ax ranges',1/tf[i+1],1/num)
            axs[i].set_ylim(1/tf[i+1],1/num)

            secax = axs[i].secondary_yaxis('right', functions=(forward, inverse))
            
            # secax.set_yticks([np.linspace(num,tf[i+1],2)])
            # secax.set_yticklabels([f"{num}",f"{tf[i+1]}"])
            
            secax.set_ylabel('period [s]')

            #axs[i].set_yscale('log')

            axs[i].set_ylabel('freq [Hz]')
            ax2 = axs[i].twiny() # Remark: twiny will create a new axes 

            n = 5

            indx = np.linspace(0,len(x)-1,n)
            # local time filtered date
            datel = dateflocal(longg, dateutc, indx)

            dateu = dateftutc(indx, dateutc)


            ax2.set_xticks(indx) # copy over the locations of the x-ticks from the first axes
            ax2.set_xticklabels(datel, rotation=20) 

            axs[i].set_xlabel('time [s]')

        # plot vertical lines for apex distance (use values in UTC at relative local time)
        for i in [1,2,3]:
            ax2.axvline(x=indx[i], color='red')
            vvalue = axs[i].get_ylim()
            # print(vvalue)
            ax2.annotate(s=f'{apexd(latt,longg, dateu.iloc[i])}', xy =((1/4 * (i+1)) -0.23,0.55), xycoords='axes fraction', verticalalignment='center', horizontalalignment='center' , rotation = 270, color='red', weight='bold')


        for i, label in enumerate(('(a)', '(b)', '(c)', '(d)')):
            axs[i].text(0.0, 0.95, label, transform=axs[i].transAxes,
            fontsize=16, va='top', color='w')


    plt.show()

filtdata(y1)

    