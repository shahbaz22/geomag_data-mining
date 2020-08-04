import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz
from crosscor_bp2 import dateftutc, dateflocal, apexd

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


def xcorr(x, y, lags1):
    # Pad shorter array if signals are different lengths
    if x.size > y.size:
        pad_amount = x.size - y.size
        y = np.append(y, np.repeat(0, pad_amount))
    elif y.size > x.size:
        pad_amount = y.size - x.size
        x = np.append(x, np.repeat(0, pad_amount))

    corr = np.correlate(x, y, mode='full')  # scale = 'none'


    corr = corr/(x.std() * y.std() * (x.size - abs(lags1)))

    return corr

#function (tranforms) for T,f duel axis

def forward(x):
    return 1 /x


def inverse(x):
    return 1 /x

md = pd.read_csv('20200121-17-21-supermag.csv')
print(md.head())

# time series to be x correlated
mdxc = pd.read_csv('20200722-17-06-supermag.csv', header = 0)

print(md.dtypes)
# index can be used as seconds count

#needs to be done
#obtain filtered data for both components
# do x cor on both data

x = md.index

y1 = md['N']
y2 = md['E']
y3 = md['Z']

y11 = mdxc['dbn_nez']
y22 = mdxc['dbe_nez']
y33 = mdxc['dbz_nez']

# label data
ld = pd.read_csv('supermag-stations.csv', delimiter=',')

longg = float(ld['RAN'].iloc[0])
latt = float(ld['RAN'].iloc[1])

dateutc= md['Date_UTC']

#Pc wave period ranges and labels
label = ['Pc2','Pc3','Pc4','Pc5']

tf = [5,10,45,150,600]


def tlxcdata( s1, s2 ):
    #function to plot frequency responce of butterworth filter and plot bandpass spectrograms
    # sample rate fs, order of butterworth filter order 
    fs = 1
    order=3

    y11 = np.zeros((len(tf),len(s1)))

    y22 = np.zeros((len(tf),len(s2)))
    fig, axs = plt.subplots(4,1, figsize=(12, 7), facecolor='w', edgecolor='k')

    fig.tight_layout(pad=3.0)

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            ts=((tf[i+1]) -num)

            #windowed time series data from butter_bandpass_filter function
            y11[i] = butter_bandpass_filter(s1, 1/(tf[i+1]),1/num, fs, order=order)

            y22[i] = butter_bandpass_filter(s2, 1/(tf[i+1]),1/num, fs, order=order)

    #----setting up tlxc

    lags = np.arange(-(s1.size - 1), s1.size)

    tlxc = np.zeros((len(tf),len(lags)))


    for i, num in enumerate(tf):
        if i<len(tf)-1:
            ts = ( (tf[i+1]) - num )

            c = ['b', 'g', 'r', 'm']

            #windowed time series data from butter_bandpass_filter function
            tlxc[i] = xcorr(y11[i],y22[i],lags)
            axs[i].plot(lags, tlxc[i], color=c[i])

            axs[i].set_ylabel('corr. coef.')
            axs[i].set_xlabel('npts')
            # axs[i].set_title(label[i], color=c[i])
            axs[i].set_ylim(-1,1)
            # axs[i].set_xlim(len(lags)/6, len(lags)/6 + ts*20)



    plt.show()



    fig, axs = plt.subplots(4,1, figsize=(12, 7), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.5)

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            #spectrogram plots
            # window nperseg must be at least 2*period due to xcor
            wsize= round(2*(tf[i+1]-tf[i]))

            f, t, sxx = signal.spectrogram(tlxc[i], 1, nperseg=wsize, window='hamming')

            print('x',round(len(x),-3))

            axs[i].pcolormesh(t, f, sxx, alpha=0.9)

            left, right = axs[i].get_xlim()

            xlim = round(len(x), 1)

            xt = np.linspace(round(left,-2),round(right,-2), 17,endpoint=True)
            # xt = np.linspace(0,2*xlim, 17,endpoint=True)

            # xl=np.linspace(-xlim, xlim, 17, endpoint=True)
            xl = np.linspace(-round(right,-2)/2,round(right,-2)/2,17,endpoint=True)


            # exeption to better tune y limits for second subplot
            if i ==1:

                axs[i].set_ylim(1/tf[i+1],(1-i*0.55)/num)
            else:

                axs[i].set_ylim(1/tf[i+1],(1-i*0.1)/num)

            secax = axs[i].secondary_yaxis('right', functions=(forward, inverse))

            axs[i].set_xticks(xt)

            axs[i].set_xticklabels(xl.astype(np.int))
         
            axs[i].set_xlabel('time lag [s]')

            secax.set_ylabel('period [s]')

            axs[i].set_ylabel('freq [Hz]')

            # plot vertical lines for apex distance (use values in UTC at relative local time)
        # for i in [1,2,3]:
        #     ax2.axvline(x=indx[i], color='red')
        #     vvalue = axs[i].get_ylim()
        #     # print(vvalue)
        #     ax2.annotate(s=f'{apexd(latt,longg, dateu.iloc[i])}', xy =((1/4 * (i+1)) -0.23,0.55), xycoords='axes fraction', verticalalignment='center', horizontalalignment='center' , rotation = 270, color='red', weight='bold')


        for i, label in enumerate(('(a)', '(b)', '(c)', '(d)')):
            axs[i].text(0.05, 0.95, label, transform=axs[i].transAxes,
            fontsize=16, va='top', color='w')

    # plt.savefig(f'tlxcorspec_{s1.name}{s2.name}_comp_Pc')

    plt.show()

tlxcdata(y1,y1)

    