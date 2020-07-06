import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz



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

    # if scale == 'biased':
    #     corr = corr / x.size
    # elif scale == 'unbiased':
    corr = corr/(x.std() * y.std()*(x.size - abs(lags1)))
    print('x.size', x.size)
    print('corr.size', corr.size)
    print(corr)
    # print(corr[-30:-1])
    # print(x.size - abs(lags))
    # divides smaller aplitdes by smaller norm factor
    return corr

#function (tranforms) for T,f duel axis

def forward(x):
    return 1 /x


def inverse(x):
    return 1 /x

md = pd.read_csv('20200121-17-21-supermag.csv')
print(md.head())

print(md.dtypes)
# index can be used as seconds count

#needs to be done
#obtain filtered data for both components
# do x cor on both data

x = md.index

y1 = md['N']
y2 = md['E']
y3 = md['Z']

plt.plot(x,y1,label='N')
plt.plot(x,y2,label='E')
plt.plot(x,y3,label='Z')

plt.legend()
plt.ylabel('nT')
plt.xlabel('seconds')
plt.show()
# print(x)



#Pc wave period ranges and labels
label = ['Pc2','Pc3','Pc4','Pc5']

tf = [5,10,45,150,600]


def tlxcdata( s1, s2 ):
    #function to plot frequency responce of butterworth filter and plot bandpass spectrograms
    # sample rate fs, order of butterworth filter order 
    fs = 1
    order=3

    y1 = np.zeros((len(tf),len(s1)))

    y2 = np.zeros((len(tf),len(s2)))
    # fig, ax = plt.subplots(nrows=1, figsize=(7,7))
    fig, axs = plt.subplots(4,1, figsize=(12, 7), facecolor='w', edgecolor='k')

    fig.tight_layout(pad=3.0)

    #[x,y ,width, hight]
    #positioning of subplots and axis for subplots
    # ax1 = fig.add_axes([0.21,0.65,0.25,0.2])
    # ax2 = fig.add_axes([0.22,0.17,0.25,0.2])

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            ts=((tf[i+1]) -num)

            #windowed time series data from butter_bandpass_filter function
            y1[i] = butter_bandpass_filter(s1, 1/(tf[i+1]),1/num, fs, order=order)

            y2[i] = butter_bandpass_filter(s2, 1/(tf[i+1]),1/num, fs, order=order)

    #----setting up tlxc

    lags = np.arange(-(s1.size - 1), s1.size)

    tlxc = np.zeros((len(tf),len(lags)))


    for i, num in enumerate(tf):
        if i<len(tf)-1:
            ts = ( (tf[i+1]) - num )

            c = ['b', 'g', 'r', 'm']

            #windowed time series data from butter_bandpass_filter function
            tlxc[i] = xcorr(y1[i],y2[i],lags)
            axs[i].plot(lags, tlxc[i], label=label[i], color=c[i])

            axs[i].set_ylabel('corr. coef.')
            axs[i].set_xlabel('npts')
            axs[i].set_title(label[i], color=c[i])
            # axs[i].set_ylim(-1,1)
            # axs[i].set_xlim(len(lags)/6, len(lags)/6 + ts*20)



    plt.show()

    # Rolling window time lagged cross correlation

    window_size = 150 #samples
    t_start = 4000
    t_end = t_start + window_size
    print(window_size)


    rss=[]
    while t_end <= len(x)-4000:
        y11 = y1[1, t_start:t_end]
        y21 = y2[1, t_start:t_end]
        lags0 = np.arange(-(y11.size - 1), y11.size)

        rs = xcorr(y11,y21,lags0)
        # not taking whole period so xcor coef is off
        print(rs, 'rs')
        # plt.plot(lags0, rs)
        # plt.show()
        print(rs.size,'rs_size')
        rss.append(rs)
        t_start = t_start + window_size
        print(t_start,'tstart')
        t_end = t_end + window_size
        print(t_end, 'tend')

    print(lags0,'lags')

    fig = plt.figure()

    ax = fig.add_subplot()

    ax.set_title(f'Rolling window (={window_size} pts) time lagged cross correlation of wave packet')
    #to get time (windows) on the x axis
    im=ax.imshow(np.transpose(rss))

    ax.set_aspect('equal')
    # ax.set_xticks(np.arange((npts/window_size)-1))
    ax.set_xlabel('windowed epochs')
    ax.set_ylabel('time lag')

    caxp = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    # get a copy of the axes rectangle 2x2 numpy array of the form [[x0, y0], [x1, y1]].
    cbar = plt.colorbar(im, cax=caxp)
    cbar.set_label('Corr. coeff.')

    plt.show()


    fig, axs = plt.subplots(4,1, figsize=(12, 7), facecolor='w', edgecolor='k')
    fig.tight_layout(pad=3.0)

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            #spectrogram plots

            f, t, sxx = signal.spectrogram(tlxc[i], 1, nperseg=round(3*tf[i+1]), window='hamming')

            axs[i].pcolormesh(t, f, sxx, alpha=0.9)
            axs[i].set_title(label[i])
            # print('ax ranges',1/tf[i+1],1/num)
            axs[i].set_ylim(1/tf[i+1],1.2/num)

            secax = axs[i].secondary_yaxis('right', functions=(forward, inverse))
            
            # secax.set_yticks([np.linspace(num,tf[i+1],2)])
            # secax.set_yticklabels([f"{num}",f"{tf[i+1]}"])
            
            secax.set_ylabel('period [s]')

            #axs[i].set_yscale('log')

            axs[i].set_ylabel('freq [Hz]')


    plt.show()

tlxcdata(y1,y2)

    