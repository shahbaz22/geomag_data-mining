import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz
from geopacktest import lvalue
from datetime import timedelta



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
    # print('x.size', x.size)
    # print('corr.size', corr.size)
    # print(corr)
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
date= md['Date_UTC']

colname = y1.name
print(date)


ml = md.iloc[0]['MLAT']
print(ml)

#start and end times for time series used later in plotting

# startt =  np.datetime64(date.iloc[0])

# endt = np.datetime64(date.iloc[-1])

# reducedt = np.arange(startt,endt,dtype='datetime64[h]')


# date = date.str.extract(f'(\d\d:\d\d:\d\d)', expand=True)

# # print(date.astype('str'))

# date = date.iloc[::3000]

# date = np.squeeze(date.values)

#code snippet for cleaning datetime to time and cleaning/sampling timestamps
# with time delta shift from UTC to local time

# creating equally spaced indices for timestamps in x

ind = np.linspace(0,len(x)-1,5)

mask = date.index.isin(ind)

date = date[mask]

# removing date from time value and storing as numpy arr

date = pd.to_datetime(date, format='%Y/%m/%d %H:%M:%S') - timedelta(hours=5)

date = date.astype(str)

date = date.str.extract(f'(\d\d:\d\d:\d\d)', expand=True)

date = np.squeeze(date.values)

# 

sarr = np.zeros(len(date))

for i, num in enumerate(date):
    print(num)
    d0 = time.strptime(str(num), "%H:%M:%S")
    sarr[i] = d0[5] + d0[4]*60 + d0[3]*3600


# plt.plot(x,y1,label='N')
# plt.plot(x,y2,label='E')
# plt.plot(x,y3,label='Z')

# plt.legend()
# plt.ylabel('nT')
# plt.xlabel('seconds')
# plt.show()
# print(x)


#function (tranforms) for Windowed epochs to time duel axis



def tlxcdata( s1, s2 ):
    # function to produce rolling window time lagged cross correleation of two singals, s1 and s2
    # plots a graph with time lags on the yaxis (2*Pc_wave_period) and windowed epochs on x axis (xmax=total_time/window size)
    # sample rate fs, order of butterworth filter order 
    fs = 1

    order = 3

    #Pc wave period ranges and labels
    label = ['Pc2','Pc3','Pc4','Pc5']

    tf = [5,10,45,150,600]

    #----setting up tlxc

    lags = np.arange(-(s1.size - 1), s1.size)

    tlxc = np.zeros((len(tf),len(lags)))

    # minimum window size 2xperiod for xcor

    fig, ax = plt.subplots(2,2, figsize=(11, 8), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.5)
    ax = ax.ravel()

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            
            window_size = int(np.round(2.3*(tf[i+1] - tf[i])))

            print('window_size',window_size)

            t_start = 0
            t_end = t_start + window_size

            rss = []

            while t_end <= len(x):

                y11 = butter_bandpass_filter(s1[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)
                y21 = butter_bandpass_filter(s2[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)


                lags0 = np.arange(-(y11.size - 1), y11.size)


                rs = xcorr(y11,y21,lags0)
                # not taking whole period so xcor coef is off
                # print(rs, 'rs')
                # plt.plot(lags0, rs)
                # plt.show()
                rss.append(rs)
                t_start = t_start + window_size
                # print(t_start,'tstart')
                t_end = t_end + window_size
                # print(t_end, 'tend')

            # print('lags0', lags0)

            # fig = plt.figure()

            # fig.tight_layout()

            # ax[i] = fig.add_subplot()

            rsst = np.transpose(rss)

            print(rsst.shape)

            print(rsst[:,0])

            ax[i].set_yticks(np.arange(0,len(lags0), tf[i]))

            ax[i].set_yticklabels(np.arange(-(y11.size - 1), y11.size, tf[i]))

            # ax[i].set_title(f'Rolling window (={window_size} pts) time lagged xcor for N,E comp. Pc{i+2} ')
            # ax[i].set_title(f'{label[i]}')
            #to get time (windows) on the x axis
            im=ax[i].imshow(rsst, aspect='auto', alpha=0.9)

            # ax.set_xticks(np.arange((npts/window_size)-1))
            ax[i].set_xlabel('windowed epochs')
            ax[i].set_ylabel('time lag')

            caxp = fig.add_axes([ax[i].get_position().x1+0.01,ax[i].get_position().y0,0.02,ax[i].get_position().height])
            # get a copy of the axes rectangle 2x2 numpy array of the form [[x0, y0], [x1, y1]].
            cbar = plt.colorbar(im, cax=caxp)
            cbar.set_label('Corr. coeff.')

            ax2 = ax[i].twiny() # Remark: twiny will create a new axes 
            # where the y-axis is shared with ax1, 
            # but the x-axis is independant - important!
            # ax2.set_xlim(ax[i].get_xlim()) # ensure the independant x-axes now span the same range
            ax2.set_xticks(ind) # copy over the locations of the x-ticks from the first axes
            ax2.set_xticklabels(date,rotation=70) # But give them a different meaning
            lvals for 

            # plot vertical lines
            for i in [1,2,3]:
                ax2.axvline(x=ind[i], color='red')
                vvalue = ax[i].get_ylim()
                print(vvalue)
                ax2.annotate(s='l=22', xy =((1/4 * (i+1)) -0.23,0.55), xycoords='axes fraction', verticalalignment='center', horizontalalignment='center' , rotation = 270, color='black', weight='bold')

            
                
    plt.savefig(f'tlxcor_{s1.name}{s2.name}_comp_Pc')

    plt.show()


tlxcdata(y2,y3)

    