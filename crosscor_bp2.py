import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz
# from geopacktest import apdvalue
from datetime import timedelta
import time



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


def xcorr(x, y, lags1, mode='full'):
    'calculate the time lagged corss covariance for two time series x and y'
    # if x.size > y.size:

    #     pad_amount = x.size - y.size

    #     y = np.append(y, np.repeat(0, pad_amount))
    # elif y.size > x.size:

    #     pad_amount = y.size - x.size

    #     x = np.append(x, np.repeat(0, pad_amount))

    corr = signal.correlate(x, y, mode=mode, method='fft')  # scale = 'none'

    print('windowedx',x)
    print('windowedy',y)
    
    y = signal.tukey(len(y),alpha=0.7)*y
    x = signal.tukey(len(x),alpha=0.7)*x

    print(len(x))

    # plt.plot(np.arange(len(x)),x,color='g')
    # plt.plot(np.arange(len(y)),y, color='b')
    # plt.show()

    lnorm1 = len(y)-abs(lags1)

    print('lagsnorm', lnorm1)

    norm1 = 1/( x.std() * y.std() * lnorm1)

    print('norm1', norm1)

    norm2 =1/(x.std() * y.std() * x.size) 

    corr = corr*norm1

    print('corr', corr)

    # corr= signal.tukey(len(corr),alpha=0.8)*corr 
    # print('x.size', x.size)
    # print('corr.size', corr.size)
    # print(corr)
    # divides smaller aplitdes by smaller norm factor
    return corr

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

# time series data
md = pd.read_csv('20200121-17-21-supermag.csv')
# header = 0 to remove title from dataset
mdxc = pd.read_csv('20200722-17-06-supermag.csv', header = 0)

# print(mdxc.head)
# label data
ld = pd.read_csv('supermag-stations.csv', delimiter=',')

# print(md.head())

longg = float(ld['RAN'].iloc[0])
latt = float(ld['RAN'].iloc[1])

longg2 = float(ld['GIM'].iloc[0])
latt2 = float(ld['GIM'].iloc[1])

# print(longg,latt,longg2,latt2)


x = md.index

y1 = md['N']
y2 = md['E']
y3 = md['Z']
dateutc= md['Date_UTC']

y11 = mdxc['dbn_nez']
y22 = mdxc['dbe_nez']
y33 = mdxc['dbz_nez']

# fig, ax = plt.subplots(constrained_layout=True)

# ax.plot(x,y11,label='N')
# ax.plot(x,y22,label='E')
# ax.plot(x,y33,label='Z')

# ax.legend()
# ax.set_ylabel('[nT]')
# ax.set_xlabel('time [s]')
# ax.set_xlim(0,len(x))

# ax2 = ax.twiny() # Remark: twiny will create a new axes 

# n = 5

# indx = np.linspace(0,len(x)-1,n)
# # local time filtered date
# datel = dateflocal(longg, dateutc, indx)

# ax2.set_xticks(indx) # copy over the locations of the x-ticks from the first axes
# ax2.set_xticklabels(datel, rotation=0)
# ax2.set_xlabel('timestamp [24 hour]')

# plt.show()


ml = md.iloc[0]['MLAT']

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


    # fig, ax = plt.subplots(4,1, figsize=(15, 8), facecolor='w', edgecolor='k')
    # fig.subplots_adjust(hspace = .5, wspace=.7)
    # ax = ax.ravel()

    # for i, num in enumerate(tf):
    #     if i<len(tf)-1:
            
    #         window_size = 8*(tf[i+1] - tf[i])

    #         print('window_size',window_size)

    #         t_start = 0
    #         t_end = t_start + window_size

    #         rss0 = []

    #         while t_end <= len(x):

    #             y11 = butter_bandpass_filter(s1[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)
    #             y21 = butter_bandpass_filter(s2[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)

    #             # mode='same' for autocor
    #             rs = xcorr(y11,y21,0, mode='valid')


    #             # where staments used to remove edge effects from data with N - lags norm as seen in test code
    #             # N - lags norm allows for better resoltuion of edges rather than only central peak
    #             # not needed if using standard pearson norm for tlcxc 

    #             # rs = np.where(rs<1, rs, 1)

    #             # rs = np.where(rs>-1, rs, -1)

    #             rss0 = np.append(rss0, rs)


    #             t_start = t_start + window_size
    #             # print(t_start,'tstart')
    #             t_end = t_end + window_size
    #             # print(t_end, 'tend')
            
    #         c=['b', 'g', 'r', 'm']
            
    #         ax[i].plot(np.arange(len(rss0)),rss0, color=c[i])
    #         ax[i].scatter(np.arange(len(rss0)),rss0,s=5, color='black')
    #         ax[i].axhline(0,color='black')

    #         # ax[i].set_ylim(-1,1)

    #         xt = np.arange(0,len(rss0),window_size)

    #         # ax[i].set_xticks(xt)

    #         # ax[i].set_xticklabels(np.arange(len(xt)))

    #         ax[i].set_xlabel('windowed epochs')

    #         ax[i].set_ylabel('Corr. coeff.')


    #         # if i <2:
    #         #     ax[i].set_xlim(0, window_size*30)
    #         # else:
    #         #     ax[i].set_xlim(0, len(rss0))
            
    #         for i, label in enumerate(('(a)', '(b)', '(c)', '(d)')):
    #             ax[i].text(0, 0.95, label, transform=ax[i].transAxes,
    #             fontsize=16, va='top')


    # plt.show()

            # print(len(rss0), np.shape(rss0))
            # print(rss0)


    
    fig, ax = plt.subplots(2,2, figsize=(11, 8), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.7)
    ax = ax.ravel()

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            # minimum window size 2xperiod for xcor

            window_size = 8*(tf[i+1] - tf[i])

            print('window_size',window_size)

            t_start = 0
            t_end = t_start + window_size

            rss = []

            while t_end <= len(x):

                y11 = butter_bandpass_filter(s1[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)
                y21 = butter_bandpass_filter(s2[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)

                print('y11',y11)

                lags0 = np.arange(-(len(y21) - 1), len(y21))

                rs = xcorr(y11,y21,lags0)


                # where staments used to remove edge effects from data with N - lags norm as seen in test code
                # N - lags norm allows for better resoltuion of edges rather than only central peak
                # not needed if using standard pearson norm for tlcxc 

                # rs = np.where(rs<1, rs, 1)

                # rs = np.where(rs>-1, rs, -1)
                # print('xcor_result',rs)


                # plt.plot(lags0,rs)
                # plt.show()

                rss.append(rs)

                t_start = t_start + window_size
                # print(t_start,'tstart')
                t_end = t_end + window_size
                # print(t_end, 'tend')

            rsst = np.transpose(rss)

            print(lags0)

            im=ax[i].imshow(rsst, aspect='auto', alpha=0.9, extent=[0,len(x)//window_size, min(lags0), max(lags0)])

            ax[i].set_xlabel('windowed epochs')

            ax[i].set_ylabel('time lag [s]')

            caxp = fig.add_axes([ax[i].get_position().x1+0.01,ax[i].get_position().y0,0.02,ax[i].get_position().height])
            # get a copy of the axes rectangle 2x2 numpy array of the form [[x0, y0], [x1, y1]].
            cbar = plt.colorbar(im, cax=caxp)
            cbar.set_label('Corr. coeff.')

            ax2 = ax[i].twiny() # Remark: twiny will create a new axes 
            # where the y-axis is shared with ax1, 
            # but the x-axis is independant - important!
            # ax2.set_xlim(ax[i].get_xlim()) # ensure the independant x-axes now span the same range
            
            # index values to select timestamps and filter time values

            n = 5

            indx = np.linspace(0,len(x)-1,n)
            # local time filtered date
            datel = dateflocal(longg, dateutc, indx)

            dateu = dateftutc(indx, dateutc)


            ax2.set_xticks(indx) # copy over the locations of the x-ticks from the first axes
            ax2.set_xticklabels(datel, rotation=70) 

            # plot vertical lines for apex distance (use values in UTC at relative local time)
    # for i in [1,2,3]:
    #     ax2.axvline(x=indx[i], color='red')
    #     vvalue = ax[i].get_ylim()
    #     # print(vvalue)
    #     # ax2.annotate(s=f'{apexd(latt,longg, dateu.iloc[i])}, {apexd(latt2,longg2, dateu.iloc[i])}', xy =((1/4 * (i+1)) -0.23,0.55), xycoords='axes fraction', verticalalignment='center', horizontalalignment='center' , rotation = 270, color='black', weight='bold')
    #     ax2.annotate(s=f'{apexd(latt,longg, dateu.iloc[i])}', xy =((1/4 * (i+1)) -0.23,0.55), xycoords='axes fraction', verticalalignment='center', horizontalalignment='center' , rotation = 270, color='black', weight='bold')

    for i, label in enumerate(('(a)', '(b)', '(c)', '(d)')):
        ax[i].text(0.05, 0.95, label, transform=ax[i].transAxes,
        fontsize=16, va='top')
            
                
    plt.savefig(f'tlxcor_{s1.name}{s2.name}_comp_Pc')

    plt.show()

# print(y1.name,'sep' ,y11.name)
# print(y2.name,'sep' ,y22.name)
tlxcdata(y1,y2)

    