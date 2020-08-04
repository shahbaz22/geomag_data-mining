import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz
from geopacktest import apdvalue
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
    
    # if norm_bm:
 
    #     y = signal.slepian(len(y),width=0.00001)*y
    #     x = signal.slepian(len(x),width=0.00001)*x
    # else:
    #     y = signal.tukey(len(y),alpha=wparr)*y
    #     x = signal.tukey(len(x),alpha=wparr)*x

    y = signal.tukey(len(y),alpha=wparr)*y
    x = signal.tukey(len(x),alpha=wparr)*x
    y = np.array(y)
    x = np.array(x)

    corr = signal.correlate(x, y, mode=mode, method='fft')

    lnorm1 = len(y)-abs(lags)

    norm1 = 1/( x.std() * y.std() * lnorm1)

    norm2 = 1/(x.std() * y.std() * x.size) 

    return corr*norm1


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

x = md.index

y1 = md['N']
y2 = md['E']
y3 = md['Z']
dateutc= md['Date_UTC']

y11 = mdxc['dbn_nez']
y22 = mdxc['dbe_nez']
y33 = mdxc['dbz_nez']

# plt.plot(signal.tukey(len(y1),alpha=1.2))
# plt.plot(signal.slepian(len(y1),width=0.0001))
# plt.show()
# y1wind = signal.tukey(len(y1),alpha=0.4)
# plt.plot(y1wind)
# plt.show()
    # fig, ax = plt.subplots(4,1, figsize=(15, 8), facecolor='w', edgecolor='k')
    # fig.subplots_adjust(hspace = .5, wspace=.7)
    # ax = ax.ravel()
def autoxcdata( s1, s2 ):
    # function to produce rolling window time lagged cross correleation of two singals, s1 and s2
    # plots a graph with time lags on the yaxis (2*Pc_wave_period) and windowed epochs on x axis (xmax=total_time/window size)
    # sample rate fs, order of butterworth filter order 
    
    fs = 1

    order = 3

    tf = [5,10,45,150,600]

    fig, ax = plt.subplots(4,1, figsize=(15, 8), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.7)
    
    for i, num in enumerate(tf):
        if i<len(tf)-1:
            
            window_size = 4*(tf[i+1] - tf[i])

            print('window_size',window_size)

            t_start = 0
            t_end = t_start + window_size

            rss0 = []

            while t_end <= len(x):

                y11 = butter_bandpass_filter(s1[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)
                y21 = butter_bandpass_filter(s2[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)

                # mode='same' for autocor
                rs = xcorr(y11,y21,0, mode='valid', wparr=0.1)

                w_start = 0.05
                window = 0.05
                w_end = w_start + window

                rs = xcorr(y11,y21,0, mode='valid',wparr=w_start)

                while np.amax(rs)>1 or np.amin(rs)<-1:

                    rs = xcorr(y11,y21,0, mode='valid',wparr=w_end)
      
                    w_start = w_start + window
                    w_end = w_end + window

                    if w_end >1:
                        break

                # where staments used to remove edge effects from data with N - lags norm as seen in test code
                # N - lags norm allows for better resoltuion of edges rather than only central peak
                # not needed if using standard pearson norm for tlcxc 


                rss0 = np.append(rss0, rs)


                t_start = t_start + window_size
                # print(t_start,'tstart')
                t_end = t_end + window_size
                # print(t_end, 'tend')

            ax = ax.ravel()
            
            c=['b', 'g', 'r', 'm']
            
            ax[i].plot(np.arange(len(rss0)),rss0, color=c[i])
            ax[i].scatter(np.arange(len(rss0)),rss0,s=5, color='black')
            ax[i].axhline(0,color='black')

            # ax[i].set_ylim(-1,1)

            xt = np.arange(0,len(rss0),window_size)

            # ax[i].set_xticks(xt)

            # ax[i].set_xticklabels(np.arange(len(xt)))

            ax[i].set_xlabel('windowed epochs')

            ax[i].set_ylabel('Corr. coeff.')


            # if i <2:
            #     ax[i].set_xlim(0, window_size*30)
            # else:
            #     ax[i].set_xlim(0, len(rss0))
            
            for i, label in enumerate(('(a)', '(b)', '(c)', '(d)')):
                ax[i].text(0, 0.95, label, transform=ax[i].transAxes,
                fontsize=16, va='top')


    plt.show()

autoxcdata(y1,y2)

def tlxcdata( s1, s2 ):
    # function to produce rolling window time lagged cross correleation of two singals, s1 and s2
    # plots a graph with time lags on the yaxis (2*Pc_wave_period) and windowed epochs on x axis (xmax=total_time/window size)
    # sample rate fs, order of butterworth filter order 
    
    fs = 1

    order = 3

    #Pc wave period ranges and labels
    label = ['Pc2','Pc3','Pc4','Pc5']

    tf = [5,10,45,150,600]
    
    # index values to select timestamps and filter time values
    n = 5

    indx = np.linspace(0,len(x)-1,n)
    # local time filtered date
    datel = dateflocal(longg, dateutc, indx)

    dateu = dateftutc(indx, dateutc)

    # loop to calculate apexd values sepatarley to speed up code
    
    apexda = np.zeros(4)
    
    for i in [1,2,3]:
        apexda[i] = apexd(latt,longg, dateu.iloc[i])
    print(apexda)

    #----setting up tlxc
    
    fig, ax = plt.subplots(2,2, figsize=(11, 8), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.7)
    ax = ax.ravel()

    for i, num in enumerate(tf):
        if i<len(tf)-1:
            # minimum window size 2xperiod for xcor

            window_size = 4*(tf[i+1] - tf[i])

            print('window_size',window_size)

            t_start = 0
            t_end = t_start + window_size

            rss = []

            while t_end <= len(x):

                y11 = butter_bandpass_filter(s1[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)
                y21 = butter_bandpass_filter(s2[t_start:t_end], 1/(tf[i+1]), 1/num, fs, order=order)

                # print('y11',y11)

                lags0 = np.arange(-(len(y21) - 1), len(y21))

                # statment to print out anomolous values
                # values above range shown with both norms 
                # hence seem like rounding errors
                # if np.amax(rs)>1 or np.amin(rs)<-1:
                #     plt.plot(rs)
                #     plt.show()

                # --- section of code to customise window size fow each collumn in tlxc plot
                # window becomes as narrow as possible for given conditions until break value
                w_start = 0.05
                window = 0.05
                w_end = w_start + window

                rs = xcorr(y11,y21,lags0,wparr=w_start)

                while np.amax(rs)>1 or np.amin(rs)<-1:

                    rs = xcorr(y11,y21,lags0,wparr=w_end)
                    
                    w_start = w_start + window
                    w_end = w_end + window

                    if w_end >1:
                        break
                # print('minrs',np.amax(rs),'maxrs',np.amin(rs))
                # print(w_end, 'alpha')
                # print(i, 'Pc index')


                # where staments used to remove edge effects from data with N - lags norm as seen in test code
                # N - lags norm allows for better resoltuion of edges rather than only central peak
                # not needed if using standard pearson norm for tlcxc 

                # rs = np.where(rs<1, rs, 50)

                # rs = np.where(rs>-1, rs, -50)


                # plt.plot(lags0,rs)
                # plt.show()

                rss.append(rs)

                t_start = t_start + window_size
                # print(t_start,'tstart')
                t_end = t_end + window_size
                # print(t_end, 'tend')
            print('max',np.amax(rss),'min',np.amin(rss))

            rsst = np.transpose(rss)

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
            ax2.set_xlim(ax[i].get_xlim()) # ensure the independant x-axes now span the same range

            ax2.set_xticks(indx) # copy over the locations of the x-ticks from the first axes
            
            ax2.set_xticklabels(datel, rotation=70) 

            # plot vertical lines for apex distance (use values in UTC at relative local time)
            for i in [1,2,3]:
                ax2.axvline(x=indx[i], color='red')
                vvalue = ax[i].get_ylim()
                # print(vvalue)
                # ax2.annotate(s=f'{apexd(latt,longg, dateu.iloc[i])}, {apexd(latt2,longg2, dateu.iloc[i])}', xy =((1/4 * (i+1)) -0.23,0.55), xycoords='axes fraction', verticalalignment='center', horizontalalignment='center' , rotation = 270, color='black', weight='bold')
                ax2.annotate(s=f'{apexda[i]}', xy =((1/4 * (i+1)) -0.23,0.55), xycoords='axes fraction', verticalalignment='center', horizontalalignment='center' , rotation = 270, color='black', weight='bold')

    
    for i, label in enumerate(('(a)', '(b)', '(c)', '(d)')):
        ax[i].text(0.05, 0.95, label, transform=ax[i].transAxes,
        fontsize=16, va='top')
            
                
    plt.savefig(f'tlxcor_{s1.name}{s2.name}_comp_Pc')

    plt.show()

tlxcdata(y1,y2)
