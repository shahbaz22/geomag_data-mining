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

md = pd.read_csv('20200121-17-21-supermag.csv')

x = md.index

y1 = md['N']
y2 = md['E']
y3 = md['Z']
dateutc= md['Date_UTC']

# npts = 200

# x = np.linspace(0, 10, npts)

# y = 1 #signal.gaussian(npts, std=50)

# y1 =  5*np.sin(3*x)*y

# y2 = 5*np.cos(3*x)*y

y11 = [ 0.          ,0.01210866,  0.08244743,  0.24974041,  0.42920564,  0.40718479,
  0.04802476, -0.51087371, -0.93236929, -0.95871827, -0.60322135, -0.076753,
  0.42775375,  0.81436987,  0.9822061,   0.80329072,  0.27768937, -0.36396404,
 -0.78582963, -0.78577094, -0.42392209,  0.04909129,  0.38389821,  0.47558552,
  0.36132142,  0.1409583,  -0.07918484, -0.21642398, -0.24056226, -0.17437194,
 -0.06618001,  0.00763182, -0.04296227, -0.20830181, -0.30867873, -0.14233138,
  0.27786262,  0.68988918,  0.81129174,  0.58091782,]

print('len(y11)',len(y11))

y22 = [ 0.,         -0.00055886, -0.00345903, -0.00899851, -0.01127104, -0.00121474,
  0.02153523,  0.04284484,  0.04600149,  0.02659502, -0.0099757,  -0.05259758,
 -0.07753518, -0.06018273, -0.00833468,  0.03771936,  0.04888929,  0.03463103,
  0.02077647,  0.01147895, -0.00817879, -0.03606534, -0.04518662, -0.0184889,
  0.02104195,  0.03342533,  0.01189114, -0.01227732, -0.01237704,  0.00312277,
  0.00840174, -0.00753681, -0.04322475, -0.08691029, -0.10281339, -0.05395201,
  0.05224448,  0.1577689,   0.20157502,  0.17056901,]

print('len(y22)',len(y22))

lags0 = np.arange(-(len(y22) - 1), len(y22))


lnorm1 = len(y11)-abs(lags0)

print(lnorm1)
def manulcor(x1,x2):
    
    if len(x1)!=len(x2):
        raise ValueError('Arrays must be the same size')

    corm=np.zeros(len(lags))

    for i in range(lags.min(),lags.max()):
        # need not mess with range as 1 added in definition of lags
        print(i,'i')
        # print('shifted array',shift5(x1,i))
        # lags1=np.roll(lags, npts)
        corm[i]=np.sum(np.multiply(shift5(x1,i),x2))#*1/(x.size-abs(lags1[i]))
        print(corm[i],'comi')

        print('shiftarr', shift5(x1,i))
    nar = x1.size - abs(rlags)
    snar = 1#np.roll(nar, npts)

    corm = corm /((x1.std() * x2.std()) * nar)

    print('nar',nar)
    print('ccor1',corm)
    print('size ccor1',len(corm))
    
    return corm

print('y11',y11)
# print('rlags',rlags)
def xcorr(x, y, lags, normlags = True, return_norm = False):

    y = signal.tukey(len(y),alpha=0.5)*y
    x = signal.tukey(len(x),alpha=0.5)*x
    y = np.array(y)
    x = np.array(x)
    corr = signal.correlate(x, y, mode='full', method='fft')

    lnorm1 = len(y)-abs(lags)

    norm1 = 1/( x.std() * y.std() * lnorm1)

    # print('norm1', norm1)

    norm2 =1/(x.std() * y.std() * x.size) 

    print('norm2' , norm2)




    if normlags:
        corr = corr*norm1
    else:
        corr = corr*norm2

    if return_norm:
        return norm1
    else:
        return corr

xcorn1 =xcorr(y11,y22, lags0)

# mask to select troublesome values, | is or for numpy array
corr_mask = (xcorn1 > 1) | (xcorn1 < -1)
print('bad cov values',xcorn1[corr_mask])
print('lags',lags0[corr_mask])
print('lnorm',lnorm1[corr_mask])
print('bad norm vals',xcorr(y11,y22,lags0,return_norm=True)[corr_mask])
# reutrns normalisation array

print('cov vals',xcorn1)
print('norm vals 1/',xcorr(y11,y22,lags0,return_norm=True))


fig, axs = plt.subplots(nrows=2,figsize=(6,6))
fig.subplots_adjust(hspace=.3)

x1 = np.arange(0, len(y11))

ax = axs[0]
ax.plot(x1, y11, 'b', label='y1')
ax.plot(x1, y22, 'r', label='y2')

# ax.set_ylim(-10, 10)
ax.legend(loc='upper right', fontsize='small', ncol=2)
# ax.set_title('Signal time series')
ax.set_xlabel('time [s]')
ax = axs[1]
# lags=lags[:-2]
# ccor=ccor[:-2]
print(xcorn1)
ax.plot(lags0, xcorn1)
# ax.plot(lags0, xcorr(y11,y22, normlags=False))

# # ax.set_ylim(-1, 1)
ax.set_ylabel('Corr. coeff.')
ax.set_xlabel('point lag')
ax.grid()
# ax.set_title('time lagged xcorrelation (1/std(y1)*std(y2)*(N-tau) norm)')
for i, label in enumerate(('(a)', '(b)')):
    axs[i].text(0.05, 0.95, label, transform=axs[i].transAxes,
    fontsize=16, va='top', color='black')
plt.show()

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

                rs = np.where((rs<1) | (rs>-1), rs, 5)

                # rs = np.where(, rs, -1)
                # print('xcor_result',rs)


                # plt.plot(lags0,rs)
                # plt.show()

                rss.append(rs)

                t_start = t_start + window_size
                # print(t_start,'tstart')
                t_end = t_end + window_size
                # print(t_end, 'tend')

            rsst = np.transpose(rss)

            im=ax[i].imshow(rsst, aspect='auto', alpha=0.9, extent=[0,len(x)//window_size, min(lags0), max(lags0)])

            ax[i].set_xlabel('windowed epochs')

            ax[i].set_ylabel('time lag [s]')

            caxp = fig.add_axes([ax[i].get_position().x1+0.01,ax[i].get_position().y0,0.02,ax[i].get_position().height])
            # get a copy of the axes rectangle 2x2 numpy array of the form [[x0, y0], [x1, y1]].
            cbar = plt.colorbar(im, cax=caxp)
            cbar.set_label('Corr. coeff.')

            # ax2 = ax[i].twiny() # Remark: twiny will create a new axes 
            # where the y-axis is shared with ax1, 
            # but the x-axis is independant - important!
            # ax2.set_xlim(ax[i].get_xlim()) # ensure the independant x-axes now span the same range
            
            # index values to select timestamps and filter time values

            # n = 5

            # indx = np.linspace(0,len(x)-1,n)
            # # local time filtered date
            # datel = dateflocal(longg, dateutc, indx)

            # dateu = dateftutc(indx, dateutc)


            # ax2.set_xticks(indx) # copy over the locations of the x-ticks from the first axes
            # ax2.set_xticklabels(datel, rotation=70) 

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
# xcorr of small section of data 
# l=60
# s=40
# y11 = y1[s:l]
# y21 = y2[s:l]
# lags0 = np.arange(-(y11.size - 1), y11.size)
# fig, axs = plt.subplots(nrows=2,figsize=(6,6))
# fig.subplots_adjust(hspace=1)
# ax = axs[0]
# ax.plot(np.arange(s,l), y11, 'b', label='y1')
# ax.plot(np.arange(s,l), y21, 'r', label='y2')
# # ax.set_ylim(-10, 10)
# ax.legend(loc='upper right', fontsize='small', ncol=2)
# ax.set_title('Signal time series')
# ax.set_xlabel('time')

# ax = axs[1]
# ax.plot(lags0, xcorr(y21,y11,lags0), 'b', label='y1')
# ax.set_title('Manual xcorrelation, 1/std(y1)*std(y2)*(N-tau) norm')
# ax.axhline(y=1)
# ax.axhline(y=-1)
# plt.show()

# ax.set_xlabel('lag of y1 relative to y2')
# # ax.set_ylim(-1.1, 1.1)

# # Rolling window time lagged cross correlation

# window_size = int(np.round(npts/60)) #samples
# t_start = 0
# t_end = t_start + window_size
# print(window_size)

# lags0 = np.arange(-(y1.size - 1), y1.size)
# print(y1,'y1')
# print(y2,'y2')
# rs = xcorr(y1,y2,lags0)


# rss=[]
# while t_end <= npts:
#     y11 = y1[t_start:t_end]
#     y21 = y2[t_start:t_end]
#     lags0 = np.arange(-(y11.size - 1), y11.size)
#     print(y1,'y1')
#     print(y2,'y2')
#     rs = xcorr(y11,y21,lags0)
#     # not taking whole period so xcor coef is off
#     print(rs, 'rs')
#     # plt.plot(lags0, rs)
#     # plt.show()
#     print(rs.size,'rs_size')
#     rss.append(rs)
#     t_start = t_start + window_size
#     print(t_start,'tstart')
#     t_end = t_end + window_size
#     print(t_end, 'tend')

# print(lags0,'lags')

# fig = plt.figure()

# ax = fig.add_subplot()

# ax.set_title(f'Rolling window (={window_size} pts) time lagged cross correlation of wave packet')
# #to get time (windows) on the x axis
# im=ax.imshow(np.transpose(rss))

# ax.set_aspect('equal')
# # ax.set_xticks(np.arange((npts/window_size)-1))
# ax.set_xlabel('windowed epochs')
# ax.set_ylabel('time lag')

# caxp = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
# # get a copy of the axes rectangle 2x2 numpy array of the form [[x0, y0], [x1, y1]].
# cbar = plt.colorbar(im, cax=caxp)
# cbar.set_label('Corr. coeff.')
# plt.show()

# #cross spectrogram of two combined signals

# windowpts = 210

# fig, axs = plt.subplots(2,1, figsize=(12, 7), facecolor='w', edgecolor='k')


# f, t, sxx = signal.spectrogram(ccor, 1, nperseg=windowpts, window='hamming')

# print(np.amax(f), 'f.max')

# fm = np.amax(f)

# axs[0].set_ylim(0,fm/5)

# print(t,'t')

# print(t.size,'t size')

# newticks = np.linspace(-np.amax(t)/2, np.amax(t)/2, t.size)

# print(newticks, 'nax')

# print(newticks.size, 'nax')

# # axs[0].set_xticks(np.linspace(-t.max/2, t.max/2, t.size))

# cax0 = axs[0].pcolormesh(newticks, f, sxx, alpha=0.7)

# axs[0].set_xlabel('npts')


# axs[0].set_ylabel('cross frequency Hz')

# axs[1].set_ylim(0,fm/5)

# f1, t1, sxx1 = signal.spectrogram(y1+y2, 1, nperseg=windowpts, window='hamming')

# cax1 = axs[1].pcolormesh(t1, f1, sxx1, alpha=0.7)

# axs[1].set_ylabel('frequency Hz')

# axs[1].set_xlabel('npts')

# # fig.colorbar(cax1)

# # fig.colorbar(cax0)

# plt.show()

