import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz

# from geopacktest import apdvalue

npts = 180

x = np.linspace(0, 7, npts)

wf = signal.tukey(npts,alpha=0.7)

randarr = np.random.rand(npts)


y = 1 #signal.gaussian(npts, std=50)

y1 =  np.sin(3*x)

y2 = np.cos(3*x +0.5)

lags = np.arange(-(len(y1) - 1), len(y1))

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

# print('rlags',rlags)
def xcorr(x, y, lags, normlags = True, return_norm = False):
    'y shifted on-top of x'

    y = signal.tukey(len(y),alpha=1)*y
    x = signal.tukey(len(x),alpha=1)*x
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

xcorn1 =xcorr(y1,y2, lags)

# mask to select troublesome values, | is or for numpy array
# corr_mask = (xcorn1 > 1) | (xcorn1 < -1)

fig, axs = plt.subplots(nrows=2,figsize=(6,6))
fig.subplots_adjust(hspace=.3)

x1 = np.arange(0, npts)

ax = axs[0]
ax.plot(x1, y1, 'b', label='y1')
ax.plot(x1, y2, 'r', label='y2')

# ax.set_ylim(-10, 10)
ax.legend(loc='upper right', fontsize='small', ncol=2)
# ax.set_title('Signal time series')
ax.set_xlabel('npts')
ax = axs[1]
# lags=lags[:-2]
# ccor=ccor[:-2]
print(xcorn1)
ax.plot(lags, xcorn1)
# ax.plot(lags0, xcorr(y11,y22, normlags=False))

# ax.set_ylim(-1, 1)
ax.set_ylabel('Corr. coeff.')
ax.set_xlabel('point lag')

ax.grid()

# ax.set_title('time lagged xcorrelation (1/std(y1)*std(y2)*(N-tau) norm)')
for i, label in enumerate(('(a)', '(b)')):
    axs[i].text(0.05, 0.95, label, transform=axs[i].transAxes,
    fontsize=16, va='top', color='black')
plt.show()

# print(y1.name,'sep' ,y11.name)
# print(y2.name,'sep' ,y22.name)


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

