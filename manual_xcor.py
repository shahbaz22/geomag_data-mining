import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fftpack import fft


npts = 2000

x = np.linspace(0, 200, npts)

y = signal.gaussian(npts, std=50)

y1 =  5*np.sin(2*x)*y

y2 = 5*np.cos(2*x)*y


# y1 =  np.sin(5*x**1.5)

# y2 = np.cos(1.2*x**2)

lags = np.arange(-(x.size - 1), x.size)

rlags=np.roll(lags,npts)
# Maybe add array sizr check into loop
# [:num] takes index values from begining up to num
# [num:] takes index values from num to end
# [:] returns the whole array
# [-num:] takes values starting from -num to end

def shift5(arr, num, fill_value=0):
    result = np.zeros(len(arr))
    if num > 0:
        result[:num] = fill_value
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = fill_value
        result[:num] = arr[-num:]
    else:
        result[:] = arr
    return result


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
    # elif scale == 'coeff':
    # corr /= np.sqrt(np.dot(x, x) * np.dot(y, y))

ccor = xcorr(y1,y2,lags)

fig, axs = plt.subplots(nrows=2,figsize=(6,6))
fig.subplots_adjust(hspace=1)

ax = axs[0]
ax.plot(x, y1, 'b', label='y1')
ax.plot(x, y2, 'r', label='y2')
# ax.set_ylim(-10, 10)
ax.legend(loc='upper right', fontsize='small', ncol=2)
ax.set_title('Signal time series')
ax.set_xlabel('time')

# ax = axs[1]
# ax.plot(rlags,manulcor(y2,y1), 'b', label='y1')
# ax.set_title('Manual xcorrelation, 1/std(y1)*std(y2)*(N-tau) norm')
# # ax.set_xlabel('lag of y1 relative to y2')
# # ax.set_ylim(-1.1, 1.1)


ax = axs[1]
# lags=lags[:-2]
# ccor=ccor[:-2]
ax.plot(lags, ccor)
ax.set_ylim(-1.1, 1.1)
ax.set_ylabel('cross-correlation')
ax.set_xlabel('lag of y1 relative to y2')
ax.set_title('time lagged xcorrelation (1/std(y1)*std(y2)*(N-tau) norm)')
plt.show()


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
# ax.set_ylim(-1.1, 1.1)

# Rolling window time lagged cross correlation

window_size = int(np.round(npts/60)) #samples
t_start = 0
t_end = t_start + window_size
print(window_size)

# lags0 = np.arange(-(y1.size - 1), y1.size)
# print(y1,'y1')
# print(y2,'y2')
# rs = xcorr(y1,y2,lags0)


rss=[]
while t_end <= npts:
    y11 = y1[t_start:t_end]
    y21 = y2[t_start:t_end]
    lags0 = np.arange(-(y11.size - 1), y11.size)
    print(y1,'y1')
    print(y2,'y2')
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

#cross spectrogram of two combined signals

windowpts = 210

fig, axs = plt.subplots(2,1, figsize=(12, 7), facecolor='w', edgecolor='k')


f, t, sxx = signal.spectrogram(ccor, 1, nperseg=windowpts, window='hamming')

print(np.amax(f), 'f.max')

fm = np.amax(f)

axs[0].set_ylim(0,fm/5)

print(t,'t')

print(t.size,'t size')

newticks = np.linspace(-np.amax(t)/2, np.amax(t)/2, t.size)

print(newticks, 'nax')

print(newticks.size, 'nax')

# axs[0].set_xticks(np.linspace(-t.max/2, t.max/2, t.size))

cax0 = axs[0].pcolormesh(newticks, f, sxx, alpha=0.7)

axs[0].set_xlabel('npts')


axs[0].set_ylabel('cross frequency Hz')

axs[1].set_ylim(0,fm/5)

f1, t1, sxx1 = signal.spectrogram(y1+y2, 1, nperseg=windowpts, window='hamming')

cax1 = axs[1].pcolormesh(t1, f1, sxx1, alpha=0.7)

axs[1].set_ylabel('frequency Hz')

axs[1].set_xlabel('npts')

# fig.colorbar(cax1)

# fig.colorbar(cax0)

plt.show()

