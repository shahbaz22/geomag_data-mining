import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal
from scipy.signal import butter, lfilter, freqz, find_peaks

# # from geopacktest import apdvalue

npts = 4000

x = np.linspace(0, 15, npts)

y = np.sin(x)

maxl = find_peaks(y, height=2, prominence = 0.1 )[0]

print(maxl)

plt.plot(x,y)
plt.scatter(x[maxl],y[maxl])
plt.show()


# yn1 = np.random.normal(loc=0.0, scale=0.3, size=len(x)) # white noise (additive guassian/normal), scale = std 

# yn2 = np.random.normal(loc=0.0, scale=0.3, size=len(x))

# y1 = np.cos(3*x) + yn1
# y2 =  np.sin(3*x) + yn2
    

# lags = np.arange(-(len(y1) - 1), len(y1))

# def xcorr(x, y, lags, normlags = True, return_norm = False):
#     'y shifted on-top of x'

#     y = signal.tukey(len(y),alpha=1)*y
#     x = signal.tukey(len(x),alpha=1)*x
#     y = np.array(y)
#     x = np.array(x)
#     corr = signal.correlate(x, y, mode='full', method='fft')

#     lnorm1 = len(y)-abs(lags)

#     print(lnorm1)

#     norm1 = 1/( x.std() * y.std() * lnorm1)
#     # corr = np.round(corr,5)
#     # lnorm1 = np.round(lnorm1,5)

#     # print('norm1', norm1)

#     norm2 =1/(x.std() * y.std() * x.size)

#     if normlags:
#         corr = corr*norm1
#     else:
#         corr = corr*norm2

#     if return_norm:
#         return norm1
#     else:
#         return corr

#     c = np.where(corr, corr>1, 1)

#     print(np.max(c),'maxc')

#     return c



# # xcorn1 =xcorr(y1,y2, lags)

# # # mask to select troublesome values, | is or for numpy array
# # # corr_mask = (xcorn1 > 1) | (xcorn1 < -1)

# # fig, axs = plt.subplots(nrows=4,figsize=(6,8))
# # fig.subplots_adjust(hspace=.8)


# # x1 = np.arange(0, npts)

# # ax = axs[0]
# # ax.plot(x1, y1, 'b', label='y1')
# # ax.plot(x1, y2, 'r', label='y2')

# # # ax.set_ylim(-10, 10)
# # ax.legend(loc='upper right', fontsize='small', ncol=2)
# # # ax.set_title('Signal time series')
# # ax.set_xlabel('npts')
# # ax.set_title(f'sin signals with white noise of std = 0.3')
# # ax = axs[1]
# # # lags=lags[:-2]
# # # ccor=ccor[:-2]
# # ax.plot(lags, xcorn1)
# # # ax.plot(lags0, xcorr(y11,y22, normlags=False))

# # # ax.set_ylim(-1, 1)
# # ax.set_ylabel('Corr. coeff.')
# # ax.set_xlabel('point lag')
# # ax.set_title('tlxcor of y1 and y2')

# # ax.grid()
# # # plt.show()

# # def XCORSNR(c):
# #     ''' function to calculate the noise to singal (max(abs(tlxc))) ratio of two sinusoidal signals with white noise
# #     for increasing noise power, scale (where the scale is standard devation, noise power)'''

# #     snr = []
# #     xcormax = []

# #     npts = 2000

# #     x = np.linspace(0, 6, npts)


# #     for i in np.linspace(0.01,3,1000):

# #         print(i, 'iii')

# #         noise_std = i

# #         yn1 = np.random.normal(loc=0.0, scale=noise_std, size=len(x))
# #         yn2 = np.random.normal(loc=0.0, scale=noise_std, size=len(x))

# #         y1n = np.cos(3*x) + yn1
# #         y2n = np.sin(3*x)+ yn2
# #         y1 = np.cos(3*x)
# #         y2 = np.sin(3*x)

# #         lags = np.arange(-(len(y1) - 1), len(y1))

# #         # plt.plot(lags,xcorr(y1n,y2n, lags))
# #         # plt.hlines([1,-1],min(lags),max(lags), 'r', '--')

# #         # plt.show()


# #         xcorn = np.max(abs(xcorr(y1n,y2n, lags, normlags=True)))

# #         xcormax.append(xcorn)


# #         snr1 = i/(np.std(y1)**2)
                       
# #         snr.append(snr1)

# #     axs[2].plot(snr, xcormax, color = c) 

# # # plt.scatter(snr,xcormax)
# # # XCORSNR('r')
# # # XCORSNR('b')
# # # XCORSNR('brown')
# # # XCORSNR('purple')
# # XCORSNR('y')
# # # XCORSNR('c')
# # # XCORSNR('g')
# # axs[2].set_xlabel('NSR')
# # axs[2].set_ylabel('max(|C(tau)|)')
# # axs[2].set_title('max(abs(tlxcor)) vs NSR with white noise')
# # axs[2].grid()
# # # axs[2].set_ylim(0,1.05)
# # # plt.show()


# def XCORSNRN():
#     ''' function to calculate the sample-size to singal (max(abs(tlxc))) ratio of two sinusoidal signals with white noise 
#     of fixed noise power, given as scale (where the scale is standard devation, noise power)'''

#     nnr = []
#     xcormax2 = []

#     t11=[]
#     t22=[]
#     t33=[]
#     t44=[]

#     for i in np.arange(10,400):

#         print(i, 'iii')

#         x = np.linspace(0,i,i*5)

#         noise_std = 0.9

#         yn1 = np.random.normal(loc=0.0, scale=noise_std, size=len(x))
#         yn2 = np.random.normal(loc=0.0, scale=noise_std, size=len(x))
        
#         lags = np.arange(-(len(x) - 1), len(x))

#         y1 = np.cos(3*x)
#         y2 = np.sin(3*x)

#         y1 = y1 +yn1
#         y2 = y2 +yn2


#         # plt.plot(lags,xcorr(y1n,y2n, lags))
#         # plt.hlines([1,-1],min(lags),max(lags), 'r', '--')

#         # plt.show()


#         # xcorn = np.max(abs(xcorr(yn1,yn2, lags, normlags=True)))


#         t1 = np.max(abs(xcorr(y1,y2, lags, normlags=True)))
        
#         t2 = np.max(abs(xcorr(y1,yn2, lags, normlags=True)))
        
#         t3 = np.max(abs(xcorr(yn1,y2, lags, normlags=True)))
        
#         t4 = np.max(abs(xcorr(yn1,yn2, lags, normlags=True)))


#         # xcormax2.append(xcorn)

#         t11.append(t1)

#         t22.append(t2)

#         t33.append(t3)

#         t44.append(t4)

                       
#         nnr.append(i*5)

#     plt.plot(nnr, t11, label ='s1 s1') 
#     plt.plot(nnr, t22, label ='s1 n2') 
#     plt.plot(nnr, t33, label ='n1 s2') 
#     plt.plot(nnr, t44, label ='n1 n2') 

#     plt.xlabel('sample size')
#     plt.ylabel('max(|C(tau)|)')
#     # axs[2].set_title('')
#     plt.grid()
#     # axs[2].set_ylim(0,1.2)

#     #   plt.plot(nnr, xcormax2) 

#     plt.title(f'tlxcor with white noise of {noise_std} std')

#     plt.legend()

#     plt.show()

# # XCORSNR('r')
# # XCORSNR('b')
# # XCORSNR('brown')
# # XCORSNR('purple')
# XCORSNRN()
# # XCORSNR('c')
# # XCORSNR('g')










