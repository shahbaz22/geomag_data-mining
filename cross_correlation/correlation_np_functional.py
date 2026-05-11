import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

npts = 500
x = np.linspace(0, 120, npts)
y1 =  np.sin(x/2) 
y = 5 * np.cos(x/2)
y3 = signal.gaussian(npts, std=30)
# + np.random.randn(npts)


lags = np.arange(-npts + 1, npts)

#lags matches the dimension of ccov
ccov = np.correlate(y1 - y1.mean(), y3 - y3.mean(), mode='full')
print('y1.mean',y1.mean())
print('y3.mean',y3.mean())
# ‘full’ returns the convolution at each point of overlap
ccor = ccov / (npts * y1.std() * y3.std())

fig, axs = plt.subplots(nrows=2)
fig.subplots_adjust(hspace=0.4)
ax = axs[0]
ax.plot(x, y1, 'b', label='y1')
ax.plot(x, y3, 'r', label='y2')
# ax.set_ylim(-10, 10)
ax.legend(loc='upper right', fontsize='small', ncol=2)
maxlag = lags[np.argmax(ccor)]
print(np.argmax(ccor))
ax.set_title(f'max correlation at lag {maxlag}')

ax = axs[1]
ax.plot(lags, ccor)
# ax.set_ylim(-1.1, 1.1)
ax.set_ylabel('cross-correlation')
ax.set_xlabel('lag of y1 relative to y2')
plt.show()


