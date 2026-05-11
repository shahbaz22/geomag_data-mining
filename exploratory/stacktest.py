import numpy as np 
import matplotlib.pyplot as plt

a=1/np.arange(1,5) #numpy array an only be used with numpy.stak
print(a)
b=1/np.arange(6,10)
print(b)
c=np.stack((a,b))
print()

wnds = 30

fig = plt.figure(figsize=(6, 3.2))

ax = fig.add_subplot(111)
ax.set_title('colorMap')
plt.imshow(c)
ax.set_aspect('equal')

# ax = fig.add_axes([0.12, 0.1, 0.78, 0.8])

plt.colorbar()
plt.show()