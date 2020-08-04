# import matplotlib as mpl
# from matplotlib import pyplot as plt
# from matplotlib.patches import Wedge, Circle
import numpy as np
from geopack import geopack, t89,t96,t01,t04
# from datetime import timedelta
# import datetime
# import pandas as pd



# def dual_half_circle(center=(0,0), radius=1, angle=90, ax=None, colors=('w','k','k'),
#                      **kwargs):
#     """
#     Add two half circles to the axes *ax* (or the current axes) with the 
#     specified facecolors *colors* rotated at *angle* (in degrees).
#     """
#     if ax is None:
#         ax = plt.gca()
#     theta1, theta2 = angle, angle + 180
#     #w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], **kwargs)
#     #w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], **kwargs)
    
#     w1 = Wedge(center, radius, theta1, theta2, fc=colors[1], **kwargs)
#     w2 = Wedge(center, radius, theta2, theta1, fc=colors[0], **kwargs)
   
#     cr = Circle(center, radius, fc=colors[2], fill=False, **kwargs)
#     for wedge in [w1, w2, cr]:
#         ax.add_artist(wedge)
#     return [w1, w2, cr]

# def setup_fig(xlim=(10,-30),ylim=(-20,20),xlabel='X GSM [Re]',ylabel='Z GSM [Re]'):

#     fig = plt.figure(figsize=(15,10))
#     ax  = fig.add_subplot(111)
#     ax.axvline(0,ls=':',color='k')
#     ax.axhline(0,ls=':',color='k')
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
    
#     ax.set_aspect('equal')
#     w1,w2,cr = dual_half_circle(ax=ax)
    
#     return ax

def apdvalue(obs_lat, obs_lon, ut, dirv=1):

    '''link to all var info https://github.com/tsssss/geopack/blob/master/geopack/geopack.py
    code to find l-number for particular magnetometer station, ut in sec 
    | par |  1   |    2    |    3    |    4    |    5    |    6    |  7   |
    | Kp  | 0,0+ | 1-,1,1+ | 2-,2,2+ | 3-,3,3+ | 4-,4,4+ | 5-,5,5+ | > 6- |'''


    
    ps = geopack.recalc(ut)

    lat_rad = np.deg2rad(obs_lat)
    lon_rad = np.deg2rad(obs_lon)

    # Convert Geodetic to Geocentric Spherical
    r, theta_rad = geopack.geodgeo(0,lat_rad,1)

    # Convert Geocentric Spherical to Geocentric Cartesian
    x_gc,y_gc,z_gc = geopack.sphcar(1,theta_rad,lon_rad,1)

    # print('GC:  ', x_gc,y_gc,z_gc,' R=',np.sqrt(x_gc**2+y_gc**2+z_gc**2))

    # Convert Geocentric Cartesian to GSM
    x_gsm,y_gsm,z_gsm = geopack.geogsm(x_gc,y_gc,z_gc,1)
    # print('GSM for station: ', x_gsm,y_gsm,z_gsm,' R=',np.sqrt(x_gsm**2+y_gsm**2+z_gsm**2))

    x,y,z,xx,yy,zz = geopack.trace(x_gsm,y_gsm,z_gsm,dir=dirv,rlim=100,r0=.99999,parmod=5,exname='t89',inname='igrf',maxloop=10000)


    # The x-axis of the GSM coordinate system is defined along the line connecting the center of the Sun to the center of the Earth.
    #  The origin is defined at the center of the Earth, and is positive towards the Sun. 
    # The y-axis is defined as the cross product of the GSM x-axis and the magnetic dipole axis; directed positive towards dusk. The z-axis is defined as the cross product of the x- and y-axes. The magnetic dipole axis lies within the xz plane.

    # ax=setup_fig() #for plotting field lines
    # ax.plot(xx,zz)
    # plt.show()

    absxx=np.abs(xx)
    l=np.amax(absxx)


    return l

apdvalue(62.82, 267.89, 1332595679)

# time range in hours
# th=24

# lvs = []
# t = []
# t1 = []

# # default values set for intial RAN data set using UTC time and shifting x axis
# for i in np.arange(0,th,1):
#     print(i)
#     lvs.append(apdvalue(62.82, 267.89, 1332584880 +i*3600))
        
#     ts = pd.to_datetime('2012-03-24 10:28:00', format='%Y/%m/%d %H:%M:%S') - timedelta(hours=6.1 - i)

#     ts = str(ts)

#     ts = ts.rsplit(' ')

#     t.append(ts[1])

#     t1.append(4.28+i)

#     print(ts)


# rt=np.arange(0,th,1)
# r=np.arange(0,th,1)


# plt.plot(t1,lvs)
# plt.xticks(t1, t, rotation=70) # copy over the locations of the x-ticks from the first axes
# plt.xlabel('time [24 hour]')
# plt.ylabel('Apex distance (Re)')
# plt.savefig("ApexdistanceRANintialdatakp26.png",bbox_inches='tight')
# plt.show()


