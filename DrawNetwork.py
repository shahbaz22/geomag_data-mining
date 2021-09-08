import numpy as np
import numpy as np
import pandas as pd
import networkx as nx
# import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from IPython.display import Image
from mpl_toolkits.basemap import Basemap as Basemap
import matplotlib.cbook
from datetime import datetime
from collections import Counter


# first want to draw network for single specific time, eventually want to run this plotting class a bunch of times (class usful to acces attributes and run methods)
# Need to creat a dictonary which maps station name to (key)conda evs
 




class DrawNetwork:
    def __init__(self, edge_list_dir: str, n: int, directed: bool):
        G = nx.read_edgelist(edge_list_dir)
        self.edge_list = sorted( list(G.edges(data=True)), key = lambda e: e[2]['attr_dict']['UTC1'])
        self.all_times = [t[2]['attr_dict']['UTC1'] for t in self.edge_list]
        # print(times_check)
        self.n = n
        self.directed = directed
        self.station_data = pd.read_csv('supermag-stations.csv')
        self.station_data.set_index('IAGA', inplace=True)



    def t_snapshot(self, target_time: str) -> None:
        date = self.edge_list[0][2]['attr_dict']['UTC1'].split(' ')[0]
        date_time = f'{date} {target_time}'
        # print(self.all_times)
        date_f = '%Y-%m-%d %H:%M:%S'
        closest_time = min(self.all_times, key=lambda x: 
            abs(datetime.strptime(x, date_f) - datetime.strptime(date_time, date_f)))
        # print(closest_time, type(closest_time))       
        t_edge_list = []
        # need to find time closest to, due to only having time windwows
        for line in self.edge_list:
            # can use indexer if needing multiple times
            if line[2]['attr_dict']['UTC1'] == date_time:
                t_edge_list.append([line[0].split('_')[0], line[1].split('_')[0]])

        if self.directed == True:
            G = nx.DiGraph(t_edge_list)
        else:
            G = nx.Graph(t_edge_list) 

        reshaped_edge_list = np.reshape(t_edge_list, len(t_edge_list)*2)

        stations = set(reshaped_edge_list)

        print('stations', stations)
        print('stations method', G.nodes())

        counts = Counter(reshaped_edge_list)

        print(self.station_data.head())

        available_stations = self.station_data[self.station_data.index.isin(stations)]

        # print('target_time',target_time)
        # print('closest_time',closest_time)
        # print(t_edge_list)

        # llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
        # are the lat/lon values of the lower left and upper right corners
        # of the map.
        # lat_ts is the latitude of true scale.
        # resolution = 'c' means use crude resolution coastlines.
        # m = Basemap(projection='merc' ,llcrnrlat=-85, urcrnrlat=85,
        #             llcrnrlon=-180, urcrnrlon=180, resolution='c')

        # maybe do a global network and a section of small pannels showing greater detail.

        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(15,15))

        lon_0 = -105; lat_0 = 40
        m1 = Basemap(ax = ax[0,0], projection='aeqd', lat_0=lat_0, lon_0=lon_0)
        # fill background.
        m1.drawmapboundary(fill_color='aqua')
        # draw coasts and fill continents.
        m1.drawcoastlines(linewidth=0.5)
        # m.fillcontinents(color='coral',lake_color='aqua')
        # 20 degree graticule.
        m1.drawparallels(np.arange(-80,81,20))
        m1.drawmeridians(np.arange(-180,180,20))
        m1.drawcountries(linewidth = 1)
        m1.drawcoastlines(linewidth = 1)

        m2 = Basemap(ax = ax[0,1],width=15.e6,height=15.e6,
            projection='gnom',lat_0=60,lon_0=-30)
        m2.drawmapboundary(fill_color='aqua')
        m2.drawcoastlines()
        # m.fillcontinents(color='coral',lake_color='aqua')
        m2.drawparallels(np.arange(10,90,20))
        m2.drawmeridians(np.arange(-180,180,30))
        m2.drawcountries(linewidth = 1)


        m3 = Basemap(ax = ax[1,0],width=15.e6,height=15.e6,
            projection='gnom',lat_0=-40.,lon_0=130)
        m3.drawcoastlines()
        m3.drawmapboundary(fill_color='aqua')
        # m.fillcontinents(color='coral',lake_color='aqua')
        m3.drawparallels(np.arange(10,90,20))
        m3.drawmeridians(np.arange(-180,180,30))
        m3.drawcountries(linewidth = 1)


        m4 = Basemap(ax = ax[1,1],width=15.e6,height=15.e6,
            projection='gnom',lat_0=60,lon_0=120)
        m4.drawcoastlines()
        m4.drawmapboundary(fill_color='aqua')
        # m.fillcontinents(color='coral',lake_color='aqua')
        m4.drawparallels(np.arange(10,90,20))
        m4.drawmeridians(np.arange(-180,180,30))
        m4.drawcountries(linewidth = 1)


        for ind1, ax in enumerate(ax.flat):
            m = [m1,m2,m3,m4]
            mx, my = m[ind1](available_stations['GEOLON'].values, available_stations['GEOLAT'].values)
            pos = {}
            print(ind1, ax)        
            for ind2, station in enumerate(available_stations.index):
                pos[station] = (mx[ind2], my[ind2])
            if ind1 ==0:
                freqs = [(1/600)*counts[station]**2.8 for station in G.nodes()]
            else:
                freqs = [0.01*counts[station]**2.8 for station in G.nodes()]
            nx.draw_networkx_nodes(G = G, pos = pos, node_list = G.nodes(), 
            node_color = 'r', alpha = 0.8, node_size = freqs, ax = ax)
            nx.draw_networkx_edges(G = G, pos = pos, edge_color='g',
            alpha=1, arrows = self.directed, ax = ax)

        fig.suptitle(f'Pc2, Network snapshot at time {closest_time}')
        plt.savefig('und_net_pc2_snapshot.png')

        plt.show()


        def net_animation(self) -> None:
            # TODO: plot animation using t_snapshot for all times
            pass

# skipps bad lines

# pd.set_option("display.max_rows", None, "display.max_columns", None)



# print(pos)

# need to think cheak weather all station coordinates data needed for each plot

drawnet = DrawNetwork('networks_data/170313_all_available_stations_networks/na0_e_window_20_peakh_0.3_num_stations_70.txt', 70, directed=False)
drawnet.t_snapshot('11:20:00')
# print(dir(DrawNetwork))
# print(drawnet.edge_list[0][2],type(drawnet.edge_list[0][2]))

# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# lat_ts is the latitude of true scale.
# resolution = 'c' means use crude resolution coastlines.

plt.show()
# print(station_data['IAGA'],len(station_data['IAGA']))



# d = {'blue':{'red': {1:[2,3,4], 2:[3,4,5,6,7]}},'orange':[1,2,3,4]}
# print(d.items())

# names = ('airline,airline_id,'
#          'source,source_id,'
#          'dest,dest_id,'
#          'codeshare,stops,equipment').split(',')
# routes = pd.read_csv(
#     'https://github.com/ipython-books/'
#     'cookbook-2nd-data/blob/master/'
#     'routes.dat?raw=true',
#     names=names,
#     header=None)

# names = ('id,name,city,country,iata,icao,lat,lon,'
#          'alt,timezone,dst,tz,type,source').split(',')
# airports = pd.read_csv(
#     'https://github.com/ipython-books/'
#     'cookbook-2nd-data/blob/master/'
#     'airports.dat?raw=true',
#     header=None,
#     names=names,
#     index_col=4,
#     na_values='\\N')
# airports_us = airports[airports['country'] ==
#                        'United States']
# print(type(airports_us), airports_us.columns)
# print(airports_us)
# # pos = {airport: (v['lon'], v['lat'])
# #        for airport, v in
# #        airports_us.to_dict().items()}

# # create new index
# airports_us_unique = airports_us.drop_duplicates(keep='first')
# print(airports_us_unique)

# for i, row in enumerate(airports_us.iteritems('iata',)):
#     print(i,row, type(row))

