import numpy as np
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
import cartopy.feature as cfeature
from datetime import datetime
from collections import Counter
from celluloid import Camera
# rework with cartopy
# first want to draw network for single specific time, eventually want to run this 
# plotting class a bunch of times (class usful to acces attributes and run methods)
# Need to creat a dictonary which maps station name to (key)conda evs


class DrawNetwork:
    def __init__(self, edge_list_dir: str, n: int, directed: bool, geomag_poles_coord:dict):
        G = nx.read_edgelist(edge_list_dir)
        self.edge_list = list(G.edges(data=True))
        all_times = [s[2]['attr_dict']['UTC1'] for s in self.edge_list]
        self.ordered_times = sorted(list(set(all_times)))
        # print(times_check)
        self.n = n
        self.directed = directed
        self.station_data = pd.read_csv('supermag-stations.csv')
        self.station_data.set_index('IAGA', inplace=True)
        self.geomag_poles_coord = geomag_poles_coord
        self.date = self.edge_list[0][2]['attr_dict']['UTC1'].split(' ')[0]


    def t_close(self, target_time: str) -> str:
        date_time = f'{self.date} {target_time}'
        date_format = '%Y-%m-%d %H:%M:%S'
        self.closest_time = min(self.all_times, key=lambda x: 
            abs(datetime.strptime(x, date_format) - datetime.strptime(date_time, date_format)))
        self.datetime_closest_time = datetime.strptime(self.closest_time,'%Y-%m-%d %H:%M:%S')
        return self.closest_time

    def cartopy_ax(self, fig: plt.figure, date_time: str, proj: str, G:nx.Graph, available_stations:pd.DataFrame, counts:dict, label:str) ->None:

        year = self.date.split('-')[0]
        gn_lat, gn_lon  = self.geomag_poles_coord[int(year)]['north']
        gs_lat, gs_lon = self.geomag_poles_coord[int(year)]['south']
        source_proj = ccrs.PlateCarree() 
        map_proj = ccrs.EckertII(central_longitude=gn_lon)
        ax = plt.axes(projection = map_proj)
        ax.stock_img()
        ax.add_feature(cfeature.GSHHSFeature(edgecolor='k'))          
        node_size_factor = 0.01
        dt = datetime.strptime(date_time,'%Y-%m-%d %H:%M:%S')
        ax.add_feature(Nightshade(dt, alpha=0.2))
        ax.scatter(x=[gn_lon,gs_lon], y=[gn_lat,gs_lat], color="orange", s=100,alpha=1, transform=source_proj)
        title = ax.text(0.5,1.05, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5}, transform=ax.transAxes, ha="center")
        title.set_text(f'{label}\n{date_time} UTC')
        pos = {}
        pos = {station: (map_proj.transform_point(v['GEOLON'], v['GEOLAT'], src_crs=source_proj))
           for station, v in
           available_stations.to_dict('index').items()}
        freqs = [node_size_factor*counts[station]**2.8 for station in G.nodes()]
        nx.draw_networkx_nodes(G = G, pos = pos, nodelist = G.nodes(), 
        node_color = 'r', alpha = 0.8, node_size = freqs, ax = ax)
        nx.draw_networkx_edges(G = G, pos = pos, edge_color='g',
        alpha=1, arrows = self.directed, ax = ax, nodelist = G.nodes(),)

    def t_snapshot(self, target_time: str, fig:plt.figure, label:str) -> None:
        
        t_edge_list = []
        # # need to find time closest to, due to only having time windwows
        for line in self.edge_list:
            # can use indexer if needing multiple times
            if line[2]['attr_dict']['UTC1'] == target_time:
                t_edge_list.append([line[0].split('_')[0], line[1].split('_')[0]])

        if self.directed == True:
            G = nx.DiGraph(t_edge_list)
        else:
            G = nx.Graph(t_edge_list)

        reshaped_edge_list = np.reshape(t_edge_list, len(t_edge_list)*2)
        stations = set(reshaped_edge_list)
        counts = Counter(reshaped_edge_list)
        avail_stations = self.station_data[self.station_data.index.isin(stations)]        
        self.cartopy_ax(fig, date_time=target_time, G=G, available_stations=avail_stations, proj='gnom', counts=counts, label = label)
        # plt.savefig(f'plots/t_snap_{label}.png')
        # plt.show()


    def animation(self, label:str) ->None:

        plt.rcParams['animation.ffmpeg_path'] = '/Users/sc/anaconda3/envs/geomag_env/bin/ffmpeg'
        fig = plt.figure(figsize=(10,6))
        camera = Camera(fig)
        for n, t in enumerate(self.ordered_times[0:10]):
            print(t,n,len(self.ordered_times))
            self.t_snapshot(t,fig,label)
            camera.snap()
        animation = camera.animate()
        animation.save('animation.mp4')

# skipps bad lines
# can show arrows, or look for conneced paths only can do a lot of analysis, need to find a way to show arrows

# data taken from http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html 
dict_geomag_poles = { 2000: {'north':   [79.6, -71.6] , 'south': [-79.6, 108.4]}, 
2005: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
2006: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
2007: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
2008: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
2009: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
2010: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
2011: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
2012: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
2013: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
2014: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
2015: {'north':   [80.4, -72.6] , 'south': [-80.4, 107.4]},
2016: {'north':   [80.4, -72.6] , 'south': [-80.4, 107.4]},
2017: {'north':   [80.5, -72.6] , 'south': [-80.5, 107.4]},
2018: {'north':   [80.5, -72.7] , 'south': [-80.5, 107.3]},
2019: {'north':   [80.6, -72.7] , 'south': [-80.6, 107.3]},
2020: {'north':   [80.7, -72.7] , 'south': [-80.7, 107.3]},
2021: {'north':   [80.7, -72.7] , 'south': [-80.7, 107.3]},
2022: {'north':   [80.7, -72.7] , 'south': [-80.7, 107.3]} }

# print(pos)

# need to think cheak weather all station coordinates data needed for each plot

drawnet = DrawNetwork('networks_data/170313_all_available_stations_networks/dna0_e_window_20_peakh_0.3_num_stations_70.txt', 70, geomag_poles_coord=dict_geomag_poles, directed=False)
# close_time = drawnet.t_close('11:20:00')
# print(close_time)
# drawnet.t_snapshot(close_time, fig = plt.figure(figsize=(15,10)))
drawnet.animation('Directed network Pc2')
# drawnet.net_animation()
# print(dir(DrawNetwork))
# print(drawnet.edge_list[0][2],type(drawnet.edge_list[0][2]))

# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# lat_ts is the latitude of true scale.
# resolution = 'c' means use crude resolution coastlines.

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

