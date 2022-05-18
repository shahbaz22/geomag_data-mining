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
import matplotlib.gridspec as gridspec
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os
from os import listdir
from os.path import isfile, join
# # rework with cartopy
# # first want to draw network for single specific time, eventually want to run this 
# # plotting class a bunch of times (class usful to acces attributes and run methods)
# # Need to creat a dictonary which maps station name to (key)conda evs


class DrawNetwork:
    def __init__(self, edge_list_dir: str, n: int, directed: bool, geomag_poles_coord:dict):
        G = nx.read_edgelist(edge_list_dir)
        self.edge_list = list(G.edges(data=True))
        all_times = [s[2]['attr_dict']['UTC'] for s in self.edge_list]
        self.ordered_times = sorted(list(set(all_times)))
        self.n = n
        self.directed = directed
        self.station_data = pd.read_csv('supermag-stations.csv')
        self.station_data.set_index('IAGA', inplace=True)
        self.geomag_poles_coord = geomag_poles_coord
        self.date = self.edge_list[0][2]['attr_dict']['UTC'].split(' ')[0]
        self.fig = plt.figure(figsize=(14,7))
        self.gs= self.fig.add_gridspec(2, 2)


    def t_close(self, target_time: str) -> str:
        date_time = f'{self.date} {target_time}'
        date_format = '%Y-%m-%d %H:%M:%S'
        self.closest_time = min(self.ordered_times, key=lambda x: 
            abs(datetime.strptime(x, date_format) - datetime.strptime(date_time, date_format)))
        self.datetime_closest_time = datetime.strptime(self.closest_time,'%Y-%m-%d %H:%M:%S')
        return self.closest_time

    def color_map_pc(self, pc_powers:list)->list:
        
        sorted_pc_powers = sorted(pc_powers)
        print(sorted_pc_powers)
        norm = mpl.colors.Normalize(vmin=min(sorted_pc_powers), vmax=max(sorted_pc_powers), clip=True)
        mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)
        c_transform = [mapper.to_rgba(i) for i in pc_powers]
        return c_transform


    def cartopy_ax(self, gridspec: plt.figure, date_time: str, G:nx.Graph, available_stations:list, proj: str, counts:dict, label:str) ->None:

        year = self.date.split('-')[0]
        gn_lat, gn_lon  = self.geomag_poles_coord[int(year)]['north']
        gs_lat, gs_lon = self.geomag_poles_coord[int(year)]['south']
        source_proj = ccrs.PlateCarree() 
        if proj == 'globe':
            map_proj = ccrs.EckertII(central_longitude=gn_lon)
            ax = self.fig.add_subplot(gridspec, projection = map_proj)
            title = ax.text(0.5,1.05, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5}, transform=ax.transAxes, ha="center")
            title.set_text(f'{label}\n{date_time} UTC')
        elif proj == 'south_pole':
            map_proj = ccrs.Orthographic(central_longitude=gs_lon, central_latitude=gs_lat)
            ax = self.fig.add_subplot(gridspec, projection = map_proj)
        elif proj == 'north_pole':
            map_proj = ccrs.Orthographic(central_longitude=gn_lon, central_latitude=gn_lat)
            ax = self.fig.add_subplot(gridspec, projection = map_proj)
        ax.stock_img()
        ax.add_feature(cfeature.GSHHSFeature(edgecolor='k'))
        dt = datetime.strptime(date_time,'%Y-%m-%d %H:%M:%S')
        ax.add_feature(Nightshade(dt, alpha=0.4))
        ax.gridlines(source_proj, xlocs=np.linspace(-180,180,11), linestyle='--', draw_labels=False)
        ax.scatter(x=[gn_lon,gs_lon], y=[gn_lat,gs_lat], color="orange", s=100,alpha=1, transform=source_proj)
        pos = {}
        pos = {station: (map_proj.transform_point(v['GEOLON'], v['GEOLAT'], src_crs=source_proj))
           for station, v in
           available_stations.to_dict('index').items()}
        # str split here to obtain pc_power
        freqs = [counts[station]**2 for station in G.nodes()]
        freqs = np.array(freqs)*0.01/2
        nx.draw_networkx_nodes(G = G, pos = pos, nodelist = G.nodes(), node_color='r', alpha = 0.8, node_size = freqs, ax = ax)
        nx.draw_networkx_edges(G = G, pos = pos, edge_color='g',
        alpha=1, arrows = self.directed, ax = ax, nodelist = G.nodes())


    def t_snapshot(self, target_time: str, label:str, index:int, save=False) -> None:
        
        t_edge_list = []
        for line in self.edge_list:
            if line[2]['attr_dict']['UTC'] == target_time:
                t_edge_list.append([line[0].split('_')[0], line[1].split('_')[0]])

        if self.directed == True:
            G = nx.DiGraph(t_edge_list)
        else:
            G = nx.Graph(t_edge_list)

        reshaped_edge_list = np.reshape(t_edge_list, len(t_edge_list)*2)
        stations = list(set(reshaped_edge_list))
        # print(stations)

        counts = Counter(reshaped_edge_list)
        avail_stations = self.station_data[self.station_data.index.isin(stations)]        
        self.cartopy_ax(self.gs[0,:], date_time=target_time, G=G, available_stations=avail_stations, 
            counts=counts, label = label, proj='globe')
        self.cartopy_ax(self.gs[1,0], date_time=target_time, G=G, available_stations=avail_stations, 
            counts=counts, label = label, proj='south_pole')
        self.cartopy_ax(self.gs[1,1], date_time=target_time, G=G, available_stations=avail_stations, 
            counts=counts, label = label, proj='north_pole')   
        if save==True:     
            plt.savefig(f'plots/movie_plots/t_snap_{label}_{index}.png')
        plt.clf()


    def split_animation(self, chunks:int, label:str, event:str) ->None:
        'for large animations, use frames less than 200'

        all_times = np.array(self.ordered_times)
        times_arr = np.array_split(all_times, chunks)
        for ind, times in enumerate(times_arr):
            print('chunk',ind)
            self.animation(label, f'{event}_section_{ind}', time=times)
    
    def animation1(self, label:str, event:str, time:list) ->None:

        # plt.rcParams['animation.ffmpeg_path'] = '/Users/sc/anaconda3/envs/geomag_env/bin/ffmpeg'
        # plt.rcParams['animation.ffmpeg_path'] = '~/anaconda3/envs/geomag/bin/ffmpeg'
        camera = Camera(self.fig)
        for n, t in enumerate(time):
            print(t,n,len(time))
            self.t_snapshot(t, label,index=None)
            camera.snap()
        animation = camera.animate()
        animation.save(f'gifs/{label}_{event}.gif', fps=1)
    
    def animation2(self, label:str, event:str, time:list) ->None:

        for n, t in enumerate(time):
            print(t,n,len(time))
            self.t_snapshot(t, label,index=n, save=True)
        
        # home settings
        # save_file_path = '/Users/sc/conda_envs/geomag/plots/movie_plots/'
        # ffmpeg_path = '/Users/sc/anaconda3/envs/geomag_env/bin/ffmpeg'
        # remote settings
        save_file_path = '/home/space/phrmfl/Shared/disk/_Users_sc_conda_envs_geomag/plots/movie_plots/'
        ffmpeg_path = '~/.conda/envs/geomag/bin/ffmpeg'        

        os.system('{} -r 4 -s 2000x2000 -f image2 -i \
            plots/movie_plots/t_snap_{}_%01d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p \
            {}{}_{}.mp4'.format(ffmpeg_path,label, save_file_path, label, event))
        os.system('rm {}*.png'.format(save_file_path))



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
file_path = 'networks_data/2015'
station_names = [f.split('_')[1] for f in listdir(file_path) if isfile(join(file_path, f))]
year='2015'
fhandel = '128_0.4_0.25_2.5_2015.txt'
path = 'networks_data/2015_nets/comp_e/'
print(station_names, len(station_names))
for comp in ['e']:
    for net_type in ['dir_net','undir_net_inphase','undir_net_antiphase']:
        for i in [0,1]:
            print(net_type, comp, year, f'Pc{i+2}')
            drawnet = DrawNetwork(f'{path}{net_type}{i}_{comp}_{fhandel}', 128, geomag_poles_coord=dict_geomag_poles, directed=False)
            # close_time = drawnet.t_close('04:46:00')
            # drawnet.t_snapshot(close_time, label=f'post-onset Pc{i+2}', index=1, save=True)           
            drawnet.animation2(label= f'{net_type}_Pc{i+2}_{comp}_{year}',event=f'{net_type}_0.4', time=drawnet.ordered_times)







