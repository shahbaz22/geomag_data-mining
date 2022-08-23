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
from NetDistr import NetDistr
import string

# station_data = pd.read_csv('supermag-stations.csv')
# station_data.set_index('IAGA', inplace=True)
# mag_coords_dict1 = {}
# mag_coords_dict2 = {}
# mag_coords_dict_diff = {}
# # A = Apex(date=2015.3)
# # for station , v in station_data.to_dict('index').items():
# #     mlat, mlon = A.convert(v['GEOLAT'], v['GEOLON'], 'geo', 'apex', height=0)
# #     mag_coords_dict1[station] = (mlon, mlat)

# for station , v in station_data.to_dict('index').items():
#     mlat, mlon = v['AACGMLAT'], v['AACGMLON']
#     mag_coords_dict2[station] = (mlon, mlat)
# print(mag_coords_dict2)
# kk

# for station , v in station_data.to_dict('index').items():
#     coords1 = mag_coords_dict1[station]
#     coords2 = mag_coords_dict2[station]
#     mag_coords_dict_diff[station]= (abs(coords1[0])-abs(coords2[0]), abs(coords1[1])-abs(coords2[1]))

# diff_list = mag_coords_dict_diff.values()
# vals = max(diff_list,key=lambda item:item[1])
# print(vals)
# print(list(mag_coords_dict_diff.keys())[list(diff_list).index(vals)])
# print(mag_coords_dict1['TAM'])

# # rework with cartopy
# # first want to draw network for single specific time, eventually want to run this 
# # plotting class a bunch of times (class usful to acces attributes and run methods)
# # Need to creat a dictonary which maps station name to (key)conda evs


class DrawNetwork:
    def __init__(self, edge_list_dir: str, n: int, comp:str, directed: bool, geomag_poles_coord:dict, cluster:bool, conn_cluster:list):
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
        self.fig = plt.figure(figsize=(13,8))
        self.gs = self.fig.add_gridspec(2, 3)
        self.cluster = cluster
        self.conn_cluster = conn_cluster
        self.dist = NetDistr(edge_list_dir, self.n)
        self.comp= comp
        # self.apex_coords()

    def apex_coords(self):
        year = int(self.date.split('-')[0])
        month = int(self.date.split('-')[1])
        self.mag_coords_dict1 = {}
        A = Apex(date=f'{year}.{month}')
        for station , v in self.station_data.to_dict('index').items():
            mlat, mlon = A.convert(v['GEOLAT'], v['GEOLON'], 'geo', 'apex', height=0)
            self.mag_coords_dict1[station] = (mlon, mlat)
        return self.mag_coords_dict1

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

    def plot_deg_dist(self, gridspec: plt.figure, target_time:str, xlims:list, ylims:list):
        deg_dist1 =  self.dist.t_snapshot_dist([target_time])
        ax = self.fig.add_subplot(gridspec)
        ax.bar(deg_dist1[target_time].keys(), deg_dist1[target_time].values(), width=1.0)
        ax.set_ylabel('freq.')
        ax.set_xlabel('degree')
        ax.grid()
        ax.set_xlim(xlims)
        ax.set_ylim(ylims) 
        return ax


    def cartopy_ax(self, gridspec: plt.figure, date_time: str, G:nx.Graph, available_stations:list, proj: str, counts:dict, label:str, geographic:bool) ->None:
        year = self.date.split('-')[0]
        gn_lat, gn_lon  = self.geomag_poles_coord[int(year)]['north']
        gs_lat, gs_lon = self.geomag_poles_coord[int(year)]['south']
        dt = datetime.strptime(date_time,'%Y-%m-%d %H:%M:%S')
        source_proj = ccrs.PlateCarree() 
        if geographic == True:
            if proj == 'globe':
                map_proj = ccrs.PlateCarree(central_longitude=gn_lon)
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
            ax.add_feature(Nightshade(dt, alpha=0.4))
            ax.gridlines(source_proj, xlocs=np.linspace(-180,180,11), linestyle='--', draw_labels=False)
            ax.scatter(x=[gn_lon,gs_lon], y=[gn_lat,gs_lat], color="orange", s=100,alpha=1, transform=source_proj)
            pos = {station: (map_proj.transform_point(v['GEOLON'], v['GEOLAT'], src_crs=source_proj))
               for station, v in
               available_stations.to_dict('index').items()}
        elif proj == 'geomag':
            map_proj = ccrs.PlateCarree() 
            ax = self.fig.add_subplot(gridspec, projection = map_proj)
            title = ax.text(0.5,1.12, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5}, transform=ax.transAxes, ha="center")
            title.set_text('Geomagnetic coordinates')
            ax.set_extent([-180, 180, -90, 90], source_proj)
            ax.gridlines(source_proj, linestyle='--', draw_labels=True)
            ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=map_proj)
            ax.set_yticks([-78.5, -60, -25.5, 25.5, 60, 80], crs=map_proj)
            ax.xformatter = LONGITUDE_FORMATTER
            ax.yformatter = LATITUDE_FORMATTER
            ax.add_feature(Nightshade(dt, alpha=0.4))
            # pos = {s: self.mag_coords_dict1[s] for s in available_stations.to_dict('index')}
            pos = {station: (map_proj.transform_point(v['AACGMLON'], v['AACGMLAT'], src_crs=source_proj))
               for station, v in
               available_stations.to_dict('index').items()}
        # str split here to obtain pc_power
        freqs = [counts[station]**2 for station in G.nodes()]
        freqs = np.array(freqs)*0.01/2
        nx.draw_networkx_nodes(G = G, pos = pos, nodelist = G.nodes(), node_color='r', alpha = 0.8, node_size=freqs, ax = ax)
        nx.draw_networkx_edges(G = G, pos = pos, edge_color='g',
        alpha=1, arrows = self.directed, ax = ax, nodelist = G.nodes())

        if self.cluster==True:
            # reason for nodes in cluster not beeing connected to each other because they may be connected to another cluster in another group
            c_nodes = [station for station, conn in counts.items() if conn <= self.conn_cluster[1] and conn >= self.conn_cluster[0]]
            c_freqs = [(counts[station]**2)*0.01 for station in c_nodes]
            sub_G = G.subgraph(c_nodes)
            nx.draw_networkx_nodes(G = sub_G, pos = pos, nodelist = sub_G.nodes(), node_color='red', alpha = 0.8, node_size = 10, ax = ax)
            nx.draw_networkx_edges(G = sub_G, pos = pos, edge_color='lime',
            alpha=1, arrows = self.directed, ax = ax, nodelist = sub_G.nodes())
        return ax


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
        # self.cartopy_ax(self.gs[0,1], date_time=target_time, G=G, available_stations=avail_stations, 
        #     counts=counts, label = label, proj='globe', geographic=True)
        self.fig.suptitle(f'{label} {self.comp}\n{target_time}\nnumber of stations = {len(stations)}')
        ax1 = self.cartopy_ax(self.gs[0,0], date_time=target_time, G=G, available_stations=avail_stations, 
            counts=counts, label = label, proj='south_pole', geographic=True)
        ax2 = self.plot_deg_dist(self.gs[0,1], target_time, xlims=[0,120], ylims=[0,35])
        ax3 = self.cartopy_ax(self.gs[0,2], date_time=target_time, G=G, available_stations=avail_stations, 
            counts=counts, label = label, proj='north_pole', geographic=True) 
        ax4 = self.cartopy_ax(self.gs[1,:], date_time=target_time, G=G, available_stations=avail_stations, 
            counts=counts, label = label, proj='geomag', geographic=False)
        plt.subplots_adjust(wspace=0.2, hspace=0.4)
        for n, ax in enumerate([ax1,ax2,ax3,ax4]):
            ax.text(-0.1, 1.1, string.ascii_lowercase[n], transform=ax.transAxes, 
                    size=16, weight='bold')
        if save==True:     
            # plt.savefig(f'plots/movie_plots/t_snap_{label}_{index}.png')
            hr_min_sec = target_time.split(' ')[1]
            plt.savefig(f'plots/t_snap_{label}_{self.comp}_{hr_min_sec}.png')

        # plt.show()
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
# print('Enter event year, 2012, 2013, 2015:')
# year = input()
# print('Enter event component:')
# comp = input()
year='2015'
stations = 128
fhandel = f'{stations}_0.3_0.25_2.5_{year}.txt'
# print(station_names, len(station_names))
# conn_rangest1 = [[10,20],[80,105]]
# conn_rangest2  = [[10,30]]
# times = ['6:25:00', '7:45:00','9:30:00']
# 11:15:00 Bz becomes negative again
times = ['04:30:00','04:47:50', '06:00:00','06:45:00','09:15:00','09:33:00' ]
for comp in ['n']:
    path = f'networks_data/{year}_nets/comp_{comp}/'
    i=0
    for t in times:
        for net_type in ['undir_net_inphase']:
                if i ==1:
                    conn_range = [30,50]
                    drawnet = DrawNetwork(f'{path}combined_net{i}.txt',
                     128, comp, geomag_poles_coord=dict_geomag_poles, directed=False, cluster=False, conn_cluster=None)
                    label = 'combined_net'
                    print(net_type, comp, year, f'Pc{i+2}')
                    close_time = drawnet.t_close(t)
                    # label used for saving and title name
                    drawnet.t_snapshot(close_time, label=label, index=1, save=True) 
                    break
                else:
                    drawnet = DrawNetwork(f'{path}{net_type}{i}_{comp}_{fhandel}',
                     128, comp, geomag_poles_coord=dict_geomag_poles, directed=False,  cluster=False, conn_cluster=None)
                    label = f'{net_type}'
                    print(net_type, comp, year, f'Pc{i+2}')
                    close_time = drawnet.t_close(t)
                    drawnet.t_snapshot(close_time, label=label, index=1, save=True) 

          
            # drawnet.animation2(label= f'{net_type}_Pc{i+2}_{comp}_{year}',event=f'{net_type}_0.4', time=drawnet.ordered_times)



