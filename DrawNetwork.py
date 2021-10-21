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
import matplotlib.gridspec as gridspec
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os
import ray
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
        self.n = n
        self.directed = directed
        self.station_data = pd.read_csv('supermag-stations.csv')
        self.station_data.set_index('IAGA', inplace=True)
        self.geomag_poles_coord = geomag_poles_coord
        self.date = self.edge_list[0][2]['attr_dict']['UTC1'].split(' ')[0]
        self.fig = plt.figure(figsize=(14,7))
        self.gs= self.fig.add_gridspec(2, 2)


    def t_close(self, target_time: str) -> str:
        date_time = f'{self.date} {target_time}'
        date_format = '%Y-%m-%d %H:%M:%S'
        self.closest_time = min(self.ordered_times, key=lambda x: 
            abs(datetime.strptime(x, date_format) - datetime.strptime(date_time, date_format)))
        self.datetime_closest_time = datetime.strptime(self.closest_time,'%Y-%m-%d %H:%M:%S')
        return self.closest_time

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
        nx.draw_networkx_nodes(G = G, pos = pos, nodelist = G.nodes(), 
        node_color = 'r', alpha = 0.8, node_size = freqs, ax = ax)
        nx.draw_networkx_edges(G = G, pos = pos, edge_color='g',
        alpha=1, arrows = self.directed, ax = ax, nodelist = G.nodes())

    def t_snapshot(self, target_time: str, label:str, index:int, save=False) -> None:
        
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
        self.cartopy_ax(self.gs[0,:], date_time=target_time, G=G, available_stations=avail_stations, counts=counts, label = label, proj='globe')
        self.cartopy_ax(self.gs[1,0], date_time=target_time, G=G, available_stations=avail_stations, counts=counts, label = label, proj='south_pole')
        self.cartopy_ax(self.gs[1,1], date_time=target_time, G=G, available_stations=avail_stations, counts=counts, label = label, proj='north_pole')   
        if save==True:     
            plt.savefig(f'plots/movie_plots/t_snap_{label}_{index}_dir_net_z.png')
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
        os.system('rm {}t_snap_{}_*.png'.format(save_file_path,label))





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

for i in [0,1]:
    drawnet = DrawNetwork(f'networks_data/170313_all_available_stations_networks/dir_net{i}_z_20_peakh_0.3_num_stations_70_170313.txt', 70, geomag_poles_coord=dict_geomag_poles, directed=False)
    close_time = drawnet.t_close('04:40:00')
    drawnet.t_snapshot(close_time, label=f'pre-onset Pc{i+2}', index=1, save=True)
    close_time = drawnet.t_close('11:20:00')
    drawnet.t_snapshot(close_time, label=f'post-onset Pc{i+2}', index=1, save=True)

# drawnet.animation2(f'Pc2',event=2, time=drawnet.ordered_times)


# need to think cheak weather all station coordinates data needed for each plot
# for l in ['undir_net_inphase','undir_net_antiphase']: 
#     comp='e'
#     for i in [0,1]:
#         print(l,i, comp)

#         if i ==0 and l=='undir_net_inphase':
#             label= f'Unirected_Pc2_{comp}_inphase'

#         elif i ==0 and l=='undir_net_antiphase':
#             label= f'Unirected_Pc2_{comp}_antiphase'

#         elif i ==1 and l=='undir_net_inphase':
#             label= f'Undirected_Pc3_{comp}_inphase' 

#         elif i ==1 and l=='undir_net_antiphase':
#             label= f'Undirected_Pc3_{comp}_antiphase'   

#         drawnet = DrawNetwork(f'networks_data/170313_all_available_stations_networks/{l}{i}_{comp}_20_peakh_0.3_num_stations_70_170313.txt', 70, geomag_poles_coord=dict_geomag_poles, directed=False)
#         # close_time = drawnet.t_close('11:20:00')
#         # print(close_time)
#         # drawnet.t_snapshot(close_time, fig = plt.figure(figsize=(15,10)))
#         drawnet.animation2(label= f'{label}', event=2, time=drawnet.ordered_times)

# for l in ['undir_net_inphase','undir_net_antiphase']: 
#     for i in [1]:
#         print(l,i)
#         # if i ==0 and l=='dir_net':
#         #     label= 'Directed_Pc2_n'

#         # elif i ==0 and l=='undir_net_inphase':
#         #     label= 'Unirected_Pc2_n_inphase'

#         # elif i ==0 and l=='undir_net_antiphase':
#         #     label= 'Unirected_Pc2_n_antiphase'

#         # elif i ==1 and l=='dir_net':
#         #     label= 'Directed_Pc3_n'

#         if i ==1 and l=='undir_net_inphase':
#             label= 'Undirected_Pc3_n_inphase' 

#         elif i ==1 and l=='undir_net_antiphase':
#             label= 'Undirected_Pc3_n_antiphase'    

#         drawnet = DrawNetwork(f'networks_data/170313_all_available_stations_networks/{l}{i}_n_20_peakh_0.3_num_stations_70_170313.txt', 70, geomag_poles_coord=dict_geomag_poles, directed=False)
#         # close_time = drawnet.t_close('11:20:00')
#         # print(close_time)
#         # drawnet.t_snapshot(close_time, fig = plt.figure(figsize=(15,10)))
#         drawnet.animation2(label= f'{label}', event=2, time=drawnet.ordered_times)

# import numpy as np
# import numpy as np
# import pandas as pd
# import networkx as nx
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# from cartopy.feature.nightshade import Nightshade
# import cartopy.feature as cfeature
# from datetime import datetime
# from collections import Counter
# from celluloid import Camera
# import matplotlib.gridspec as gridspec
# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# import os
# import matplotlib as mpl

# # rework with cartopy
# # first want to draw network for single specific time, eventually want to run this 
# # plotting class a bunch of times (class usful to acces attributes and run methods)
# # Need to creat a dictonary which maps station name to (key)conda evs


# class DrawNetwork:
#     def __init__(self, edge_list_dir: str, n: int, directed: bool, geomag_poles_coord:dict):
#         G = nx.read_edgelist(edge_list_dir)
#         self.edge_list = list(G.edges(data=True))
#         all_times = [s[2]['attr_dict']['UTC1'] for s in self.edge_list]
#         self.ordered_times = sorted(list(set(all_times)))
#         self.n = n
#         self.directed = directed
#         self.station_data = pd.read_csv('supermag-stations.csv')
#         self.station_data.set_index('IAGA', inplace=True)
#         self.geomag_poles_coord = geomag_poles_coord
#         self.date = self.edge_list[0][2]['attr_dict']['UTC1'].split(' ')[0]
#         self.fig = plt.figure(figsize=(14,7))
#         self.gs= self.fig.add_gridspec(2, 2)


#     def t_close(self, target_time: str) -> str:
#         date_time = f'{self.date} {target_time}'
#         date_format = '%Y-%m-%d %H:%M:%S'
#         self.closest_time = min(self.ordered_times, key=lambda x: 
#             abs(datetime.strptime(x, date_format) - datetime.strptime(date_time, date_format)))
#         self.datetime_closest_time = datetime.strptime(self.closest_time,'%Y-%m-%d %H:%M:%S')
#         return self.closest_time

#     def color_map_pc(self, pc_powers:list)->list:
        
#         sorted_pc_powers = sorted(pc_powers)
#         print(sorted_pc_powers)
#         norm = mpl.colors.Normalize(vmin=min(sorted_pc_powers), vmax=max(sorted_pc_powers), clip=True)
#         mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)
#         c_transform = [mapper.to_rgba(i) for i in pc_powers]
#         return c_transform


#     def cartopy_ax(self, gridspec: plt.figure, date_time: str, G:nx.Graph, available_stations:list, proj: str, counts:dict, pc_powers:list, label:str) ->None:

#         year = self.date.split('-')[0]
#         gn_lat, gn_lon  = self.geomag_poles_coord[int(year)]['north']
#         gs_lat, gs_lon = self.geomag_poles_coord[int(year)]['south']
#         source_proj = ccrs.PlateCarree() 
#         if proj == 'globe':
#             map_proj = ccrs.EckertII(central_longitude=gn_lon)
#             ax = self.fig.add_subplot(gridspec, projection = map_proj)
#             title = ax.text(0.5,1.05, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5}, transform=ax.transAxes, ha="center")
#             title.set_text(f'{label}\n{date_time} UTC')
#         elif proj == 'south_pole':
#             map_proj = ccrs.Orthographic(central_longitude=gs_lon, central_latitude=gs_lat)
#             ax = self.fig.add_subplot(gridspec, projection = map_proj)
#         elif proj == 'north_pole':
#             map_proj = ccrs.Orthographic(central_longitude=gn_lon, central_latitude=gn_lat)
#             ax = self.fig.add_subplot(gridspec, projection = map_proj)
#         ax.stock_img()
#         ax.add_feature(cfeature.GSHHSFeature(edgecolor='k'))          
#         dt = datetime.strptime(date_time,'%Y-%m-%d %H:%M:%S')
#         ax.add_feature(Nightshade(dt, alpha=0.4))
#         ax.gridlines(source_proj, xlocs=np.linspace(-180,180,11), linestyle='--', draw_labels=False)
#         ax.scatter(x=[gn_lon,gs_lon], y=[gn_lat,gs_lat], color="orange", s=100,alpha=1, transform=source_proj)
#         pos = {}
#         pos = {station: (map_proj.transform_point(v['GEOLON'], v['GEOLAT'], src_crs=source_proj))
#            for station, v in
#            available_stations.to_dict('index').items()}
#         # str split here to obtain pc_power
#         freqs = [counts[station]**2 for station in G.nodes()]
#         freqs = np.array(freqs)*0.01/2
#         nx.draw_networkx_nodes(G = G, pos = pos, nodelist = G.nodes(), 
#         node_color = self.color_map_pc(pc_powers), alpha = 0.8, node_size = freqs, ax = ax)
#         nx.draw_networkx_edges(G = G, pos = pos, edge_color='g',
#         alpha=1, arrows = self.directed, ax = ax, nodelist = G.nodes())


#     def t_snapshot(self, target_time: str, label:str, index:int, pc_ind:int, save=False) -> None:
        
#         t_edge_list = []
#         for line in self.edge_list:
#             if line[2]['attr_dict']['UTC1'] == target_time:
#                 t_edge_list.append([line[0], line[1]])

#         if self.directed == True:
#             G = nx.DiGraph(t_edge_list)
#         else:
#             G = nx.Graph(t_edge_list)

#         reshaped_edge_list = np.reshape(t_edge_list, len(t_edge_list)*2)
#         node_list = list(set(reshaped_edge_list))
#         stations = [n.split('_')[0] for n in node_list]
#         pc_power_list = [float(n.split('_')[pc_ind+2]) for n in node_list]
#         counts = Counter(reshaped_edge_list)
#         avail_stations = self.station_data[self.station_data.index.isin(stations)]        
#         self.cartopy_ax(self.gs[0,:], date_time=target_time, G=G, available_stations=avail_stations, 
#             counts=counts, label = label, proj='globe', pc_powers=pc_power_list)
#         self.cartopy_ax(self.gs[1,0], date_time=target_time, G=G, available_stations=avail_stations, 
#             counts=counts, label = label, proj='south_pole', pc_powers=pc_power_list)
#         self.cartopy_ax(self.gs[1,1], date_time=target_time, G=G, available_stations=avail_stations, 
#             counts=counts, label = label, proj='north_pole', pc_powers=pc_power_list)   
#         if save==True:     
#             plt.savefig(f'plots/movie_plots/t_snap_{label}_{index}.png')
#         plt.show()
#         plt.clf()


#     def split_animation(self, chunks:int, label:str, event:str) ->None:
#         'for large animations, use frames less than 200'

#         all_times = np.array(self.ordered_times)
#         times_arr = np.array_split(all_times, chunks)
#         for ind, times in enumerate(times_arr):
#             print('chunk',ind)
#             self.animation(label, f'{event}_section_{ind}', time=times)
    
#     def animation1(self, label:str, event:str, time:list) ->None:

#         # plt.rcParams['animation.ffmpeg_path'] = '/Users/sc/anaconda3/envs/geomag_env/bin/ffmpeg'
#         # plt.rcParams['animation.ffmpeg_path'] = '~/anaconda3/envs/geomag/bin/ffmpeg'
#         camera = Camera(self.fig)
#         for n, t in enumerate(time):
#             print(t,n,len(time))
#             self.t_snapshot(t, label,index=None)
#             camera.snap()
#         animation = camera.animate()
#         animation.save(f'gifs/{label}_{event}.gif', fps=1)
    
#     def animation2(self, label:str, event:str, time:list) ->None:

#         for n, t in enumerate(time):
#             print(t,n,len(time))
#             self.t_snapshot(t, label,index=n, save=True)
        
#         # home settings
#         # save_file_path = '/Users/sc/conda_envs/geomag/plots/movie_plots/'
#         # ffmpeg_path = '/Users/sc/anaconda3/envs/geomag_env/bin/ffmpeg'
#         # remote settings
#         save_file_path = '/home/space/phrmfl/Shared/disk/_Users_sc_conda_envs_geomag/plots/movie_plots/'
#         ffmpeg_path = '~/.conda/envs/geomag/bin/ffmpeg'        

#         os.system('{} -r 4 -s 2000x2000 -f image2 -i \
#             plots/movie_plots/t_snap_{}_%01d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p \
#             {}{}_{}.mp4'.format(ffmpeg_path,label, save_file_path, label, event))
#         os.system('rm {}*.png'.format(save_file_path))





# # skipps bad lines
# # can show arrows, or look for conneced paths only can do a lot of analysis, need to find a way to show arrows

# # data taken from http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html 
# dict_geomag_poles = { 2000: {'north':   [79.6, -71.6] , 'south': [-79.6, 108.4]}, 
# 2005: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
# 2006: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
# 2007: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
# 2008: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
# 2009: {'north':   [79.8, -71.8] , 'south': [-79.8, 108.2]},
# 2010: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
# 2011: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
# 2012: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
# 2013: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
# 2014: {'north':   [80.1, -72.2] , 'south': [-80.1, 107.8]},
# 2015: {'north':   [80.4, -72.6] , 'south': [-80.4, 107.4]},
# 2016: {'north':   [80.4, -72.6] , 'south': [-80.4, 107.4]},
# 2017: {'north':   [80.5, -72.6] , 'south': [-80.5, 107.4]},
# 2018: {'north':   [80.5, -72.7] , 'south': [-80.5, 107.3]},
# 2019: {'north':   [80.6, -72.7] , 'south': [-80.6, 107.3]},
# 2020: {'north':   [80.7, -72.7] , 'south': [-80.7, 107.3]},
# 2021: {'north':   [80.7, -72.7] , 'south': [-80.7, 107.3]},
# 2022: {'north':   [80.7, -72.7] , 'south': [-80.7, 107.3]} }

# drawnet = DrawNetwork(f'networks_data/170313_all_available_stations_networks/dir_net0_e_20_peakh_0.3_num_stations_70_170313.txt', 70, geomag_poles_coord=dict_geomag_poles, directed=False)
# close_time = drawnet.t_close('11:20:00')
# # # print(close_time)
# drawnet.t_snapshot(close_time, label='test',index=None, pc_ind=0)
# # drawnet.animation2(f'Pc2',event=2, time=drawnet.ordered_times)


# # # need to think cheak weather all station coordinates data needed for each plot
# # for l in ['dna','na']: 
# #     for i in [0,1]:
# #         print(l,i)
# #         if i ==0 and l=='dna':
# #             label= 'Directed_Pc2_e'
# #         elif i ==0 and l=='na':
# #             label= 'Unirected_Pc2_e'
# #         elif i ==1 and l=='dna':
# #             label= 'Directed_Pc3_e'
# #         elif i ==1 and l=='na':
# #             label= 'Undirected_Pc3_e'    

# #         drawnet = DrawNetwork(f'networks_data/170313_all_available_stations_networks/{l}{i}_e_window_20_peakh_0.3_num_stations_70.txt', 70, geomag_poles_coord=dict_geomag_poles, directed=False)
# #         # close_time = drawnet.t_close('11:20:00')
# #         # print(close_time)
# #         # drawnet.t_snapshot(close_time, fig = plt.figure(figsize=(15,10)))
# #         if i ==0:
# #             drawnet.split_animation(chunks=20, label = f'{label}',event=2)
# #         if i ==1:
# #             drawnet.animation(label= f'{label}', event=2, time=drawnet.ordered_times)






