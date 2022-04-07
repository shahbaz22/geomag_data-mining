import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
# from aacgmv2 import get_aacgm_coord
from datetime import datetime
from tqdm import tqdm
import matplotlib.pyplot as plt
class NetDistr:
    def __init__(self, edge_list_dir: str, n: int):
        G = nx.read_edgelist(edge_list_dir)
        self.edge_list = list(G.edges(data=True))
        self.n = n
        station_data = pd.read_csv('supermag-stations.csv')
        # 0 index is of key is GEOLON, 1 is GEOLAT, 2 MLON, 3 MLAT
        self.stations_dict =  station_data.set_index('IAGA').T.to_dict('list')
        all_times = [s[2]['attr_dict']['UTC'] for s in self.edge_list]
        self.ordered_times = sorted(list(set(all_times)))
        self.date = self.edge_list[0][2]['attr_dict']['UTC'].split(' ')[0]


    def t_close(self, target_time: str) -> str:
        date_time = f'{self.date} {target_time}'
        date_format = '%Y-%m-%d %H:%M:%S'
        self.closest_time = min(self.ordered_times, key=lambda x: 
            abs(datetime.strptime(x, date_format) - datetime.strptime(date_time, date_format)))
        self.datetime_closest_time = datetime.strptime(self.closest_time,'%Y-%m-%d %H:%M:%S')
        return self.closest_time

    def coverage(self) -> pd.DataFrame:

        node_counter_dict = dict()
        cov_dict = {}
        for ind,edge in tqdm(enumerate(self.edge_list)):
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = Counter()
            # specifying/setting dict key and asinging count to counter 
            node_counter_dict[edge_time][edge[0]] += 1
            node_counter_dict[edge_time][edge[1]] += 1
        
        for time, vals in node_counter_dict.items():
            # print(time,vals)
            # break
            stations  = [station for station in vals.keys()]
            print(time,len(stations), len(set(stations)))
            cov_dict[time] = np.round(len(stations)/self.n, 2)

        df = pd.DataFrame.from_dict(cov_dict, orient='index')
        df = df.sort_index()
        print(df)

    def create_pivottable(self) -> pd.DataFrame:

        node_counter_dict = dict()
        for ind,edge in tqdm(enumerate(self.edge_list)):
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = Counter()
            # specifying/setting dict key and asinging count to counter 
            node_counter_dict[edge_time][edge[0]] += 1
            node_counter_dict[edge_time][edge[1]] += 1
        # print(node_counter_dict.keys())
        # for v in node_counter_dict.values():
        #     a = [i for i in v.keys()]

            
        pivot_table = pd.DataFrame([])
        for time, node_counter in node_counter_dict.items():
            # print(time,node_counter)
            new_row = pd.DataFrame(index=[time], columns=[np.arange(1,self.n)])
            pivot_table = pivot_table.append(new_row)

            degree_counts = Counter(list(node_counter.values()))
            # print(degree_counts)
            for degree_key, degree_count in degree_counts.items():
                # print(degree_key, degree_count)
                pivot_table.loc[time, degree_key] = degree_count

                
        pivot_table = pivot_table.sort_index()
        pivot_table.fillna(value=np.nan, inplace=True)
        return pivot_table


    def create_pivottable_by_mlt(self, mlt_lower_lim:float) -> pd.DataFrame:

        node_counter_dict = dict()
        for edge in tqdm(self.edge_list):

            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = Counter()

           # needs to be writen out because calling function many times slows down code
            n1_name = edge[0].split('_')[0]
            n2_name = edge[1].split('_')[0]
            glon1 = self.stations_dict[n1_name][2]
            glon2 = self.stations_dict[n2_name][2]
            dt = datetime.strptime(edge_time,'%Y-%m-%d %H:%M:%S')
            utc_hours = np.round(dt.hour + dt.minute/60, 2)
            glon1 = float(glon1)
            glon2 = float(glon2)
            geo_north_lon = -72
            mlt1 = (24 + (utc_hours + (glon1 + geo_north_lon)/15))
            mlt2 = (24 + (utc_hours + (glon2 + geo_north_lon)/15))
            mltrange = abs(mlt1-mlt2)

            if mlt_lower_lim < mltrange:
                node_counter_dict[edge_time][edge[0]] += 1
                node_counter_dict[edge_time][edge[1]] += 1
            
        pivot_table = pd.DataFrame([])
        for time, node_counter in node_counter_dict.items():
            new_row = pd.DataFrame(index=[time], columns=[np.arange(1,self.n)])
            pivot_table = pivot_table.append(new_row)

            degree_counts = Counter(list(node_counter.values()))
            for degree_key, degree_count in degree_counts.items():
                pivot_table.loc[time, degree_key] = degree_count

                
        pivot_table = pivot_table.sort_index()
        pivot_table.fillna(value=np.nan, inplace=True)
        return pivot_table

    def create_t_snapshot_distrobution(self, time:str, mlt_lower_lim:float) -> dict:

        node_counter_dict_mlt1 = Counter()
        node_counter_dict_mlt2 = Counter()
        date, *_ = self.edge_list[0][2]['attr_dict']['UTC'].split(' ')
        date_time = f'{date} {time}'

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time==date_time:
                node_counter_dict_mlt1[edge[0]] += 1
                node_counter_dict_mlt1[edge[1]] += 1

                n1_name = edge[0].split('_')[0]
                n2_name = edge[1].split('_')[0]
                glon1 = self.stations_dict[n1_name][2]
                glon2 = self.stations_dict[n2_name][2]
                dt = datetime.strptime(edge_time,'%Y-%m-%d %H:%M:%S')
                utc_hours = np.round(dt.hour + dt.minute/60, 2)
                glon1 = float(glon1)
                glon2 = float(glon2)
                geo_north_lon = -72
                mlt1 = (24 + (utc_hours + (glon1 + geo_north_lon)/15))
                mlt2 = (24 + (utc_hours + (glon2 + geo_north_lon)/15))
                mltrange = abs(mlt1-mlt2)

                if mlt_lower_lim < mltrange:
                    node_counter_dict_mlt2[edge[0]] += 1
                    node_counter_dict_mlt2[edge[1]] += 1

        degree_counts_mlt1 = Counter(list(node_counter_dict_mlt1.values()))
        degree_counts_mlt2 = Counter(list(node_counter_dict_mlt2.values()))

        return degree_counts_mlt1, degree_counts_mlt2

    def create_t_snapshot_connection_lengths(self, time:str) -> list:

        mlt_ranges = []
        date, *_ = self.edge_list[0][2]['attr_dict']['UTC'].split(' ')
        date_time = f'{date} {time}'

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time==date_time:
                n1_name = edge[0].split('_')[0]
                n2_name = edge[1].split('_')[0]
                # actually magnetic longitude
                glon1 = self.stations_dict[n1_name][2]
                glon2 = self.stations_dict[n2_name][2]
                dt = datetime.strptime(edge_time,'%Y-%m-%d %H:%M:%S')
                utc_hours = np.round(dt.hour + dt.minute/60, 2)
                glon1 = float(glon1)
                glon2 = float(glon2)
                geo_north_lon = -72
                mlt1 = (24 + (utc_hours + (glon1 + geo_north_lon)/15))
                mlt2 = (24 + (utc_hours + (glon2 + geo_north_lon)/15))
                mltrange = abs(mlt1-mlt2)
                mlt_ranges.append(mltrange)
                # create bins by mlt range

        return mlt_ranges
   
    def create_t_snapshot_lags(self, time:str) -> dict:

        connections_lag = []
        date, *_ = self.edge_list[0][2]['attr_dict']['UTC'].split(' ')
        date_time = f'{date} {time}'

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time==date_time:
                edge_lag = edge[2]['attr_dict']['lag']
                connections_lag.append(edge_lag)

        return connections_lag

    def station_coverage(self) -> pd.DataFrame:

        node_counter_dict = dict()
        for ind,edge in enumerate(self.edge_list):
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = []
            # specifying/setting dict key and asinging count to counter 
            n1 = edge[0].split('_')[0]
            n2 = edge[1].split('_')[0]
            node_counter_dict[edge_time].append([n1,n2])
            
        df_coverage = pd.DataFrame([])
        for time, nodes in node_counter_dict.items():
            unique_nodes = set(np.reshape(nodes,2*len(nodes)))
            new_row = pd.DataFrame(index=[time], columns=['num_stations'])
            df_coverage = df_coverage.append(new_row)
            df_coverage.loc[time,'num_stations'] = len(unique_nodes)
                
        df_coverage = df_coverage.sort_index()
        return df_coverage
   
    def ts_con_length_ratio_by_mlt(self, mlt_lower_lim:float) -> pd.DataFrame:

        edge_counter_dict = dict()
        for edge in tqdm(self.edge_list):

            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = {'long':0, 'short':0}

           # needs to be writen out because calling function many times slows down code
            n1_name = edge[0].split('_')[0]
            n2_name = edge[1].split('_')[0]
            glon1 = self.stations_dict[n1_name][2]
            glon2 = self.stations_dict[n2_name][2]
            dt = datetime.strptime(edge_time,'%Y-%m-%d %H:%M:%S')
            utc_hours = np.round(dt.hour + dt.minute/60, 2)
            glon1 = float(glon1)
            glon2 = float(glon2)
            geo_north_lon = -72
            mlt1 = (24 + (utc_hours + (glon1 + geo_north_lon)/15))
            mlt2 = (24 + (utc_hours + (glon2 + geo_north_lon)/15))
            mltrange = abs(mlt1-mlt2)

            if mlt_lower_lim < mltrange:
                edge_counter_dict[edge_time]['long'] += 1
            else:
                edge_counter_dict[edge_time]['short'] += 1

        edge_ratio_dict = {}

        for time, v in edge_counter_dict.items():
            edge_ratio_dict[time] = v['long']/v['short']
        
        edge_ratio = pd.DataFrame(edge_ratio_dict.items(), columns=['times','long/short'])
        edge_ratio['times'] =  pd.to_datetime(edge_ratio['times'], format='%Y-%m-%d %H:%M:%S')
        edge_ratio = edge_ratio.set_index('times')
        edge_ratio = edge_ratio.sort_index()

        return edge_ratio

    def ts_con_length_ratio_by_long(self, mlon_range_lim) -> pd.DataFrame:

        edge_counter_dict = dict()
        for edge in self.edge_list:

            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = {'long':0, 'short':0}

           # needs to be writen out because calling function many times slows down code
            n1_name = edge[0].split('_')[0]
            n2_name = edge[1].split('_')[0]
            # check longitiude
            mlon1 = self.stations_dict[n1_name][2]
            mlon2 = self.stations_dict[n2_name][2]

            mlon_range = abs(mlon1 - mlon2)
            if mlon_range_lim < mlon_range:
                edge_counter_dict[edge_time]['long'] += 1
            else:
                edge_counter_dict[edge_time]['short'] += 1

        edge_ratio_dict = {}

        for time, v in edge_counter_dict.items():
            edge_ratio_dict[time] = v['long']/v['short']
        edge_ratio = pd.DataFrame(edge_ratio_dict.items(), columns=['times','long/short'])
        edge_ratio['times'] =  pd.to_datetime(edge_ratio['times'], format='%Y-%m-%d %H:%M:%S')
        edge_ratio = edge_ratio.set_index('times')
        edge_ratio = edge_ratio.sort_index()

        return edge_ratio


    def t_snapshot_dist(self, times:list) -> dict:

        deg_dists = {}
        for t in times:
            nearest_t = self.t_close(t)
            deg_dists[nearest_t] = Counter()

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time in deg_dists.keys():
                deg_dists[edge_time][edge[0]] += 1
                deg_dists[edge_time][edge[1]] += 1

        for t in deg_dists.keys():
            deg_dists[t] = Counter(list(deg_dists[t].values()))

        return deg_dists

    def t_snapshot_dist_power(self, times:list) -> dict:

        power_dist = {}

        for t in times:
            nearest_t = self.t_close(t)
            power_dist[nearest_t] = []

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            n1_label = edge[0]
            n2_label = edge[1]
            if edge_time in power_dist.keys():
                power_dist[edge_time].append(n1_label)
                power_dist[edge_time].append(n2_label)

        for t in power_dist.keys():
            powers = [float(v.split('_')[2]) for v in list(set(power_dist[t]))]
            power_dist[t] = np.round(np.log10(powers),1)

        return power_dist

    def average_pc_ts(self) -> list:

        node_counter_dict = dict()
        for ind,edge in tqdm(enumerate(self.edge_list)):
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = []
            # specifying/setting dict key and asinging count to counter 
            node_counter_dict[edge_time].add(edge[0])
            node_counter_dict[edge_time].add(edge[1])

        # loop through times and average
            
        pivot_table = pd.DataFrame([])
        for time, node_counter in node_counter_dict.items():
            # print(time,node_counter)
            new_row = pd.DataFrame(index=[time], columns=[np.arange(1,self.n)])
            pivot_table = pivot_table.append(new_row)

            degree_counts = Counter(list(node_counter.values()))
            # print(degree_counts)
            for degree_key, degree_count in degree_counts.items():
                # print(degree_key, degree_count)
                pivot_table.loc[time, degree_key] = degree_count

                
        pivot_table = pivot_table.sort_index()
        pivot_table.fillna(value=np.nan, inplace=True)
        return pivot_table

def bounds_arrs(arrs:dict)->int:
    vals = [[max(i), min(i)] for i in arrs.values()]
    vals = np.array(vals).flatten()
    l_val = np.ceil(min(vals))
    u_val = np.ceil(max(vals))
    return l_val, u_val

def bounds_arrs_deg(arrs:dict)->int:
    # bnch of counter in vals
    vals = [i for i in arrs.values()]
    vals = [[max(i), min(i)] for i in arrs.values()]
    vals = np.array(vals).flatten()
    l_val = min(vals)
    u_val = max(vals)
    return l_val, u_val
# create a function to plot all histograms for given snapshots
path ='networks_data/2015_nets/comp_e'
file_handel = '_e_20_peakh_0.3_n_128_2015_0.2_2.0.txt'
file_names = ['dir_net','undir_net_inphase','undir_net_antiphase']
for i in [0,1]:
    for f in file_names:
        file_short = f'{f}{i}'
        file_full =f'{file_short}{file_handel}'
        dist = NetDistr(f'{path}/{file_full}', 128)
        # # pivottable.edge_list
        # can create a function which plots all the distrobutions with all the file names
        # and specified file path
        times = ['04:28:20', '04:45:50', '06:00:00' ]
        power_dist = dist.t_snapshot_dist_power(times)
        deg_dist =  dist.t_snapshot_dist(times)
        fig, ax = plt.subplots(nrows=2,ncols=3, facecolor='w', edgecolor='k',figsize=(10,8))
        fig.suptitle(" ".join(file_short.split('_')[0:3]))
        lb1, ub1 = bounds_arrs(power_dist)
        lb2, ub2 = bounds_arrs_deg(deg_dist)
        bins1 = np.linspace(lb1,ub1,20)
        labs = ['pre-onset','onset','post-onset']
        for ind, val in enumerate(power_dist.keys()):
            ax[0,ind].set_title(f'{val}\n{labs[ind]}')
            ax[0,ind].hist(power_dist[val], bins=bins1, edgecolor="red")
            ax[0,ind].set_xlabel('log(nT^2)')
            ax[1,ind].bar(deg_dist[val].keys(),deg_dist[val].values(), edgecolor="red", width=1.0 )
            ax[1,ind].set_xlim(lb2,ub2)
            ax[1,ind].set_xlabel('degree')
        ax[0,0].set_ylabel('freq')
        ax[1,0].set_ylabel('freq')
        plt.savefig(f'plots/dist_{file_short}.png')
        # plt.show()
# # er = dist.ts_con_length_ratio_by_mlt(4)
# # er.plot()
# # plt.show()
# fig, ax = plt.subplots(3,1)
# for ind, t in enumerate(['04:28:20', '04:45:50', '06:00:00' ]):
#     tc = dist.t_close(t)
#     print('tc', tc)
#     lags_arr = dist.create_t_snapshot_lags(t)
#     ax[ind].axvline(0, color='r')
#     if ind ==0:
#         bins = np.linspace(min(lags_arr),max(lags_arr),2)
 
#     # bins = np.linspace(min(lags_arr),max(lags_arr),10)
#     ax[ind].hist(lags_arr, bins=bins, edgecolor="k")
#     ax[ind].set_xticks(bins)

# plt.show()


# 04:45:50 UTC (SSC) both short and long range, directed/undirected After SSC the short 
# range populations bifurcate. One population 
# becomes highly connected at during SME enhancement around 06:00:00 UTC
# pt = pivottable.create_pivottable()

# pt = pivottable.create_pivottable_by_mlt(10)
# # print(g[1][0].split('_')[0])
# # pt = pivottable.create_pivottable()
# print(pt, type(pt))
# print(list(pt.index))
# ax = sb.heatmap(pt.T,cmap = sb.cm.rocket)
# ax.invert_yaxis()
# plt.grid()
# plt.ylabel('degrees')
# plt.show()
