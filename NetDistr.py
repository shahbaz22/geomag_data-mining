import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
# from aacgmv2 import get_aacgm_coord
from datetime import datetime
from tqdm import tqdm
from math import radians, cos, sin, asin, sqrt

class NetDistr:
    def __init__(self, edge_list_dir: str, n:int):
        G = nx.read_edgelist(edge_list_dir)
        self.edge_list = list(G.edges(data=True))
        self.n = n
        station_data = pd.read_csv('supermag-stations.csv')
        self.station_data = station_data.set_index('IAGA')
                # 0 index is of key is GEOLON, 1 is GEOLAT, 2 MLON, 3 MLAT
        self.num_edges = len(self.edge_list)
        if self.num_edges>2:
            self.stations_dict = self.station_data.T.to_dict('list')
            all_times = [s[2]['attr_dict']['UTC'] for s in self.edge_list]
            self.ordered_times = sorted(list(set(all_times)))

    def geo_coords_dict(self):
        geo_coords_dict = {}
        for station , v in self.station_data.to_dict('index').items():
            glat, glon = v['GEOLAT'], v['GEOLON']
            geo_coords_dict[station] = (glon, glat)
        return geo_coords_dict    

    def mag_coords_dict(self):
        mag_coords_dict = {}
        for station , v in self.station_data.to_dict('index').items():
            mlat, mlon = v['AACGMLAT'], v['AACGMLON']
            mag_coords_dict[station] = (mlon, mlat)
        return mag_coords_dict

    def t_close(self, target_time: str) -> str:
        date_time = f'{self.date} {target_time}'
        date_format = '%Y-%m-%d %H:%M:%S'
        self.closest_time = min(self.ordered_times, key=lambda x: 
            abs(datetime.strptime(x, date_format) - datetime.strptime(date_time, date_format)))
        self.datetime_closest_time = datetime.strptime(self.closest_time,'%Y-%m-%d %H:%M:%S')
        return self.closest_time
    
    def complete_t_arr(self, dt_start:str, dt_end:str, window_multiple:int, pc:int) ->list:
        t_arr = pd.date_range(f'{dt_start}', f'{dt_end}', freq='s')
        pc_period = [5,10,45,150,600]
        # custom window sizes for each band using the lower limit for each pc band period range
        window_size_a = np.multiply(np.diff(pc_period), window_multiple)
        #custom step_size for different Pc bands multiplied by amount of window overlap // ensures whole number eg. *3//4 means 25% overlap
        step_size_a = np.zeros(len(window_size_a))
        step_size_a[0:2] = (window_size_a[0:2] * 1) // 2
        # ensure stepsize is an integer
        step_size_a = step_size_a.astype(np.int32)
        # first window
        t_start = 0
        t_end = t_start + window_size_a[pc]
        t_mid = t_end//2
        step_size = step_size_a[pc]
        # indices for time values
        t_mid_arr = []
        t_mid_arr.append(t_mid)
        while t_end <= len(t_arr):
            t_start = t_start + step_size                
            t_end = t_end + step_size
            t_mid = t_mid +step_size
            t_mid_arr.append(t_mid)
        return t_arr[t_mid_arr]


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
            stations  = [station for station in vals.keys()]
            # print(time,len(stations), len(set(stations)))
            cov_dict[time] = np.round(len(stations)/self.n, 2)

        df = pd.DataFrame.from_dict(cov_dict, orient='index')
        df = df.sort_index()

    def create_pivottable(self, start_dt:str, end_dt:str, pc:str) -> pd.DataFrame:

        node_counter_dict = dict()
        for ind,edge in tqdm(enumerate(self.edge_list)):
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = Counter()
            # specifying/setting dict key and asinging count to counter 
            node_counter_dict[edge_time][edge[0]] += 1
            node_counter_dict[edge_time][edge[1]] += 1

        pivot_table = pd.DataFrame([])
        # add start and end-times used to create equally spaced tracer for 
        # pivot table and heatmap most useful for large window-sizes i.e Pc3
        # if wndow size and step size is the same than I can create a uniform time array
        if pc == 'pc3':
            times = self.complete_t_arr(start_dt, end_dt,20, pc=1)
        else:
            times = self.complete_t_arr(start_dt, end_dt,20, pc=0)
        for t in times:
            t = datetime.strftime(t,'%Y-%m-%d %H:%M:%S')
            new_row = pd.DataFrame(index=[t], columns=[np.arange(1,self.n)])
            pivot_table = pivot_table.append(new_row)
            pivot_table.loc[t, 1] = 1

        for time, node_counter in node_counter_dict.items():
            new_row = pd.DataFrame(index=[time], columns=[np.arange(1,self.n)])
            
            degree_counts = Counter(list(node_counter.values()))
            for degree_key, degree_count in degree_counts.items():
                pivot_table.loc[time, degree_key] = degree_count

        pivot_table = pivot_table.sort_index()
        pivot_table.fillna(value=np.nan, inplace=True)
        return pivot_table

    def create_pivottable_lags(self, start_dt:str, end_dt:str, pc:str, period_mult:int, maxlags:list) -> pd.DataFrame:
        
        if pc == 'pc3':
            times = self.complete_t_arr(start_dt, end_dt,20, pc=1)
            period = 45*period_mult
            max_lag = maxlags[1]

        else:
            times = self.complete_t_arr(start_dt, end_dt,20, pc=0)
            period = 10*period_mult
            max_lag = maxlags[0]

        lag_counter_dict = dict()
        for ind,edge in tqdm(enumerate(self.edge_list)):
            edge_time = edge[2]['attr_dict']['UTC']
            edge_lag = edge[2]['attr_dict']['lag']

            if edge_time not in lag_counter_dict:
                lag_counter_dict[edge_time] = Counter()
            # specifying/setting dict key and asinging count to counter 
            lag_counter_dict[edge_time][edge_lag] += 1

        pivot_table = pd.DataFrame([])
        # add start and end-times used to create equally spaced tracer for 
        # pivot table and heatmap most useful for large window-sizes i.e Pc3
        # if wndow size and step size is the same than I can create a uniform time array

        for t in times:
            t = datetime.strftime(t,'%Y-%m-%d %H:%M:%S')
            new_row = pd.DataFrame(index=[t], columns=[np.arange(-period,period)])
            pivot_table = pivot_table.append(new_row)
            pivot_table.loc[t, 2] = 0
        for time, lag_counter in lag_counter_dict.items():
            new_row = pd.DataFrame(index=[time], columns=[np.arange(-period,period)])

            for lag_key, lag_count in lag_counter.items():
                if abs(lag_key) ==1 or lag_key ==0 or abs(lag_key) ==2:
                    if lag_count > max_lag:
                        pivot_table.loc[time, lag_key] = max_lag
                        continue
                    else:
                        pivot_table.loc[time, lag_key] =lag_count
                        continue
    
                pivot_table.loc[time, lag_key] = lag_count

        pivot_table = pivot_table.sort_index()
        pivot_table.fillna(value=np.nan, inplace=True)
        return pivot_table
    
    def create_pivottable_lags_mlt(self, mlt_low_lim:float, start_dt:str, end_dt:str, pc:str, period_mult:int) -> pd.DataFrame:
        
        if pc == 'pc3':
            times = self.complete_t_arr(start_dt, end_dt,20, pc=1)
            period = 45*period_mult
        else:
            times = self.complete_t_arr(start_dt, end_dt,20, pc=0)
            period = 10*period_mult

        lag_counter_dict = dict()
        for ind,edge in tqdm(enumerate(self.edge_list)):

            edge_time = edge[2]['attr_dict']['UTC']
            edge_lag = edge[2]['attr_dict']['lag']
            if edge_time not in lag_counter_dict:
                lag_counter_dict[edge_time] = Counter()            
            
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

            if mlt_low_lim < mltrange:
                lag_counter_dict[edge_time][edge_lag] += 1

        pivot_table = pd.DataFrame([])

        for t in times:
            t = datetime.strftime(t,'%Y-%m-%d %H:%M:%S')
            new_row = pd.DataFrame(index=[t], columns=[np.arange(-period,period)])
            pivot_table = pivot_table.append(new_row)
            pivot_table.loc[t, 2] = 0
        for time, lag_counter in lag_counter_dict.items():
            new_row = pd.DataFrame(index=[time], columns=[np.arange(-period,period)])

            for lag_key, lag_count in lag_counter.items():
                pivot_table.loc[time, lag_key] = lag_count

        pivot_table = pivot_table.sort_index()
        pivot_table.fillna(value=np.nan, inplace=True)
        return pivot_table

    def create_pivottable_by_mlt(self, mlt_low_lim:float, start_dt:str, end_dt:str, pc:str) -> pd.DataFrame:
        geo_cords = self.geo_coords_dict()
        node_counter_dict = dict()
        for edge in tqdm(self.edge_list):

            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = Counter()
           # needs to be writen out because calling function many times slows down code
            n1_name = edge[0].split('_')[0]
            n2_name = edge[1].split('_')[0]
            glon1 = geo_cords[n1_name][0]
            glon2 = geo_cords[n2_name][0]
            dt = datetime.strptime(edge_time,'%Y-%m-%d %H:%M:%S')
            utc_hours = np.round(dt.hour + dt.minute/60, 2)
            geo_north_lon = -72
            mlt1 = (24 + (utc_hours + (glon1 + geo_north_lon)/15))
            mlt2 = (24 + (utc_hours + (glon2 + geo_north_lon)/15))
            mltrange = abs(mlt1-mlt2)

            if mlt_low_lim < mltrange:
                node_counter_dict[edge_time][edge[0]] += 1
                node_counter_dict[edge_time][edge[1]] += 1

        pivot_table = pd.DataFrame([])
        if pc == 'pc3':
            times = self.complete_t_arr(start_dt, end_dt,20, pc=1)
        else:
            times = self.complete_t_arr(start_dt, end_dt,20, pc=0)
        for t in times:
            t = datetime.strftime(t,'%Y-%m-%d %H:%M:%S')
            new_row = pd.DataFrame(index=[t], columns=[np.arange(1,self.n)])
            pivot_table = pivot_table.append(new_row)
            pivot_table.loc[t, 1] = 1

        for time, node_counter in node_counter_dict.items():
            new_row = pd.DataFrame(index=[time], columns=[np.arange(1,self.n)])

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

    def dict_to_pd_df(self, pairwise_dict:dict, key:str)->pd.DataFrame:
        df = pd.DataFrame(pairwise_dict.items(), columns=['times', key])
        df['times'] =  pd.to_datetime(df['times'], format='%Y-%m-%d %H:%M:%S')
        df = df.set_index('times')
        df = df.sort_index()
        return df

    def ts_con_length_ratio_by_length(self) -> pd.DataFrame:
        geo_cords = self.geo_coords_dict()
        edge_counter_dict = dict()
        for edge in tqdm(self.edge_list):
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = {'long':0, 'short':0}
           # needs to be writen out because calling function many times slows down code
            n1_name = edge[0].split('_')[0]
            n2_name = edge[1].split('_')[0]
            glon1 = geo_cords[n1_name][0]
            glon2 = geo_cords[n2_name][0]
            glat1 = geo_cords[n1_name][1]
            glat2 = geo_cords[n2_name][1]
            lon1, lat1, lon2, lat2 = map(radians, [glon1, glat1, glon2, glat2])
            # haversine formula 
            dlon = lon2 - lon1 
            dlat = lat2 - lat1 
            a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
            c = 2 * asin(sqrt(a)) 
            # Radius of earth in kilometers is 6371
            km = 6371* c
            # rounded up length of the US canada border including alaska
            if km>9000:
                edge_counter_dict[edge_time]['long'] += 1
            else:
                edge_counter_dict[edge_time]['short'] += 1

        ratio_dict = {}
        long_conn = {}
        short_conn ={}
        for time, v in edge_counter_dict.items():
            long_conn[time] = v['long']
            short_conn[time] = v['short']
            if v['short'] ==0:
                ratio_dict[time] = 0            
            else:
                ratio_dict[time] = v['long']/v['short']

        ratio_dict = self.dict_to_pd_df(ratio_dict,'total')                       
        long_conn = self.dict_to_pd_df(long_conn,'total')
        short_conn = self.dict_to_pd_df(short_conn,'total')
        return ratio_dict, long_conn, short_conn
    
    def ts_con_length_ratio_by_mlt(self, mlt_lower_lim:float) -> pd.DataFrame:
        geo_cords = self.geo_coords_dict()
        edge_counter_dict = dict()
        
        for edge in tqdm(self.edge_list):
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = {'long':0, 'short':0}
           # needs to be writen out because calling function many times slows down code
            n1_name = edge[0].split('_')[0]
            n2_name = edge[1].split('_')[0]
            glon1 = geo_cords[n1_name][0]
            glon2 = geo_cords[n2_name][0]
            dt = datetime.strptime(edge_time,'%Y-%m-%d %H:%M:%S')
            utc_hours = np.round(dt.hour + dt.minute/60, 2)
            geo_north_lon = -72
            mlt1 = (24 + (utc_hours + (glon1 + geo_north_lon)/15))
            mlt2 = (24 + (utc_hours + (glon2 + geo_north_lon)/15))
            mltrange = abs(mlt1-mlt2)
            if mlt_lower_lim < mltrange:
                edge_counter_dict[edge_time]['long'] += 1
            else:
                edge_counter_dict[edge_time]['short'] += 1

        ratio_dict = {}
        long_conn = {}
        short_conn ={}
        for time, v in edge_counter_dict.items():
            long_conn[time] = v['long']
            short_conn[time] = v['short']
            if v['short'] ==0:
                ratio_dict[time] = 0            
            else:
                ratio_dict[time] = v['long']/v['short']

        ratio_dict = self.dict_to_pd_df(ratio_dict,'total')                       
        long_conn = self.dict_to_pd_df(long_conn,'total')
        short_conn = self.dict_to_pd_df(short_conn,'total')
        return ratio_dict, long_conn, short_conn

    def ts_all_conn(self) -> pd.DataFrame:

        edge_counter_dict = dict()
        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = 0
            edge_counter_dict[edge_time] += 1
       
        return self.dict_to_pd_df(edge_counter_dict,'total')

    def ts_conj_div_north_south_conn(self) -> pd.DataFrame:
        m_coords_dict = self.mag_coords_dict()
        edge_counter_dict = dict()
        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = {'conj':0, 'n-s':0}
            n1_label = edge[0].split('_')[0]
            n2_label = edge[1].split('_')[0]
            mlon1 = m_coords_dict[n1_label][0]
            mlon2 = m_coords_dict[n2_label][0]
            mlat1 = m_coords_dict[n1_label][1]
            mlat2 = m_coords_dict[n2_label][1]
            
            if abs(mlon1 - mlon2) <= 2:
                edge_counter_dict[edge_time]['conj'] += 1

            if mlat1>0 and mlat2<0:
                edge_counter_dict[edge_time]['n-s'] += 1
            elif mlat2>0 and mlat1<0:
                edge_counter_dict[edge_time]['n-s'] += 1
        
        conj_north_south = {}
        conj = {}
        north_south = {}
        for time, v in edge_counter_dict.items():
            conj[time] =  v['conj']
            north_south[time] = v['n-s']
            if v['n-s']==0:
                conj_north_south[time] = 0
            else:
                conj_north_south[time] = v['conj']/v['n-s']

        conj_north_south = self.dict_to_pd_df(conj_north_south,'total')                       
        conj = self.dict_to_pd_df(conj,'total')
        north_south = self.dict_to_pd_df(north_south,'total')
        return conj_north_south, conj, north_south


    def ts_north_south_div_north_conn(self) -> pd.DataFrame:
        m_coords_dict = self.mag_coords_dict()
        edge_counter_dict = dict()
        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']

            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = {'n':0, 'n-s':0}
            n1_label = edge[0].split('_')[0]
            n2_label = edge[1].split('_')[0]
            mlat1 = m_coords_dict[n1_label][1]
            mlat2 = m_coords_dict[n2_label][1]

            if mlat1>0 and mlat2>0:
                edge_counter_dict[edge_time]['n'] += 1
            elif mlat1>0 and mlat2<0:
                edge_counter_dict[edge_time]['n-s'] += 1
            elif mlat2>0 and mlat1<0:
                edge_counter_dict[edge_time]['n-s'] += 1

        north_south_north = {}
        north_south = {}
        north = {}
        for time, v in edge_counter_dict.items():
            north_south[time] = v['n-s']
            north[time] = v['n']
            if v['n']==0:
                north_south_north[time] = 0
            else:
                north_south_north[time] = v['n-s']/v['n']
        
        north_south_north = self.dict_to_pd_df(north_south_north,'total')                       
        north_south = self.dict_to_pd_df(north_south,'total')
        north = self.dict_to_pd_df(north,'total')
        return north_south_north, north_south, north

    def ts_ew_div_we_conn(self) -> pd.DataFrame:
        m_coords_dict = self.mag_coords_dict()
        edge_counter_dict = dict()
        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']

            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = {'ew':0, 'we':0}
            n1_label = edge[0].split('_')[0]
            n2_label = edge[1].split('_')[0]
            mlon1 = m_coords_dict[n1_label][0]
            mlon2 = m_coords_dict[n2_label][0]

            if mlon1>mlon2:
                    edge_counter_dict[edge_time]['we'] += 1
            elif mlon2> mlon1:
                edge_counter_dict[edge_time]['ew'] += 1
        
        ratio_dict = {}
        we_dict = {}
        ew_dict = {}
        for time, v in edge_counter_dict.items():
            ew_dict[time] = v['ew']
            we_dict[time] = v['we']
            if v['we']==0:
                ratio_dict[time]= 0
            else:
                ratio_dict[time] = v['ew']/v['we']

        ratio = self.dict_to_pd_df(ratio_dict,'total')                       
        ew = self.dict_to_pd_df(ew_dict,'total')
        we = self.dict_to_pd_df(we_dict,'total')
        return ratio, ew, we
        
    def ts_ns_div_sn_conn(self) -> pd.DataFrame:
        m_coords_dict = self.mag_coords_dict()
        edge_counter_dict = dict()
        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']

            if edge_time not in edge_counter_dict:
                edge_counter_dict[edge_time] = {'ns':0, 'sn':0}
            n1_label = edge[0].split('_')[0]
            n2_label = edge[1].split('_')[0]
            mlat1 = m_coords_dict[n1_label][1]
            mlat2 = m_coords_dict[n2_label][1]

            if mlat1>mlat2:
                    edge_counter_dict[edge_time]['ns'] += 1
            elif mlat2> mlat1:
                    edge_counter_dict[edge_time]['sn'] += 1

        ratio_dict = {}
        ns_dict = {}
        sn_dict = {}
        for time, v in edge_counter_dict.items():
            ns_dict[time] = v['ns']
            sn_dict[time] = v['sn']
            if v['sn']==0:
                ratio_dict[time]= 0
            else:
                ratio_dict[time] = v['ns']/v['sn']

        ratio = self.dict_to_pd_df(ratio_dict,'total')                       
        ns = self.dict_to_pd_df(ns_dict,'total')
        sn = self.dict_to_pd_df(sn_dict,'total')
        return ratio, ns, sn


    def t_snapshot_dist(self, times:list) -> dict:

        deg_dists = {}
        for t in times:
            deg_dists[t] = Counter()

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time in deg_dists.keys():
                deg_dists[edge_time][edge[0]] += 1
                deg_dists[edge_time][edge[1]] += 1

        for t in deg_dists.keys():
            deg_dists[t] = Counter(list(deg_dists[t].values()))

        return deg_dists

    def t_snapshot_dist_mlt(self, times:list, mlt_lower_lim) -> dict:
        geo_cords = self.geo_coords_dict()
        m_coords_dict = self.mag_coords_dict()

        deg_dists = {}
        deg_dists2 = {}
        for t in times:
            deg_dists[t] = Counter()
            deg_dists2[t] = Counter()

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time in deg_dists.keys():
                n1_name = edge[0].split('_')[0]
                n2_name = edge[1].split('_')[0]
                glon1 = geo_cords[n1_name][0]
                glon2 = geo_cords[n2_name][0]
                mlat1 = m_coords_dict[n1_name][1]
                mlat2 = m_coords_dict[n2_name][1]
                dt = datetime.strptime(edge_time,'%Y-%m-%d %H:%M:%S')
                utc_hours = np.round(dt.hour + dt.minute/60, 2)
                geo_north_lon = -72
                mlt1 = (24 + (utc_hours + (glon1 + geo_north_lon)/15))
                mlt2 = (24 + (utc_hours + (glon2 + geo_north_lon)/15))
                mltrange = abs(mlt1-mlt2)
                if mlt_lower_lim < mltrange:
                    deg_dists[edge_time][edge[0]] += 1
                    deg_dists[edge_time][edge[1]] += 1
                    if mlat1>0 and mlat2<0:
                        deg_dists2[edge_time][edge[0]] += 1
                        deg_dists2[edge_time][edge[1]] += 1
                    elif mlat1<0 and mlat2>0:
                        deg_dists2[edge_time][edge[0]] += 1
                        deg_dists2[edge_time][edge[1]] += 1        
        for t in deg_dists.keys():
            deg_dists[t] = Counter(list(deg_dists[t].values()))
        for t in deg_dists2.keys():
            deg_dists2[t] = Counter(list(deg_dists2[t].values()))
        # short mlt and ns_short_mlt
        return deg_dists, deg_dists2

    def t_snapshot_dist_chain(self, times:list) -> dict:
        m_coords_dict = self.mag_coords_dict()
        deg_dists = {}
        for t in times:
            deg_dists[t] = Counter()

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC']
            if edge_time in deg_dists.keys():
                n1_name = edge[0].split('_')[0]
                n2_name = edge[1].split('_')[0]
                mlon1 = m_coords_dict[n1_name][0]
                mlon2 = m_coords_dict[n2_name][0]
                if abs(mlon1 - mlon2) <= 2:
                    deg_dists[edge_time][n1_name] += 1
                    deg_dists[edge_time][n2_name] += 1

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


    def pc_power(self, ulf_power_filename:str, times:list, band:int):
        '''mean ulf powercs for pc single band for corespondig times'''
        pc_powers = pd.read_csv(ulf_power_filename)
        pc_powers['Date_UTC'] = pd.to_datetime(pc_powers['Date_UTC'])
        dict_pcp = {}
        for t in times:
            t = self.t_close(t)
            ts = pc_powers[pc_powers['Date_UTC'] == t][['Date_UTC', 'IAGA','PC2_IPOW','PC3_IPOW','PC4_IPOW','PC5_IPOW']]
            powers = ['PC2_IPOW', 'PC3_IPOW', 'PC4_IPOW', 'PC5_IPOW']
            dict_pcp[t] = np.array(ts[powers[band]].values)
        return dict_pcp


