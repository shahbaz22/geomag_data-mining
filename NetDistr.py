import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sb
import os


class NetDistr:
    def __init__(self, edge_list_dir: str, n: int):
        G = nx.read_edgelist(edge_list_dir)
        self.edge_list = list(G.edges(data=True))
        self.n = n

    def create_pivottable(self) -> pd.DataFrame:

        node_counter_dict = dict()
        for edge in self.edge_list:

            edge_time = edge[2]['attr_dict']['UTC1']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = Counter()
            # specifying/setting dict key and asinging count to counter 
            node_counter_dict[edge_time][edge[0]] += 1
            node_counter_dict[edge_time][edge[1]] += 1
            
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
        for edge in self.edge_list:

            edge_time = edge[2]['attr_dict']['UTC1']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = Counter()
 
            mltrange = abs(edge[2]['attr_dict']['MLT1'] - edge[2]['attr_dict']['MLT2'])
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
        date, *_ = self.edge_list[0][2]['attr_dict']['UTC1'].split(' ')
        date_time = f'{date} {time}'

        for edge in self.edge_list:
            edge_time = edge[2]['attr_dict']['UTC1']
            if edge_time==date_time:
                node_counter_dict_mlt1[edge[0]] += 1
                node_counter_dict_mlt1[edge[1]] += 1
                mltrange = abs(edge[2]['attr_dict']['MLT1'] - edge[2]['attr_dict']['MLT2'])
                if mlt_lower_lim < mltrange:
                    node_counter_dict_mlt2[edge[0]] += 1
                    node_counter_dict_mlt2[edge[1]] += 1

        degree_counts_mlt1 = Counter(list(node_counter_dict_mlt1.values()))
        degree_counts_mlt2 = Counter(list(node_counter_dict_mlt2.values()))

        return degree_counts_mlt1, degree_counts_mlt2



                


