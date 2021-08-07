import pandas as pd
import pickle
import networkx as nx
import ast
import numpy as np


class PivotTable:
    def __init__(self, edge_list_dir: str, n: int):
        G = nx.read_edgelist('test_e_170313_15_0.3.txt')
        self.edge_list = list(G.edges(data=True))
        self.n = n
        
    def create_PivotTable(self) -> None:
        curr_index = 0
        curr_time = edge_list[curr_index][2]['attr_dict']['UTC1']
        
        self.pivot_table = pd.DataFrame([])
        
        while curr_index < len(edge_list):
            node_counter = Counter()
            while curr_index < len(edge_list) and edge_list[curr_index][2]['attr_dict']['UTC1'] == curr_time:
                curr_edge = edge_list[curr_index]
                node_counter[curr_edge[0]] += 1
                node_counter[curr_edge[1]] += 1
                curr_index +=1
            
            if curr_index >= len(edge_list):
                break
                
                
            new_row = pd.DataFrame(index=[curr_time],columns=[np.arange(1,self.n)])
            self.pivot_table = self.pivot_table.append(new_row)
            
            degree_counts = Counter(list(node_counter.values()))
            for degree_key, degree_count in degree_counts.items():
                self.pivot_table.loc[curr_time, degree_key] = degree_count
                
            curr_time = edge_list[curr_index][2]['attr_dict']['UTC1']
            
    def plot_pivot_table(self) -> None:
        # TODO: PLot Heatmap
        pass



