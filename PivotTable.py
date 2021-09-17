import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter

class PivotTable:
    def __init__(self, edge_list_dir: str, n: int):
        G = nx.read_edgelist(edge_list_dir)
        self.edge_list = list(G.edges(data=True))
        self.n = n
        # initializing function runs create_PivotTable() function  creats self.pivot table object
        self.create_PivotTable()
        
    def create_PivotTable(self) -> None:

        node_counter_dict = dict()
        for edge in self.edge_list:

            edge_time = edge[2]['attr_dict']['UTC1']
            if edge_time not in node_counter_dict:
                node_counter_dict[edge_time] = Counter()
            node_counter_dict[edge_time][edge[0]] += 1
            node_counter_dict[edge_time][edge[1]] += 1
            
          
        self.pivot_table = pd.DataFrame([])
        for time, node_counter in node_counter_dict.items():
            # print(time,node_counter)
            new_row = pd.DataFrame(index=[time], columns=[np.arange(1,self.n)])
            self.pivot_table = self.pivot_table.append(new_row)

            degree_counts = Counter(list(node_counter.values()))
            # print(degree_counts)
            for degree_key, degree_count in degree_counts.items():
                # print(degree_key, degree_count)
                self.pivot_table.loc[time, degree_key] = degree_count

                
        self.pivot_table = self.pivot_table.sort_index()
