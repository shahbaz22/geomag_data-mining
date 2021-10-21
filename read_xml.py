import pandas as pd
import pandas_read_xml as pdx
import struct
import scipy
from scipy import io
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
from datetime import datetime
import time 
import numpy as np

def compare_long_dist(new_data_folder_path:str, original_data_file:str, times_original_data:int, year:str)->None:

	station_names = [f.split('_')[1] for f in listdir(new_data_folder_path) if isfile(join(new_data_folder_path, f))]
	all_stations_data = pd.read_csv('supermag-stations.csv')
	longitudes = all_stations_data[all_stations_data['IAGA'].isin(station_names)][['IAGA','GEOLON']]

	orig_stations_data = pd.read_csv(original_data_file)
	# data set with only station names to be used as a label
	s_labels = orig_stations_data.drop_duplicates(subset=['IAGA'])[['IAGA','GEOLON']]
	# data set with only station names to be used as a label
	# checking all stations for incomplete data 4*3600 in the lenght of the data in seconds, keeping only complete data sets
	plt.hist(longitudes['GEOLON'], alpha=0.5, label='new_data', bins=15)
	plt.hist(s_labels['GEOLON'],alpha=0.5, label='original_data',bins=15)
	plt.legend()
	plt.title(f'comparison of longitudation station coverage {year} event')
	plt.xlabel('latt')
	plt.ylabel('freq')
	plt.savefig(f'plots/long_coverage_comp_{year}_data.png')
	plt.show()

# compare_long_dist('networks_data/2012','networks_data/spd_ulfpower.csv', 4*3600, 2012)
# compare_long_dist('networks_data/2015','networks_data/17march2013_4am.csv', 4*3600,2015)
# filename='networks_data/2012/2012_A08_1s_final.xdr'
# # df = pdx.read_xml(filename)
# fn = scipy.io.readsav(filename, idict=None, python_dict=True, uncompressed_file_name=None, verbose=False)
# print(fn)
# # solutions for dealing with weird time formatting stime given as time in seconds from the begining of the year
# # could add number of seconds from epoch time to current year then add seconds from year
# # time.ctime
# year = '2012'
# unix_year = time.mktime(datetime.strptime(year, "%Y").timetuple())
# print(unix_year)
# unix_start_time = unix_year + fn['stime'] 
# unix_end_time = unix_year + fn['etime'] 

# print(datetime.utcfromtimestamp(unix_start_time).strftime('%Y-%m-%d %H:%M:%S'))
# print(datetime.utcfromtimestamp(unix_end_time).strftime('%Y-%m-%d %H:%M:%S')
stations_data = pd.read_csv('supermag-stations.csv')
path = 'networks_data/2012'
station_names = [f.split('_')[1] for f in listdir(path) if isfile(join(path, f))]
counter_dict = {}
for i in ['e1','n1','z1']:
	counter_dict[i] = 0
for station in station_names:
	filename = f'{path}/2012_{station}_1s_final.xdr'
	fn = scipy.io.readsav(filename, idict=None, python_dict=True, uncompressed_file_name=None, verbose=False)
	for i in ['e1','n1','z1']:
		num_gaps = np.count_nonzero(fn[i] == 999999.0)
		pcent_amount = 100* num_gaps/len(fn[i])
		print(pcent_amount, i, station,'long',float(stations_data[stations_data['IAGA']==station]['GEOLON']))
		if pcent_amount <20:
			counter_dict[i] += 1
	# if pcent_amount>30:
	# 	filtarr = np.where(fn[i]==max(fn[i]),np.nan,fn[i])
	# 	# print(max(fn[i]))
	# 	for i in ['e1','n1','z1']:
	# 		plt.plot(np.arange(len(filtarr)),filtarr, label=i)
	# 	plt.legend()
	# 	plt.title(f'{station}')
	# 	plt.show()
	# 	plt.cla()
for i in ['e1','n1','z1']:
	print('num of complete datasts',counter_dict[i],i)
 





