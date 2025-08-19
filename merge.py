import pandas as pd
import os
import sys
from itertools import groupby

#change directory to where files are located
#os.chdir(sys.argv[1])

#tell interpreter where files are located
directory = os.fsencode(sys.argv[1])

#make list of files
file_list = [os.fsdecode(file) for file in os.listdir(directory) if os.fsdecode(file).endswith('.txt')]

#if there are other files in the directory other than the counts files,
#iterate through file list and keep only counts files
file_list_2 = []
for file in file_list:
	#if '27986' in file: #change quoted text to something specific to all counts files
	file_list_2.append(file)

#need to get ID column for usage later
id_col = pd.read_csv(file_list_2[0], sep = '\t', usecols = [0])

# sort list
# essential for grouping
file_list_2.sort()

for i in file_list:
	print('_'.join(i.split('_')[:-1]))

# using lambda + itertools.groupby() + split()
# group similar substrings
res = [list(i) for j, i in groupby(file_list_2, lambda a: '_'.join(a.split('_')[:-1]))]

#generate list of basenames for output files
out_files = ['_'.join(x[0].split('_')[:-1]) for x in res]

#define function to merge columns with same names together
def same_merge(x): return ','.join(x[x.notnull()].astype(str))

#len(res) = # of groups, which should be equal to len(out_files)
for group, file_name in zip(range(len(res)), range(len(out_files))):
	#read in all files using only count column
	file_map = map(lambda x: pd.read_csv(x, usecols = [1], sep = '\t'), res[group])
	#merge all dataframes by row
	df = pd.concat(file_map, axis = 1, ignore_index = True)
	#sum numeric columns (count) and output new df with ID column
	new_df = pd.DataFrame({'ID': id_col['ID'], 'count': df.sum(axis = 1)})
	#new_df.to_csv('/mnt/c/Users/jdlab/Desktop/tiling_mpra/MGI_051523_2/data/' + out_files[file_name] + '.txt', index = False, sep = '\t')
	#new_df.to_csv(out_files[file_name] + '.txt', index = False, sep = '\t')
	new_df.to_csv(sys.argv[2] + '/' + out_files[file_name] + '.txt', index = False, sep = '\t')
