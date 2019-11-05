#convert fastq files to csv files
#nmake your next program lead this csv file into a
# pandas dataframe with 3 columns: cluster_identity, sequence, Q_score

import pandas as pd
import subprocess
import numpy as np

f = open("lib2_09212019_low_clusters/practice_fastq", "r")
fastq = f.read().splitlines()
f.close()

seq = pd.DataFrame(fastq, columns = ["row"]) # columns = ['cluster_identity', 'sequence', 'Q_score'])
#subprocess.call()
seq2 = []
pos = 0

#for index, row in seq.iterrows():
 #   if seq.row != '+\n':
  #      seq2.append(pos)
   #     pos = pos + 1
    #else:
     #   pos = pos +1
column_names = ['cluster_identity', 'sequence', 'Q_score']

seq3 = pd.DataFrame(columns = ['cluster_identity'])

col1 = 0
col2 = 1
col3 = 3
for index, row in seq.iterrows():
    seq3.append({pos:[seq[col1],seq[col2],seq[col3]]})
    col1 = col1 + 3
    col2 = col2 + 3
    col3 = col3 + 3
#for index, row in seq.iterrows():
 #   if seq.contains('+\n'):
  #      pass
   # else:
    #    list.append(row)

#seq2 = pd.DataFrame(row)


print("This is seq 3", seq3.head())



#df = pd.read_csv('lib2_09212019_low_clusters/tecto_S1_L001_R1_001.fastq', sep="/n", engine = 'python')


#f = open("lib2_09212019_low_clusters/tecto_S1_L001_R1_001.fastq", "r")
#f.readlines()

#f.close()

