
import numpy as np
import pandas as pd


seq_df = pd.read_csv("lib2_09212019_low_clusters/R1_low_den_09272019_tecto.csv")

seq_array = np.array(seq_df)

#RNAP_5end = "TATGCTATAATTATTTC"
#RNAP_3end = "ATACGATATTAATAAAG"

RNAP_5end_short = "GCTATAATTAT"
RNAP_5end_shorter = "GCTATAAT"
RNAP_5end = "TATGCTATAATTATT"
RNAP_3end = "CGATATTAATA"
rvs_RNAP_5end = "ATAATTATAGC"
rvs_RNAP_3end = "TATTAATATCG"
#print(seq_array)
#print(seq_array.shape)

#tecto_array = np.empty(0,3)
pos = 0
tecto_numbers = 0
for row in seq_array:
    if RNAP_5end_shorter in seq_array[pos][2]:
        
        #line_array = np.array(seq_array[pos])
        #tecto_array = np.vstack((tecto_array, line_array))
        pos += 1
        tecto_numbers += 1
    else:
        pos += 1


print(tecto_numbers)

#tecto_df = seq_df.loc
#RNAP 5 PRIME IS IN THE R1 TECTO SEQUENCES!!!!!!
#IN 424674 READS (If there's ~55,000 unique tectos, if they were evenly distributed there would be N=~7.7 for each!)