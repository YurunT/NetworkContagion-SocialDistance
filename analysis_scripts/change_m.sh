#!/bin/bash
source activate ytian

# -modelname: {'mask', 'mutation'}
# -itemname:  {'es', 'pe'}
# -mdl: if need to specify a certian mean_degree_list, or please comment this option

# for p_mask in 0.45
# do
# 	echo $p_mask
#     time python analysis.py -modelname mask -itemname pe -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 10 -maxd 10 -ns 1 -nc 40 -change 0 -msg changemdl_8to10 -mdl  8 9 10
# done

for T in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
	echo $T
    time python analysis.py -modelname mask -itemname pe -m 0.45 -T $T -tm1 0.3 -tm2 0.7 -mind 5 -maxd 5 -ns 1 -nc 40 -change 1 -msg varyT # -mdl  8 9 10
done