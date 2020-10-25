#!/bin/bash
source activate ytian

# -modelname: {'mask', 'mutation'}
# -itemname:  {'es', 'pe'}
# -mdl: if need to specify a certian mean_degree_list, or please comment this option

for p_mask in 0.45
do
	echo $p_mask
    time python analysis.py -modelname mask -itemname pe -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 1 -maxd 1 -ns 1 -nc 40 -change 0 -msg testtnn #-mdl  8 9 10
done
