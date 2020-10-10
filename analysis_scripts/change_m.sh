#!/bin/bash
source activate ytian

# -modelname: {'mask', 'mutation'}
# -itemname:  {'es', 'pe'}
# -mdl: if need to specify a certian mean_degree_list, or please comment this option


for p_mask in  0.45
do
	echo $p_mask
    time python analysis.py -modelname mutation -itemname es -m $p_mask -T 0.6 -tm1 0.3  -tm2 0.7 -mind 0 -maxd 10 -ns 50 -nc 40 -change 0 -msg test -mdl 0.2 0.3 0.4 
done
