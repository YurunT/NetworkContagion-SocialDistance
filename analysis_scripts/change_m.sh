#!/bin/bash
source activate ytian
for p_mask in 0.45
do
	echo $p_mask
    time python analysis.py -modelname mask -itemname es -m $p_mask 0.3 0.25 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.7 0.2 1 -mind 0 -maxd 10 -ns 100 -nc 40 -change 0 -msg 0to10_3_types #-mdl  8 9 10
done
