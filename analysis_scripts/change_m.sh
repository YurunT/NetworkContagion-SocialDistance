#!/bin/bash
source activate ytian
for p_mask in 0.45
do
	echo $p_mask
    time python analysis.py -modelname mask -itemname es -m $p_mask 0.55 -T 0.6 -tm1 0.3 1 -tm2 0.7 1 -mind 0 -maxd 10 -ns 50 -nc 40 -change 0 -msg change_A #-mdl  8 9 10
done
