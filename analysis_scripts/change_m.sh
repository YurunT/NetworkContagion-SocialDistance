#!/bin/bash
source activate ytian

for p_mask in  0.6 
do
	echo $p_mask
    time python analysis.py -modelname mask -itemname es -m $p_mask -T 0.7 -tm1 0.4  -tm2 0.6 -mind 3 -maxd 3 -ns 1 -nc 40 -change 0 -msg testmutation -mdl 0.2 0.3 0.4 
done
