#!/bin/bash
source activate ytian
# for T_mask1 in  0.3   0.5  0.8
for T_mask1 in  0.3  
do
# for T_mask2 in  0.4  0.6  0.9 
for T_mask2 in  0.4  
do
	echo $T_mask1
	echo $T_mask2
	time python Mask_simulation.py -n 5000 -e 100 -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2  -md 10 -ns 50 -nc 40 -cp 500 -change 2
done
done