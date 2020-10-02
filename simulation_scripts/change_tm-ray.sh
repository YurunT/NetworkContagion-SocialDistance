#!/bin/bash
source activate ytian
# for T_mask1 in  0.3   0.5  0.8
for T_mask1 in  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
do
# for T_mask2 in  0.4  0.6  0.9 

	echo $T_mask1
	
	time python Mask-2-Tmask12_newpath-Ray-limitCore-100save-iidmask.py -n 5000000 -e 200 -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 0.7  -mind 5 -maxd 5 -ns 1 -nc 40 -cp 100 -change 2
done

for T_mask2 in  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
do
echo $T_mask2

time python Mask-2-Tmask12_newpath-Ray-limitCore-100save-iidmask.py -n 5000000 -e 200 -th 0.01 -m 0.45 -T 0.6 -tm1 0.3 -tm2 $T_mask2  -mind 5 -maxd 5 -ns 1 -nc 40 -cp 100 -change 2

done
