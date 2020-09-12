#!/bin/bash
source activate ytian
for T_mask1 in  0.3   0.5  0.8
do
for T_mask2 in  0.4  0.6  0.9 
do
	echo $T_mask1
	echo $T_mask2
#	python Mask-2-Tmask12_newpath-Ray.py -n 500000 -e 10000 -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50
    time python anaEvProb-Parellel.py  -n 2000000  -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2
    time python MaskModel_PE_top_bottom-startfrom.py  -n 2000000  -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2
done

done
