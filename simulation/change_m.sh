#!/bin/bash
source activate ytian
# for p_mask in 0.2 0.45 0.6 0.9
for p_mask in 0.6
do
	echo $p_mask
    time python Mask_simulation.py  -n 500 -e 10 -th 0.01 -m $p_mask -T 0.6 -tm1 1 -tm2 0 -md 10 -ns 50 -nc 40 -change 0
done