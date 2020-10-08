#!/bin/bash
source activate ytian
# for p_mask in  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
for p_mask in 0.6 
# for p_mask in 0.45
do
echo $p_mask

time python ani-ray.py  -n 50 -e 10 -th 0.3 -m $p_mask -T 0.7 -tm1 0.4 -tm2 0.6 -mind 5 -maxd 5 -ns 1 -nc 40 -cp 2 -change 0

done
