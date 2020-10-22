#!/bin/bash
source activate ytian
# for p_mask in  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
for p_mask in 0.45
do
echo change m: $p_mask
time python sim-ray.py -modelname mask  -n 5000 -e 10000 -th 0.05 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 11 -ns 50 -nc 40 -cp 1000 -change 0 -msg 0to11 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0
done
