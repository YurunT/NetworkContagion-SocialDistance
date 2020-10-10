#!/bin/bash
source activate ytian
# for p_mask in  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
for p_mask in 0.45
do
echo change m: $p_mask
time python sim-ray.py -modelname mask  -n 50 -e 10 -th 0.05 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 5 -change 0 -msg test -mdl 0.2 0.3 0.4 
done
