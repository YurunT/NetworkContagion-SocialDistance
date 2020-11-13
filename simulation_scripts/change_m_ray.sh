#!/bin/bash
source activate ytian
for p_mask in 0.45
do
echo change m: $p_mask
time python sim-ray.py -modelname mask  -n 5000000 -e 1000 -th 0.05 -m $p_mask 0.3 0.25 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.7 0.2 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg 3tps_n5M_e1k_5to10 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0
done
