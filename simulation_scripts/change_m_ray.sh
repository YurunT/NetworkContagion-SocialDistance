#!/bin/bash
source activate ytian
# for p_mask in  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
for p_mask in 0.45
do
echo change m: $p_mask
time python sim-ray.py -modelname mask  -n 5000 -e 1000 -th 0.05 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 3 -ns 50 -nc 40 -cp 100 -change 0 -msg 0to3 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0
done

for p_mask in 0.45
do
echo change m: $p_mask
time python sim-ray.py -modelname mask  -n 5000 -e 1000 -th 0.05 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 5 -ns 50 -nc 40 -cp 100 -change 0 -msg 0to5 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0
done

for p_mask in 0.45
do
echo change m: $p_mask
time python sim-ray.py -modelname mask  -n 5000 -e 1000 -th 0.05 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 7 -ns 50 -nc 40 -cp 100 -change 0 -msg 0to7 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0
done

for p_mask in 0.45
do
echo change m: $p_mask
time python sim-ray.py -modelname mask  -n 5000 -e 1000 -th 0.05 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 9 -ns 50 -nc 40 -cp 100 -change 0 -msg 0to9 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0
done