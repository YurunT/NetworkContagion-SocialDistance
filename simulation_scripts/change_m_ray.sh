#!/bin/bash
source activate ytian

# echo strategy0
# time python sim-ray.py -modelname mask  -n 5000 -e 10000 -th 0.05 -m 0 0 1 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 1000 -change 0 -msg strategy0 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

# echo strategy1
# time python sim-ray.py -modelname mask  -n 5000 -e 10000 -th 0.05 -m 1 0 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 1000 -change 0 -msg strategy1 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo strategy2
time python sim-ray.py -modelname mask  -n 500000 -e 10000 -th 0.05 -m 0 1 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 1000 -change 0 -msg strategy2 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo strategy3
time python sim-ray.py -modelname mask  -n 500000 -e 10000 -th 0.05 -m 0.5 0.5 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 1000 -change 0 -msg strategy3 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0
