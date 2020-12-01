#!/bin/bash
source activate ytian

echo es_strategy0
time python sim-ray.py -modelname mask -itemname es -n 500000 -e 1000 -th 0.05 -m 1 -T 0.6 -tm1 1 -tm2 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg strategy0_single #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo es_strategy1
time python sim-ray.py -modelname mask -itemname es -n 500000 -e 1000 -th 0.05 -m 1 -T 0.6 -tm1 0.3 -tm2 0.2 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg strategy1_single #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo es_strategy2
time python sim-ray.py -modelname mask -itemname es -n 500000 -e 1000 -th 0.05 -m 1 -T 0.6 -tm1 0.5 -tm2 0.5 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg strategy2_single #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo es_strategy3
time python sim-ray.py -modelname mask -itemname es -n 500000 -e 1000 -th 0.05 -m 0.5 0.5 -T 0.6 -tm1 0.3 0.5 -tm2 0.2 0.5 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg strategy3_single #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0