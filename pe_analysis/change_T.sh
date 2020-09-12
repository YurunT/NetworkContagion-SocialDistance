#!/bin/bash
source activate ytian
# for T in   0.4  0.7  0.9 0.6
for T in   0.8
do
	echo $T
    time python Mask_PE_Analysis.py     -n 2000000  -th 0.01 -m 0.45 -T $T -tm1 0.4 -tm2 0.6 -md 10 -ns 50 -nc 40 -change 1
    time python Mutation_PE_Analysis.py -n 2000000  -th 0.01 -m 0.45 -T $T -tm1 0.4 -tm2 0.6 -md 10 -ns 50 -nc 40 -change 1
done