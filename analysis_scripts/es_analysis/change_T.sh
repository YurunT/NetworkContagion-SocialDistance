#!/bin/bash
source activate ytian
for T in   0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
# for T in   0.4  0.7  0.9 0.6
# for T in   0.5
do
	echo $T
    time python Mask_ES_Analysis.py     -m 0.45 -T $T -tm1 0.3  -tm2 0.7 -mind 5 -maxd 5 -ns 1 -nc 40 -change 1
#     time python Mutation_ES_Analysis.py -n 2000000  -th 0.01 -m 0.6 -T $T -tm1 0.3 -tm2 0.7 -md 10 -ns 50 -nc 40 -change 1
    
#     time python Mask_ES_Analysis.py     -n 2000000  -th 0.01 -m 0.6 -T $T -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 1
#     time python Mutation_ES_Analysis.py -n 2000000  -th 0.01 -m 0.6 -T $T -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 1
done