#!/bin/bash
source activate ytian

for p_mask in  0.6 
do
	echo $p_mask
     
    time python Mask_ES_Analysis.py     -m $p_mask -T 0.7 -tm1 0.4  -tm2 0.6 -mind 3 -maxd 3 -ns 1 -nc 40 -change 0
done






# for p_mask in  0.45
# do
# 	echo $p_mask
     
#     time python Mask_ES_Analysis.py     -m $p_mask -T 0.6 -tm1 0.3  -tm2 0.7 -mind 2 -maxd 3 -ns 10 -nc 40 -change 0
# #     time python Mutation_ES_Analysis.py -m $p_mask -T 0.6 -tm1 0.3  -tm2 0.7 -mind 0 -maxd 10 -ns 50 -nc 40 -change 0
# done