#!/bin/bash
source activate ytian
# for p_mask in 0.2 0.45 0.6 0.9
# for p_mask in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9

for p_mask in 0.45 
do
	echo $p_mask
     
       time python Mask_PE_Analysis.py     -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 2 -maxd 3 -ns 10 -nc 60 -change 0
       
#        time python Mutation_PE_Analysis.py -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 5 -maxd 5 -ns 1 -nc 60 -change 0
done

for p_mask in 0.6 
do
	echo $p_mask
     
       time python Mask_PE_Analysis.py     -n 2000000  -th 0.01 -m $p_mask -T 0.7 -tm1 0.4 -tm2 0.6 -mind 2 -maxd 3 -ns 10 -nc 60 -change 0
       
#        time python Mutation_PE_Analysis.py -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 5 -maxd 5 -ns 1 -nc 60 -change 0
done