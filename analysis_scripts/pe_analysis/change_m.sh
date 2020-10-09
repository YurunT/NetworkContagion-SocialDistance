#!/bin/bash
source activate ytian

for p_mask in 0.45 
do
	echo $p_mask
     
       time python Mutation_PE_Analysis.py     -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 3 -maxd 3 -ns 1 -nc 60 -change 0
done