#!/bin/bash
source activate ytian
for p_mask in 0.1 0.3 0.5 0.7 0.9 
do
	echo $T
        python Mask-2-Tmask12_newpath.py -n 500000 -e 10000 -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 50
done
