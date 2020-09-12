#!/bin/bash
source activate ytian
for p_mask in 0.2 0.45 0.6 0.9
#  for p_mask in 0.45
do
	echo $p_mask
        #python Mask-2-Tmask12_newpath-Ray.py -n 500000 -e 10000 -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 50
	time python Mask-2-Tmask12_newpath-Ray-limitCore-100save.py -n 5000 -e 100 -th 0.01 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.4 -md 10 -ns 50 -nc 40 -cp 100 -change 0
done
