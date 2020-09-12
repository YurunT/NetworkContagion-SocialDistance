#!/bin/bash
source activate ytian
# for p_mask in 0.2 0.45 0.6 0.9
for p_mask in 0.6
do
	echo $p_mask

    time python Mask-2-Save_parameters.py  -n 5000 -e 100 -th 0.01 -m $p_mask -T 1 -tm1 1 -tm2 0 -md 10 -ns 50 -nc 40 -change 0
    time python Mask-2-Save_parameters.py  -n 5000 -e 100 -th 0.01 -m $p_mask -T 1 -tm1 0 -tm2 1 -md 10 -ns 50 -nc 40 -change 0
#     time python Mask-2-Save_parameters-randomstart.py  -n 5000 -e 100 -th 0.01 -m $p_mask -T 0.6 -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 0
done
