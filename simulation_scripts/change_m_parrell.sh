#!/bin/bash
source activate ytian
for p_mask in 0.45
do
echo $p_mask

    time python Mask_simulation-wrappednoglobal-iidmask.py -n 5000000 -e 10000 -th 0.05 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 10  -ns 50 -nc 40 -change 0
done
