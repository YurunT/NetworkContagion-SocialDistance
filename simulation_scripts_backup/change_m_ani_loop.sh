#!/bin/bash
source activate ytian
for p_mask in 0.2 0.3 0.4 0.5 0.6 0.8
do
echo $p_mask

    time python mask_model_sim.py  -n 5000000 -e 10000 -th 0.05 -m $p_mask -T 0.8 -tm1 0.4 -tm2 0.5  -nc 40 -change 0
done
