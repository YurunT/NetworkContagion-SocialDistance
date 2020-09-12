#!/bin/bash
source activate ytian
for p_mask in 0.2 0.45 0.6 0.9
# for p_mask in 0.45
do
	echo $p_mask
        #python Mask-2-Tmask12_newpath-Ray.py -n 500000 -e 10000 -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 50
    time python anaEvProb-Parellel.py                 -n 2000000  -th 0.01 -m $p_mask -T 0.4 -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 0
    time python MaskModel_PE_top_bottom-startfrom.py  -n 2000000  -th 0.01 -m $p_mask -T 0.4 -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 0
done
