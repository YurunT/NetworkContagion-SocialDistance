#!/bin/bash
source activate ytian
for T in   0.4  0.6  0.8 
do
	echo $T
#	python Mask-2-Tmask12_newpath-Ray.py -n 500000 -e 10000 -th 0.01 -m 0.45 -T $T -tm1 0.4 -tm2 0.6 -md 10 -ns 50
    time python anaEvProb-Parellel.py  -n 2000000  -th 0.01 -m 0.45 -T $T -tm1 0.4 -tm2 0.9 -md 10 -ns 50 -nc 40 -change 1
    time python MaskModel_PE_top_bottom-startfrom.py  -n 2000000  -th 0.01 -m 0.45 -T $T -tm1 0.4 -tm2 0.9 -md 10 -ns 50 -nc 40 -change 1
done
