#!/bin/bash
source activate ytian
for T in   0.4  0.7  0.9 0.6
# for T in   0.5
do
	echo $T

# 	time python Mask-2-Tmask12_newpath-Ray-limitCore-100save.py -n 500000 -e 10000 -th 0.01 -m 0.45 -T $T -tm1 0.4 -tm2 0.6 -md 10 -ns 50 -nc 40 -cp 500 -change 1
#    time python analysisEvolutionFinal-Parellel.py -n 2000000  -th 0.01 -m 0.6 -T $T -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 1
#     time python MaskModelDerivation-vR-startfrom0-Tmask12-recover.py  -n 2000000  -th 0.01 -m 0.45 -T $T -tm1 0.4 -tm2 0.6 -md 10 -ns 50 -nc 40 -change 1

#     time python MaskModelDerivation-vR-startfrom0-Tmask12-recover.py -n 2000000  -th 0.01 -m 0.6 -T $T -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 1
#     time python MaskModelDerivation-vR-startfrom0-play.py -n 2000000  -th 0.01 -m 0.6 -T $T -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 1
#        time python MaskModelDerivation-vR-startfrom0-Tmask12-recover-Ray.py  -n 2000000  -th 0.01 -m 0.6 -T $T -tm1 0.4 -tm2 1 -md 10 -ns 50 -nc 40 -change 1
#        time python analysisEvolutionFinal-Parellel-Ray.py                    -n 2000000  -th 0.01 -m 0.6 -T $T -tm1 0.4 -tm2 1 -md 10 -ns 50 -nc 40 -change 1

#'''
#Sep 10th correct version
#'''

       time python MaskModelDerivation-vR-startfrom0-play-correct.py           -n 2000000  -th 0.01 -m 0.8 -T $T -tm1 0.3 -tm2 0.7 -md 10 -ns 50 -nc 40 -change 1
       
       time python analysisEvolutionFinal-Parellel-COrrect-Rashad.py -n 2000000  -th 0.01 -m 0.8 -T $T -tm1 0.3 -tm2 0.7 -md 10 -ns 50 -nc 40 -change 1
       
#'''
#Sep 10th correct version
#'''
done
