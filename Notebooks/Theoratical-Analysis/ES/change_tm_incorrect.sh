#!/bin/bash
source activate ytian
for T_mask1 in  0.3  0.5  0.8
# for T_mask1 in   0.5  
do
for T_mask2 in  0.4  0.6  0.9 
# for T_mask2 in  0.8
do
	echo $T_mask1
	echo $T_mask2

# 	time python Mask-2-Tmask12_newpath-Ray-limitCore-100save.py -n 500000 -e 10000 -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2  -md 10 -ns 50 -nc 40 -cp 500 -change 2
#     time python analysisEvolutionFinal-Parellel.py  -n 2000000  -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2

#     time python MaskModelDerivation-vR-startfrom0-Tmask12-with_prabh.py                -n 2000000  -th 0.01 -m 0.6 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2
    
#     time python MaskModelDerivation-vR-startfrom0-play.py  -n 2000000  -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2

#        time python MaskModelDerivation-vR-startfrom0-Tmask12-recover-Ray.py  -n 2000000  -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2
#        time python analysisEvolutionFinal-Parellel-correct-Rashad.py                    -n 2000000  -th 0.01 -m 0.45 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2

'''
Sep 10th incorrect version
'''

       time python MaskModelDerivation-vR-startfrom0-play.py            -n 2000000  -th 0.01 -m 0.6 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2
       
       time python analysisEvolutionFinal-Parellel-INcorrect-Rashad.py  -n 2000000  -th 0.01 -m 0.6 -T 0.6 -tm1 $T_mask1 -tm2 $T_mask2 -md 10 -ns 50 -nc 40 -change 2

'''
Sep 10th incorrect version
'''
    
done

done
