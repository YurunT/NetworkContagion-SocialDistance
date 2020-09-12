#!/bin/bash
source activate ytian
# for p_mask in 0.2 0.45 0.6 0.9

for p_mask in 0.6
do
	echo $p_mask
     
# 	time python Mask-2-Tmask12_newpath-Ray-limitCore-100save.py -n 500000 -e 10000 -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 50 -nc 40 -cp 500 -change 0
#     time python analysisEvolutionFinal-Parellel.py  -n 2000000          -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 10 -nc 40 -change 0
#     time python MaskModelDerivation-vR-startfrom0-Tmask12-recover.py  -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 50 -nc 40 -change 0
    
     
#      time python MaskModelDerivation-vR-startfrom0-play.py  -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 10 -nc 40 -change 0
#        time python analysisEvolutionFinal-Parellel-Prabh.py  -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 10 -nc 40 -change 0


#'''
#Sep 10th correct version
#'''


#        time python Mask-basedonRashad-Parellel-COrrect.py                      -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 0

       time python MaskModelDerivation-vR-startfrom0-play-correct.py            -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 1   -tm2 0 -md 10 -ns 50 -nc 40 -change 0
       
       time python analysisEvolutionFinal-Parellel-COrrect-Rashad.py            -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 1   -tm2 0 -md 10 -ns 50 -nc 40 -change 0
       
       time python MaskModelDerivation-vR-startfrom0-play-correct.py            -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0   -tm2 1 -md 10 -ns 50 -nc 40 -change 0
       
       time python analysisEvolutionFinal-Parellel-COrrect-Rashad.py            -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0   -tm2 1 -md 10 -ns 50 -nc 40 -change 0



       
       
#'''
#Sep 10th correct version
#'''       
       
# Sep 11th Osman
#        time python MaskModelDerivation-vR-startfrom0-play-correct-Osman.py          -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 0
       
       

# Sep 11th Osman
       
       
#        time python MaskModelDerivation-vR-startfrom0-Tmask12-recover-Ray.py  -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -md 10 -ns 50 -nc 40 -change 0
#        time python analysisEvolutionFinal-Parellel-Ray.py  -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 0
       
       
#        time python analysisEvolutionFinal-Parellel-correct-Rashad.py         -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.5 -tm2 0.5 -md 10 -ns 50 -nc 40 -change 0
#        time python analysisEvolutionFinal-Parellel-correct-Rashad.py         -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.4 -tm2 0.6 -md 10 -ns 50 -nc 40 -change 0
#        time python analysisEvolutionFinal-Parellel-correct-Rashad.py         -n 2000000  -th 0.01 -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -md 10 -ns 50 -nc 40 -change 0

     
done
