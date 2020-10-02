
#!/bin/bash
source activate ytian
for T in  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
do
	echo $T

    time python Mask-2-Tmask12_newpath-Ray-limitCore-100save-iidmask.py -n 50000 -e 1000 -th 0.3 -m 0.45 -T $T -tm1 0.3 -tm2 0.7 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 1
done
