#!/bin/bash
source activate ytian

# -modelname: {'mask', 'mutation'}
# -itemname:  {'es', 'pe'}
# -mdl: if need to specify a certian mean_degree_list, or please comment this option


for p_mask in  0.45
do
	echo $p_mask
    time python analysis.py -modelname mask -itemname pe -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 3 -ns 50 -nc 40 -change 0 -msg 0to3 #-mdl 2.44897959  2.65306122  2.85714286  3.06122449  3.26530612  3.46938776
done

for p_mask in  0.45
do
	echo $p_mask
    time python analysis.py -modelname mask -itemname pe -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 5 -ns 50 -nc 40 -change 0 -msg 0to5 #-mdl 2.44897959  2.65306122  2.85714286  3.06122449  3.26530612  3.46938776
done

for p_mask in  0.45
do
	echo $p_mask
    time python analysis.py -modelname mask -itemname pe -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 7 -ns 50 -nc 40 -change 0 -msg 0to7 #-mdl 2.44897959  2.65306122  2.85714286  3.06122449  3.26530612  3.46938776
done

for p_mask in  0.45
do
	echo $p_mask
    time python analysis.py -modelname mask -itemname pe -m $p_mask -T 0.6 -tm1 0.3 -tm2 0.7 -mind 0 -maxd 9 -ns 50 -nc 40 -change 0 -msg 0to9 #-mdl 2.44897959  2.65306122  2.85714286  3.06122449  3.26530612  3.46938776
done
