#!/bin/bash
source activate ytian
for m0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
m1=$(echo "1-$m0"|bc)

echo ----------------EXP STARTS----------------
echo m0: $m0
echo m1: $m1

time python analysis.py -modelname mask -itemname es -m $m0 $m1 -T 0.6 -tm1 0.2 0.5 -tm2 0.3 0.5 -mind 20 -maxd 20 -ns 1 -nc 40 -change 0 -msg es_check_m_impact_md20_$m0 #-mdl  8 9 10

time python analysis.py -modelname mask -itemname pe -m $m0 $m1 -T 0.6 -tm1 0.2 0.5 -tm2 0.3 0.5 -mind 20 -maxd 20 -ns 1 -nc 40 -change 0 -msg pe_check_m_impact_md20_$m0 #-mdl  8 9 10
echo ----------------EXP ENDS----------------
done



