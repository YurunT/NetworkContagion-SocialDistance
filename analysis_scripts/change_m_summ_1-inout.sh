#!/bin/bash
source activate ytian


# for m0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
# do

# m1=$(echo "1-$m0"|bc)

# echo ----------------EXP STARTS----------------
# echo m0: $m0
# echo m1: $m1


# # time python analysis.py -modelname mask -itemname es -m $m0 $m1 -T 0.6 -tm1 0.6 0.9 -tm2 0.9 0.6 -mind 10 -maxd 10 -ns 1 -nc 40 -change 0 -msg es_check_inout3_impact_md10_m2_0_$m0 #-mdl  8 9 10

# time python analysis.py -modelname mask -itemname pe -m $m0 $m1 -T 0.6 -tm1 0.3 0.4 -tm2 0.4 0.3 -mind 10 -maxd 10 -ns 1 -nc 40 -change 0 -msg pe_check_inout2_impact_md10_m2_0_$m0 #-mdl  8 9 10
# echo ----------------EXP ENDS----------------
# done



for m0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7
do
for m2 in 0.2
do
m1=$(echo "1-$m0-$m2"|bc)

echo ----------------EXP STARTS----------------
echo m0: $m0
echo m1: $m1
echo m2: $m2


time python analysis.py -modelname mask -itemname es -m $m0 $m1 $m2 -T 0.6 -tm1 0.3 0.7 1 -tm2 0.7 0.3 1 -mind 10 -maxd 10 -ns 1 -nc 40 -change 0 -msg es_check_inout_impact_md10_m2_0.2_$m0 #-mdl  8 9 10

time python analysis.py -modelname mask -itemname pe -m $m0 $m1 $m2 -T 0.6 -tm1 0.3 0.7 1 -tm2 0.7 0.3 1 -mind 10 -maxd 10 -ns 1 -nc 40 -change 0 -msg pe_check_inout_impact_md10_m2_0.2_$m0 #-mdl  8 9 10
echo ----------------EXP ENDS----------------
done
done

for m0 in 0.1 0.2 0.3 0.4 0.5 
do
for m2 in 0.4
do
m1=$(echo "1-$m0-$m2"|bc)

echo ----------------EXP STARTS----------------
echo m0: $m0
echo m1: $m1
echo m2: $m2


time python analysis.py -modelname mask -itemname es -m $m0 $m1 $m2 -T 0.6 -tm1 0.3 0.7 1 -tm2 0.7 0.3 1 -mind 10 -maxd 10 -ns 1 -nc 40 -change 0 -msg es_check_inout_impact_md10_m2_0.4_$m0 #-mdl  8 9 10

time python analysis.py -modelname mask -itemname pe -m $m0 $m1 $m2 -T 0.6 -tm1 0.3 0.7 1 -tm2 0.7 0.3 1 -mind 10 -maxd 10 -ns 1 -nc 40 -change 0 -msg pe_check_inout_impact_md10_m2_0.4_$m0 #-mdl  8 9 10
echo ----------------EXP ENDS----------------
done
done