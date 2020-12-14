#!/bin/bash
source activate ytian


# for m0 in 0.7 0.8
# do
# for m2 in 0.1
# do
# m1=$(echo "1-$m0-$m2"|bc)

# echo m0: $m0
# echo m1: $m1
# echo m2: $m2

#     time python sim-ray.py -modelname mask -itemname es -n 50000 -e 100 -m $m0 $m1 $m2 -T 0.6 -tm1 0.2 0.5 1 -tm2 0.3 0.5 1 -mind 15 -maxd 15 -ns 1 -nc 40 -change 0 -cp 50 -msg es_check_m_impact_md15_m2_0.1_$m0 #-mdl  8 9 10

# done
# done



for m0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
m1=$(echo "1-$m0"|bc)

echo m0: $m0
echo m1: $m1


    time python sim-ray.py -modelname mask -itemname es -n 50000 -e 100 -m $m0 $m1 -T 0.6 -tm1 0.2 0.5 -tm2 0.3 0.5 -mind 10 -maxd 10 -ns 1 -nc 40 -change 0 -cp 50 -msg es_check_m_impact_md10_m2_0_$m0 #-mdl  8 9 10

done


# for m0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8
# do
# for m2 in 0.1
# do
# m1=$(echo "1-$m0-$m2"|bc)

# echo m0: $m0
# echo m1: $m1
# echo m2: $m2

#     time python sim-ray.py -modelname mask -itemname es -n 50000 -e 100 -m $m0 $m1 $m2 -T 0.6 -tm1 0.2 0.3 1 -tm2 0.3 0.2 1 -mind 15 -maxd 15 -ns 1 -nc 40 -change 0 -cp 50 -msg es_check_inout_impact_md15_m2_0.1_$m0 #-mdl  8 9 10

# done
# done

# for m0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7
# do
# for m2 in 0.2
# do
# m1=$(echo "1-$m0-$m2"|bc)

# echo m0: $m0
# echo m1: $m1
# echo m2: $m2

#     time python sim-ray.py -modelname mask -itemname es -n 50000 -e 100 -m $m0 $m1 $m2 -T 0.6 -tm1 0.2 0.3 1 -tm2 0.3 0.2 1 -mind 15 -maxd 15 -ns 1 -nc 40 -change 0 -cp 50 -msg es_check_inout_impact_md15_m2_0.2_$m0 #-mdl  8 9 10

# done
# done