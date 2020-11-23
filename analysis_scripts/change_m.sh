#!/bin/bash
source activate ytian
# for p_mask in 0.45
# do
# 	echo $p_mask
#     time python analysis.py -modelname mask -itemname es -m $p_mask 0.3 0.25 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.7 0.2 1 -mind 0 -maxd 10 -ns 100 -nc 40 -change 0 -msg es_0_to_10_baseline #-mdl  8 9 10
# done

# for p_mask in 0.45
# do
# 	echo $p_mask
#     time python analysis.py -modelname mask -itemname es -m $p_mask  0.55 -T 0.6 -tm1  0.3 1 -tm2 0.7 1 -mind 0 -maxd 10 -ns 50 -nc 40 -change 0 -msg change_vec_es #-mdl  8 9 10
# done


# for p_mask in 0.45
# do
# 	echo $p_mask
#     time python analysis.py -modelname mask -itemname pe -m $p_mask  0.55 -T 0.6 -tm1  0.3 1 -tm2 0.7 1 -mind 0 -maxd 10 -ns 50 -nc 40 -change 0 -msg change_vec_pe #-mdl  8 9 10
# done


echo es_strategy0
time python analysis.py -modelname mask -itemname es -n 500000 -e 1000 -th 0.05 -m 0 0 1 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg es_strategy0 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo es_strategy1
time python analysis.py -modelname mask -itemname es -n 500000 -e 1000 -th 0.05 -m 1 0 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg es_strategy1 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo es_strategy2
time python analysis.py -modelname mask -itemname es -n 500000 -e 1000 -th 0.05 -m 0 1 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg es_strategy2 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo es_strategy3
time python analysis.py -modelname mask -itemname es -n 500000 -e 1000 -th 0.05 -m 0.5 0.5 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg es_strategy3 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0


echo pe_strategy0
time python analysis.py -modelname mask -itemname pe -n 500000 -e 1000 -th 0.05 -m 0 0 1 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg pe_strategy0 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo pe_strategy1
time python analysis.py -modelname mask -itemname pe -n 500000 -e 1000 -th 0.05 -m 1 0 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg pe_strategy1 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo pe_strategy2
time python analysis.py -modelname mask -itemname pe -n 500000 -e 1000 -th 0.05 -m 0 1 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg pe_strategy2 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0

echo pe_strategy3
time python analysis.py -modelname mask -itemname pe -n 500000 -e 1000 -th 0.05 -m 0.5 0.5 0 -T 0.6 -tm1 0.3 0.5 1 -tm2 0.2 0.5 1 -mind 0 -maxd 10 -ns 50 -nc 40 -cp 100 -change 0 -msg pe_strategy3 #-mdl 2.6 2.64 2.68 2.72 2.76 2.8 2.84 2.88 2.92 2.96 3.0
