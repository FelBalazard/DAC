#!/bin/bash
touch result_auc
for sh in 0.18 0.4 0.66;do
for N in 500 1500 3000 5000;do
for shsim in `seq 0 9`; do
qsub sim2sl.sh $sh $N $shsim
done
done
done
