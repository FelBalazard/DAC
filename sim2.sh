#!/bin/bash
touch result_auc
for sh in 0.18 0.4 0.66;do
for N in 500 1500 3000 5000;do
qsub sim2sl.sh $sh $N
done
done
