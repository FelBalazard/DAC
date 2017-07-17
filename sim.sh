#!/bin/bash
touch result
for K in 0.002 0.006 0.01;do
for N in 500 1500 3000 5000;do
for shsim in `seq 0 9`; do
qsub simsl.sh $K $N $shsim
done
done
done
