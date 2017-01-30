#!/bin/bash
touch result_power
for name in "X219" "X238" "X239" "X260" "X342" "X418" "X423";do
qsub sim3sl.sh $name 
done
