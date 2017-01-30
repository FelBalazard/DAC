#!/bin/bash
#$ -N sim
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -V
#This script is called by sim3.sh

 
Rscript ./sim_power.R $1
