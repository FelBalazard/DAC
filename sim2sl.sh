#!/bin/bash
#$ -N sim
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -V
#This script is called by sim.sh

 
Rscript ./sim_sh_N.R $1 $2
