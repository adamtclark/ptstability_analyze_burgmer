#!/bin/bash
#------------------------------------------
#$ -S /bin/bash
#$ -V
#$ -j yes
#$ -q all.q
#$ -l h_vmem=4G
#$ -N pts_all
#$ -m e
#$ -cwd
#$ -l h_rt=14:00:00 
#$ -pe smp 1
#$ -t 9-16

date
module load R

cd /cl_tmp/clarka/ptstability_analyze_burgmer/
./analyze_burgmer_casestudy.R $SGE_TASK_ID
date
