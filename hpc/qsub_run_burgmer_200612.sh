#!/bin/bash
#------------------------------------------
#$ -N pts_all
#$ -o iman.out$JOB_ID
#$ -j y
#$ -S /bin/bash
#$ -l h_rt=14:00:00 
#$ -l h_vmem=6G
#$ -pe smp 1
#$ -t 1-152

date
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load R/3.5.1

cd /home/clarka/ptstability_analyze_burgmer/
./analyze_burgmer.R $SGE_TASK_ID
date
