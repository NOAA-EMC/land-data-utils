#!/bin/bash

#BATCH --job-name=cmc_C96
#SBATCH -t 01:00:00
#SBATCH -A bigmem
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o cmc_run.out
#SBATCH -e cmc_run.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

./run_cmc_snow.sh



