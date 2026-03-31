#!/bin/bash

#BATCH --job-name=snodas_run
#SBATCH -t 03:55:00
#SBATCH -A bigmem
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o snodas_run.out
#SBATCH -e snodas_run.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

./run_snodas.sh
