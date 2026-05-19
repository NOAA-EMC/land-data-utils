#!/bin/bash

#BATCH --job-name=modis_run
#SBATCH -t 07:55:00
#SBATCH -A bigmem
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o modis_run.out
#SBATCH -e modis_run.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

./run_modis_albedo.sh



