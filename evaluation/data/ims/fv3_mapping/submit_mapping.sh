#!/bin/bash

#BATCH --job-name=fv3_mapping
#SBATCH -t 03:55:00
#SBATCH -A bigmem
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o mapping.out
#SBATCH -e mapping.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

source mods_bash 

./create_fv3_mapping.exe



