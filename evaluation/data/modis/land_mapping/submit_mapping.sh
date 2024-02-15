#!/bin/bash

#BATCH --job-name=fv3_mapping_albedo
#SBATCH -t 05:30:00
#SBATCH -A bigmem
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o mapping_albedo.out
#SBATCH -e mapping_albedo.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 

./create_land_mapping.exe



