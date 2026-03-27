#!/bin/bash

#BATCH --job-name=IMS1km_mapping
#SBATCH -t 07:55:00
#SBATCH -A bigmem
#SBATCH -A da-cpu
#SBATCH --qos=batch
#SBATCH -o IMS1km_mapping.out
#SBATCH -e IMS1km_mapping.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

module load ncl

ncl create_pt_ims1km_mapping_prototype.ncl



