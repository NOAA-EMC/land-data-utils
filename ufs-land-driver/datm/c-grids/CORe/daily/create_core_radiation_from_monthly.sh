#!/bin/sh -l
#
# -- Specify queue
#SBATCH -q batch
#
# -- Specify under which account a job should run
#SBATCH --account=fv3-cpu
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=create_core_radiation
#SBATCH -o create_core_radiation.out
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.
#
# -- Request tasks, this should correspond to the number of lines in your regrid-tasks file
#SBATCH --ntasks=93
#
#
#SBATCH --time=1:30:00

module purge
module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
module load stack-oneapi/2024.2.1 
module load stack-intel-oneapi-mpi/2021.13 
module load netcdf-fortran/4.6.1

srun -l --multi-prog radiation-tasks.CORe

