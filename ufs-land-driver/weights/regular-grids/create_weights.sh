#!/bin/sh -l
#
# -- Request n cores
#SBATCH --ntasks=2
#
# -- Specify queue
#SBATCH -q debug
#
# -- Specify a maximum wallclock
#SBATCH --time=0:05:00
#
# -- Specify under which account a job should run
#SBATCH --account=fv3-cpu
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=ufs-land-create_weights
#SBATCH -o ufs-land-create_weights.out
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.

module purge
module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
module load stack-oneapi/2024.2.1
module load stack-intel-oneapi-mpi/2021.13
module load esmf/8.8.0
module load ncl/6.6.2

# set parameters for weights generation
#
# grid_string            : regular grid name, e.g., global_0.1
# datm_source            : ERA5 or CORe or GDAS or CDAS
# datm_source_file       : a datm source file to extract info for SCRIP file
# destination_scrip_path : location of the destination SCRIP file

grid_string="global_0.1"
datm_source="ERA5"
datm_source_file="/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/original/1980/ERA5_forcing_1980-01-01.nc"
destination_scrip_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

# create scrip file for data atmosphere sources, will replace existing file

datm_scrip_file=$datm_source"_SCRIP.nc"

# create the strings for the ncl parameter file

echo "datm_scrip_file = $datm_scrip_file" > regrid_parameter_assignment
echo "datm_source_file = $datm_source_file" >> regrid_parameter_assignment

eval "time ncl create_datm_scrip.ncl"

output_path=$grid_string"/"

destination_scrip_file=$destination_scrip_path$output_path"ufs-land_"$grid_string"_SCRIP.nc"

if [ -d $output_path ]; then 
  echo "BEWARE: output_path directory exists and overwriting is allowed"
else
  mkdir -p $output_path
fi

# create weights file for bilinear interpolation

interpolation_method="bilinear"
weights_filename=$datm_source"-"$grid_string"_"$interpolation_method"_wts.nc"
echo "Creating weights file: "$weights_filename

srun -n $SLURM_NTASKS time ESMF_RegridWeightGen --netcdf4 --ignore_degenerate \
       --source $datm_scrip_file \
       --destination $destination_scrip_file \
       --weight $output_path$weights_filename --method $interpolation_method

# create weights file for neareststod interpolation

interpolation_method="neareststod"
weights_filename=$datm_source"-"$grid_string"_"$interpolation_method"_wts.nc"
echo "Creating weights file: "$weights_filename

srun -n $SLURM_NTASKS time ESMF_RegridWeightGen --netcdf4 --ignore_degenerate \
       --source $datm_scrip_file \
       --destination $destination_scrip_file \
       --weight $output_path$weights_filename --method $interpolation_method

rm regrid_parameter_assignment
rm PET*
