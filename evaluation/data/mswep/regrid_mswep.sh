#!/bin/sh -l
#
# -- Request n cores
#SBATCH --partition=u1-service
#SBATCH --ntasks=24
#
# -- Specify queue
#SBATCH -q debug
#
# -- Specify a maximum wallclock
#SBATCH --time=0:10:00
#
# -- Specify under which account a job should run
#SBATCH --account=fv3-cpu
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=regrid_mswep
#SBATCH -o regrid_mswep.out
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
# atm_res      : fv3 grid resolution
# ocn_res      : ocean resolution, not used for AQM or ARC regional grids
# grid_version : 20231027 - append directory date string
#                AQM - AQM regional grid
#                ARC - UFS-Arctic regional grid
# grid_extent  : total - use all grids (e.g., global or entire regional)
#                subset - regional cutout, limits below
# subset_name  : if subset, name for subset, e.g., conus
# data_source  : mswep
# data_source_file       : a datm source file to extract info for SCRIP file
# interpolation_method   : ESMF method, e.g., bilinear,neareststod
# destination_scrip_path : location of the destination SCRIP file

yyyy_begin="2010"
yyyy_end="2010"
atm_res="C192"
ocn_res="mx025"
grid_version="20231027"
grid_extent="total"
subset_name="conus"
data_source="MSWEP"
data_source_directory="/scratch4/NCEPDEV/land/data/evaluation/MSWEP/original/monthly/"
data_destination_directory="/scratch4/NCEPDEV/land/data/evaluation/MSWEP/"
interpolation_method="conserve"
destination_scrip_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

if [[ $grid_version == "20231027" ]] || [[ $grid_version == "20240917" ]] ; then 
  grid_string=$atm_res.$ocn_res
  if [[ $grid_extent == "subset" ]]; then
    grid_string=$grid_string.$subset_name
  fi
elif [[ $grid_version == "AQM" ]] || [[ $grid_version == "ARC" ]]; then 
  grid_string=$atm_res.$grid_extent
else
  echo "ERROR: unknown combination"
  echo "ERROR: grid_version = $grid_version"
  echo "ERROR: grid_extent = $grid_extent"
  echo "NOTE:  subset not currently supported for regional grids"
  exit 1
fi

output_path=$grid_string"/"

destination_scrip_file=$destination_scrip_path$output_path"ufs-land_"$grid_string"_SCRIP.nc"

if [ -d $output_path ]; then 
  echo "BEWARE: output_path directory exists and overwriting is allowed"
else
  mkdir -p $output_path
fi

# create source scrip file

data_scrip_file=$data_source"_SCRIP.nc"
data_source_file=${data_source_directory}"/201001.nc"
weights_filename=$data_source"-"$grid_string"_"$interpolation_method"_wts.nc"

if [[ -e $output_path$weights_filename ]]; then

  echo "$weights_filename exists so no need to create it"

else

# create the ncl parameter file

  echo "data_scrip_file = $data_scrip_file" > regrid_parameter_assignment
  echo "data_source_file = $data_source_file" >> regrid_parameter_assignment

  eval "time ncl create_source_scrip.ncl"

  rm regrid_parameter_assignment

# create weights file

  echo "Creating weights file: "$weights_filename

  srun -n $SLURM_NTASKS time ESMF_RegridWeightGen --netcdf4 --ignore_degenerate \
       --source $data_scrip_file \
       --destination $destination_scrip_file \
       --weight $output_path$weights_filename --method $interpolation_method \
       --ignore_unmapped

  rm $data_scrip_file
  rm PET*

fi   # weights file exists

destination_directory=$data_destination_directory$grid_string"/"

echo "data_source_directory = $data_source_directory" > regrid_parameter_assignment
echo "destination_directory = $destination_directory" >> regrid_parameter_assignment
echo "weights_filename = $output_path$weights_filename" >> regrid_parameter_assignment
echo "grid_string = $grid_string" >> regrid_parameter_assignment
echo "yyyy_begin = $yyyy_begin" >> regrid_parameter_assignment
echo "yyyy_end = $yyyy_end" >> regrid_parameter_assignment

if [ -d $destination_directory ]; then 
  echo "BEWARE: destination_directory directory exists and overwriting is allowed"
else
  mkdir -p $destination_directory
fi

eval "time ncl regrid_mswep.ncl"

rm regrid_parameter_assignment
rm -Rf $output_path

