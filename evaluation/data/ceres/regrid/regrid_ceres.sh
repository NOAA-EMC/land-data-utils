#!/bin/sh -l
#
# -- Request n cores
#SBATCH --ntasks=1
#
# -- Specify queue
#SBATCH -q debug
#
# -- Specify under which account a job should run
#SBATCH --account=fv3-cpu
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=regrid_ceres
#SBATCH -o regrid_ceres.out
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.
#
# -- Specify a maximum wallclock
# -- C96  : ~5 minutes  C96 conus : ~5 minutes
# -- C192 : ~30 minutes
# -- C384 : ~1 hours
# -- C768 : ~1.5 hours
# -- C1152: ~3 hours
#
#SBATCH --time=0:05:00

module purge
module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
module load stack-oneapi/2024.2.1
module load stack-intel-oneapi-mpi/2021.13
module load ncl/6.6.2

yyyy_begin="2011"
yyyy_end="2025"
atm_res="C96"
ocn_res="mx100"
grid_version="20231027"
grid_extent="total"
subset_name="conus"
data_source="CERES"
data_source_directory="/scratch4/NCEPDEV/land/data/evaluation/CERES/yearly/"
data_destination_directory="/scratch4/NCEPDEV/land/data/evaluation/CERES/"
weights_directory="/scratch4/NCEPDEV/land/data/evaluation/CERES/weights/"
vector_directory="/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/"
interpolation_method="neareststod"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

if [[ $grid_version == "20231027" ]] ; then 
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

destination_directory=$data_destination_directory$grid_string"/"
weights_filename=$weights_directory$grid_string"/"$data_source"-"$grid_string"_"$interpolation_method"_wts.nc"
vector_filename=$vector_directory$grid_string"/ufs-land_"$grid_string"_corners.nc"

echo "data_source_directory = $data_source_directory" > regrid_parameter_assignment
echo "destination_directory = $destination_directory" >> regrid_parameter_assignment
echo "weights_filename = $weights_filename" >> regrid_parameter_assignment
echo "vector_filename = $vector_filename" >> regrid_parameter_assignment
echo "grid_string = $grid_string" >> regrid_parameter_assignment
echo "yyyy_begin = $yyyy_begin" >> regrid_parameter_assignment
echo "yyyy_end = $yyyy_end" >> regrid_parameter_assignment

if [ -d $destination_directory ]; then 
  echo "BEWARE: destination_directory directory exists and overwriting is allowed"
else
  mkdir -p $destination_directory
fi

cmdparm="'variable_set="\"radiation"\"' "
eval "time ncl regrid_ceres.ncl $cmdparm"

cmdparm="'variable_set="\"clouds"\"' "
eval "time ncl regrid_ceres.ncl $cmdparm"

cmdparm="'variable_set="\"surface"\"' "
eval "time ncl regrid_ceres.ncl $cmdparm"

rm -f regrid_parameter_assignment
