#!/bin/sh -l
#
# -- Request n cores
#SBATCH --ntasks=2
#
# -- Specify queue
#SBATCH -q debug
#
# -- Specify under which account a job should run
#SBATCH --partition=u1-service
#SBATCH --account=fv3-cpu
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=ufs-land-create_inputs
#SBATCH -o ufs-land-create_inputs.out
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.
#
#SBATCH --time=0:10:00

module purge
module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
module load stack-oneapi/2024.2.1
module load stack-intel-oneapi-mpi/2021.13
module load esmf/8.8.0
module load ncl/6.6.2

# set parameters for grid generation
#
# atm_res      : fv3 grid resolution
# ocn_res      : ocean resolution, not used for AQM or ARC regional grids
# grid_version : 20231027 - append directory date string
#                AQM - AQM regional grid
#                ARC - UFS-Arctic regional grid
# fixfile_path : top level path for fix files
# grid_extent  : total - use all grids (e.g., global or entire regional)
#                subset - regional cutout, limits below
# subset_name  : if subset, name for subset, e.g., conus
# subset_maxlat: cutout maximum latitude
# subset_minlat: cutout minimum latitude
# subset_maxlon: cutout maximum longitude
# subset_minlon: cutout minimum longitude

grid_string="global_0.1"
atm_res="C1152"
ocn_res="mx025"
grid_version="20231027"
grid_extent="total"
subset_name=""
source_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

if [[ $grid_version == "20231027" ]] ; then 
  source_grid_string=$atm_res.$ocn_res
  if [[ $grid_extent == "subset" ]]; then
    source_grid_string=$source_grid_string.$subset_name
  fi
elif [[ $grid_version == "AQM" ]] || [[ $grid_version == "ARC" ]]; then 
  source_grid_string=$atm_res.$grid_extent
else
  echo "ERROR: unknown combination"
  echo "ERROR: grid_version = $grid_version"
  echo "ERROR: grid_extent = $grid_extent"
  echo "NOTE:  subset not currently supported for regional grids"
  exit 1
fi

output_path=$grid_string"/"

if [ -d $output_path ]; then 
  echo "ERROR: directory $output_path exists and overwriting is prevented"
  echo "ERROR: remove $output_path and resubmit"
  exit 2
else
  mkdir -p $output_path
fi

# create the strings for the ncl parameter file

echo "grid_string = $grid_string" > regrid_parameter_assignment
echo "output_path = $output_path" >> regrid_parameter_assignment

# create the grid corners file, this is used in follow-on ncl scripts and for other tools

eval "time ncl create_corners.ncl"

# create the SCRIP file, this is used for ESMF regridding

eval "time ncl create_scrip.ncl"

# create weights file for neareststod interpolation of static file

destination_scrip_file=$output_path"ufs-land_"$grid_string"_SCRIP.nc"
source_scrip_file=$source_path$source_grid_string"/ufs-land_"$source_grid_string"_SCRIP.nc"
interpolation_method="neareststod"
weights_filename=$source_grid_string"-"$grid_string"_"$interpolation_method"_wts.nc"
echo "Creating weights file: "$weights_filename
echo "using source file: "$source_scrip_file
echo "and destination file: "$destination_scrip_file

srun -n $SLURM_NTASKS time ESMF_RegridWeightGen --netcdf4 --ignore_degenerate \
       --source $source_scrip_file \
       --destination $destination_scrip_file \
       --weight $output_path$weights_filename --method $interpolation_method

rm PET*

source_static_file=$source_path$source_grid_string"/ufs-land_"$source_grid_string"_static_fields.nc"

echo "grid_string = $grid_string" > regrid_parameter_assignment
echo "output_path = $output_path" >> regrid_parameter_assignment
echo "source_static_file = $source_static_file" >> regrid_parameter_assignment
echo "weights_filename = $output_path$weights_filename" >> regrid_parameter_assignment

# create the static fields file, this is used to create the inputs to the driver

eval "time ncl extract_static.ncl"

rm regrid_parameter_assignment
rm $output_path$weights_filename
