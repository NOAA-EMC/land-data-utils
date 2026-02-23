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
#SBATCH --job-name=ufs-land-create_inputs
#SBATCH -o ufs-land-create_inputs.out
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.
#
# -- Specify a maximum wallclock
# -- C96  : ~1 minute
# -- C192 : ~1 minute
# -- C384 : ~2 minutes
# -- C768 : ~6 minutes
# -- C1152: ~15 minutes
#
#SBATCH --time=0:01:00

module purge
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

atm_res="C96"
ocn_res="mx100"
grid_version="20231027"
fixfile_path="/scratch3/NCEPDEV/global/role.glopara/fix/orog/"
grid_extent="subset"
subset_name="conus"
subset_maxlat="53.0"
subset_minlat="25.0"
subset_maxlon="293.0"
subset_minlon="235.0"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

# set full fix file based on grid version

if [[ $grid_version == "20231027" ]] ; then 
  fixfile_path=$fixfile_path$grid_version"/"
  is_global="True"
  grid_string=$atm_res.$ocn_res
  if [[ $grid_extent == "subset" ]]; then
    grid_string=$grid_string.$subset_name
  fi
elif [[ $grid_version == "AQM" ]] || [[ $grid_version == "ARC" ]]; then 
  grid_string=$atm_res.$grid_extent
  is_global="False"
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

echo "atm_res = $atm_res" > regrid_parameter_assignment
echo "ocn_res = $ocn_res" >> regrid_parameter_assignment
echo "grid_string = $grid_string" >> regrid_parameter_assignment
echo "output_path = $output_path" >> regrid_parameter_assignment
echo "fixfile_path = $fixfile_path" >> regrid_parameter_assignment
echo "grid_extent = $grid_extent" >> regrid_parameter_assignment
echo "is_global = $is_global" >> regrid_parameter_assignment
echo "subset_maxlat = $subset_maxlat" >> regrid_parameter_assignment
echo "subset_minlat = $subset_minlat" >> regrid_parameter_assignment
echo "subset_maxlon = $subset_maxlon" >> regrid_parameter_assignment
echo "subset_minlon = $subset_minlon" >> regrid_parameter_assignment

# create the grid corners file, this is used in follow-on ncl scripts and for other tools

eval "time ncl extract_corners.ncl"

# create the static fields file, this is used to create the inputs to the driver

eval "time ncl extract_static.ncl"

# create the SCRIP file, this is used for ESMF regridding

eval "time ncl create_scrip.ncl"

rm regrid_parameter_assignment
