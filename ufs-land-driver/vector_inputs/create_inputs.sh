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
# grid_version : hr3 - append directory date string (not supporting other global options for now)
#                AQM - AQM regional grid
#                ARC - UFS-Arctic regional grid
# fixfile_path : top level path for fix files
# grid_extent  : global or conus for standard global FV3 or cutout regional
#                AQM for AQM regional grid
#                ARC for UFS-Arctic regional grid

atm_res="C918"
ocn_res="mx100"
grid_version="ARC"
fixfile_path="/scratch4/BMC/ufs-artic/Kristin.Barton/files/ufs_arctic_development/fix_files/mesh_files/C918/"
grid_extent="ARC"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

# set full fix file based on grid version

if [ $grid_version = "hr3" ]; then 
  fixfile_path=$fixfile_path"20231027/"
elif [ $grid_version = "AQM" ]; then 
  fixfile_path=$fixfile_path
elif [ $grid_version = "ARC" ]; then 
  fixfile_path=$fixfile_path
else
  echo "ERROR: unknown fixfile_path $fixfile_path"
  exit 1
fi

# the default location for output files is $atm_res.$ocn_res

if [ $grid_extent = "global" ]; then 
  output_path=$atm_res.$ocn_res"/"
elif [ $grid_extent = "AQM" ]; then 
  output_path=$atm_res.$grid_extent"/"
elif [ $grid_extent = "ARC" ]; then 
  output_path=$atm_res.$grid_extent"/"
elif [ $grid_extent = "conus" ]; then 
  output_path=$atm_res.$ocn_res.$grid_extent"/"
else
  echo "ERROR: unknown grid_extent $grid_extent"
  exit 2
fi

if [ -d $output_path ]; then 
  echo "ERROR: directory $output_path exists and overwriting is prevented"
  echo "ERROR: remove $output_path and resubmit"
  exit 3
else
  mkdir -p $output_path
fi

# create the strings for the ncl parameter file

echo "atm_res = $atm_res" > regrid_parameter_assignment
echo "ocn_res = $ocn_res" >> regrid_parameter_assignment
echo "grid_version = $grid_version" >> regrid_parameter_assignment
echo "output_path = $output_path" >> regrid_parameter_assignment
echo "fixfile_path = $fixfile_path" >> regrid_parameter_assignment
echo "grid_extent = $grid_extent" >> regrid_parameter_assignment

# create the grid corners file, this is used in follow-on ncl scripts and for other tools

eval "time ncl extract_corners.ncl"

# create the static fields file, this is used to create the inputs to the driver

eval "time ncl extract_static.ncl"

# create the SCRIP file, this is used for ESMF regridding

eval "time ncl create_scrip.ncl"


