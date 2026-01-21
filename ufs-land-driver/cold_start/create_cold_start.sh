#!/bin/sh -l
#
# -- Request n cores
#SBATCH --ntasks=1
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
#SBATCH --job-name=create_cold_start
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.

module purge
module load ncl/6.6.2

atm_res="C96"
ocn_res="mx100"
grid_version="hr3"
grid_extent="global"
datm_source="ERA5"
datm_source_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/"
static_file_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/"
yyyy=1980
mm=1
dd=1
hh=0
timestep=60

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

# the default location for output files is $atm_res.$ocn_res
# create a variable $res, e.g., C96.mx100 or C96.mx100.conus

if [ $grid_extent = "global" ]; then 
  res=$atm_res.$ocn_res
else
  res=$atm_res.$ocn_res.$grid_extent
fi

# create a variable $grid, e.g., C96.mx100_hr3 or C96.mx100.conus_hr3
# create a variable $output_path, e.g., C96.mx100/ or C96.mx100.conus/

grid=$res"_"$grid_version
output_path=$res"/"

if [ -d $output_path ]; then 
  echo "BEWARE: output_path directory exists and overwriting is allowed"
else
  echo "creating directory: "$output_path
  mkdir -p $output_path
fi

# create static filename and check if it exists

static_filename=$static_file_path$output_path"ufs-land_"$grid"_static_fields.nc"

if [ -e $static_filename ]; then 
  echo "using static_filename:"$static_filename
else
  echo "ERROR: static_filename does not exist: "$static_filename
  exit 3
fi

# create the cold start initial conditions file

echo "Creating cold start IC file"

echo "yyyy = $yyyy" > parameter_assignment
echo "mm = $mm" >> parameter_assignment
echo "dd = $dd" >> parameter_assignment
echo "hh = $hh" >> parameter_assignment
echo "timestep = $timestep" >> parameter_assignment
echo "ic_preamble = "$output_path$datm_source"-"$grid >> parameter_assignment
echo "datm_source = "$datm_source >> parameter_assignment
echo "datm_source_path = "$datm_source_path$output_path$datm_source"-"$grid >> parameter_assignment
echo "static_filename = "$static_filename >> parameter_assignment

eval "time ncl create_cold_start.ncl"

rm -f parameter_assignment
