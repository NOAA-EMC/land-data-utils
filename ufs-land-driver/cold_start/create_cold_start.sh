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
#SBATCH --job-name=ufs-land-create_cold_start
#SBATCH -o ufs-land-create_cold_start.out
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

if [ $grid_version = "hr3" ]; then 
  grid=$atm_res.$ocn_res"_hr3"
elif [ $grid_version = "AQM" ]; then 
  grid=$atm_res.$grid_version
elif [ $grid_version = "ARC" ]; then 
  grid=$atm_res.$grid_version
else
  echo "ERROR: unknown grid_version $grid_version"
  exit 2
fi

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
  exit 3
fi

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
