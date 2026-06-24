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

# set parameters for grid generation
#
# grid_string            : regular grid name, e.g., global_0.1
# datm_source            : ERA5 or CORe or GDAS or CDAS
# datm_source_path       : path to datm source files
# static_file_path       : location of the destination static file
# yyyy,mm,dd,hh,timestep : create a cold start file timestep before yyyy-mm-dd_hh

grid_string="global_0.1"
datm_source="ERA5"
datm_source_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/"
static_file_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/"
yyyy=2021
mm=1
dd=1
hh=0
timestep=60

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

output_path=$grid_string"/"

if [ -d $output_path ]; then 
  echo "ERROR: directory $output_path exists and overwriting is prevented"
  echo "ERROR: remove $output_path and resubmit"
  exit 2
else
  mkdir -p $output_path
fi

# create static filename and check if it exists

static_filename=$static_file_path$output_path"ufs-land_"$grid_string"_static_fields.nc"

if [[ -e $static_filename ]]; then 
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
echo "ic_preamble = "$output_path$datm_source"-"$grid_string >> parameter_assignment
echo "datm_source = "$datm_source >> parameter_assignment
echo "datm_source_path = "$datm_source_path$output_path$datm_source"-"$grid_string >> parameter_assignment
echo "static_filename = "$static_filename >> parameter_assignment

eval "time ncl create_cold_start.ncl"

rm -f parameter_assignment
