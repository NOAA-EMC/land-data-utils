#!/bin/sh -l
#
# -- Specify queue
#SBATCH -q batch
#
# -- Specify under which account a job should run
#SBATCH --account=fv3-cpu
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=ufs-land-create_datm
#SBATCH -o ufs-land-create_datm.out
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.
#
# -- Request tasks, this should correspond to the number of lines in your regrid-tasks file
#SBATCH --ntasks=20
#
# -- Specify a maximum wallclock
# -- global_0.1  : ~1 hours
#
#SBATCH --time=1:00:00

module purge
module load ncl/6.6.2

# set parameters for datm generation
#
# grid_string                 : regular grid name, e.g., global_0.1
# datm_source                 : ERA5 or CORe or GDAS or CDAS
# datm_source_path            : path to datm source files
# elevation_source_filename   : path to datm elevation file
# static_file_path            : location of the destination static file
# weights_path                : location of the ESMF regrid weights file
# interpolation_method1       : primary interpolation method
# interpolation_method2       : secondary interpolation method
# preicp_interpolation_method : "method1","method2","all"
# regrid_tasks_file           : file to split the regrid task across processors

grid_string="global_0.1"
datm_source="ERA5"
datm_source_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/original/"
elevation_source_filename="/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/original/elevation/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc"
static_file_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/"
weights_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/weights/"
interpolation_method1="bilinear"
interpolation_method2="neareststod"
precip_interpolation_method="method1"
regrid_tasks_file="regrid-tasks.ERA5"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

# check if elevation file exists

if [[ -e $elevation_source_filename ]]; then 
  echo "using elevation_source_filename:"$elevation_source_filename
else
  echo "ERROR: elevation_source_filename does not exist"
  exit 1
fi

output_path=$grid_string"/"

if [ -d $output_path ]; then 
  echo "BEWARE: output_path directory exists and overwriting is allowed"
else
  echo "creating directory: "$output_path
  mkdir -p $output_path
fi

# create weights filename for method 1 and check if it exists

weights_method1_filename=$weights_path$output_path$datm_source"-"$grid_string"_"$interpolation_method1"_wts.nc"

if [[ -e $weights_method1_filename ]]; then 
  echo "using weights_method1_filename:"$weights_method1_filename
else
  echo "ERROR: weights_method1_filename does not exist: "$weights_method1_filename
  exit 4
fi

# create weights filename for method 2 and check if it exists

weights_method2_filename=$weights_path$output_path$datm_source"-"$grid_string"_"$interpolation_method2"_wts.nc"

if [[ -e $weights_method2_filename ]]; then 
  echo "using weights_method2_filename:"$weights_method2_filename
else
  echo "ERROR: weights_method2_filename does not exist: "$weights_method2_filename
  exit 5
fi

# create static filename and check if it exists

static_filename=$static_file_path$output_path"ufs-land_"$grid_string"_static_fields.nc"

if [[ -e $static_filename ]]; then 
  echo "using static_filename:"$static_filename
else
  echo "ERROR: static_filename does not exist: "$static_filename
  exit 6
fi

if [[ $precip_interpolation_method = "all"     ]] ||  
   [[ $precip_interpolation_method = "method1" ]] ||
   [[ $precip_interpolation_method = "method2" ]]; then 
  echo "using precip_interpolation_method:"$precip_interpolation_method
else
  echo "ERROR: precip_interpolation_method not set correctly: "$precip_interpolation_method
  exit 7
fi

# create elevation filename

elevation_filename="elevation_"$datm_source"_"$grid_string".nc"

echo "creating elevation_filename:"$elevation_filename

# create elevation difference file

echo "static_filename = $static_filename" > regrid_parameter_assignment
echo "elevation_source_filename = $elevation_source_filename" >> regrid_parameter_assignment
echo "datm_source = $datm_source" >> regrid_parameter_assignment
echo "weights_filename = $weights_method1_filename" >> regrid_parameter_assignment
echo "elevation_filename = $elevation_filename" >> regrid_parameter_assignment

eval "/usr/bin/time ncl create_vector_elevation.ncl"

# regrid the source data atmosphere

echo "Creating datm files"

echo "elevation_filename = $elevation_filename" > regrid_parameter_assignment
echo "weights_method1_filename = $weights_method1_filename" >> regrid_parameter_assignment
echo "weights_method2_filename = $weights_method2_filename" >> regrid_parameter_assignment
echo "precip_interpolation_method = $precip_interpolation_method" >> regrid_parameter_assignment
echo "interpolation_method1 = $interpolation_method1" >> regrid_parameter_assignment
echo "interpolation_method2 = $interpolation_method2" >> regrid_parameter_assignment
echo "datm_source_path = $datm_source_path" >> regrid_parameter_assignment
echo "output_preamble = "$output_path$datm_source"-"$grid_string >> regrid_parameter_assignment

srun -l --multi-prog $regrid_tasks_file

rm -f $elevation_filename
rm -f regrid_parameter_assignment

