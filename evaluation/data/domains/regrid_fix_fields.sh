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
#SBATCH --job-name=regrid_fix_fields
#SBATCH -o regrid_fix_fields.out
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
# atm_res               : fv3 grid resolution
# ocn_res               : ocean resolution, not used for AQM or ARC regional grids
# out_res               : output resolution of ncep grid: 1p00, 0p50, 0p25, sfc
# variable2regrid       : variable name, both in the filename and name in file
# destination_directory : location to put the resulting file
# grid_version          : 20231027 - directory from fix archive
# interpolation_method  : ESMF method, e.g., bilinear,neareststod

atm_res="C1152"
ocn_res="mx025"
out_res="0p25"
variable2regrid="vegetation_type"
destination_directory="/scratch4/NCEPDEV/land/data/evaluation/domains/output_grids/"
grid_version="20231027"
interpolation_method="neareststod"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

mosaic_file="/scratch4/NCEPDEV/land/data/fix/mosaic_files/${atm_res}_mosaic.nc"

if [[ $out_res == "1p00" ]] ; then 
  destination_file="/scratch4/NCEPDEV/land/data/evaluation/domains/output_grids/grid_specs/lat-lon.land.1p00.nc"
elif [[ $out_res == "0p50" ]] ; then 
  destination_file="/scratch4/NCEPDEV/land/data/evaluation/domains/output_grids/grid_specs/lat-lon.land.0p50.nc"
elif [[ $out_res == "0p25" ]] ; then 
  destination_file="/scratch4/NCEPDEV/land/data/evaluation/domains/output_grids/grid_specs/lat-lon.land.0p25.nc"
elif [[ $out_res == "sfc" ]] ; then 
  destination_file="/scratch4/NCEPDEV/land/data/evaluation/domains/output_grids/grid_specs/gaussian.land.3052x1536.nc"
  out_res="gaussian.3052x1536"
else
  echo "ERROR: unknown out_res = $out_res"
  exit 1
fi

destination_scrip_file="GFS_output_grid_${out_res}_SCRIP.nc"

# create the ncl parameter file

echo "destination_scrip_file = $destination_scrip_file" > regrid_parameter_assignment
echo "destination_file = $destination_file" >> regrid_parameter_assignment

eval "time ncl create_output_scrip.ncl"

rm regrid_parameter_assignment

# create weights file

weights_filename=$atm_res.$ocn_res"-"$out_res"_"$interpolation_method"_wts.nc"

echo "Creating weights file: "$weights_filename
echo "using mosaic_file: "$mosaic_file
echo "and destination_scrip_file: "$destination_scrip_file

srun -n $SLURM_NTASKS time ESMF_RegridWeightGen --netcdf4 --ignore_degenerate \
       --source $mosaic_file \
       --destination $destination_scrip_file \
       --weight $weights_filename --method $interpolation_method \
       --ignore_unmapped

rm $destination_scrip_file
rm PET*

output_filename=$destination_directory$atm_res.$ocn_res"-"$out_res"."$variable2regrid".nc"
if [[ $variable2regrid == "vegetation_greenness" ]] ; then 
  output_filename=$destination_directory$atm_res.$ocn_res"-"$out_res".vegetation_fraction.nc"
fi

echo "Creating output_filename: $output_filename"
echo "atm_res = $atm_res"
echo "ocn_res = $ocn_res"
echo "grid_version = $grid_version"
echo "weights_filename = $weights_filename"
echo "variable2regrid = $variable2regrid"

echo "atm_res = $atm_res" > regrid_parameter_assignment
echo "ocn_res = $ocn_res" >> regrid_parameter_assignment
echo "grid_version = $grid_version" >> regrid_parameter_assignment
echo "weights_filename = $weights_filename" >> regrid_parameter_assignment
echo "variable2regrid = $variable2regrid" >> regrid_parameter_assignment
echo "output_filename = $output_filename" >> regrid_parameter_assignment

eval "time ncl regrid_fix_fields.ncl"

rm regrid_parameter_assignment
rm $weights_filename

