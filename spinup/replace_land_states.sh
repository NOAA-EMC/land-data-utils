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
#SBATCH --job-name=replace_land_states
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.
#SBATCH -o replace_land_states.out

module purge
module load stack-oneapi/2024.2.1 stack-intel-oneapi-mpi/2021.13 py-numpy py-netcdf4

# Set the start and end times and time increment for sfc_data files

start_yyyy=1993
start_mm=03
start_dd=01
start_hh=00
end_yyyy=1993
end_mm=03
end_dd=01
end_hh=00
time_increment="1 month"

# set the location of the original non-spinup sfc_data files, these will be copied to a new location

original_sfc_data_path="/path-to-original/"

# set the location where you want the sfc_data files with the spinup states to be stored

updated_sfc_data_path="/path-to-updated/"

# set the path to the vector files from the spinup

spinup_path="/path-to-spinup/"

# set the spinup file format ["restart","output"]

spinup_file_format="restart"

#################################################################################
#  shouldn't need to modify anything below
#################################################################################

current_yyyy=$start_yyyy
current_mm=$start_mm
current_dd=$start_dd
current_hh=$start_hh
current_yyyymmdd="$current_yyyy$current_mm$current_dd"
current_yyyymmddhh="$current_yyyy$current_mm$current_dd$current_hh"

while [ $current_yyyymmddhh -le $end_yyyy$end_mm$end_dd$end_hh ]; do

  echo "Starting state replacment for $current_yyyymmddhh"

  mkdir -p $updated_sfc_data_path/${current_yyyymmdd}
  cp $original_sfc_data_path/${current_yyyymmdd}/sfc_data* $updated_sfc_data_path/${current_yyyymmdd}

  sfc_data_path="$updated_sfc_data_path/${current_yyyymmdd}/"
  if [ $spinup_file_format = "restart" ]; then 
    spinup_file=$spinup_path"ufs_land_restart.${current_yyyy}-${current_mm}-${current_dd}_${current_hh}-00-00.nc"
  elif [ $spinup_file_format = "output" ]; then 
    spinup_file=$spinup_path"ufs_land_output.${current_yyyy}-${current_mm}-${current_dd}_${current_hh}-00-00.nc"
  else
    echo "ERROR: spinup format unknown: "$spinup_file_format
    exit 1
  fi

  python replace_land_states.py $sfc_data_path $spinup_file

# Increment the time

  current_yyyy=$(date -d "${current_yyyymmdd} ${current_hh} +${time_increment}" +"%Y")
  current_mm=$(date -d "${current_yyyymmdd} ${current_hh} +${time_increment}" +"%m")
  current_dd=$(date -d "${current_yyyymmdd} ${current_hh} +${time_increment}" +"%d")
  current_hh=$(date -d "${current_yyyymmdd} ${current_hh} +${time_increment}" +"%H")
  current_yyyymmdd="$current_yyyy$current_mm$current_dd"
  current_yyyymmddhh="$current_yyyy$current_mm$current_dd$current_hh"

done
