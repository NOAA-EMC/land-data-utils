#!/bin/sh -l
#
# -- Request n cores
#SBATCH --partition=u1-service
#SBATCH --ntasks=1
#
# -- Specify queue
#SBATCH -q batch
#
# -- Specify a maximum wallclock
#SBATCH --time=0:05:00
#
# -- Specify under which account a job should run
#SBATCH --account=fv3-cpu
#
# -- Set the name of the job, or Slurm will default to the name of the script
#SBATCH --job-name=create_climo_mswep
#SBATCH -o create_climo_mswep.out
#
# -- Tell the batch system to set the working directory to the current working directory
#SBATCH --chdir=.

module purge
module load nco

data_directory="/scratch4/NCEPDEV/land/data/evaluation/MSWEP/original/monthly/"
monthly_file=$data_directory/"MSWEP.v3.monthly_climatology.2015-2024.nc"
annual_file=$data_directory/"MSWEP.v3.climatology.2015-2024.nc"

if [ -e $monthly_file ]; then 
  echo "BEWARE: remove $monthly_file"
  rm -f $monthly_file
fi
if [ -e $annual_file ]; then 
  echo "BEWARE: remove $annual_file"
  rm -f $annual_file
fi

for mm in {01..12}
do

  imm=$((10#$mm-1))

#create a monthly mean file

  echo "creating: ${mm}.mean.nc"
  ncra -h $data_directory/201[5-9]*${mm}.nc $data_directory/202[0-4]*${mm}.nc ${mm}.mean.nc

done  # end of mm loop

ncrcat -h *mean.nc $monthly_file
ncatted -h -a description,global,c,c,'MSWEP monthly climatology 2015-2024' $monthly_file

ncra -h -y ttl $monthly_file $annual_file
ncatted -h -a description,global,c,c,'MSWEP annual climatology 2015-2024' $annual_file
ncatted -h -a units,precipitation,m,c,'mm' $annual_file

rm *mean.nc

