#!/bin/sh -l
#
#SBATCH --partition=u1-service
#SBATCH --ntasks=1
#SBATCH -q batch
#SBATCH -t 24:00:00
#SBATCH --account=fv3-cpu
#SBATCH --job-name=CORE-download1
#SBATCH --chdir=.
#SBATCH -o 2025.out

ens_member=12
source_dir="/NCEPDEV/cpc-om/Permanent/Leigh.Zhang/core/flux"

for yyyy in {2025..2025}
do

for mm in {07..12}
do
   echo "starting $source_dir/flux_${yyyy}${mm}.tar"
   htar -tvf $source_dir/flux_${yyyy}${mm}.tar > tempfile
   grep mem0$ens_member tempfile | cut -c 68-92 > filelist.${yyyy}${mm}
   htar -xvf $source_dir/flux_${yyyy}${mm}.tar -L filelist.${yyyy}${mm}
   rm tempfile filelist.${yyyy}${mm}

done

done
