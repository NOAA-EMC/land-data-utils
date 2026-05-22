#!/bin/sh -l
#

ens_member=12
source_dir="/NCEPDEV/cpc-om/Permanent/Leigh.Zhang/core/flux"

for yyyy in {2025..2025}
do

for mm in {06..06}
do
   echo "starting $source_dir/flux_${yyyy}${mm}.tar"
   htar -tvf $source_dir/flux_${yyyy}${mm}.tar > tempfile.${yyyy}${mm}
   grep mem0$ens_member tempfile.${yyyy}${mm} | cut -c 68-92 > filelist.${yyyy}${mm}
   htar -xvf $source_dir/flux_${yyyy}${mm}.tar -L filelist.${yyyy}${mm}
#   rm tempfile.${yyyy}${mm} filelist.${yyyy}${mm}

done

done
