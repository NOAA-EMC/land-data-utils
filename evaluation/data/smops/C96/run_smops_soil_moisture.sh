#! /bin/sh -l

export scripts=/scratch2/NCEPDEV/land/data/evaluation/ush/
export EXEC_DIR=./

# for the whole period run
sdate=20180101
edate=20221231

# for a test run
#sdate=20180701
#edate=20180701

while [ $sdate -le $edate ]; do

\rm date_input.txt
  
echo $sdate  > date_input.txt

${EXEC_DIR}/regrid_smopsSM_cgrid.exe

sdate=`$scripts/finddate.sh $sdate d+1`

done
