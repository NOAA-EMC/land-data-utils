#! /bin/sh -l

export scripts=/scratch2/NCEPDEV/land/data/evaluation/ush
export EXEC_DIR=../sorc

sdate=20240101
edate=20241231

#sdate=20150101
#edate=20150101

while [ $sdate -le $edate ]; do
   year=`echo $sdate |cut -c1-4`
   mon=`echo $sdate |cut -c5-6`
   day=`echo $sdate |cut -c7-8`

\rm date_input.txt
  
echo $year >> date_input.txt
echo $mon >> date_input.txt
echo $day >> date_input.txt

${EXEC_DIR}/regrid_snodas_to_gaussian.exe
echo $sdate

sdate=`$scripts/finddate.sh $sdate d+1`

done
