#! /bin/sh -l

syear=2012
eyear=2020

while [ $syear -le $eyear ]; do

# check if input file exists

export input_file=year_days.txt

if [ -f $input_file ]; then
\rm $input_file
fi

b=`expr $syear % 4`
if [ $b -eq 0 ] 
then
ndays=366
else
ndays=365
fi

echo $syear > year_days.txt
echo $ndays >>year_days.txt

echo $syear, $ndays

../sorc/regrid_CMCsnow_to1Dvector.exe

((syear++))

done

