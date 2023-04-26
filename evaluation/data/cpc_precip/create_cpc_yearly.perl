#!/usr/bin/perl

$yyyy_beg   = 2000;
$yyyy_end   = 2021;

$datadir = "/scratch2/NCEPDEV/land/data/evaluation/cpc_precip/orig";

@nums = ("00","01","02","03","04","05","06","07","08","09", 
         "10","11","12","13","14","15","16","17","18","19", 
         "20","21","22","23","24","25","26","27","28","29", 
         "30","31","32","33","34","35","36","37","38","39",
         "40","41","42","43","44","45","46","47","48","49",
         "50","51","52","53","54","55","56","57","58","59",
         "60","61","62","63","64","65","66","67","68","69",
         "70","71","72","73","74","75","76","77","78","79",
         "80","81","82","83","84","85","86","87","88","89",
         "90","91","92","93","94","95","96","97","98","99");

for($yyyy=$yyyy_beg; $yyyy<=$yyyy_end; $yyyy++)
 {

  print("  starting $yyyy \n");
  
  print("  starting precip.$yyyy.ttl.nc \n");
  system("ncap2 -h -O -s 'where(precip<0) precip=-9999' $datadir/daily/precip.$yyyy.nc blah.nc");
  system("ncatted -h -O -a _FillValue,precip,a,f,-9999 blah.nc");
  system("ncatted -h -O -a missing_value,precip,m,f,-9999 blah.nc");
  system("ncra -h -O -y ttl -v precip blah.nc $datadir/yearly/precip.$yyyy.ttl.nc");
  system("ncatted -h -O -a valid_range,precip,d,, $datadir/yearly/precip.$yyyy.ttl.nc");
  system("ncatted -h -O -a actual_range,precip,d,, $datadir/yearly/precip.$yyyy.ttl.nc");
  system("rm blah.nc");

 }  #yyyy loop
