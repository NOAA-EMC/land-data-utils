#!/usr/bin/perl

# sbatch -A fv3-cpu --time=2:00:00 -n 1 --wrap 'perl extract_site_forcing.perl >& extract_site_forcing.out'

@nums = ("00","01","02","03","04","05","06","07","08","09", 
         "10","11","12","13","14","15","16","17","18","19", 
         "20","21","22","23","24","25","26","27","28","29", 
         "30","31");

$site_name    = "wb_test";
@locs2extract = ( 697060, 790117, 817975, 1049556, 1242140, 1242523, 1384173 );
$static_file  = "/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/global_0.1/ufs-land_global_0.1_static_fields.nc";
$init_file    = "/scratch4/NCEPDEV/land/data/ufs-land-driver/cold_start/global_0.1/ERA5-global_0.1_cold_start_2020-12-31_23:00:00.nc";
$forcing_path = "/scratch3/NCEPDEV/stmp/Michael.Barlage/ufs-land-driver/datm/global_0.1/ERA5-global_0.1_datm_";
$subset_path  = "/scratch3/NCEPDEV/stmp/Michael.Barlage/ufs-land-driver/";

$num_sites = @locs2extract;
print("Number of sites: $num_sites \n");

$site_string = "";
for($isite=0; $isite<$num_sites; $isite++)
 {
  $site_string = "$site_string -d location,$locs2extract[$isite]";
 }

# create a cold start file

system("mkdir -p $subset_path/cold_start/${site_name}");
$outname = "$subset_path/cold_start/${site_name}/ERA5-${site_name}_cold_start_2020-12-31_23:00:00.nc";
system("ncks -h $site_string $init_file $outname");

# create a static file

system("mkdir -p $subset_path/vector_inputs/${site_name}");
$outname = "$subset_path/vector_inputs/${site_name}/ufs-land_${site_name}_static_fields.nc";
system("ncks -h $site_string $static_file $outname");

@leap_days    = (0,31,29,31,30,31,30,31,31,30,31,30,31);
@nonleap_days = (0,31,28,31,30,31,30,31,31,30,31,30,31);

system("mkdir -p $subset_path/datm/${site_name}");

for($yyyy=2021; $yyyy<=2021; $yyyy++)
 {
  @days = @nonleap_days;
  if($yyyy%4 == 0 ) {@days = @leap_days}

for($mm=1; $mm<=12; $mm++)
 {

for($dd=1; $dd<=$days[$mm]; $dd++)
 {

  $datestring = "$yyyy-$nums[$mm]-$nums[$dd]";
  $inname = "$forcing_path${datestring}.nc";
  $outname = "$subset_path/datm/${site_name}/ERA5-${site_name}_datm_${datestring}.nc";
  
  system("ncks -h $site_string $inname $outname");
  
  print("Done extracting: $datestring \n");

  }
  }
  }
