#!/usr/bin/perl

# module load wgrib2 nco
# sbatch --partition=u1-service -A fv3-cpu -n 1 --time=4:00:00 -q batch --wrap "perl create_core_monthly.perl"

# CORe grib variables

# 7:884182:d=1980010100:WIND:10 m above ground:anl:ENS=+12
# 8:1086323:d=1980010100:PRES:surface:anl:ENS=+12
# 20:2999030:d=1980010100:SPFH:2 m above ground:anl:ENS=+12
# 21:3162866:d=1980010100:TMP:2 m above ground:anl:ENS=+12
# 58:5891323:d=1980010100:DLWRF:surface:0-3 hour ave fcst:ENS=+12
# 60:6224366:d=1980010100:DSWRF:surface:0-3 hour ave fcst:ENS=+12
# 114:12345817:d=1980010100:PRATE:surface:0-3 hour ave fcst:ENS=+12

# 11:1573741:d=1980010100:TMP:1 hybrid pressure level:anl:ENS=+12
# 12:1730289:d=1980010100:SPFH:1 hybrid pressure level:anl:ENS=+12
# 13:1894322:d=1980010100:UGRD:1 hybrid pressure level:anl:ENS=+12
# 14:2079711:d=1980010100:VGRD:1 hybrid pressure level:anl:ENS=+12

# Utilities
@nums = ("00","01","02","03","04","05","06","07","08","09","10",
              "11","12","13","14","15","16","17","18","19","20",
              "21","22","23","24","25","26","27","28","29","30",
              "31","32","33","34","35","36","37","38","39","40",
              "41","42","43","44","45","46","47","48","49","50",
              "51","52","53","54","55","56","57","58","59","60",
              "61","62","63","64","65","66","67","68","69","70",
              "71","72","73","74","75","76","77","78","79","80",
              "81","82","83","84","85","86","87","88","89","90",
              "91","92","93","94","95","96","97","98","99");

@ddinmm = (31,28,31,30,31,30,31,31,30,31,30,31);

$yyyy = 2026;
$mm = 5;
$final_day = 2; # this day only has 00Z time
$pathroot="/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/CORe/original";

for($dd=1; $dd<=$final_day; $dd++)
 {

 $netcdf_datestring = "$yyyy-$nums[$mm]-$nums[$dd]";
 $daily_filename = "${pathroot}/netcdf_monthly/CORe_forcing_${netcdf_datestring}.nc";
 $monthly_filename = "${pathroot}/netcdf_monthly/CORe_forcing_$yyyy-$nums[$mm].nc";

 $hh_end = 21;
 if ($dd == $final_day) {$hh_end = 1}

for($hh=0; $hh<=$hh_end; $hh=$hh+3)
 {

 $core_datestring = "$yyyy$nums[$mm]$nums[$dd]$nums[$hh]";
 
 print("Starting $core_datestring \n");
 
 $core_filename = "${pathroot}/mem012/flx_${core_datestring}_mem012.grb";

 print("Starting $core_filename \n");

# use wgrib2 to extract into state and flux files because they have different times
 system("wgrib2 $core_filename | egrep -i 'WIND|TMP:2 m above ground:anl|SPFH:2 m above ground:anl|PRES:surface' | wgrib2 -i $core_filename -netcdf state.nc");
 system("wgrib2 $core_filename | egrep -i 'DLWRF:surface:0|DSWRF:surface:0|PRATE:surface' | wgrib2 -i $core_filename -netcdf flux.nc");

# copy state vars to flux file to preserve time of state, flux fields will be forward averages
 system("cp flux.nc $daily_filename.$nums[$hh]");
 system("ncks -h -A -v WIND_10maboveground,PRES_surface,SPFH_2maboveground,TMP_2maboveground state.nc $daily_filename.$nums[$hh]");

 system("rm state.nc flux.nc");

 }

# create an intermediate daily file
 system("ncrcat -h $daily_filename.* $daily_filename");
 system("rm $daily_filename.*");
  
 }

# create the final monthly file
 system("ncrcat -h ${pathroot}/netcdf_monthly/CORe_forcing_$yyyy-$nums[$mm]-*.nc $monthly_filename");
 system("rm ${pathroot}/netcdf_monthly/CORe_forcing_$yyyy-$nums[$mm]-*.nc");

 system("ncrename -h -v DLWRF_surface,downward_longwave $monthly_filename");
 system("ncrename -h -v DSWRF_surface,downward_solar $monthly_filename");
 system("ncrename -h -v PRATE_surface,precipitation $monthly_filename");
 system("ncrename -h -v PRES_surface,surface_pressure $monthly_filename");
 system("ncrename -h -v SPFH_2maboveground,specific_humidity $monthly_filename");
 system("ncrename -h -v TMP_2maboveground,temperature $monthly_filename");
 system("ncrename -h -v WIND_10maboveground,wind_speed $monthly_filename");

 system("ncks -h -O -4 -L 1 $monthly_filename $monthly_filename");
