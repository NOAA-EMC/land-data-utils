
begin

yyyy             = keyword_values("parameter_assignment","yyyy","integer")
  mm             = keyword_values("parameter_assignment","mm","integer")
  dd             = keyword_values("parameter_assignment","dd","integer")
  hh             = keyword_values("parameter_assignment","hh","integer")
  timestep       = keyword_values("parameter_assignment","timestep","integer")
ic_preamble      = keyword_values("parameter_assignment","ic_preamble","string")
datm_source_path = keyword_values("parameter_assignment","datm_source_path","string")
datm_source      = keyword_values("parameter_assignment","datm_source","string")
static_filename  = keyword_values("parameter_assignment","static_filename","string")

; access the static file for lat/lon later

static_file = addfile(static_filename,"r")

; get temperature from forcing as estimate for soil temperature

if(datm_source .eq. "ERA5") then
  temperature_filename = datm_source_path+"_datm_"+yyyy+"-"+tostring_with_format(mm,"%2.2d")+"-"+tostring_with_format(dd,"%2.2d")+".nc"
elseif(datm_source .eq. "CORe") then
  temperature_filename = datm_source_path+"_datm_"+yyyy+"-"+tostring_with_format(mm,"%2.2d")+".nc"
else
  print("Unknown datm: "+datm_source)
  exit
end if

temperature_file = addfile(temperature_filename,"r")
temperature_in = temperature_file->temperature(0,:)

; set dimensions from forcing data

num_locations = dimsizes(temperature_in)
soil_levels = 4

; figure out the time
;  input_time = time sent into the script
;  utc_date = date of input_time - timestep
;  ic_time = sec since 1970 format of utc_date
;  ic_date = filename date string

input_time = cd_inv_calendar(yyyy,mm,dd,hh,0,0,"seconds since 1970-01-01 00:00:00", 0)

ic_time = input_time
ic_time = (/ input_time - timestep*60 /)

utc_date = cd_calendar(ic_time,0)

ic_yyyy = tointeger(utc_date(0,0))
ic_mm   = tointeger(utc_date(0,1))
ic_dd   = tointeger(utc_date(0,2))
ic_hh   = tointeger(utc_date(0,3))
ic_nn   = tointeger(utc_date(0,4))
ic_date = stringtocharacter(ic_yyyy+"-"+tostring_with_format(ic_mm,"%2.2d")+"-"+tostring_with_format(ic_dd,"%2.2d")+"_"+tostring_with_format(ic_hh,"%2.2d")+":"+tostring_with_format(ic_nn,"%2.2d")+":00")
  
print("Creating cold start IC for date: "+ic_date)

ic_filename = ic_preamble+"_cold_start_"+ic_date+".nc"
  system("if [ -e "+ic_filename+" ]; then rm -f "+ic_filename+ ";fi")
ic_file = addfile(ic_filename,"c")

time = new(1, double)
 time@long_name = "time" 
 time@units = "seconds since 1970-01-01 00:00:00" 
 time!0 = "ncl_scalar"

date = new(19, character)
 date@ong_name = "UTC date"
 date!0 = "date_length"

latitude = new(num_locations,float)
 latitude@description = "latitude"
 latitude@units = "degrees_north"
 latitude!0 = "location"

longitude = new(num_locations,float)
 longitude@description = "longitude"
 longitude@units = "degrees_east"
 longitude!0 = "location"

snow_water_equivalent = new(num_locations,float)
 snow_water_equivalent@description = "snow water equivalent"
 snow_water_equivalent@units = "mm"
 snow_water_equivalent!0 = "location"

snow_depth = new(num_locations,float)
 snow_depth@description = "snow depth"
 snow_depth@units = "m"
 snow_depth!0 = "location"

canopy_water = new(num_locations,float)
 canopy_water@description = "canopy surface water"
 canopy_water@units = "mm"
 canopy_water!0 = "location"

skin_temperature = new(num_locations,float)
 skin_temperature@description = "surface skin temperature"
 skin_temperature@units = "K"
 skin_temperature!0 = "location"

soil_temperature = new((/soil_levels,num_locations/),float)
 soil_temperature@description = "soil temperature"
 soil_temperature@units = "K"
 soil_temperature!0 = "soil_levels"
 soil_temperature!1 = "location"

soil_moisture = new((/soil_levels,num_locations/),float)
 soil_moisture@description = "volumetric soil moisture"
 soil_moisture@units = "m3/m3"
 soil_moisture!0 = "soil_levels"
 soil_moisture!1 = "location"

soil_liquid = new((/soil_levels,num_locations/),float)
 soil_liquid@description = "volumetric soil liquid" ;
 soil_liquid@units = "m3/m3" ;
 soil_liquid!0 = "soil_levels"
 soil_liquid!1 = "location"

soil_level_thickness = new(soil_levels,float)
 soil_level_thickness@description = "soil level thicknesses"
 soil_level_thickness@units = "m"
 soil_level_thickness!0 = "soil_levels"

soil_level_nodes = new(soil_levels,float)
 soil_level_nodes@description = "soil node depths"
 soil_level_nodes@units = "m"
 soil_level_nodes!0 = "soil_levels"

time = (/ ic_time /)
date = (/ ic_date(0:18) /)  ; stringtochar puts a null at end

; grab some grid information from static file

latitude  = (/ static_file->latitude  /)
longitude = (/ static_file->longitude /)
vegetation_type = (/ static_file->vegetation_category /)

; snow and canopy water set to zero

snow_water_equivalent = 0.0
snow_water_equivalent = where(vegetation_type.eq.15,600.0,snow_water_equivalent)
snow_depth            = 0.0
snow_depth            = where(vegetation_type.eq.15,2.0,snow_depth)
canopy_water          = 0.0

; all temperatures set to 2m temperature, limit glacier to freezing

skin_temperature = where(vegetation_type.eq.15.and.temperature_in.gt.273.16,273.16,temperature_in)
soil_temperature(0,:) = where(vegetation_type.eq.15.and.temperature_in.gt.273.16,273.16,temperature_in)
soil_temperature(1,:) = where(vegetation_type.eq.15.and.temperature_in.gt.273.16,273.16,temperature_in)
soil_temperature(2,:) = where(vegetation_type.eq.15.and.temperature_in.gt.273.16,273.16,temperature_in)
soil_temperature(3,:) = where(vegetation_type.eq.15.and.temperature_in.gt.273.16,273.16,temperature_in)

; arbitrarily set soil moisture to 0.3

soil_moisture = 0.3
soil_liquid   = 0.3

soil_level_thickness = (/ 0.1, 0.3, 0.6, 1.0 /)
soil_level_nodes = (/ 0.05, 0.25, 0.7, 1.5 /)

ic_file->time                  = time
ic_file->date                  = date
ic_file->latitude              = latitude
ic_file->longitude             = longitude
ic_file->snow_water_equivalent = snow_water_equivalent
ic_file->snow_depth            = snow_depth
ic_file->canopy_water          = canopy_water
ic_file->skin_temperature      = skin_temperature
ic_file->soil_temperature      = soil_temperature
ic_file->soil_moisture         = soil_moisture
ic_file->soil_liquid           = soil_liquid
ic_file->soil_level_thickness  = soil_level_thickness
ic_file->soil_level_nodes      = soil_level_nodes

end

