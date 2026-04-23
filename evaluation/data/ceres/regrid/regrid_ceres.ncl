load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

data_source_directory = keyword_values("regrid_parameter_assignment","data_source_directory","string")
destination_directory = keyword_values("regrid_parameter_assignment","destination_directory","string")
     weights_filename = keyword_values("regrid_parameter_assignment","weights_filename","string")
      vector_filename = keyword_values("regrid_parameter_assignment","vector_filename","string")
          grid_string = keyword_values("regrid_parameter_assignment","grid_string","string")
           yyyy_begin = keyword_values("regrid_parameter_assignment","yyyy_begin","integer")
             yyyy_end = keyword_values("regrid_parameter_assignment","yyyy_end","integer")

if(variable_set.eq."radiation") then
  vars2process = (/"ini_sfc_sw_down_all_daily","ini_sfc_sw_up_all_daily","ini_sfc_lw_down_all_daily","ini_sfc_lw_up_all_daily"/)
  output_varnames = (/"shortwave_down","shortwave_up","longwave_down","longwave_up" /)
  filename_start = data_source_directory+"/CERES_radiation_"
  outname_start  = destination_directory+"/CERES_radiation_"
end if
if(variable_set.eq."surface") then
  vars2process = (/"aux_ocean_daily","aux_snow_daily"/)
  output_varnames = (/"water_percent","snow_percent" /)
  filename_start = data_source_directory+"/CERES_surface_"
  outname_start  = destination_directory+"/CERES_surface_"
end if
if(variable_set.eq."clouds") then
  vars2process = (/"cldarea_high_daily","cldarea_low_daily","cldarea_mid_high_daily","cldarea_mid_low_daily","cldarea_total_daily"/)
  output_varnames = (/"cloud_area_high","cloud_area_low","cloud_area_midhigh","cloud_area_midlow","cloud_area_total" /)
  filename_start = data_source_directory+"/CERES_clouds_"
  outname_start  = destination_directory+"/CERES_clouds_"
end if

latlon_file = addfile(vector_filename, "r")
  latlon_lat = latlon_file->latitude
  latlon_lon = latlon_file->longitude

do iyyyy = yyyy_begin, yyyy_end

  filename = filename_start+iyyyy+".nc"
      infile = addfile(filename,"r")

  ceres_time = infile->time

;;;;;;;;;;;;;;;;;;;;
  outname  = outname_start+iyyyy+"."+grid_string+".nc"

  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname,"c")
  outfile->time = ceres_time
  outfile->ceres_regrid_latitude = latlon_lat
  outfile->ceres_regrid_longitude = latlon_lon

do ivar = 0, dimsizes(vars2process)-1

  var2process = vars2process(ivar)
  output_varname = output_varnames(ivar)

  print("Starting: "+output_varname)
  
  ceres = infile->$var2process$
  
  Opt                = True
  ceres_regrid = ESMF_regrid_with_weights(ceres,weights_filename,Opt)

  ceres_regrid!0 = "time"
  ceres_regrid!1 = "location"
  delete(ceres_regrid@lat1d)
  delete(ceres_regrid@lon1d)
  ceres_regrid&time = ceres_time
  ceres_regrid@long_name = var2process
  
  outfile->$output_varname$ = ceres_regrid

  delete(ceres)
  delete(ceres_regrid)

end do

  delete(ceres_time)

end do
  
end

