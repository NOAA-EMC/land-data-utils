load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

         elevation_filename = keyword_values("regrid_parameter_assignment",         "elevation_filename","string")
precip_interpolation_method = keyword_values("regrid_parameter_assignment","precip_interpolation_method","string")
      interpolation_method1 = keyword_values("regrid_parameter_assignment",      "interpolation_method1","string")
      interpolation_method2 = keyword_values("regrid_parameter_assignment",      "interpolation_method2","string")
   weights_method1_filename = keyword_values("regrid_parameter_assignment",   "weights_method1_filename","string")
   weights_method2_filename = keyword_values("regrid_parameter_assignment",   "weights_method2_filename","string")
           datm_source_path = keyword_values("regrid_parameter_assignment",           "datm_source_path","string")
            output_preamble = keyword_values("regrid_parameter_assignment",            "output_preamble","string")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; some regrid options that are not used in production
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Opt = True
;Opt@Debug          = True
;Opt@PrintTimings   = True

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read elevation adjustment file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

do_elevation_adjustment = True
elevation_difference_limit = 2000.0 ; an arbitrary limit on elevation difference
if(do_elevation_adjustment) then
  infile = addfile(elevation_filename,"r")
  elevation_difference = infile->elevation_difference
  elevation_difference = where(elevation_difference .gt.      elevation_difference_limit,      elevation_difference_limit, elevation_difference)
  elevation_difference = where(elevation_difference .lt. -1.0*elevation_difference_limit, -1.0*elevation_difference_limit, elevation_difference)
end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the vector size from the weights file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

infile = addfile(weights_method1_filename,"r")
numlocs = dimsizes(infile->frac_b)

do iyyyy = yyyy_beg, yyyy_end
do imm = mm_beg, mm_end

  numdays = days_in_month(iyyyy,imm)

  time_start = cd_inv_calendar(iyyyy,imm,1,0,0,0,"seconds since 1970-01-01 00:00:00", 0)

  time = new(8*numdays, double)
    time@units = "seconds since 1970-01-01 00:00:00"
    time!0 = "time"

  do ihh = 0, 8*numdays-1
    time(ihh) = time_start + ihh * todouble(3600.0*3)
  end do

  datestring = ""+iyyyy+"-"+tostring_with_format(imm,"%2.2d")

  infile = addfile(datm_source_path+"/CORe_forcing_"+datestring+".nc","r")
  print("Starting: "+datestring)

  solar_radiation = new((/8*numdays,numlocs/), float)
    solar_radiation@units = "W/m2"
    solar_radiation!0 = "time"
    solar_radiation!1 = "location"
    solar_radiation&time = time
    solar_radiation@interpolation_weights = weights_method1_filename
  longwave_radiation = solar_radiation
  temperature = solar_radiation
    temperature@units = "K"
  specific_humidity = solar_radiation
    specific_humidity@units = "kg/kg"
  speed = solar_radiation
    speed@units = "m/s"
  surface_pressure = solar_radiation
    surface_pressure@units = "Pa"
  if(precip_interpolation_method .eq. "all" .or. precip_interpolation_method .eq. "method1") then
    precipitation_method1_varname = "precipitation_"+interpolation_method1
    precipitation_method1 = solar_radiation
    precipitation_method1@units = "mm/s"
    precipitation_method1@interpolation_weights = weights_method1_filename
  end if
  if(precip_interpolation_method .eq. "all" .or. precip_interpolation_method .eq. "method2") then
    precipitation_method2_varname = "precipitation_"+interpolation_method2
    precipitation_method2 = solar_radiation
    precipitation_method2@units = "mm/s"
    precipitation_method2@interpolation_weights = weights_method2_filename
  end if
  
;;;;;;;;;;;;;;;;;;;;

  data = infile->downward_solar

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  solar_radiation = (/ var_regrid /)
  
;;;;;;;;;;;;;;;;;;;;

  data = infile->downward_longwave

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  longwave_radiation = (/ var_regrid /)

;;;;;;;;;;;;;;;;;;;;

  data = infile->temperature

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)

  if(do_elevation_adjustment) then
    do itime = 0, 23
      var_regrid(itime,:) = var_regrid(itime,:) + elevation_difference * 0.0065
    end do
  end if
  
  temperature = (/ var_regrid /)

;;;;;;;;;;;;;;;;;;;;

  data = infile->specific_humidity

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  specific_humidity = (/ var_regrid /)

;;;;;;;;;;;;;;;;;;;;

  data = infile->wind_speed

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  speed = (/ var_regrid /)

;;;;;;;;;;;;;;;;;;;;

  data = infile->surface_pressure

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  surface_pressure = (/ var_regrid /)

;;;;;;;;;;;;;;;;;;;;

  data = infile->precipitation

  if(precip_interpolation_method .eq. "all" .or. precip_interpolation_method .eq. "method1") then

    var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
    precipitation_method1 = (/ var_regrid /)

  end if
  
  if(precip_interpolation_method .eq. "all" .or. precip_interpolation_method .eq. "method2") then

    var_regrid = ESMF_regrid_with_weights(data,weights_method2_filename,Opt)
  
    precipitation_method2 = (/ var_regrid /)

  end if
  
;;;;;;;;;;;;;;;;;;;;

  outname = output_preamble+"_datm_"+datestring+".nc"
   system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname,"c")
  outfile->solar_radiation = solar_radiation
  outfile->longwave_radiation = longwave_radiation
  outfile->temperature = temperature
  outfile->specific_humidity = specific_humidity
  outfile->wind_speed = speed
  outfile->surface_pressure = surface_pressure
  if(precip_interpolation_method .eq. "all" .or. precip_interpolation_method .eq. "method1") then

    outfile->$precipitation_method1_varname$ = precipitation_method1
    delete(precipitation_method1)

  end if
  
  if(precip_interpolation_method .eq. "all" .or. precip_interpolation_method .eq. "method2") then

    outfile->$precipitation_method2_varname$ = precipitation_method2
    delete(precipitation_method2)

  end if
  
;  system("ncks -h -O -4 -L 1 "+outname+" "+outname)
  
  delete(data)
  delete(var_regrid)
  delete(time)
  delete(solar_radiation)
  delete(longwave_radiation)
  delete(temperature)
  delete(specific_humidity)
  delete(speed)
  delete(surface_pressure)

end do
end do

end

