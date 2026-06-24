load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

         elevation_filename = keyword_values("regrid_parameter_assignment",         "elevation_filename","string")
            static_filename = keyword_values("regrid_parameter_assignment",            "static_filename","string")
                grid_extent = keyword_values("regrid_parameter_assignment",                "grid_extent","string")
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
; get some grid information from the static file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

infile = addfile(static_filename,"r")
atm_res = infile@atm_res
ocn_res = infile@ocn_res
grid_string = infile@grid_string
fixfile_path = infile@fixfile_path
grid_extent = infile@grid_extent
is_global = infile@is_global
if(grid_extent.eq."subset") then
  subset_maxlat = infile@subset_maxlat
  subset_minlat = infile@subset_minlat
  subset_maxlon = infile@subset_maxlon
  subset_minlon = infile@subset_minlon
end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the vector size from the weights file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

infile = addfile(weights_method1_filename,"r")
numlocs = dimsizes(infile->frac_b)

do iyyyy = yyyy_beg, yyyy_end
do imm = mm_beg, mm_end

  numdays = days_in_month(iyyyy,imm)

do idd = 1, numdays
  
  time_start = cd_inv_calendar(iyyyy,imm,idd,0,0,0,"seconds since 1970-01-01 00:00:00", 0)

  time = new(24, double)
    time@units = "seconds since 1970-01-01 00:00:00"
    time!0 = "time"

  do ihh = 0, 23
    time(ihh) = time_start + ihh * todouble(3600.0)
  end do

  datestring = ""+iyyyy+"-"+tostring_with_format(imm,"%2.2d")+"-"+tostring_with_format(idd,"%2.2d")

  infile = addfile(datm_source_path+iyyyy+"/ERA5_forcing_"+datestring+".nc","r")
  print("Starting: "+datestring)

  outname = output_preamble+"_datm_"+datestring+".nc"
   system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname,"c")
  
;;;;;;;;;;;;;;;;;;;;

  solar_radiation = new((/24,numlocs/), float)
    solar_radiation@units = "W/m2"
    solar_radiation!0 = "time"
    solar_radiation!1 = "location"
    solar_radiation&time = time
    solar_radiation@interpolation_weights = weights_method1_filename

  data = infile->downward_solar

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  solar_radiation = (/ var_regrid /)

  outfile->solar_radiation = solar_radiation
  
  delete(data)
  delete(var_regrid)
  delete(solar_radiation)

;;;;;;;;;;;;;;;;;;;;

  longwave_radiation = new((/24,numlocs/), float)
    longwave_radiation@units = "W/m2"
    longwave_radiation!0 = "time"
    longwave_radiation!1 = "location"
    longwave_radiation&time = time
    longwave_radiation@interpolation_weights = weights_method1_filename

  data = infile->downward_longwave

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  longwave_radiation = (/ var_regrid /)

  outfile->longwave_radiation = longwave_radiation

  delete(data)
  delete(var_regrid)
  delete(longwave_radiation)

;;;;;;;;;;;;;;;;;;;;

  temperature = new((/24,numlocs/), float)
    temperature@units = "K"
    temperature!0 = "time"
    temperature!1 = "location"
    temperature&time = time
    temperature@interpolation_weights = weights_method1_filename

  data = infile->temperature

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)

  if(do_elevation_adjustment) then
    do itime = 0, 23
      var_regrid(itime,:) = var_regrid(itime,:) + elevation_difference * 0.0065
    end do
  end if
  
  temperature = (/ var_regrid /)

  outfile->temperature = temperature

  delete(data)
  delete(var_regrid)
  delete(temperature)

;;;;;;;;;;;;;;;;;;;;

  specific_humidity = new((/24,numlocs/), float)
    specific_humidity@units = "kg/kg"
    specific_humidity!0 = "time"
    specific_humidity!1 = "location"
    specific_humidity&time = time
    specific_humidity@interpolation_weights = weights_method1_filename

  data = infile->specific_humidity

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  specific_humidity = (/ var_regrid /)

  outfile->specific_humidity = specific_humidity

  delete(data)
  delete(var_regrid)
  delete(specific_humidity)

;;;;;;;;;;;;;;;;;;;;

  speed = new((/24,numlocs/), float)
    speed@units = "m/s"
    speed!0 = "time"
    speed!1 = "location"
    speed&time = time
    speed@interpolation_weights = weights_method1_filename

  data = infile->wind_speed

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  speed = (/ var_regrid /)

  outfile->wind_speed = speed

  delete(data)
  delete(var_regrid)
  delete(speed)

;;;;;;;;;;;;;;;;;;;;

  surface_pressure = new((/24,numlocs/), float)
    surface_pressure@units = "Pa"
    surface_pressure!0 = "time"
    surface_pressure!1 = "location"
    surface_pressure&time = time
    surface_pressure@interpolation_weights = weights_method1_filename

  data = infile->surface_pressure

  var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
  surface_pressure = (/ var_regrid /)

  outfile->surface_pressure = surface_pressure

  delete(data)
  delete(var_regrid)
  delete(surface_pressure)

;;;;;;;;;;;;;;;;;;;;

  data = infile->precipitation

  if(precip_interpolation_method .eq. "all" .or. precip_interpolation_method .eq. "method1") then

    precipitation_method1 = new((/24,numlocs/), float)
     precipitation_method1@units = "mm/s"
     precipitation_method1!0 = "time"
     precipitation_method1!1 = "location"
     precipitation_method1&time = time
     precipitation_method1@interpolation_weights = weights_method1_filename

    var_regrid = ESMF_regrid_with_weights(data,weights_method1_filename,Opt)
  
    precipitation_method1 = (/ var_regrid /)

    precipitation_method1_varname = "precipitation_"+interpolation_method1
    outfile->$precipitation_method1_varname$ = precipitation_method1

    delete(precipitation_method1)
    delete(var_regrid)

  end if
  
  if(precip_interpolation_method .eq. "all" .or. precip_interpolation_method .eq. "method2") then

    precipitation_method2 = new((/24,numlocs/), float)
     precipitation_method2@units = "mm/s"
     precipitation_method2!0 = "time"
     precipitation_method2!1 = "location"
     precipitation_method2&time = time
     precipitation_method2@interpolation_weights = weights_method2_filename

    var_regrid = ESMF_regrid_with_weights(data,weights_method2_filename,Opt)
  
    precipitation_method2 = (/ var_regrid /)

    precipitation_method2_varname = "precipitation_"+interpolation_method2
    outfile->$precipitation_method2_varname$ = precipitation_method2

    delete(precipitation_method2)
    delete(var_regrid)

  end if
  
  delete(data)

;;;;;;;;;;;;;;;;;;;;

  outfile@atm_res       = atm_res
  outfile@ocn_res       = ocn_res
  outfile@grid_string   = grid_string
  outfile@fixfile_path  = fixfile_path
  outfile@grid_extent   = grid_extent
  outfile@is_global     = is_global
  if(grid_extent.eq."subset") then
    outfile@subset_maxlat = subset_maxlat
    outfile@subset_minlat = subset_minlat
    outfile@subset_maxlon = subset_maxlon
    outfile@subset_minlon = subset_minlon
  end if
  
;  system("ncks -h -O -4 -L 1 "+outname+" "+outname)
  
  delete(time)

end do
end do
end do

end

