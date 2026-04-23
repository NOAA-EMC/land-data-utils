
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

data_source_file = keyword_values("regrid_parameter_assignment","data_source_file","string")
 data_scrip_file = keyword_values("regrid_parameter_assignment","data_scrip_file" ,"string")

;************************************************
; get bukovsky fields
;************************************************

  bukovsky_file = addfile(data_source_file, "r")
  
  bukovsky_lat = bukovsky_file->lat
;  bukovsky_lat@units = "degrees_north"
  bukovsky_lon = bukovsky_file->lon
;  bukovsky_lon@units = "degrees_east"

; Create src grid
;----------------------------------------------------------------------
  Opt                = True
  Opt@ForceOverwrite = True
  Opt@PrintTimings   = True
  Opt@Title          = "Bukovsky Grid"

  rectilinear_to_SCRIP(data_scrip_file,bukovsky_lat,bukovsky_lon,Opt)

end

