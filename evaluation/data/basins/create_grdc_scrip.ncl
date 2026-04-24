
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

data_source_file = keyword_values("regrid_parameter_assignment","data_source_file","string")
 data_scrip_file = keyword_values("regrid_parameter_assignment","data_scrip_file" ,"string")

;************************************************
; get grdc fields
;************************************************

  grdc_file = addfile(data_source_file, "r")
  
  grdc_var = (/ grdc_file->Band1 /)
;  grdc_mask = where(ismissing(grdc_var),0,1)   ; if mask added, then resulting mask will find a source for all destination

  
  grdc_lat = grdc_file->lat
;  grdc_lat@units = "degrees_north"
  grdc_lon = grdc_file->lon
;  grdc_lon@units = "degrees_east"

; Create src grid
;----------------------------------------------------------------------
  Opt                = True
  Opt@ForceOverwrite = True
  Opt@PrintTimings   = True
  Opt@Title          = "grdc Grid"
;  Opt@GridMask       = grdc_mask

  rectilinear_to_SCRIP(data_scrip_file,grdc_lat,grdc_lon,Opt)

end

