
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

      destination_file = keyword_values("regrid_parameter_assignment","destination_file"  ,"string")
destination_scrip_file = keyword_values("regrid_parameter_assignment","destination_scrip_file"   ,"string")

;************************************************
; get datm fields
;************************************************

  datm_file = addfile(destination_file, "r")

  datm_lat = datm_file->latitude
  datm_lat@units = "degrees_north"
  datm_lon = datm_file->longitude
  datm_lon@units = "degrees_east"

; Create scrip file
;----------------------------------------------------------------------
  Opt                = True
  Opt@ForceOverwrite = True
  Opt@PrintTimings   = True
  Opt@Title          = "DATM scrip file"

  rectilinear_to_SCRIP(destination_scrip_file,datm_lat,datm_lon,Opt)

end

