
; create SCRIP file for data atmosphere source
;   ncl create_datm_scrip.ncl 'datm_source_file="$datm_source_file"' \
;                             'datm_scrip_file="$datm_scrip_file"'

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

datm_source_file = keyword_values("regrid_parameter_assignment","datm_source_file","string")
 datm_scrip_file = keyword_values("regrid_parameter_assignment","datm_scrip_file" ,"string")

;************************************************
; get datm fields
;************************************************

  datm_file = addfile(datm_source_file, "r")

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
;  Opt@GridMask       = datm_mask  ; assuming no mask exists

  rectilinear_to_SCRIP(datm_scrip_file,datm_lat,datm_lon,Opt)

end

