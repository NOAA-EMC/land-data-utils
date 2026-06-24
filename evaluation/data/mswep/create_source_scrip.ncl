
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  data_source_file = keyword_values("regrid_parameter_assignment","data_source_file"  ,"string")
   data_scrip_file = keyword_values("regrid_parameter_assignment","data_scrip_file"   ,"string")

;************************************************
; get mswep fields, use mswep missing as mask for land grids
;************************************************

  mswep_file = addfile(data_source_file, "r")
  mswep_var = mswep_file->precipitation(0,:,:)
  delete(mswep_var@_FillValue)
  mswep_mask = where(mswep_var.ge.0,1,0)   ; if mask added, then resulting mask will find a source for all destination

  mswep_lat = mswep_file->lat
;  mswep_lat@units = "degrees_north"
  mswep_lon = mswep_file->lon
;  mswep_lon@units = "degrees_east"

; Create src grid
;----------------------------------------------------------------------
  Opt                = True
  Opt@ForceOverwrite = True
  Opt@PrintTimings   = True
  Opt@Title          = "mswep Grid"
  Opt@GridMask       = mswep_mask

  rectilinear_to_SCRIP(data_scrip_file,mswep_lat,mswep_lon,Opt)

end

