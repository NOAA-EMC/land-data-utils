load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

verbose = False
  cmc_latlon_filename = "CMC_polar_stereo_grid.nc"
  cmc_scrip_filename  = "CMC_polar_stereo_SCRIP.nc"

;************************************************
; get cmc fields
;************************************************

  cmc_file = addfile(cmc_latlon_filename, "r")
 
  cmc_lat = cmc_file->latitude_center
  cmc_lat@units = "degrees_north"
  cmc_lon = cmc_file->longitude_center
  cmc_lon@units = "degrees_east"

  cmc_lat_corner = cmc_file->latitude_corner
  cmc_lat_corner@units = "degrees_north"
  cmc_lon_corner = cmc_file->longitude_corner
  cmc_lon_corner@units = "degrees_east"

; Create src grid
;----------------------------------------------------------------------
  Opt                = True
  Opt@ForceOverwrite = True
  Opt@PrintTimings   = True
  Opt@Title          = "CMC Snow Grid"
  Opt@GridCornerLat  = cmc_lat_corner
  Opt@GridCornerLon  = cmc_lon_corner

  curvilinear_to_SCRIP(cmc_scrip_filename,cmc_lat,cmc_lon,Opt)

end
