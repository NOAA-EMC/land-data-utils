; this script produces SCRIP files used for ESMF regridding: full grid, and separate for bare/snow/veg

begin

atm_res      = keyword_values("regrid_parameter_assignment",     "atm_res","string")
ocn_res      = keyword_values("regrid_parameter_assignment",     "ocn_res","string")
grid_version = keyword_values("regrid_parameter_assignment","grid_version","string")
output_path  = keyword_values("regrid_parameter_assignment", "output_path","string")
fixfile_path = keyword_values("regrid_parameter_assignment","fixfile_path","string")
grid_extent  = keyword_values("regrid_parameter_assignment", "grid_extent","string")

if(grid_extent.eq."global") then
    corners_filename      = output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_corners.nc"
  pts_scrip_filename      = output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_SCRIP.nc"
  pts_scrip_filename_bare = output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_SCRIP_bare.nc"
  pts_scrip_filename_snow = output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_SCRIP_snow.nc"
  pts_scrip_filename_veg  = output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_SCRIP_veg.nc"
     static_filename      = output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_static_fields.nc"
elseif(grid_extent.eq."conus") then
    corners_filename      = output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_corners.nc"
  pts_scrip_filename      = output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_SCRIP.nc"
  pts_scrip_filename_bare = output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_SCRIP_bare.nc"
  pts_scrip_filename_snow = output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_SCRIP_snow.nc"
  pts_scrip_filename_veg  = output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_SCRIP_veg.nc"
     static_filename      = output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_static_fields.nc"
elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
    corners_filename      = output_path+"ufs-land_"+atm_res+"."+grid_extent+"_corners.nc"
  pts_scrip_filename      = output_path+"ufs-land_"+atm_res+"."+grid_extent+"_SCRIP.nc"
  pts_scrip_filename_bare = output_path+"ufs-land_"+atm_res+"."+grid_extent+"_SCRIP_bare.nc"
  pts_scrip_filename_snow = output_path+"ufs-land_"+atm_res+"."+grid_extent+"_SCRIP_snow.nc"
  pts_scrip_filename_veg  = output_path+"ufs-land_"+atm_res+"."+grid_extent+"_SCRIP_veg.nc"
     static_filename      = output_path+"ufs-land_"+atm_res+"."+grid_extent+"_static_fields.nc"
else
  print("problem with grid_extent in create_scrip")
  exit
end if

; get grid information

  pts_file  = addfile(corners_filename, "r")

  lat1d = pts_file->latitude
  lon1d = pts_file->longitude
  corner_lats = pts_file->latitude_corners
  corner_lons = pts_file->longitude_corners

  veg_file  = addfile(static_filename, "r")

  vegtype = veg_file->vegetation_category

; create destination grid

  system("if [ -e "+pts_scrip_filename+" ]; then rm -f "+pts_scrip_filename+ ";fi")
  pts_scrip = addfile(pts_scrip_filename,"c")

  system("if [ -e "+pts_scrip_filename_bare+" ]; then rm -f "+pts_scrip_filename_bare+ ";fi")
  pts_scrip_bare = addfile(pts_scrip_filename_bare,"c")

  system("if [ -e "+pts_scrip_filename_snow+" ]; then rm -f "+pts_scrip_filename_snow+ ";fi")
  pts_scrip_snow = addfile(pts_scrip_filename_snow,"c")

  system("if [ -e "+pts_scrip_filename_veg+" ]; then rm -f "+pts_scrip_filename_veg+ ";fi")
  pts_scrip_veg = addfile(pts_scrip_filename_veg,"c")

  grid_size = dimsizes(lat1d)
  print(grid_size)
  grid_corners = 4
  grid_rank = 1

  grid_dims = new((/grid_rank/),integer)
  grid_dims!0 = "grid_rank"
  grid_dims(0) = grid_size

  grid_center_lat = new((/grid_size/),double)
  grid_center_lat!0 = "grid_size"
  grid_center_lat@units = "degrees"
  grid_center_lat = (/ lat1d /)

  grid_center_lon = grid_center_lat
  grid_center_lon = (/ lon1d /)

  grid_corner_lat = new((/grid_size,grid_corners/),double)
  grid_corner_lat!0 = "grid_size"
  grid_corner_lat!1 = "grid_corners"
  grid_corner_lat@units = "degrees"
  grid_corner_lat = (/ corner_lats /)
  
  grid_corner_lon = grid_corner_lat
  grid_corner_lon = (/ corner_lons /)

; full grid contains all land grids so imask=1 everywhere
  
  grid_imask = new((/grid_size/),integer)
  grid_imask!0 = "grid_size"
  grid_imask = 1
  
  pts_scrip->grid_dims = grid_dims
  pts_scrip->grid_center_lat = grid_center_lat
  pts_scrip->grid_center_lon = grid_center_lon
  pts_scrip->grid_imask = grid_imask
  pts_scrip->grid_corner_lat = grid_corner_lat
  pts_scrip->grid_corner_lon = grid_corner_lon

; bare grid contains only bare grids so imask=1 only when vegtype=16
  
  grid_imask = (/ where(vegtype.eq.16,1,0) /)
  
  pts_scrip_bare->grid_dims = grid_dims
  pts_scrip_bare->grid_center_lat = grid_center_lat
  pts_scrip_bare->grid_center_lon = grid_center_lon
  pts_scrip_bare->grid_imask = grid_imask
  pts_scrip_bare->grid_corner_lat = grid_corner_lat
  pts_scrip_bare->grid_corner_lon = grid_corner_lon

; snow grid contains only snow grids so imask=1 only when vegtype=15
  
  grid_imask = (/ where(vegtype.eq.15,1,0) /)
  
  pts_scrip_snow->grid_dims = grid_dims
  pts_scrip_snow->grid_center_lat = grid_center_lat
  pts_scrip_snow->grid_center_lon = grid_center_lon
  pts_scrip_snow->grid_imask = grid_imask
  pts_scrip_snow->grid_corner_lat = grid_corner_lat
  pts_scrip_snow->grid_corner_lon = grid_corner_lon

; veg grid contains only veg grids so imask=1 only when vegtype/=15 and vegtype/=16
  
  grid_imask = (/ where(vegtype.eq.16.or.vegtype.eq.15,0,1) /)
  
  pts_scrip_veg->grid_dims = grid_dims
  pts_scrip_veg->grid_center_lat = grid_center_lat
  pts_scrip_veg->grid_center_lon = grid_center_lon
  pts_scrip_veg->grid_imask = grid_imask
  pts_scrip_veg->grid_corner_lat = grid_corner_lat
  pts_scrip_veg->grid_corner_lon = grid_corner_lon

end

