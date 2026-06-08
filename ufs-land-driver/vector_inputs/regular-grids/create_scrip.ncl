; this script produces SCRIP files used for ESMF regridding: full grid, and separate for bare/snow/veg

begin

grid_string   = keyword_values("regrid_parameter_assignment",  "grid_string","string")
output_path   = keyword_values("regrid_parameter_assignment",  "output_path","string")

    corners_filename      = output_path+"ufs-land_"+grid_string+"_corners.nc"
  pts_scrip_filename      = output_path+"ufs-land_"+grid_string+"_SCRIP.nc"

; get grid information

  pts_file  = addfile(corners_filename, "r")

  lat1d = pts_file->latitude
  lon1d = pts_file->longitude
  corner_lats = pts_file->latitude_corners
  corner_lons = pts_file->longitude_corners

; create destination grid

  system("if [ -e "+pts_scrip_filename+" ]; then rm -f "+pts_scrip_filename+ ";fi")
  pts_scrip = addfile(pts_scrip_filename,"c")

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

end

