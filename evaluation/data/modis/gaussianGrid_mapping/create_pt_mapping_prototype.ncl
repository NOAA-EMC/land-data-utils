
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
  
  outname = "/scratch2/NCEPDEV/land/data/evaluation/MODIS/fix_20231027/gaussianGrid_mapping/modis_prototype_mapping.nc"
   system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname, "c")

  corners_file = addfile("/scratch2/NCEPDEV/land/data/evaluation/domains/prototype/prototype_corners.nc", "r")
  prototype_lat = corners_file->latitude
  prototype_lon = corners_file->longitude
   
  modis_file = addfile("/scratch2/NCEPDEV/land/data/evaluation/MODIS/albedo/orig/MCD43C3.A2011001.061.nc", "r")
  modis_lat = modis_file->lat
  modis_lon = modis_file->lon
  modis_lon = where(modis_lon.lt.0,modis_lon+360,modis_lon)

  pt_i_modis = new((/7200/),short)
   pt_i_modis!0 = "lon"
  pt_j_modis = new((/3600/),short)
  pt_j_modis!0 = "lat"
 
  pt_i_modis = toshort(ind_nearest_coord(modis_lon,prototype_lon,0))
  pt_j_modis = toshort(ind_nearest_coord(modis_lat,prototype_lat,0))
    
  pt_j_modis = pt_j_modis + toshort(1)
  pt_i_modis = pt_i_modis + toshort(1)
  
  outfile->pt_i_modis = pt_i_modis
  outfile->pt_j_modis = pt_j_modis
  
end

