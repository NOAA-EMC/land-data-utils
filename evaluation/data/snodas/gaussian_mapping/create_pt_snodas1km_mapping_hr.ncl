
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
  
  outname = "/scratch2/NCEPDEV/land/data/evaluation/SNODAS/fix_20231027/gaussianGrid_mapping/snodas1km_hr_mapping.nc"
   system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname, "c")

  corners_file = addfile("/scratch2/NCEPDEV/land/data/evaluation/domains/hr/hr_corners.nc", "r")
  hr_lat = corners_file->latitude
  hr_lon = corners_file->longitude
  lat_edges = corners_file->latitude_corners
  lon_edges = corners_file->longitude_corners
   
  ims_file = addfile("/scratch2/NCEPDEV/land/data/evaluation/SNODAS/orig/SNODAS_unmasked_20131001.nc", "r")
  snodas_lat = ims_file->lat
  snodas_lon = ims_file->lon
  snodas_lon = where(snodas_lon.lt.0,snodas_lon+360,snodas_lon)

  pt_i_snodas = new((/8192/),short)
   pt_i_snodas!0 = "lon"
  pt_j_snodas = new((/4096/),short)
  pt_j_snodas!0 = "lat"
 
  pt_i_snodas = toshort(ind_nearest_coord(snodas_lon,hr_lon,0))
  pt_j_snodas = toshort(ind_nearest_coord(snodas_lat,hr_lat,0))
    
  pt_j_snodas = pt_j_snodas + toshort(1)
  pt_i_snodas = pt_i_snodas + toshort(1)
  
  outfile->pt_i_snodas = pt_i_snodas
  outfile->pt_j_snodas = pt_j_snodas
  
end 
