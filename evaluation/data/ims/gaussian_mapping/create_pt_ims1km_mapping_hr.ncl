
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
  
  outname = "/scratch2/NCEPDEV/land/data/evaluation/IMS/fix_20231027/gaussianGrid_mapping/IMS1km_hr_mapping.nc"
   system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname, "c")

  corners_file = addfile("/scratch2/NCEPDEV/land/data/evaluation/domains/hr/hr_corners.nc", "r")
  prototype_lat = corners_file->latitude
  prototype_lon = corners_file->longitude
  lat_edges = corners_file->latitude_corners
  lon_edges = corners_file->longitude_corners
   
  ims_file = addfile("/scratch2/NCEPDEV/land/data/evaluation/IMS/fix_coords/IMS1km_latlon.nc", "r")
  ims_lat = ims_file->ims_lat
  ims_lon = ims_file->ims_lon
  
  pt_i_ims = new((/24576,24576/),short)
   pt_i_ims!0 = "idim"
   pt_i_ims!0 = "jdim"
  pt_j_ims = pt_i_ims
  
  do imsi = 0, 24575

   if(imsi%10.eq.0) then
     print("Starting: "+ imsi)
   end if

    lat2find = ims_lat(imsi,:)
    pt_j_ims(imsi,:) = toshort(ind_nearest_coord(lat2find,prototype_lat,0))

    lon2find = ims_lon(imsi,:)
    pt_i_ims(imsi,:) = toshort(ind_nearest_coord(lon2find,prototype_lon,0))
    
  end do
  
  pt_j_ims = pt_j_ims + toshort(1)
  pt_i_ims = pt_i_ims + toshort(1)
  
  outfile->pt_i_ims = pt_i_ims
  outfile->pt_j_ims = pt_j_ims
  
end
