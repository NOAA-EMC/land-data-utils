
begin

atm_res      = keyword_values("regrid_parameter_assignment",     "atm_res","string")
ocn_res      = keyword_values("regrid_parameter_assignment",     "ocn_res","string")
grid_version = keyword_values("regrid_parameter_assignment","grid_version","string")
output_path  = keyword_values("regrid_parameter_assignment", "output_path","string")
fixfile_path = keyword_values("regrid_parameter_assignment","fixfile_path","string")
grid_extent  = keyword_values("regrid_parameter_assignment", "grid_extent","string")

; create filename and overwrite if present

if(grid_extent.eq."global") then
  outname = output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_corners.nc"
elseif(grid_extent.eq."conus") then
  outname = output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_corners.nc"
elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
  outname = output_path+"ufs-land_"+atm_res+"."+grid_extent+"_corners.nc"
else
  print("problem with grid_extent in extract_corners")
  exit
end if
  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
newfile = addfile(outname,"c")

; set the maximum dimension

if(grid_extent.eq."global" .or. grid_extent.eq."conus") then
  maskfile = addfile(fixfile_path+atm_res+"/"+atm_res+"."+ocn_res+"_oro_data.tile1.nc","r")
elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
  maskfile = addfile(fixfile_path+atm_res+"_oro_data.tile7.halo0.nc","r")
end if

maskdims = getfiledimsizes(maskfile)

maxdim = maskdims(0)*maskdims(1)
if(grid_extent.eq."global" .or. grid_extent.eq."conus") then
  numtiles = 6
  maxdim = maxdim * 6
elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
  numtiles = 1
end if

latitude           = new(maxdim,double)
latitude!0         = "location"
longitude          = latitude
latitude_corners   = new((/maxdim,4/),double)
latitude_corners!0 = "location"
latitude_corners!1 = "corners"
longitude_corners  = latitude_corners

nloc = -1

do itile = 1, numtiles

print("Starting tile: "+itile)

 if(grid_extent.eq."global" .or. grid_extent.eq."conus") then
   orography_file = addfile(fixfile_path+atm_res+"/"+atm_res+"."+ocn_res+"_oro_data.tile"+itile+".nc","r")
   grid_file = addfile(fixfile_path+atm_res+"/"+atm_res+"_grid.tile"+itile+".nc","r")
 elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
   orography_file = addfile(fixfile_path+"/"+atm_res+"_oro_data.tile7.halo0.nc","r")
   grid_file = addfile(fixfile_path+atm_res+"_grid.tile7.halo0.nc","r")
 end if

 inmask = orography_file->land_frac 
 inlat  = orography_file->geolat
 inlon  = orography_file->geolon
 inx    = grid_file->x
 iny    = grid_file->y

 ndims = dimsizes(inmask)
 
 if(grid_extent.eq."conus") then

   inmask = where(inmask.gt.0 .and. \
                  inlat.gt.25  .and.  inlat.lt.53 .and. \
                  inlon.gt.235 .and.  inlon.lt.293, 1, 0)
 
 end if

 do idim0 = 0, ndims(0)-1
 do idim1 = 0, ndims(1)-1

   if(inmask(idim0,idim1).gt.0) then
     nloc = nloc + 1
     latitude(nloc)            = (/ iny   (idim0*2+1,idim1*2+1) /)
     longitude(nloc)           = (/ inx   (idim0*2+1,idim1*2+1) /)
     latitude_corners(nloc,0)  = (/ iny   (idim0*2  ,idim1*2) /)
     latitude_corners(nloc,1)  = (/ iny   (idim0*2  ,idim1*2+2) /)
     latitude_corners(nloc,2)  = (/ iny   (idim0*2+2,idim1*2+2) /)
     latitude_corners(nloc,3)  = (/ iny   (idim0*2+2,idim1*2) /)
     longitude_corners(nloc,0) = (/ inx   (idim0*2,  idim1*2) /)
     longitude_corners(nloc,1) = (/ inx   (idim0*2  ,idim1*2+2) /)
     longitude_corners(nloc,2) = (/ inx   (idim0*2+2,idim1*2+2) /)
     longitude_corners(nloc,3) = (/ inx   (idim0*2+2,idim1*2) /)
   end if
   
 end do
 end do
 
 print("number of cumulative locs (1-based): "+(nloc+1))

end do

newfile->latitude = latitude(0:nloc)
newfile->longitude = longitude(0:nloc)
newfile->latitude_corners = latitude_corners(0:nloc,:)
newfile->longitude_corners = longitude_corners(0:nloc,:)

end

