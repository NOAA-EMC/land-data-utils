
begin

atm_res       = keyword_values("regrid_parameter_assignment",      "atm_res","string")
ocn_res       = keyword_values("regrid_parameter_assignment",      "ocn_res","string")
grid_string   = keyword_values("regrid_parameter_assignment",  "grid_string","string")
output_path   = keyword_values("regrid_parameter_assignment",  "output_path","string")
fixfile_path  = keyword_values("regrid_parameter_assignment", "fixfile_path","string")
grid_extent   = keyword_values("regrid_parameter_assignment",  "grid_extent","string")
is_global     = keyword_values("regrid_parameter_assignment",    "is_global","string")
subset_maxlat = keyword_values("regrid_parameter_assignment","subset_maxlat", "float")
subset_minlat = keyword_values("regrid_parameter_assignment","subset_minlat", "float")
subset_maxlon = keyword_values("regrid_parameter_assignment","subset_maxlon", "float")
subset_minlon = keyword_values("regrid_parameter_assignment","subset_minlon", "float")

; create filename and overwrite if present

outname = output_path+"ufs-land_"+grid_string+"_corners.nc"
  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
newfile = addfile(outname,"c")

; set the maximum dimension

if(is_global.eq."True") then
  maskfile = addfile(fixfile_path+atm_res+"/"+atm_res+"."+ocn_res+"_oro_data.tile1.nc","r")
else
  maskfile = addfile(fixfile_path+atm_res+"_oro_data.tile7.halo0.nc","r")
end if

maskdims = getfiledimsizes(maskfile)

maxdim = maskdims(0)*maskdims(1)
if(is_global.eq."True") then
  numtiles = 6
  maxdim = maxdim * 6
else
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

 if(is_global.eq."True") then
   orography_file = addfile(fixfile_path+atm_res+"/"+atm_res+"."+ocn_res+"_oro_data.tile"+itile+".nc","r")
   grid_file = addfile(fixfile_path+atm_res+"/"+atm_res+"_grid.tile"+itile+".nc","r")
 else
   orography_file = addfile(fixfile_path+"/"+atm_res+"_oro_data.tile7.halo0.nc","r")
   grid_file = addfile(fixfile_path+atm_res+"_grid.tile7.halo0.nc","r")
 end if

 inmask = orography_file->land_frac 
 inlat  = orography_file->geolat
 inlon  = orography_file->geolon
 inx    = grid_file->x
 iny    = grid_file->y

 ndims = dimsizes(inmask)
 
 if(grid_extent.eq."subset") then

   inmask = where(inmask.gt.0 .and. \
                  inlat.gt.subset_minlat .and. inlat.lt.subset_maxlat .and. \
                  inlon.gt.subset_minlon .and. inlon.lt.subset_maxlon, 1, 0)
 
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

newfile@atm_res       = atm_res
newfile@ocn_res       = ocn_res
newfile@grid_string   = grid_string
newfile@fixfile_path  = fixfile_path
newfile@grid_extent   = grid_extent
newfile@is_global     = is_global
if(grid_extent.eq."subset") then
  newfile@subset_maxlat = subset_maxlat
  newfile@subset_minlat = subset_minlat
  newfile@subset_maxlon = subset_maxlon
  newfile@subset_minlon = subset_minlon
end if

end

