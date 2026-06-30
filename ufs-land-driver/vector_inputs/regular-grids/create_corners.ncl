
begin

grid_string   = keyword_values("regrid_parameter_assignment",  "grid_string","string")
output_path   = keyword_values("regrid_parameter_assignment",  "output_path","string")

grid_filename = "/scratch4/NCEPDEV/land/data/evaluation/GSCD/original/GSCD_v2.0_postprocess_0.1deg.nc"
grid_variable       = "QMEAN"
grid_latitude_name  = "latitude"
grid_longitude_name = "longitude"
grid_resolution     = 0.1

; user defined section for grid specification

grid_file = addfile(grid_filename,"r")

 latitude_in = grid_file->$grid_latitude_name$
longitude_in = grid_file->$grid_longitude_name$
     data_in = grid_file->$grid_variable$
   data_mask = .not.ismissing(data_in)

; make longitude conform to [0:360]

longitude_in = where(longitude_in.lt.0,longitude_in+360.0,longitude_in)

; create filename and overwrite if present

outname = output_path+"ufs-land_"+grid_string+"_corners.nc"
  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
newfile = addfile(outname,"c")

; set the maximum dimension

filedims = dimsizes(data_in)
maxdim = num(data_mask)

latitude           = new(maxdim,float)
latitude!0         = "location"
longitude          = latitude
latitude_corners   = new((/maxdim,4/),float)
latitude_corners!0 = "location"
latitude_corners!1 = "corners"
longitude_corners  = latitude_corners

nloc = -1

do idim0 = 0, filedims(0)-1
do idim1 = 0, filedims(1)-1

  if(data_mask(idim0,idim1)) then
    nloc = nloc + 1
    latitude(nloc)            = (/  latitude_in(idim0) /)
    longitude(nloc)           = (/ longitude_in(idim1) /)

    latitude_corners(nloc,0)  = (/  latitude_in(idim0) + grid_resolution/2 /)
    latitude_corners(nloc,1)  = (/  latitude_in(idim0) + grid_resolution/2 /)
    latitude_corners(nloc,2)  = (/  latitude_in(idim0) - grid_resolution/2 /)
    latitude_corners(nloc,3)  = (/  latitude_in(idim0) - grid_resolution/2 /)
    longitude_corners(nloc,0) = (/ longitude_in(idim1) - grid_resolution/2 /)
    longitude_corners(nloc,1) = (/ longitude_in(idim1) + grid_resolution/2 /)
    longitude_corners(nloc,2) = (/ longitude_in(idim1) + grid_resolution/2 /)
    longitude_corners(nloc,3) = (/ longitude_in(idim1) - grid_resolution/2 /)

  end if
   
end do
end do
 
 print("number of cumulative locs (1-based): "+(nloc+1))

newfile->latitude = latitude(0:nloc)
newfile->longitude = longitude(0:nloc)
newfile->latitude_corners = latitude_corners(0:nloc,:)
newfile->longitude_corners = longitude_corners(0:nloc,:)

newfile@grid_string   = grid_string
newfile@grid_filename = grid_filename

end

