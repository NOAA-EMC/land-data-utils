
begin

atm_res          = keyword_values("regrid_parameter_assignment","atm_res","string")
ocn_res          = keyword_values("regrid_parameter_assignment","ocn_res","string")
grid_version     = keyword_values("regrid_parameter_assignment","grid_version","string")
grid_path        = keyword_values("regrid_parameter_assignment","grid_path","string")
weights_filename = keyword_values("regrid_parameter_assignment","weights_filename","string")
variable2regrid  = keyword_values("regrid_parameter_assignment","variable2regrid","string")
output_filename  = keyword_values("regrid_parameter_assignment","output_filename","string")

remap_file = addfile(weights_filename,"r")

source_lookup = remap_file->col
destination_lookup = remap_file->row
source_dims = remap_file->src_grid_dims
destination_dims = remap_file->dst_grid_dims
destination_lon = remap_file->yc_b
destination_lat = remap_file->xc_b

latitude  = destination_lon(0::destination_dims(0))
 latitude@units = "degrees_north"
longitude = destination_lat(0:destination_dims(0)-1)
 longitude@units = "degrees_east"

tile_size = toint(sqrt(source_dims(0)/6))
tile_area = tile_size*tile_size

if(variable2regrid.eq."vegetation_greenness") then
  var2regrid = new((/tile_size,tile_size,6,12/),float)
else
  var2regrid = new((/tile_size,tile_size,6/),float)
end if

do itile = 1, 6
  tile_filename = grid_path+grid_version+"/"+atm_res+"/sfc/"+atm_res+"."+ocn_res+"."+variable2regrid+".tile"+itile+".nc"
  tile_file = addfile(tile_filename,"r")
  print("accessing file: "+tile_filename)
  if(variable2regrid.eq."vegetation_greenness") then
    do imonth = 0, 11
      var2regrid(:,:,itile-1,imonth) = tile_file->$variable2regrid$(imonth,:,:)
    end do
  else
    var2regrid(:,:,itile-1) = tile_file->$variable2regrid$(0,:,:)
  end if
end do

var_mapped = new((/destination_dims(1),destination_dims(0)/),float)
  var_mapped!0 = "latitude"
  var_mapped!1 = "longitude"
  var_mapped&latitude = latitude
  var_mapped&longitude = longitude

do location = 0, dimsizes(source_lookup)-1

;print(source_lookup(location))

  tile = (source_lookup(location)-1)/tile_area  ; 0 - 5
  this_tile_loc = source_lookup(location) - tile_area*tile
  source_y = (this_tile_loc-1)/tile_size  ; 0 - tile_size-1
  source_x = this_tile_loc - source_y*tile_size - 1

  destination_y = (destination_lookup(location)-1)/destination_dims(0) ; start at 0
  destination_x = destination_lookup(location) - destination_y*destination_dims(0) - 1

;print("sy,sx,t,dy,dx: "+source_y+" "+source_x+" "+tile+" "+destination_y+" "+destination_x)

  if(variable2regrid.eq."vegetation_greenness") then
    var_mapped(destination_y,destination_x) = (/ dim_max(var2regrid(source_y,source_x,tile,:)) /)  
  else
    var_mapped(destination_y,destination_x) = (/ var2regrid(source_y,source_x,tile) /)
  end if

end do    


outname = output_filename
 system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
 outfile = addfile(outname,"c")

; outfile->latitude         = latitude
; outfile->longitude        = longitude
if(variable2regrid.eq."vegetation_greenness") then
 outfile->vegetation_fraction = var_mapped  
else
 outfile->$variable2regrid$ = var_mapped  
end if

outfile@atm_res          = atm_res
outfile@ocn_res          = ocn_res
outfile@grid_version     = grid_version
outfile@grid_path        = grid_path
outfile@weights_filename = weights_filename


end


