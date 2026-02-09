
begin

atm_res      = keyword_values("regrid_parameter_assignment",     "atm_res","string")
ocn_res      = keyword_values("regrid_parameter_assignment",     "ocn_res","string")
grid_version = keyword_values("regrid_parameter_assignment","grid_version","string")
output_path  = keyword_values("regrid_parameter_assignment", "output_path","string")
fixfile_path = keyword_values("regrid_parameter_assignment","fixfile_path","string")
grid_extent  = keyword_values("regrid_parameter_assignment", "grid_extent","string")

; create filename and overwrite if present

if(grid_extent.eq."global") then
  outname = output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_static_fields.nc"
elseif(grid_extent.eq."conus") then
  outname = output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_static_fields.nc"
elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
  outname = output_path+"ufs-land_"+atm_res+"."+grid_extent+"_static_fields.nc"
else
  print("problem with grid_extent in extract_static")
  exit
end if
  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
newfile = addfile(outname,"c")

; set the vector dimension size

if(grid_extent.eq."global") then
  infile = addfile(output_path+"ufs-land_"+atm_res+"."+ocn_res+"_"+grid_version+"_corners.nc","r")
elseif(grid_extent.eq."conus") then
  infile = addfile(output_path+"ufs-land_"+atm_res+"."+ocn_res+"."+grid_extent+"_"+grid_version+"_corners.nc","r")
elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
  infile = addfile(output_path+"ufs-land_"+atm_res+"."+grid_extent+"_corners.nc","r")
end if

vectordims = getfiledimsizes(infile)

maxdim = vectordims(0)

latitude = new(maxdim,float)
  latitude!0 = "location"
  latitude = -9999.
latitude@_FillValue = -9999.
  latitude@units = "degrees_north"
  longitude             = latitude
  longitude@units = "degrees_east"
vegetation_type       = latitude
  vegetation_type@units = "categorical"
  vegetation_type@long_name = "IGBP Vegetation Category"
soil_type             = latitude
  soil_type@units = "categorical"
  soil_type@long_name = "Soil Texture Class"
soil_color             = latitude
  soil_color@units = "categorical"
  soil_color@long_name = "Soil Color Class"
elevation             = latitude
  elevation@units = "m"
  elevation@long_name = "Elevation"
slope_type            = latitude
  slope_type@units = "categorical"
  slope_type@long_name = "Slope Category"
substrate_temperature = latitude
  substrate_temperature@units = "K"
  substrate_temperature@long_name = "Deep Soil Temperature"
maximum_snow_albedo   = latitude
  maximum_snow_albedo@units = "-"
  maximum_snow_albedo@long_name = "Maximum Albedo Snow-Covered Land"
vegetation_greenness = new((/12,maxdim/),float)
  vegetation_greenness!0 = "time"
  vegetation_greenness!1 = "location"
  vegetation_greenness@units = "%"
  vegetation_greenness@long_name = "Monthly Green Vegetation"
visible_black_sky_albedo = vegetation_greenness
  visible_black_sky_albedo@units = "-"
  visible_black_sky_albedo@long_name = "Visible Black Sky Albedo"
visible_white_sky_albedo = vegetation_greenness
  visible_white_sky_albedo@units = "-"
  visible_white_sky_albedo@long_name = "Visible White Sky Albedo"
near_IR_black_sky_albedo = vegetation_greenness
  near_IR_black_sky_albedo@units = "-"
  near_IR_black_sky_albedo@long_name = "Near IR Black Sky Albedo"
near_IR_white_sky_albedo = vegetation_greenness
  near_IR_white_sky_albedo@units = "-"
  near_IR_white_sky_albedo@long_name = "Near IR White Sky Albedo"

cube_tile = new(maxdim,integer)
  cube_tile!0 = "location"
cube_i = cube_tile
cube_j = cube_tile

slope_category = cube_i
  slope_category@units = "categorical"
  slope_category@long_name = "Slope Category"
soil_category = cube_i
  soil_category@units = "categorical"
  soil_category@long_name = "Soil Texture Class"
soil_color_category = cube_i
  soil_color_category@units = "categorical"
  soil_color_category@long_name = "Soil Color Class"
vegetation_category = cube_i
  vegetation_category@units = "categorical"
  vegetation_category@long_name = "IGBP Vegetation Category"
land_mask = cube_i
  land_mask@units = "categorical"
  land_mask@long_name = "1=land"
  land_mask = 1

emissivity = latitude
  emissivity@units = "-"
  emissivity@long_name = "Emissivity"
  emissivity = 0.95

z0_monthly = vegetation_greenness
  z0_monthly@units = "-"
  z0_monthly@long_name = "Roughness Length"
  z0_monthly = 0.1
lai_monthly = vegetation_greenness
  lai_monthly@units = "m2/m2"
  lai_monthly@long_name = "Leaf Area Index"
  lai_monthly = 4.0
albedo_monthly = vegetation_greenness
  albedo_monthly@units = "-"
  albedo_monthly@long_name = "Total Albedo"
  albedo_monthly = 0.25

total_loc = -1

if(grid_extent.eq."global" .or. grid_extent.eq."conus") then
  numtiles = 6
elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
  numtiles = 1
end if

do itile = 1, numtiles

print("Starting tile: "+itile)

 if(grid_extent.eq."global" .or. grid_extent.eq."conus") then
  maskfile = addfile(fixfile_path+atm_res+"/"+atm_res+"."+ocn_res+"_oro_data.tile"+itile+".nc","r")
  typefile = addfile(fixfile_path+atm_res+"/sfc/"+atm_res+"."+ocn_res+".vegetation_type.tile"+itile+".nc","r")
  mxalfile = addfile(fixfile_path+atm_res+"/sfc/"+atm_res+"."+ocn_res+".maximum_snow_albedo.tile"+itile+".nc","r")
  slopfile = addfile(fixfile_path+atm_res+"/sfc/"+atm_res+"."+ocn_res+".slope_type.tile"+itile+".nc","r")
  albdfile = addfile(fixfile_path+atm_res+"/sfc/"+atm_res+"."+ocn_res+".snowfree_albedo.tile"+itile+".nc","r")
  subtfile = addfile(fixfile_path+atm_res+"/sfc/"+atm_res+"."+ocn_res+".substrate_temperature.tile"+itile+".nc","r")
  grenfile = addfile(fixfile_path+atm_res+"/sfc/"+atm_res+"."+ocn_res+".vegetation_greenness.tile"+itile+".nc","r")
  soilfile = addfile(fixfile_path+atm_res+"/sfc/"+atm_res+"."+ocn_res+".soil_type.tile"+itile+".nc","r")
  colrfile = addfile(fixfile_path+atm_res+"/sfc/"+atm_res+"."+ocn_res+".soil_color.tile"+itile+".nc","r")
 elseif(grid_extent.eq."AQM" .or. grid_extent.eq."ARC") then
  maskfile = addfile(fixfile_path+"/"+atm_res+"_oro_data.tile7.halo0.nc","r")
  typefile = addfile(fixfile_path+"/sfc/"+atm_res+".vegetation_type.tile7.halo0.nc","r")
  mxalfile = addfile(fixfile_path+"/sfc/"+atm_res+".maximum_snow_albedo.tile7.halo0.nc","r")
  slopfile = addfile(fixfile_path+"/sfc/"+atm_res+".slope_type.tile7.halo0.nc","r")
  albdfile = addfile(fixfile_path+"/sfc/"+atm_res+".snowfree_albedo.tile7.halo0.nc","r")
  subtfile = addfile(fixfile_path+"/sfc/"+atm_res+".substrate_temperature.tile7.halo0.nc","r")
  grenfile = addfile(fixfile_path+"/sfc/"+atm_res+".vegetation_greenness.tile7.halo0.nc","r")
  soilfile = addfile(fixfile_path+"/sfc/"+atm_res+".soil_type.tile7.halo0.nc","r")
  colrfile = addfile(fixfile_path+"/sfc/"+atm_res+".soil_color.tile7.halo0.nc","r")
 end if

 inmask = maskfile->land_frac 
 inlat  = maskfile->geolat
 inlon  = maskfile->geolon

 ndims = dimsizes(inmask)
 
 if(grid_extent.eq."conus") then

   inmask = where(inmask.gt.0 .and. \
                  inlat.gt.25  .and.  inlat.lt.53 .and. \
                  inlon.gt.235 .and.  inlon.lt.293, 1, 0)
 
 end if

nloc = total_loc

 do idim0 = 0, ndims(0)-1
 do idim1 = 0, ndims(1)-1
 
   if(inmask(idim0,idim1).gt.0) then
     nloc = nloc + 1
     cube_tile(nloc) = itile
     cube_j(nloc)    = idim0 + 1
     cube_i(nloc)    = idim1 + 1
     latitude(nloc)               = (/ inlat (idim0,idim1) /)
     longitude(nloc)              = (/ inlon (idim0,idim1) /)
   end if
   
 end do
 end do
 
delete(inlat)
delete(inlon)

print("done with latitude and longitude, step 1 of 7")

 inelev = maskfile->orog_filt
 inslop = slopfile->slope_type(0,:,:)

nloc = total_loc

 do idim0 = 0, ndims(0)-1
 do idim1 = 0, ndims(1)-1
 
   if(inmask(idim0,idim1).gt.0) then
     nloc = nloc + 1
     elevation(nloc)              = (/ inelev(idim0,idim1) /)
     slope_type(nloc)             = (/ inslop(idim0,idim1) /)
   end if
   
 end do
 end do

delete(inelev)
delete(inslop)

print("done with elevation and slope, step 2 of 7")

 intype = typefile->vegetation_type(0,:,:)
 insoil = soilfile->soil_type(0,:,:)
 incolr = colrfile->soil_color(0,:,:)

nloc = total_loc

 do idim0 = 0, ndims(0)-1
 do idim1 = 0, ndims(1)-1
 
   if(inmask(idim0,idim1).gt.0) then
     nloc = nloc + 1
     vegetation_type(nloc)        = (/ intype(idim0,idim1) /)
     soil_type(nloc)              = (/ insoil(idim0,idim1) /)
     soil_color(nloc)             = (/ incolr(idim0,idim1) /)
   end if
   
 end do
 end do

delete(intype)
delete(insoil)
delete(incolr)

print("done with vegetation, color and soil types, step 3 of 7")

 inmxal = mxalfile->maximum_snow_albedo(0,:,:)
 insubt = subtfile->substrate_temperature(0,:,:)

nloc = total_loc

 do idim0 = 0, ndims(0)-1
 do idim1 = 0, ndims(1)-1
 
   if(inmask(idim0,idim1).gt.0) then
     nloc = nloc + 1
     substrate_temperature(nloc)  = (/ insubt(idim0,idim1) /)
     maximum_snow_albedo(nloc)    = (/ inmxal(idim0,idim1) /)
   end if
   
 end do
 end do
 
delete(insubt)
delete(inmxal)

print("done with deep temp and max snow, step 4 of 7")

do imo = 0, 11
 ingren = grenfile->vegetation_greenness(imo,:,:)
 
nloc = total_loc

 do idim0 = 0, ndims(0)-1
 do idim1 = 0, ndims(1)-1
 
   if(inmask(idim0,idim1).gt.0) then
     nloc = nloc + 1
     vegetation_greenness(imo,nloc)     = (/ ingren(idim0,idim1) /)
   end if
   
 end do
 end do
 
print("done with greenness month: "+(imo+1))

end do

delete(ingren)

print("done with greenness, step 5 of 7")

do imo = 0, 11
 inavdr = albdfile->visible_black_sky_albedo(imo,:,:)
 inavdf = albdfile->visible_white_sky_albedo(imo,:,:)

nloc = total_loc

 do idim0 = 0, ndims(0)-1
 do idim1 = 0, ndims(1)-1
 
   if(inmask(idim0,idim1).gt.0) then
     nloc = nloc + 1
     visible_black_sky_albedo(imo,nloc) = (/ inavdr(idim0,idim1) /)
     visible_white_sky_albedo(imo,nloc) = (/ inavdf(idim0,idim1) /)
   end if
   
 end do
 end do
 
print("done with visible month: "+(imo+1))

end do

delete(inavdr)
delete(inavdf)

print("done with visible albedo, step 6 of 7")

do imo = 0, 11
 inandr = albdfile->near_IR_black_sky_albedo(imo,:,:)
 inandf = albdfile->near_IR_white_sky_albedo(imo,:,:)

nloc = total_loc

 do idim0 = 0, ndims(0)-1
 do idim1 = 0, ndims(1)-1
 
   if(inmask(idim0,idim1).gt.0) then
     nloc = nloc + 1
     near_IR_black_sky_albedo(imo,nloc) = (/ inandr(idim0,idim1) /)
     near_IR_white_sky_albedo(imo,nloc) = (/ inandf(idim0,idim1) /)
   end if
   
 end do
 end do
 
print("done with NIR month: "+(imo+1))

end do

delete(inandr)
delete(inandf)

print("done with NIR albedo, step 7 of 7")

total_loc = nloc

print("number of cumulative locs (0-based): " + total_loc)

end do

; do a quick check

if(any(vegetation_type(0:nloc).le.0) .or. any(vegetation_type(0:nloc).eq.17)) then
  print("problem with veg type")
  exit
end if

if(any(soil_type(0:nloc).le.0) .or. any(soil_type(0:nloc).eq.14)) then
  print("problem with soil type")
  exit
end if

if(any(soil_color(0:nloc).le.0) .or. any(soil_color(0:nloc).ge.21)) then
  print("problem with soil color")
  exit
end if

vegetation_category = tointeger(vegetation_type)
soil_category = tointeger(soil_type)
slope_category = tointeger(slope_type)
soil_color_category = tointeger(soil_color)

newfile->cube_tile             = cube_tile(0:nloc)
newfile->cube_i                = cube_i(0:nloc)
newfile->cube_j                = cube_j(0:nloc)
newfile->latitude              = latitude(0:nloc)
newfile->longitude             = longitude(0:nloc)
newfile->land_mask             = land_mask(0:nloc)
newfile->vegetation_category   = vegetation_category(0:nloc)
newfile->soil_category         = soil_category(0:nloc)
newfile->soil_color_category   = soil_color_category(0:nloc)
newfile->elevation             = elevation(0:nloc)
newfile->slope_category        = slope_category(0:nloc)
newfile->deep_soil_temperature = substrate_temperature(0:nloc)
newfile->max_snow_albedo       = maximum_snow_albedo(0:nloc)
newfile->emissivity            = emissivity(0:nloc)

newfile->gvf_monthly              = vegetation_greenness(:,0:nloc)
newfile->visible_black_sky_albedo = visible_black_sky_albedo(:,0:nloc)
newfile->visible_white_sky_albedo = visible_white_sky_albedo(:,0:nloc)
newfile->near_IR_black_sky_albedo = near_IR_black_sky_albedo(:,0:nloc)
newfile->near_IR_white_sky_albedo = near_IR_white_sky_albedo(:,0:nloc)

newfile->z0_monthly               = z0_monthly(:,0:nloc)
newfile->lai_monthly              = lai_monthly(:,0:nloc)
newfile->albedo_monthly           = albedo_monthly(:,0:nloc)


  soil_level_thickness = new(4,float)
   soil_level_thickness@description = "soil level thicknesses"
   soil_level_thickness@units = "m"
   soil_level_thickness!0 = "soil_levels"
  soil_level_nodes = new(4,float)
   soil_level_nodes@description = "soil node depths"
   soil_level_nodes@units = "m"
   soil_level_nodes!0 = "soil_levels"

  soil_level_thickness = (/ 0.1, 0.3, 0.6, 1.0 /)
  soil_level_nodes = (/ 0.05, 0.25, 0.7, 1.5 /)

newfile->soil_level_thickness = soil_level_thickness
newfile->soil_level_nodes = soil_level_nodes

newfile@soil_class_source = "STATSGO"
newfile@land_cover_source = "IGBP"
newfile@isice = 15
newfile@isurban = 13
newfile@iswater = 17

end
