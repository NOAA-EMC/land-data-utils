
begin

grid_string        = keyword_values("regrid_parameter_assignment",  "grid_string","string")
output_path        = keyword_values("regrid_parameter_assignment",  "output_path","string")
source_static_file = keyword_values("regrid_parameter_assignment",  "source_static_file","string")
weights_filename   = keyword_values("regrid_parameter_assignment",  "weights_filename","string")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; create filename and overwrite if present
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

outname = output_path+"ufs-land_"+grid_string+"_static_fields.nc"
  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
outfile = addfile(outname,"c")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set the vector dimension size and create variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cornerfile = addfile(output_path+"ufs-land_"+grid_string+"_corners.nc","r")

vectordims = getfiledimsizes(cornerfile)
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
deep_soil_temperature = latitude
  deep_soil_temperature@units = "K"
  deep_soil_temperature@long_name = "Deep Soil Temperature"
max_snow_albedo   = latitude
  max_snow_albedo@units = "-"
  max_snow_albedo@long_name = "Maximum Albedo Snow-Covered Land"
gvf_monthly = new((/12,maxdim/),float)
  gvf_monthly!0 = "time"
  gvf_monthly!1 = "location"
  gvf_monthly@units = "%"
  gvf_monthly@long_name = "Monthly Green Vegetation"
visible_black_sky_albedo = gvf_monthly
  visible_black_sky_albedo@units = "-"
  visible_black_sky_albedo@long_name = "Visible Black Sky Albedo"
visible_white_sky_albedo = gvf_monthly
  visible_white_sky_albedo@units = "-"
  visible_white_sky_albedo@long_name = "Visible White Sky Albedo"
near_IR_black_sky_albedo = gvf_monthly
  near_IR_black_sky_albedo@units = "-"
  near_IR_black_sky_albedo@long_name = "Near IR Black Sky Albedo"
near_IR_white_sky_albedo = gvf_monthly
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

z0_monthly = gvf_monthly
  z0_monthly@units = "-"
  z0_monthly@long_name = "Roughness Length"
  z0_monthly = 0.1
lai_monthly = gvf_monthly
  lai_monthly@units = "m2/m2"
  lai_monthly@long_name = "Leaf Area Index"
  lai_monthly = 4.0
albedo_monthly = gvf_monthly
  albedo_monthly@units = "-"
  albedo_monthly@long_name = "Total Albedo"
  albedo_monthly = 0.25

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; some regrid options that are not used in production
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Opt = True
;Opt@Debug          = True
;Opt@PrintTimings   = True

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read corner file for latitude and longitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cornerfile = addfile(output_path+"ufs-land_"+grid_string+"_corners.nc","r")

  latitude = (/ cornerfile->latitude /)
  outfile->latitude = latitude

  longitude = (/ cornerfile->longitude /)
  outfile->longitude = longitude

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; individually regrid each variable from the source static file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  infile = addfile(source_static_file,"r")

  data = infile->cube_tile
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  cube_tile = (/ toint(var_regrid) /)
  outfile->cube_tile = cube_tile
  delete(data)
  delete(var_regrid)
  delete(cube_tile)

  data = infile->cube_i
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  cube_i = (/ toint(var_regrid) /)
  outfile->cube_i = cube_i
  delete(data)
  delete(var_regrid)
  delete(cube_i)

  data = infile->cube_j
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  cube_j = (/ toint(var_regrid) /)
  outfile->cube_j = cube_j
  delete(data)
  delete(var_regrid)
  delete(cube_j)

  data = infile->land_mask
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  land_mask = (/ toint(var_regrid) /)
  outfile->land_mask = land_mask
  delete(data)
  delete(var_regrid)
  delete(land_mask)

  data = infile->vegetation_category
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  vegetation_category = (/ toint(var_regrid) /)
  outfile->vegetation_category = vegetation_category
  delete(data)
  delete(var_regrid)
  delete(vegetation_category)

  data = infile->soil_category
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  soil_category = (/ toint(var_regrid) /)
  outfile->soil_category = soil_category
  delete(data)
  delete(var_regrid)
  delete(soil_category)

  data = infile->soil_color_category
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  soil_color_category = (/ toint(var_regrid) /)
  outfile->soil_color_category = soil_color_category
  delete(data)
  delete(var_regrid)
  delete(soil_color_category)

  data = infile->elevation
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  elevation = (/ var_regrid /)
  outfile->elevation = elevation
  delete(data)
  delete(var_regrid)
  delete(elevation)

  data = infile->slope_category
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  slope_category = (/ toint(var_regrid) /)
  outfile->slope_category = slope_category
  delete(data)
  delete(var_regrid)
  delete(slope_category)

  data = infile->deep_soil_temperature
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  deep_soil_temperature = (/ var_regrid /)
  outfile->deep_soil_temperature = deep_soil_temperature
  delete(data)
  delete(var_regrid)
  delete(deep_soil_temperature)

  data = infile->max_snow_albedo
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  max_snow_albedo = (/ var_regrid /)
  outfile->max_snow_albedo = max_snow_albedo
  delete(data)
  delete(var_regrid)
  delete(max_snow_albedo)

  data = infile->emissivity
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  emissivity = (/ var_regrid /)
  outfile->emissivity = emissivity
  delete(data)
  delete(var_regrid)
  delete(emissivity)

  data = infile->gvf_monthly
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  gvf_monthly = (/ var_regrid /)
  outfile->gvf_monthly = gvf_monthly
  delete(data)
  delete(var_regrid)
  delete(gvf_monthly)

  data = infile->visible_black_sky_albedo
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  visible_black_sky_albedo = (/ var_regrid /)
  outfile->visible_black_sky_albedo = visible_black_sky_albedo
  delete(data)
  delete(var_regrid)
  delete(visible_black_sky_albedo)

  data = infile->visible_white_sky_albedo
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  visible_white_sky_albedo = (/ var_regrid /)
  outfile->visible_white_sky_albedo = visible_white_sky_albedo
  delete(data)
  delete(var_regrid)
  delete(visible_white_sky_albedo)

  data = infile->near_IR_black_sky_albedo
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  near_IR_black_sky_albedo = (/ var_regrid /)
  outfile->near_IR_black_sky_albedo = near_IR_black_sky_albedo
  delete(data)
  delete(var_regrid)
  delete(near_IR_black_sky_albedo)

  data = infile->near_IR_white_sky_albedo
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  near_IR_white_sky_albedo = (/ var_regrid /)
  outfile->near_IR_white_sky_albedo = near_IR_white_sky_albedo
  delete(data)
  delete(var_regrid)
  delete(near_IR_white_sky_albedo)

  data = infile->z0_monthly
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  z0_monthly = (/ var_regrid /)
  outfile->z0_monthly = z0_monthly
  delete(data)
  delete(var_regrid)
  delete(z0_monthly)

  data = infile->lai_monthly
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  lai_monthly = (/ var_regrid /)
  outfile->lai_monthly = lai_monthly
  delete(data)
  delete(var_regrid)
  delete(lai_monthly)

  data = infile->albedo_monthly
  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  albedo_monthly = (/ var_regrid /)
  outfile->albedo_monthly = albedo_monthly
  delete(data)
  delete(var_regrid)
  delete(albedo_monthly)

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

  outfile->soil_level_thickness = soil_level_thickness
  outfile->soil_level_nodes = soil_level_nodes

  outfile@soil_class_source = "STATSGO"
  outfile@land_cover_source = "IGBP"
  outfile@isice = 15
  outfile@isurban = 13
  outfile@iswater = 17

  outfile@source_static_file = source_static_file
  outfile@weights_filename   = weights_filename

end
