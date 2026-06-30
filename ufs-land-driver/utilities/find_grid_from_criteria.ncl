begin

; this script will find the first location in the ufs_land_driver vector that
;  satisfies a set of criteria
;
; the result is 1-based meaning you can use that location in the ufs_land_driver namelist
;    to run a single point simulation

; set the vector static field file and mask file of the desired resolution

static_file_name = "/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/C192.mx025/ufs-land_C192.mx025_static_fields.nc"
mask_file_name   = "/scratch4/NCEPDEV/land/data/evaluation/masks/global_testing/U1D_C192.mx025_VegSoil_masks.nc"

; criteria

vegetation_type_search = True          ; True or False
vegetation2find        = 16            ; vegetation type to find

soil_type_search       = True          ; True or False
soil2find              = 7             ; soil type to find

mask_search            = True          ; True or False
region2find            = "GreatBasin"  ; region name from the mask_file

latlon_search          = True          ; True or False
latbox_max             = 50.0          ; maximum latitude of lat/lon region
latbox_min             = 30.0          ; minimum latitude of lat/lon region
lonbox_max             = 270.0         ; maximum longitude (0:360) of lat/lon region
lonbox_min             = 230.0        ; minimum longitude (0:360) of lat/lon region

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; shouldn't need to modify below this line
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

num_locations = 0
if(vegetation_type_search .or. soil_type_search .or. latlon_search) then
  infile = addfile(static_file_name,"r")
  vegetation_type = infile->vegetation_category
  soil_type = infile->soil_category
  latitude = infile->latitude
  longitude = infile->longitude
  num_locations = dimsizes(vegetation_type)
end if

if(mask_search) then
  infile = addfile(mask_file_name,"r")
  region_mask = infile->$region2find$
  if(num_locations .eq. 0) then
    num_locations = dimsizes(region_mask)
  end if
end if

if(num_locations .eq. 0) then 
  print("No search defined")
  exit
end if

found = False

do iloc = 0, num_locations-1

  if(vegetation_type_search) then
    if(vegetation_type(iloc) .eq. vegetation2find) then
      found = True
    else
      found = False
    end if
  end if

  if(soil_type_search .and. found) then
    if(soil_type(iloc) .eq. soil2find) then
      found = True
    else
      found = False
    end if
  end if

  if(mask_search .and. found) then
    if(region_mask(iloc) .eq. 1) then
      found = True
    else
      found = False
    end if
  end if

  if(latlon_search .and. found) then
    if( latitude(iloc) .ge. latbox_min .and.  latitude(iloc) .le. latbox_max .and. \
       longitude(iloc) .ge. lonbox_min .and. longitude(iloc) .le. lonbox_max) then
      found = True
    else
      found = False
    end if
  end if

  if(found) then
    print("Found first location with the criteria: ")
    if(vegetation_type_search) then
      print(" Vegetation type: "+vegetation2find)
    end if

    if(soil_type_search) then
      print(" Soil type: "+soil2find)
    end if

    if(mask_search) then
      print(" Region: "+region2find)
    end if

    if(latlon_search) then
      print(" Lat/Lon Box: ["+latbox_min+":"+latbox_max+","+lonbox_min+":"+lonbox_max+"]")
    end if

    print("at 1-based location: "+(iloc+1))
    print(" Vegetation type: "+vegetation_type(iloc))
    print(" Soil type: "+soil_type(iloc))
    print(" Region: "+region2find)
    print(" Latitude: "+latitude(iloc))
    print(" Longitude: "+longitude(iloc))
    exit
  end if

end do

print("No location with the criteria: ")

if(vegetation_type_search) then
 print("Vegetation type: "+vegetation2find)
end if

if(soil_type_search) then
  print("Soil type: "+soil2find)
end if

if(mask_search) then
  print("Region: "+region2find)
end if

if(latlon_search) then
  print(" Lat/Lon Box: ["+latbox_min+":"+latbox_max+","+lonbox_min+":"+lonbox_max+"]")
end if

end
