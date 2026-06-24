;
; script to create an elevation difference file to adjust datm forcing for elevation differences
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

static_filename           = keyword_values("regrid_parameter_assignment","static_filename"          ,"string")
elevation_source_filename = keyword_values("regrid_parameter_assignment","elevation_source_filename","string")
datm_source               = keyword_values("regrid_parameter_assignment","datm_source"              ,"string")
weights_filename          = keyword_values("regrid_parameter_assignment","weights_filename"         ,"string")
elevation_filename        = keyword_values("regrid_parameter_assignment","elevation_filename"       ,"string")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; model elevation comes from vector static file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  infile = addfile(static_filename,"r")
  vector_elevation = infile->elevation
  numlocs = dimsizes(vector_elevation)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; create elevation variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  elevation = new(numlocs, float)
    elevation@units = "m"
    elevation!0 = "locations"
  elevation_difference = elevation
    elevation_difference@description = "elevation difference between source and model grid"
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read elevation source
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  infile = addfile(elevation_source_filename,"r")

  if(datm_source .eq. "ERA5") then

    data = 1.0*infile->Z(0,:,:)  ; 
    data = data / 9.81           ; convert geopotential to elevation

  elseif(datm_source .eq. "CDAS") then

    data = 1.0*infile->elevation(::-1,:)  ; reverse latitude to conform with weights file

  elseif(datm_source .eq. "GDAS") then

    data = 1.0*infile->orog(0,::-1,:)     ; reverse latitude to conform with weights file

  elseif(datm_source .eq. "CORe") then

    data = 1.0*infile->HGT_surface(0,:,:)

  else

    print("Unknown datm_source: "+datm_source)
    exit

  end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; regrid elevation to vector
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  Opt                = True

  var_regrid = ESMF_regrid_with_weights(data,weights_filename,Opt)
  
  elevation = (/ var_regrid /)
  
  elevation_difference = elevation - vector_elevation

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; write elevation file and clobber if it exists
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  system("if [ -e "+elevation_filename+" ]; then rm -f "+elevation_filename+ ";fi")
  outfile = addfile(elevation_filename,"c")
  outfile->elevation = elevation
  outfile->elevation_difference = elevation_difference

end

