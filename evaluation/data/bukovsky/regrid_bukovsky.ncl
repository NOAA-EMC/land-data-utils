load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

     data_source_file = keyword_values("regrid_parameter_assignment","data_source_file","string")
destination_directory = keyword_values("regrid_parameter_assignment","destination_directory","string")
     weights_filename = keyword_values("regrid_parameter_assignment","weights_filename","string")
          grid_string = keyword_values("regrid_parameter_assignment","grid_string","string")

infile = addfile(data_source_file,"r")

bukovsky = (/ infile->region /)
region_names = infile->region_names
  
Opt                = True

var_regrid = ESMF_regrid_with_weights(bukovsky,weights_filename,Opt)
  var_regrid!0 = "location"
  delete(var_regrid@lat1d)
  delete(var_regrid@lon1d)

  outname = destination_directory+"/BUKOVSKY."+grid_string+".nc"
   system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname,"c")
  outfile->regions = var_regrid
  outfile->region_names = region_names
  
end

