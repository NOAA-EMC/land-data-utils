load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  data_source_directory = keyword_values("regrid_parameter_assignment","data_source_directory","string")
  destination_directory = keyword_values("regrid_parameter_assignment","destination_directory","string")
       weights_filename = keyword_values("regrid_parameter_assignment","weights_filename","string")
            grid_string = keyword_values("regrid_parameter_assignment","grid_string","string")
             yyyy_begin = keyword_values("regrid_parameter_assignment","yyyy_begin","integer")
               yyyy_end = keyword_values("regrid_parameter_assignment","yyyy_end","integer")

vars2process = (/"precipitation"/)

do iyyyy = yyyy_begin, yyyy_end
do imm = 1, 12

outname  = destination_directory+"/MSWEP.v3."+grid_string+"."+iyyyy+"-"+tostring_with_format(imm,"%2.2d")+".nc"
  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
outfile = addfile(outname,"c")
  
do ivar = 0, dimsizes(vars2process)-1

  var2process = vars2process(ivar)

  print("Starting: "+var2process)
  
  inname = data_source_directory+"/"+iyyyy+tostring_with_format(imm,"%2.2d")+".nc"
  print("Reading: "+inname)
  infile = addfile(inname,"r")

  indata = infile->$var2process$(0,:,:)
  
  Opt                = True
  mswep_regrid = ESMF_regrid_with_weights(indata,weights_filename,Opt)

  mswep_regrid!0 = "location"
  delete(mswep_regrid@lat1d)
  delete(mswep_regrid@lon1d)
  
  outfile->$var2process$ = mswep_regrid

  delete(indata)
  delete(mswep_regrid)

end do

end do

end do

end
