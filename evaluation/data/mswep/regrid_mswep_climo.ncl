load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  data_source_directory = keyword_values("regrid_parameter_assignment","data_source_directory","string")
  destination_directory = keyword_values("regrid_parameter_assignment","destination_directory","string")
       weights_filename = keyword_values("regrid_parameter_assignment","weights_filename","string")
            grid_string = keyword_values("regrid_parameter_assignment","grid_string","string")

vars2process = (/"precipitation"/)

outname  = destination_directory+"/MSWEP.v3."+grid_string+".climatology.2015-2024.nc"
  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
outfile = addfile(outname,"c")
  
do ivar = 0, dimsizes(vars2process)-1

  var2process = vars2process(ivar)

  print("Starting: "+var2process)
  
  inname = data_source_directory+"/MSWEP.v3.climatology.2015-2024.nc"
  print("Reading: "+inname)
  infile = addfile(inname,"r")

  indata = infile->$var2process$
  
  Opt                = True
  mswep_regrid = ESMF_regrid_with_weights(indata,weights_filename,Opt)

  mswep_regrid!0 = "time"
  mswep_regrid!1 = "location"
  delete(mswep_regrid@lat1d)
  delete(mswep_regrid@lon1d)
  
  outfile->$var2process$ = mswep_regrid

  delete(indata)
  delete(mswep_regrid)

end do

outname  = destination_directory+"/MSWEP.v3."+grid_string+".monthly_climatology.2015-2024.nc"
  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
outfile = addfile(outname,"c")
  
do ivar = 0, dimsizes(vars2process)-1

  var2process = vars2process(ivar)

  print("Starting: "+var2process)
  
  inname = data_source_directory+"/MSWEP.v3.monthly_climatology.2015-2024.nc"
  print("Reading: "+inname)
  infile = addfile(inname,"r")

  indata = infile->$var2process$
  
  Opt                = True
  mswep_regrid = ESMF_regrid_with_weights(indata,weights_filename,Opt)

  mswep_regrid!0 = "time"
  mswep_regrid!1 = "location"
  delete(mswep_regrid@lat1d)
  delete(mswep_regrid@lon1d)
  
  outfile->$var2process$ = mswep_regrid

  delete(indata)
  delete(mswep_regrid)

end do

end
