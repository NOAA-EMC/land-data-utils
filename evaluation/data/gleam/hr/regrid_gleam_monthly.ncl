load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

vars2process = (/"E","Eb","Ei","Es","Et","S","SMroot","SMsurf"/)

varlongnames = (/"Actual ET","E bare soil","E canopy interception","E snow", \
                "E transpiration","Evap stress","SM root zone","SM surface" /)

hr_filename = "/scratch2/NCEPDEV/land/data/evaluation/domains/hr/hr1_vegsoil.nc"
hr_file = addfile(hr_filename, "r")
  hr_lat = hr_file->lat
  hr_lon = hr_file->lon

do ivar = 0, 7

var2process = vars2process(ivar)
varlongname = varlongnames(ivar)

print("Starting: "+varlongname)

;filename = "/scratch2/NCEPDEV/land/data/evaluation/GLEAM/orig/"+var2process+"/v3.7b/"+var2process+"_2003-2022_GLEAM_v3.7b_MO.nc"
;outname  = "/scratch2/NCEPDEV/land/data/evaluation/GLEAM/hr/monthly/GLEAM_v3.7b_"+var2process+"_hr.nc"
;tstart = 192  ; jan 2019

filename = "/scratch2/NCEPDEV/land/data/evaluation/GLEAM/orig/"+var2process+"/v3.8a/"+var2process+"_1980-2022_GLEAM_v3.8a_MO.nc"
outname  = "/scratch2/NCEPDEV/land/data/evaluation/GLEAM/hr/monthly/GLEAM_v3.8a_"+var2process+"_hr.nc"
tstart = 468  ; jan 2019

weights = "/scratch2/NCEPDEV/land/data/evaluation/GLEAM/hr/GLEAM-hr_bilinear_wts.nc"

Opt                = True
    
  infile = addfile(filename,"r")

  in_time = infile->time(tstart:)
  gleam_time = (/ toint(in_time) /)     ; cd_calendar cant handle int64
  gleam_time@units = in_time@units
  print("gleam_time: "+gleam_time(0))

  utc_date = cd_calendar(gleam_time,0)
  yyyy = tointeger(utc_date(0,0))
  mm = tointeger(utc_date(0,1))
  dd = tointeger(utc_date(0,2))
  print("date: "+yyyy+" "+mm+" "+dd)

  time = new(dimsizes(gleam_time), double)
;   time!0 = "time"
;   time@units = "seconds since 1970-01-01 00:00:00"
  time = cd_inv_calendar(utc_date(:,0),utc_date(:,1),utc_date(:,2),utc_date(:,3),utc_date(:,4),utc_date(:,5),"seconds since 1970-01-01 00:00:00", 0)
  
;  print("time: "+time(0))
  
  gleam = infile->$var2process$(tstart:,:,:)
  print(dimsizes(gleam))
  
  gleam_regrid = ESMF_regrid_with_weights(gleam,weights,Opt)

  gleam_regrid!0 = "time"
  gleam_regrid!1 = "latitude"
  gleam_regrid!2 = "longitude"
  delete(gleam_regrid&time)
  gleam_regrid&time = time
  gleam_regrid&latitude = hr_lat
  gleam_regrid&longitude = hr_lon
  gleam_regrid@long_name = varlongname
  
;;;;;;;;;;;;;;;;;;;;

  system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
  outfile = addfile(outname,"c")
  outfile->$var2process$ = gleam_regrid
;  outfile->time = gleam_time
;  outfile->latitude = hr_lat
;  outfile->longitude = hr_lon

delete(gleam)
delete(gleam_regrid)

end do
  
end

