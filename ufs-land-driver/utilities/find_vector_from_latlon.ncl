begin

; this script will find the location in the ufs_land_driver vector given latitude and longitude
;
; the result is 1-based meaning you can use that location in the ufs_land_driver namelist
;    to run a single point simulation

; set the vector static field file of the desired resolution

static_file_name = "/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/global_0.1/ufs-land_global_0.1_static_fields.nc"

; create vectors of locations to find (lon ranges -180:180)

lats2find = (/  35.25,  42.75,  37.05, -2.55,  -2.55, -19.35, 19.05 /)
lons2find = (/ -82.95, -89.95, -77.55, 21.75, -68.85, 133.05, 80.85 /)

; print_level: 0 = only print final location
;              1 = also print search results

print_level = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; shouldn't need to modify below this line
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; longitude in static files is 0:360

lons2find = where(lons2find.lt.0,lons2find+360.0,lons2find)

infile = addfile(static_file_name,"r")
latitude = infile->latitude
longitude = infile->longitude

do iloc = 0, dimsizes(lats2find)-1

max_dist = 9999999.0
vector_lookup = -1

do icheck = 0, dimsizes(latitude)-1
  
  if(abs(lats2find(iloc) - latitude (icheck)) .lt. 1.0 .and. \
     abs(lons2find(iloc) - longitude(icheck)) .lt. 1.0) then    ; try to speed things up

    this_dist = gc_latlon(lats2find(iloc),lons2find(iloc),latitude(icheck),longitude(icheck),20,4)
    
    if(this_dist .lt. max_dist) then
      max_dist = this_dist
      vector_lookup = icheck
      if(print_level .gt. 0) then
        print("icheck: "+icheck+" max_dist: "+max_dist)
      end if
    end if

  end if

end do

if(vector_lookup .ge. 0) then
  print("found lat,lon: "+lats2find(iloc)+","+lons2find(iloc)+" at 1-based location: "+(vector_lookup+1)+" with distance: "+max_dist)
else
  print("did not find lat,lon: "+lats2find(iloc)+","+lons2find(iloc))
end if

end do

end
