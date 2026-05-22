begin

lat2find = 40.0
lon2find = 75.0
var2check = "downward_solar"

daily_file = addfile("/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/CORe/original/netcdf_daily/CORe_forcing_1980-01-01.nc","r")
monthly_file = addfile("/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/CORe/original/netcdf_monthly/CORe_forcing_1980-01.nc","r")

latitude = daily_file->latitude
longitude = daily_file->longitude

ifind = ind_nearest_coord(lon2find,longitude, 0)
jfind = ind_nearest_coord(lat2find,latitude, 0)

daily = daily_file->$var2check$(:,jfind,ifind)
monthly = monthly_file->$var2check$(0:8,jfind,ifind)

print(daily_file->$var2check$(:,jfind,ifind))
print(monthly_file->$var2check$(0:8,jfind,ifind))

wks  = gsn_open_wks("png", "check_daily_interp")
plot = new(1,graphic)

time = ispan(0,24,3)
time2 = ispan(0,23,1)

res = True
res@xyLineColors  = (/"red","blue","red","blue","red","blue","black","black"/)
res@xyMarkLineModes = (/"MarkLines","Lines","Lines","Lines","Lines","Lines","Lines","Lines","Lines"/)
res@xyLineThicknessF = 2.0
res@xyMonoDashPattern = False
res@xyDashPatterns = (/0,0,1,1,2,2,0,0/)
res@xyMonoMarker = True
res@xyMarker = 16
res@xyMarkerSizeF = 0.01
res@trYMinF = min((/min(daily),min(monthly)/))
res@trYMaxF = max((/max(daily),max(monthly)/))

res2 = res
res2@xyMarkLineModes = (/"Markers","Lines","Lines","Lines","Lines","Lines","Lines","Lines","Lines"/)
res2@xyMarkerSizeF = 0.015

plot = gsn_csm_x2y2(wks,time2,time,daily,monthly,res,res2)

end
