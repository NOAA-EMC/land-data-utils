program create_shell_daily_files

use netcdf
implicit none

integer, parameter :: yyyy_beg = 1980
integer, parameter :: yyyy_end = 2026
integer, parameter ::   mm_beg = 1
integer, parameter ::   mm_end = 12

character*256 :: monthly_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/CORe/original/netcdf_monthly/"
character*256 ::   daily_path="/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/CORe/original/netcdf_daily/"

integer :: iyyyy,imm,itime,number_days,next_mm,next_yyyy,iday,ihour,next_time,prev_time
integer :: ncid, dimid_lat, dimid_lon, dimid_time, varid, status        ! netcdf identifiers
integer :: number_times                      ! number of times in each month
integer :: number_lats                       ! number latitudes
integer :: number_lons                       ! number longitudes
integer :: time_period_begin

character*256 :: filename
logical       :: file_exists

real*8, allocatable :: time_monthly(:)
real*8, allocatable :: time_daily(:)
real  , allocatable :: latitude(:)
real  , allocatable :: longitude(:)

! loop through yyyy and mm

yyyy_loop: do iyyyy = yyyy_beg, yyyy_end
  mm_loop: do imm   = mm_beg, mm_end

! read the current month and next month (first time only)

  write(filename,'(a13,i4,a1,i2.2,a3)') "CORe_forcing_", iyyyy, "-", imm, ".nc"
  filename = trim(monthly_path)//trim(filename)

  inquire(file=filename, exist=file_exists)
  if(.not.file_exists) cycle mm_loop

  print *, "Reading: ", trim(filename)

  status = nf90_open(filename, NF90_NOWRITE, ncid)
    if (status /= nf90_noerr) call handle_err(status,"problem opening file")

  status = nf90_inq_dimid(ncid, "time", dimid_time)
  status = nf90_inquire_dimension(ncid, dimid_time, len=number_times)
    if(status /= nf90_noerr) call handle_err(status, "problem getting time dimension length")

  status = nf90_inq_dimid(ncid, "latitude", dimid_lat)
  status = nf90_inquire_dimension(ncid, dimid_lat, len=number_lats)
    if(status /= nf90_noerr) call handle_err(status, "problem getting lat dimension length")

  status = nf90_inq_dimid(ncid, "longitude", dimid_lon)
  status = nf90_inquire_dimension(ncid, dimid_lon, len=number_lons)
    if(status /= nf90_noerr) call handle_err(status, "problem getting lon dimension length")

  allocate(time_monthly(number_times))
  allocate(latitude(number_lats))
  allocate(longitude(number_lons))

  status = nf90_inq_varid(ncid, "time", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire time")
  status = nf90_get_var(ncid, varid , time_monthly)
    if (status /= nf90_noerr) call handle_err(status, "problem getting time")

  status = nf90_inq_varid(ncid, "latitude", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire latitude")
  status = nf90_get_var(ncid, varid , latitude)
    if (status /= nf90_noerr) call handle_err(status, "problem getting latitude")

  status = nf90_inq_varid(ncid, "longitude", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire longitude")
  status = nf90_get_var(ncid, varid , longitude)
    if (status /= nf90_noerr) call handle_err(status, "problem getting longitude")

  status = nf90_close(ncid)
  
! now loop through the days to create hourly values

  allocate(time_daily(24))

  number_days = number_times / 8

  do iday  = 1, number_days
    do ihour = 1, 8

      prev_time = (iday-1)*8+ihour

      time_period_begin = time_monthly(prev_time)
      time_daily((ihour-1)*3+1) = time_period_begin
      time_daily((ihour-1)*3+2) = time_period_begin + 3600.d0
      time_daily((ihour-1)*3+3) = time_period_begin + 7200.d0

    end do

! create the output filename and netcdf file (overwrite old)

  write(filename,'(a13,i4,a1,i2.2,a1,i2.2,a3)') "CORe_forcing_", iyyyy, "-", imm, "-", iday, ".nc"
  filename = trim(daily_path)//trim(filename)

  print*, "Creating ", trim(filename)

  status = nf90_create(filename, NF90_NETCDF4, ncid)
    if (status /= nf90_noerr) call handle_err(status,"problem creating output file")

! Define dimensions in the file.

  status = nf90_def_dim(ncid, "latitude"   , number_lats    , dimid_lat)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "longitude"  , number_lons    , dimid_lon)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "time"       , NF90_UNLIMITED , dimid_time)
    if (status /= nf90_noerr) call handle_err(status)
  
! Define variables in the file.

  status = nf90_def_var(ncid, "time", NF90_DOUBLE, (/dimid_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "verification time")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "seconds since 1970-01-01 00:00:00.0 0:00")
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "latitude", NF90_FLOAT, (/dimid_lat/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "latitude")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "degrees_north")
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "longitude", NF90_FLOAT, (/dimid_lon/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "longitude")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "degrees_east")
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "temperature", NF90_FLOAT, (/dimid_lon,dimid_lat,dimid_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Temperature")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "K")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "level", "2 m above ground")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "_FillVale", 9.999e20)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "specific_humidity", NF90_FLOAT, (/dimid_lon,dimid_lat,dimid_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Specific Humidity")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "kg/kg")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "level", "2 m above ground")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "_FillValue", 9.999e20)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "surface_pressure", NF90_FLOAT, (/dimid_lon,dimid_lat,dimid_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Pressure")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "Pa")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "level", "surface")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "_FillValue", 9.999e20)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "wind_speed", NF90_FLOAT, (/dimid_lon,dimid_lat,dimid_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Wind Speed")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "m/s")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "level", "10 m above ground")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "_FillValue", 9.999e20)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "precipitation", NF90_FLOAT, (/dimid_lon,dimid_lat,dimid_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Precipitation Rate")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "kg/m^2/s")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "level", "surface")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "_FillValue", 9.999e20)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "downward_longwave", NF90_FLOAT, (/dimid_lon,dimid_lat,dimid_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Downward Long-Wave Rad. Flux")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "W/m^2")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "level", "surface")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "_FillValue", 9.999e20)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "downward_solar", NF90_FLOAT, (/dimid_lon,dimid_lat,dimid_time/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Downward Short-Wave Radiation Flux")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "units", "W/m^2")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "level", "surface")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "_FillValue", 9.999e20)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_enddef(ncid)

  status = nf90_inq_varid(ncid, "time", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire time")
  status = nf90_put_var(ncid, varid , time_daily)
    if (status /= nf90_noerr) call handle_err(status, "problem putting time")
  
  status = nf90_inq_varid(ncid, "latitude", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire latitude")
  status = nf90_put_var(ncid, varid , latitude)
    if (status /= nf90_noerr) call handle_err(status, "problem putting latitude")
  
  status = nf90_inq_varid(ncid, "longitude", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire longitude")
  status = nf90_put_var(ncid, varid , longitude)
    if (status /= nf90_noerr) call handle_err(status, "problem putting longitude")
  
  status = nf90_close(ncid)  

  end do

  deallocate(time_daily)
  deallocate(time_monthly)
  deallocate(latitude)
  deallocate(longitude)

end do   mm_loop
end do yyyy_loop

end program
  
subroutine handle_err(status, text)
  use netcdf
  integer, intent ( in) :: status
  character(len=*), optional :: text
 
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    if(present(text)) print *, text
    stop "Stopped"
  end if
end subroutine handle_err
