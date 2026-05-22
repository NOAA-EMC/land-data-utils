program create_core_hourly_from_monthly

! create hourly T,q,ps,wind speed,LWdn,precip from monthly three-hourly files
!   T,q,ps,wind are instantaneous analysis so 00Z = 00Z, 03Z = 03Z, etc
!     01Z is linearly interpolated from 00Z and 03Z
!   LWdn, precip are average over the three hours so
!     00Z = 01Z = 02Z = 3hr value from 00Z file
! sbatch --partition=u1-service -A fv3-cpu -n 1 --mem=3g --time=3:00:00 -q batch --wrap "./create_core_hourly_from_monthly.exe"

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
real  , allocatable :: temperature_monthly(:,:,:)
real  , allocatable :: temperature_daily(:,:,:)
real  , allocatable :: specific_humidity_monthly(:,:,:)
real  , allocatable :: specific_humidity_daily(:,:,:)
real  , allocatable :: surface_pressure_monthly(:,:,:)
real  , allocatable :: surface_pressure_daily(:,:,:)
real  , allocatable :: wind_speed_monthly(:,:,:)
real  , allocatable :: wind_speed_daily(:,:,:)
real  , allocatable :: downward_longwave_monthly(:,:,:)
real  , allocatable :: downward_longwave_daily(:,:,:)
real  , allocatable :: precipitation_monthly(:,:,:)
real  , allocatable :: precipitation_daily(:,:,:)


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
  allocate(temperature_monthly(number_lons,number_lats,number_times+1))  ! add one to time for the beginning of next month
  allocate(specific_humidity_monthly(number_lons,number_lats,number_times+1))  ! add one to time for the beginning of next month
  allocate(surface_pressure_monthly(number_lons,number_lats,number_times+1))  ! add one to time for the beginning of next month
  allocate(wind_speed_monthly(number_lons,number_lats,number_times+1))  ! add one to time for the beginning of next month
  allocate(downward_longwave_monthly(number_lons,number_lats,number_times+1))  ! add one to time for the beginning of next month
  allocate(precipitation_monthly(number_lons,number_lats,number_times+1))  ! add one to time for the beginning of next month
  
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

  status = nf90_inq_varid(ncid, "temperature", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire temperature")
  status = nf90_get_var(ncid, varid , temperature_monthly(:,:,1:number_times))
    if (status /= nf90_noerr) call handle_err(status, "problem getting temperature")

  status = nf90_inq_varid(ncid, "specific_humidity", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire specific_humidity")
  status = nf90_get_var(ncid, varid , specific_humidity_monthly(:,:,1:number_times))
    if (status /= nf90_noerr) call handle_err(status, "problem getting specific_humidity")

  status = nf90_inq_varid(ncid, "surface_pressure", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire surface_pressure")
  status = nf90_get_var(ncid, varid , surface_pressure_monthly(:,:,1:number_times))
    if (status /= nf90_noerr) call handle_err(status, "problem getting surface_pressure")

  status = nf90_inq_varid(ncid, "wind_speed", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire wind_speed")
  status = nf90_get_var(ncid, varid , wind_speed_monthly(:,:,1:number_times))
    if (status /= nf90_noerr) call handle_err(status, "problem getting wind_speed")

  status = nf90_inq_varid(ncid, "downward_longwave", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire downward_longwave")
  status = nf90_get_var(ncid, varid , downward_longwave_monthly(:,:,1:number_times))
    if (status /= nf90_noerr) call handle_err(status, "problem getting downward_longwave")

  status = nf90_inq_varid(ncid, "precipitation", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire precipitation")
  status = nf90_get_var(ncid, varid , precipitation_monthly(:,:,1:number_times))
    if (status /= nf90_noerr) call handle_err(status, "problem getting precipitation")

  status = nf90_close(ncid)
  
  if(imm == 12) then
    next_mm   = 1
    next_yyyy = iyyyy + 1
  else
    next_mm   = imm + 1
    next_yyyy = iyyyy
  end if

  write(filename,'(a13,i4,a1,i2.2,a3)') "CORe_forcing_", next_yyyy, "-", next_mm, ".nc"
  filename = trim(monthly_path)//trim(filename)

  inquire(file=filename, exist=file_exists)

! this file will not exist for the last partial month

  if(file_exists) then

    print *, "Reading: ", trim(filename)

    status = nf90_open(filename, NF90_NOWRITE, ncid)
      if (status /= nf90_noerr) call handle_err(status,"problem opening file")
  
    status = nf90_inq_varid(ncid, "temperature", varid)
      if(status /= nf90_noerr) call handle_err(status, "problem inquire temperature")
    status = nf90_get_var(ncid, varid , temperature_monthly(:,:,number_times+1), start=(/1,1,1/), count=(/number_lons,number_lats,1/))
      if (status /= nf90_noerr) call handle_err(status, "problem getting temperature")

    status = nf90_inq_varid(ncid, "specific_humidity", varid)
      if(status /= nf90_noerr) call handle_err(status, "problem inquire specific_humidity")
    status = nf90_get_var(ncid, varid , specific_humidity_monthly(:,:,number_times+1), start=(/1,1,1/), count=(/number_lons,number_lats,1/))
      if (status /= nf90_noerr) call handle_err(status, "problem getting specific_humidity")

    status = nf90_inq_varid(ncid, "surface_pressure", varid)
      if(status /= nf90_noerr) call handle_err(status, "problem inquire surface_pressure")
    status = nf90_get_var(ncid, varid , surface_pressure_monthly(:,:,number_times+1), start=(/1,1,1/), count=(/number_lons,number_lats,1/))
      if (status /= nf90_noerr) call handle_err(status, "problem getting surface_pressure")

    status = nf90_inq_varid(ncid, "wind_speed", varid)
      if(status /= nf90_noerr) call handle_err(status, "problem inquire wind_speed")
    status = nf90_get_var(ncid, varid , wind_speed_monthly(:,:,number_times+1), start=(/1,1,1/), count=(/number_lons,number_lats,1/))
      if (status /= nf90_noerr) call handle_err(status, "problem getting wind_speed")

    status = nf90_inq_varid(ncid, "downward_longwave", varid)
      if(status /= nf90_noerr) call handle_err(status, "problem inquire downward_longwave")
    status = nf90_get_var(ncid, varid , downward_longwave_monthly(:,:,number_times+1), start=(/1,1,1/), count=(/number_lons,number_lats,1/))
      if (status /= nf90_noerr) call handle_err(status, "problem getting downward_longwave")

    status = nf90_inq_varid(ncid, "precipitation", varid)
      if(status /= nf90_noerr) call handle_err(status, "problem inquire precipitation")
    status = nf90_get_var(ncid, varid , precipitation_monthly(:,:,number_times+1), start=(/1,1,1/), count=(/number_lons,number_lats,1/))
      if (status /= nf90_noerr) call handle_err(status, "problem getting precipitation")

    status = nf90_close(ncid)

  end if  ! next month exists

! now loop through the days to create hourly values

  allocate(             time_daily(                        24))
  allocate(      temperature_daily(number_lons,number_lats,24))
  allocate(specific_humidity_daily(number_lons,number_lats,24))
  allocate( surface_pressure_daily(number_lons,number_lats,24))
  allocate(       wind_speed_daily(number_lons,number_lats,24))
  allocate(downward_longwave_daily(number_lons,number_lats,24))
  allocate(    precipitation_daily(number_lons,number_lats,24))

  number_days = number_times / 8

  do iday  = 1, number_days
  do ihour = 1, 8

    prev_time = (iday-1)*8+ihour
    next_time = (iday-1)*8+ihour+1

    time_period_begin = time_monthly(prev_time)
    time_daily((ihour-1)*3+1) = time_period_begin
    time_daily((ihour-1)*3+2) = time_period_begin + 3600.d0
    time_daily((ihour-1)*3+3) = time_period_begin + 7200.d0

    temperature_daily(:,:,(ihour-1)*3+1) = temperature_monthly(:,:,prev_time)
    temperature_daily(:,:,(ihour-1)*3+2) = (2.0 * temperature_monthly(:,:,prev_time) + 1.0 * temperature_monthly(:,:,next_time))/3.0
    temperature_daily(:,:,(ihour-1)*3+3) = (1.0 * temperature_monthly(:,:,prev_time) + 2.0 * temperature_monthly(:,:,next_time))/3.0

    specific_humidity_daily(:,:,(ihour-1)*3+1) = specific_humidity_monthly(:,:,prev_time)
    specific_humidity_daily(:,:,(ihour-1)*3+2) = (2.0 * specific_humidity_monthly(:,:,prev_time) + 1.0 * specific_humidity_monthly(:,:,next_time))/3.0
    specific_humidity_daily(:,:,(ihour-1)*3+3) = (1.0 * specific_humidity_monthly(:,:,prev_time) + 2.0 * specific_humidity_monthly(:,:,next_time))/3.0

    surface_pressure_daily(:,:,(ihour-1)*3+1) = surface_pressure_monthly(:,:,prev_time)
    surface_pressure_daily(:,:,(ihour-1)*3+2) = (2.0 * surface_pressure_monthly(:,:,prev_time) + 1.0 * surface_pressure_monthly(:,:,next_time))/3.0
    surface_pressure_daily(:,:,(ihour-1)*3+3) = (1.0 * surface_pressure_monthly(:,:,prev_time) + 2.0 * surface_pressure_monthly(:,:,next_time))/3.0

    wind_speed_daily(:,:,(ihour-1)*3+1) = wind_speed_monthly(:,:,prev_time)
    wind_speed_daily(:,:,(ihour-1)*3+2) = (2.0 * wind_speed_monthly(:,:,prev_time) + 1.0 * wind_speed_monthly(:,:,next_time))/3.0
    wind_speed_daily(:,:,(ihour-1)*3+3) = (1.0 * wind_speed_monthly(:,:,prev_time) + 2.0 * wind_speed_monthly(:,:,next_time))/3.0

    downward_longwave_daily(:,:,(ihour-1)*3+1) = downward_longwave_monthly(:,:,prev_time)
    downward_longwave_daily(:,:,(ihour-1)*3+2) = downward_longwave_monthly(:,:,prev_time)
    downward_longwave_daily(:,:,(ihour-1)*3+3) = downward_longwave_monthly(:,:,prev_time)

    precipitation_daily(:,:,(ihour-1)*3+1) = precipitation_monthly(:,:,prev_time)
    precipitation_daily(:,:,(ihour-1)*3+2) = precipitation_monthly(:,:,prev_time)
    precipitation_daily(:,:,(ihour-1)*3+3) = precipitation_monthly(:,:,prev_time)

  end do

! create the output filename and netcdf file (overwrite old)

  write(filename,'(a13,i4,a1,i2.2,a1,i2.2,a3)') "CORe_forcing_", iyyyy, "-", imm, "-", iday, ".nc"
  filename = trim(daily_path)//trim(filename)

  status = nf90_open(filename, NF90_WRITE, ncid)
    if (status /= nf90_noerr) call handle_err(status,"problem opening output file")

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
  
  status = nf90_inq_varid(ncid, "temperature", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire temperature")
  status = nf90_put_var(ncid, varid , temperature_daily)
    if (status /= nf90_noerr) call handle_err(status, "problem putting temperature")
  
  status = nf90_inq_varid(ncid, "specific_humidity", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire specific_humidity")
  status = nf90_put_var(ncid, varid , specific_humidity_daily)
    if (status /= nf90_noerr) call handle_err(status, "problem putting specific_humidity")
  
  status = nf90_inq_varid(ncid, "surface_pressure", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire surface_pressure")
  status = nf90_put_var(ncid, varid , surface_pressure_daily)
    if (status /= nf90_noerr) call handle_err(status, "problem putting surface_pressure")
  
  status = nf90_inq_varid(ncid, "wind_speed", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire wind_speed")
  status = nf90_put_var(ncid, varid , wind_speed_daily)
    if (status /= nf90_noerr) call handle_err(status, "problem putting wind_speed")
  
  status = nf90_inq_varid(ncid, "precipitation", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire precipitation")
  status = nf90_put_var(ncid, varid , precipitation_daily)
    if (status /= nf90_noerr) call handle_err(status, "problem putting precipitation")
  
  status = nf90_inq_varid(ncid, "downward_longwave", varid)
    if(status /= nf90_noerr) call handle_err(status, "problem inquire downward_longwave")
  status = nf90_put_var(ncid, varid , downward_longwave_daily)
    if (status /= nf90_noerr) call handle_err(status, "problem putting downward_longwave")
  
  status = nf90_close(ncid)
  
  end do

  deallocate(               latitude)
  deallocate(              longitude)
  deallocate(             time_daily)
  deallocate(      temperature_daily)
  deallocate(specific_humidity_daily)
  deallocate( surface_pressure_daily)
  deallocate(       wind_speed_daily)
  deallocate(downward_longwave_daily)
  deallocate(    precipitation_daily)
  deallocate(             time_monthly)
  deallocate(      temperature_monthly)
  deallocate(specific_humidity_monthly)
  deallocate( surface_pressure_monthly)
  deallocate(       wind_speed_monthly)
  deallocate(downward_longwave_monthly)
  deallocate(    precipitation_monthly)

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

subroutine date_from_since(since_date, sec_since, current_date)

! create date string from since_date date string and number of seconds since
! offset_ss can be used to adjust for time zone

implicit none
character*19 :: since_date, current_date  ! format: yyyy-mm-dd hh:nn:ss
double precision      :: sec_since
integer      :: count_sec, count_sav
integer      :: since_yyyy, since_mm, since_dd, since_hh, since_nn, since_ss
integer      :: current_yyyy, current_mm, current_dd, current_hh, current_nn, current_ss
logical      :: leap_year = .false.
integer      :: num_leap = 0
integer      :: iyyyy, imm, limit
integer, dimension(12), parameter :: days_in_mm = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  current_date = "xxxx-xx-xx xx:xx:xx"

  limit = huge(1)
  if(sec_since > limit .or. sec_since < 0) stop "not capable of dealing with sec_since"

  read(since_date( 1: 4),  '(i4)') since_yyyy
  read(since_date( 6: 7),  '(i2)') since_mm
  read(since_date( 9:10),  '(i2)') since_dd
  read(since_date(12:13),  '(i2)') since_hh
  read(since_date(15:16),  '(i2)') since_nn
  read(since_date(18:19),  '(i2)') since_ss
  
! reset to beginning of since_year

  count_sec = 0
  if(mod(since_yyyy,4) == 0) leap_year = .true.
  
  do imm = 1,since_mm-1
    count_sec = count_sec - days_in_mm(imm)*86400
  end do

  if(leap_year .and. since_mm > 2) count_sec = count_sec - 86400
  count_sec = count_sec - (since_dd - 1) * 86400
  count_sec = count_sec - (since_hh) * 3600
  count_sec = count_sec - (since_nn) * 60
  count_sec = count_sec - (since_ss)
  
! not worrying about the complexity of non-recent leap years 

! find the year

  current_yyyy = since_yyyy

  do while (count_sec <= sec_since)
    count_sav = count_sec
    if(mod(current_yyyy,4) == 0) then
      count_sec = count_sec + 366*86400
    else
      count_sec = count_sec + 365*86400
    end if
    current_yyyy = current_yyyy + 1
  end do
  
  count_sec = count_sav
  current_yyyy = current_yyyy - 1
  
! find the month

  current_mm = 0!since_mm

  leap_year = .false.
  if(mod(current_yyyy,4) == 0) leap_year = .true.

  do while (count_sec <= sec_since)
    count_sav = count_sec
    current_mm = current_mm + 1
    if(leap_year .and. current_mm == 2) then
      count_sec = count_sec + 29*86400
    else
      count_sec = count_sec + days_in_mm(current_mm)*86400
    end if
  end do
  
  count_sec = count_sav
  
! find the day

  current_dd = 0!since_dd

  do while (count_sec <= sec_since)
    current_dd = current_dd + 1
    count_sav = count_sec
    count_sec = count_sec + 86400
  end do
  
  count_sec = count_sav
  
! find the hour

  current_hh = -1!since_hh

  do while (count_sec <= sec_since)
    current_hh = current_hh + 1
    count_sav = count_sec
    count_sec = count_sec + 3600
  end do
  
  count_sec = count_sav
  
! find the minute

  current_nn = -1!since_nn

  do while (count_sec <= sec_since)
    current_nn = current_nn + 1
    count_sav = count_sec
    count_sec = count_sec + 60
  end do
  
  count_sec = count_sav
  
  current_ss = sec_since - count_sec

  write(current_date( 1: 4),  '(i4.4)') current_yyyy
  write(current_date( 6: 7),  '(i2.2)') current_mm
  write(current_date( 9:10),  '(i2.2)') current_dd
  write(current_date(12:13),  '(i2.2)') current_hh
  write(current_date(15:16),  '(i2.2)') current_nn
  write(current_date(18:19),  '(i2.2)') current_ss

end subroutine date_from_since
  
subroutine calc_sec_since(since_date, current_date, offset_ss, sec_since)

! calculate number of seconds between since_date and current_date

double precision      :: sec_since
character*19 :: since_date, current_date  ! format: yyyy-mm-dd hh:nn:ss
integer      :: offset_ss
integer      :: since_yyyy, since_mm, since_dd, since_hh, since_nn, since_ss
integer      :: current_yyyy, current_mm, current_dd, current_hh, current_nn, current_ss
logical      :: leap_year = .false.
integer      :: iyyyy, imm
integer, dimension(12), parameter :: days_in_mm = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  sec_since = 0

  read(since_date( 1: 4),  '(i4)') since_yyyy
  read(since_date( 6: 7),  '(i2)') since_mm
  read(since_date( 9:10),  '(i2)') since_dd
  read(since_date(12:13),  '(i2)') since_hh
  read(since_date(15:16),  '(i2)') since_nn
  read(since_date(18:19),  '(i2)') since_ss

  read(current_date( 1: 4),  '(i4)') current_yyyy
  read(current_date( 6: 7),  '(i2)') current_mm
  read(current_date( 9:10),  '(i2)') current_dd
  read(current_date(12:13),  '(i2)') current_hh
  read(current_date(15:16),  '(i2)') current_nn
  read(current_date(18:19),  '(i2)') current_ss

! not worrying about the complexity of non-recent leap years 

! calculate number of seconds in all years  
  do iyyyy = since_yyyy, current_yyyy
    if(mod(iyyyy,4) == 0) then
      sec_since = sec_since + 366*86400
    else
      sec_since = sec_since + 365*86400
    end if
  end do
  
! remove seconds from since_year 
  if(mod(since_yyyy,4) == 0) leap_year = .true.
  
  do imm = 1,since_mm-1
    sec_since = sec_since - days_in_mm(imm)*86400
  end do

  if(leap_year .and. since_mm > 2) sec_since = sec_since - 86400
  
  sec_since = sec_since - (since_dd - 1) * 86400
  
  sec_since = sec_since - (since_hh) * 3600

  sec_since = sec_since - (since_nn) * 60

  sec_since = sec_since - (since_ss)
  
! remove seconds in current_year 
  leap_year = .false.
  if(mod(current_yyyy,4) == 0) leap_year = .true.
  
  do imm = current_mm+1, 12
    sec_since = sec_since - days_in_mm(imm)*86400
  end do
  if(leap_year .and. current_mm < 3) sec_since = sec_since - 86400
  
  sec_since = sec_since - (days_in_mm(current_mm) - current_dd) * 86400
  
  sec_since = sec_since - (23 - current_hh) * 3600

  sec_since = sec_since - (59 - current_nn) * 60

  sec_since = sec_since - (60 - current_ss)
  
  sec_since = sec_since + offset_ss
  
end subroutine calc_sec_since

subroutine calc_cosine_zenith(now_date,latitude,longitude,cosine_zenith)

implicit none

character*19        :: now_date    ! UTC time
real                :: latitude
real                :: longitude
real                :: cosine_zenith
logical             :: print_debug = .false.

double precision    :: sec_since
double precision    :: julian_day, julian_century, geom_mean_long_sun, geom_mean_anom_sun
double precision    :: sun_eq_of_ctr, sun_true_long, sun_app_long, obliq_corr
double precision    :: vary, eccent_earth_orbit,mean_obliq_ecliptic
double precision    :: sun_declin, eq_of_time, true_solar_time, hour_angle
double precision    :: solar_zenith_angle, solar_elevation_angle
double precision    :: approx_atmo_refraction, solar_elevation_angle_corrected
double precision    :: solar_zenith_angle_corrected
double precision    :: cos_solar_zenith_angle, sin_solar_elevation_angle
double precision    :: cos_solar_zenith_angle_corrected, sin_solar_elevation_angle_corrected
double precision    :: hour, minute, second
double precision    :: degrad 

degrad = 4.d0 * datan(1.d0) / 180.d0

call calc_sec_since("1970-01-01 00:00:00", now_date, 0, sec_since)

read(now_date(12:13), *) hour
read(now_date(15:16), *) minute
read(now_date(18:19), *) second

! According to Astronomical Algorithms (Meeus), Julian Day is based on year -4712 at noon UTC
! JD of 2000 12Z should be 2451545.0
! The equations below are taken from https://gml.noaa.gov/grad/solcalc/calcdetails.html

julian_day = sec_since/86400.0 + 2440587.5  ! this adjustment assumes 1970 reference date

julian_century = (julian_day - 2451545.0) / 36525.0

geom_mean_long_sun = mod(280.46646 + julian_century * (36000.76983 + julian_century * 0.0003032),360.0)

geom_mean_anom_sun = 357.52911 + julian_century * (35999.05029 - 0.0001537 * julian_century)

sun_eq_of_ctr = sin(degrad*geom_mean_anom_sun) * &
                (1.914602 - julian_century * (0.004817 + 0.000014 * julian_century)) + &
                sin(degrad*2.0*geom_mean_anom_sun) * &
                (0.019993 - 0.000101 * julian_century) + &
                sin(degrad*3.0*geom_mean_anom_sun)*0.000289

sun_true_long = geom_mean_long_sun + sun_eq_of_ctr

sun_app_long = sun_true_long - 0.00569 - 0.00478 * sin(degrad*(125.04-1934.136*julian_century))

mean_obliq_ecliptic = 23.0 + &
                       (26.0 + ( 21.448 - julian_century * &
                       (46.815 + julian_century * (0.00059 - julian_century * 0.001813)))/60)/60

obliq_corr = mean_obliq_ecliptic + 0.00256*cos(degrad*(125.04 - 1934.136 * julian_century))

sun_declin = asin(sin(degrad*obliq_corr)*sin(degrad*sun_app_long)) / degrad

vary = tan(degrad*obliq_corr/2.0)*tan(degrad*obliq_corr/2.0)

eccent_earth_orbit = 0.016708634 - julian_century * (0.000042037 + 0.0000001267 * julian_century)

eq_of_time = 4.d0 * (1.d0/degrad) * ( vary *  sin(2.d0*degrad*geom_mean_long_sun) - &
             2.d0 * eccent_earth_orbit * sin(degrad*geom_mean_anom_sun) + &
             4.d0 * eccent_earth_orbit * vary * sin(degrad*geom_mean_anom_sun) * &
               cos(2.d0*degrad*geom_mean_long_sun) - &
             0.5d0 * vary * vary * sin(4.d0*degrad*geom_mean_long_sun) - &
             1.25d0 * eccent_earth_orbit * eccent_earth_orbit * sin(2.d0*degrad*geom_mean_anom_sun) )

true_solar_time = modulo(hour*60.0 + minute + second/60.0 + eq_of_time + 4.0 * longitude, 1440.d0)

if(true_solar_time < 0 ) then  ! using modulo above, this shouldn't happen
  hour_angle = true_solar_time/4.0 + 180
else
  hour_angle = true_solar_time/4.0 - 180
end if

solar_zenith_angle = (1.d0/degrad)*(acos(sin(degrad*latitude)*sin(degrad*sun_declin) + &
                                  cos(degrad*latitude)*cos(degrad*sun_declin)*cos(degrad*hour_angle)))

cos_solar_zenith_angle = cos(degrad*solar_zenith_angle)

solar_elevation_angle = 90.0 - solar_zenith_angle

sin_solar_elevation_angle = sin(degrad*solar_elevation_angle)

if(solar_elevation_angle > 85.0) then

  approx_atmo_refraction = 0.0
  
elseif(solar_elevation_angle > 5.0) then

  approx_atmo_refraction = 58.1 / tan(degrad*solar_elevation_angle) - &
             0.07 / (tan(degrad*solar_elevation_angle))**3 + &
             0.000086 / (tan(degrad*solar_elevation_angle))**5

elseif(solar_elevation_angle > -0.575) then

  approx_atmo_refraction = 1735.0 + solar_elevation_angle * &
                              ( -518.2 + solar_elevation_angle * &
                              (  103.4 + solar_elevation_angle * &
                              ( -12.79 + solar_elevation_angle * 0.711)))
  
else

  approx_atmo_refraction = -20.772 / tan(degrad*solar_elevation_angle)
  
end if

approx_atmo_refraction = approx_atmo_refraction / 3600.0

solar_elevation_angle_corrected = solar_elevation_angle + approx_atmo_refraction

solar_zenith_angle_corrected = solar_zenith_angle - approx_atmo_refraction

cos_solar_zenith_angle_corrected = cos(degrad*solar_zenith_angle_corrected)

sin_solar_elevation_angle_corrected = sin(degrad*solar_elevation_angle_corrected)

cosine_zenith = real(cos_solar_zenith_angle_corrected)

if(print_debug) then

print *, "                     julian_day:", julian_day
print *, "                 julian_century:", julian_century
print *, "             geom_mean_long_sun:", geom_mean_long_sun
print *, "             geom_mean_anom_sun:", geom_mean_anom_sun
print *, "                  sun_eq_of_ctr:", sun_eq_of_ctr
print *, "                  sun_true_long:", sun_true_long
print *, "                   sun_app_long:", sun_app_long
print *, "            mean_obliq_ecliptic:", mean_obliq_ecliptic
print *, "                     obliq_corr:", obliq_corr
print *, "                     sun_declin:", sun_declin
print *, "                           vary:", vary
print *, "             eccent_earth_orbit:", eccent_earth_orbit
print *, "                     eq_of_time:", eq_of_time
print *, "                true_solar_time:", true_solar_time
print *, "                     hour_angle:", hour_angle
print *, "             solar_zenith_angle:", solar_zenith_angle
print *, "         cos_solar_zenith_angle:", cos_solar_zenith_angle
print *, "          solar_elevation_angle:", solar_elevation_angle
print *, "      sin_solar_elevation_angle:", sin_solar_elevation_angle
print *, "         approx_atmo_refraction:", approx_atmo_refraction
print *, "   solar_zenith_angle_corrected:", solar_zenith_angle_corrected
print *, "solar_elevation_angle_corrected:", solar_elevation_angle_corrected
print *, "          cos_zenith_angle_corr:", cos_solar_zenith_angle_corrected
print *, "            sin_elev_angle_corr:", sin_solar_elevation_angle_corrected

end if

end subroutine calc_cosine_zenith

