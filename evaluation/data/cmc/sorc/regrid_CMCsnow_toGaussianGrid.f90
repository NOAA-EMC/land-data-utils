program data_regrid_gaussian

  use netcdf
  implicit none

  double precision      :: sec_since

  character*256         :: source_path
  character*256         :: destination_path
  
  character*256         :: source_filename
  character*256         :: destination_filename
  character*256         :: weights_filename

  character*19     :: current_date  ! current date
  character*19     :: since_date = "1970-01-01 00:00:00"
  character*4      :: yyyy, year
  character*2      :: mm, dd
  character*10     :: idate
  character*10     :: fv3_grid 

  integer, parameter           :: source_lats           = 706
  integer, parameter           :: source_lons           = 706
  integer :: destination_locs, weight_locs, destination_lats, destination_lons

  integer :: latloc, lonloc, iwt, ndays
  integer :: i, j, itotal, id, iyear, im, io, ierr
  integer :: offset_ss

  integer, dimension(3)                       :: start, count
  real   , dimension(source_lons,source_lats) :: snowd
  real   , dimension(source_lons,source_lats) :: source_input

  integer,  allocatable :: source_lookup(:), destination_lookup(:)
  real*8 ,  allocatable :: weights(:)
  real   ,  allocatable :: cmcSnowDepth(:)
  real   ,  allocatable :: cmcSnowDepthWT(:)

  real   ,  allocatable :: xcb(:)    ! destination lon
  real   ,  allocatable :: ycb(:)    ! destination lat

! Define 2D array 
  real   ,  allocatable :: cmcSnow(:,:)
  real   ,  allocatable :: lon2d(:,:)
  real   ,  allocatable :: lat2d(:,:)
! Define 1D lat and lon
  real   , allocatable :: lon1d(:)
  real   , allocatable :: lat1d(:)

  integer :: ncid, dimid, varid, status   ! netcdf identifiers
  integer :: dim_id_i, dim_id_j,dim_id_t  ! netcdf dimension identifiers

  logical :: file_exists

  real*4  fillVal;
  fillVal = -9999.0

  offset_ss=0

  namelist/regrid_cmc_nml/ destination_locs, weight_locs, destination_lons, destination_lats, fv3_grid, source_path, weights_filename, destination_path

! read namelist

  inquire(file='regrid_cmc.nml', exist=file_exists)

  if (.not. file_exists) then
      print *, 'namelistfile does not exist, exiting'
      stop 10
  endif

  open (action='read', file='regrid_cmc.nml', iostat=ierr, newunit=io)
  read (nml=regrid_cmc_nml, iostat=ierr, unit=io)
  close (io)

  allocate(source_lookup(weight_locs))
  allocate(destination_lookup(weight_locs))
  allocate(weights(weight_locs))
  allocate(cmcSnowDepth(destination_locs))
  allocate(cmcSnowDepthWT(destination_locs))
  allocate(xcb(destination_locs))
  allocate(ycb(destination_locs))

  allocate(cmcSnow(destination_lons, destination_lats))
  allocate(lon2d(destination_lons, destination_lats))
  allocate(lat2d(destination_lons, destination_lats))

  allocate(lon1d(destination_lons))
  allocate(lat1d(destination_lats))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read weights file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  status = nf90_open(trim(weights_filename), NF90_NOWRITE, ncid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inq_varid(ncid, "col", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , source_lookup)
    if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_inq_varid(ncid, "row", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , destination_lookup)
    if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_inq_varid(ncid, "S", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , weights)
    if (status /= nf90_noerr) call handle_err(status)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read lattitude and longitude data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   status = nf90_inq_varid(ncid, "yc_b", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , ycb)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "xc_b", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_get_var(ncid, varid , xcb)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_close(ncid)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read year and numbers of days for a given year
  open(20, file='year_days.txt', status='old')
  read(20,*) year, ndays
  close(20)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read source data in ascii format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  source_filename=trim(source_path)//'cmc_sdepth_dly_'//year//'_v01.2.txt'
  open(10, file=trim(source_filename), status='old')
  do id = 1, ndays 

  read(10,*) idate ! date in format YYYYMMDDHH where HH is always 00 
  destination_filename=trim(destination_path)//trim(fv3_grid)//'_cmc_'//idate//'.nc'
  
  yyyy=idate(1:4)
  mm=idate(5:6)
  dd=idate(7:8)

  write(*,*) yyyy, mm, dd
  current_date=yyyy//'-'//mm//'-'//dd//' 00:00:00'
  call calc_sec_since(since_date, current_date, offset_ss, sec_since)

  do i=1,706
  read(10,'(706(1x,f5.1))') (snowd(i,j),j=1,706) !cm 
  enddo
!!!!!!!!!!!! convert cm into mm !!!!!!!!!!!!!!!!!!!!!!!!!!
  do i=1,706
   do j=1,706
   source_input(i,j)=snowd(i,j)*10.0
   end do
  end do
  !write(*,*) source_input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regrid the data
  cmcSnowDepth = 0.0
  cmcSnowDepthWT = 0.0

  do iwt = 1, weight_locs

    latloc = source_lookup(iwt)/source_lons + 1
    lonloc = source_lookup(iwt) - (latloc-1)*source_lons
       
    if(source_input(lonloc,latloc) >= 0.0 ) then
    cmcSnowDepth(destination_lookup(iwt)) = cmcSnowDepth(destination_lookup(iwt)) + weights(iwt) * source_input(lonloc,latloc)
    cmcSnowDepthWT(destination_lookup(iwt)) = cmcSnowDepthWT(destination_lookup(iwt)) + weights(iwt)
    end if

  end do
  
  do iwt = 1, destination_locs

    if(cmcSnowDepthWT(iwt) > 0.0 ) then
    cmcSnowDepth(iwt) = cmcSnowDepth(iwt)/cmcSnowDepthWT(iwt)
    else
    cmcSnowDepth(iwt) = -9999.0
    end if

  end do

  where(ycb < 0.0) cmcSnowDepth=-9999.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert 1D vector into 2D array for cmc snow depth
! Keep lats and lons as 1D to save storage space 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  itotal=0
  do j = 1, destination_lats
  do i = 1, destination_lons
    itotal=itotal+1
    cmcSnow(i,j) = cmcSnowDepth(itotal)
    lat2d(i,j) = ycb(itotal)
    lon2d(i,j) = xcb(itotal)
  end do
  end do
  
  do j = 1,destination_lats
    lat1d(j) = lat2d(1,j)
  end do

  do i = 1, destination_lons
    lon1d(i) = lon2d(i,1)
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the destination file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  status = nf90_create(destination_filename, NF90_CLOBBER, ncid)
    if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.

  status = nf90_def_dim(ncid, "lat"   , destination_lats , dim_id_j)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "lon"   , destination_lons , dim_id_i)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "time"   , NF90_UNLIMITED , dim_id_t)
    if (status /= nf90_noerr) call handle_err(status)
  
! Define variables in the file.

  status = nf90_def_var(ncid, "time", NF90_DOUBLE, dim_id_t, varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "time")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "units", "seconds since "//since_date)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "lat", NF90_FLOAT, (/dim_id_j/), varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "latitude")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "units", "degrees_north")
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "lon", NF90_FLOAT, (/dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "longitude")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "unit", "degrees_east")
    if (status /= nf90_noerr) call handle_err(status) 

  status = nf90_def_var(ncid, "snowDepth", NF90_FLOAT, (/dim_id_i, dim_id_j, dim_id_t/), varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "long_name", "cmc snow depth")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "units", "mm")
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "_FillValue", fillVal)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call handle_err(status)

! Write variables in the file.
  
  status = nf90_inq_varid(ncid, "lat", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lat1d, start = (/1/), count = (/destination_lats/))
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "lon", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lon1d, start = (/1/), count = (/destination_lons/))
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "time", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , sec_since)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_inq_varid(ncid, "snowDepth", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , cmcSnow, start = (/1,1,1/), count = (/destination_lons, destination_lats,1/))
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, varid, "_FillValue", fillVal)
    if (status /= nf90_noerr) call handle_err(status)
 
  status = nf90_close(ncid)
  
  enddo                ! end of ndays loop

end program

  subroutine handle_err(status)
    use netcdf
    integer, intent ( in) :: status
 
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine handle_err

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate time in seconds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
