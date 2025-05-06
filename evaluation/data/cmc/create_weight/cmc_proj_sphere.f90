program cmc_proj

use ps_routines
use netcdf

implicit none

! CMC values
!            Lat (deg. N)      Lon (deg. E)
! LL Corner  1.665461E-01     -125.000000                   
! LR Corner  8.516040E-02      144.918700                   
! UR Corner  4.005048E-03       55.000000 
! UL Corner  8.516040E-02      -34.918730 

TYPE(proj_info)  :: proj
  
character*8        :: startloc = 'CENTER  '
!   character*8 :: startloc = 'SWCORNER'
integer, parameter :: nx = 706
integer, parameter :: ny = 706
real               :: i,j,lat,lon
integer            :: ix,iy

real, dimension(ny,nx)   :: lat_center, lon_center
real, dimension(4,ny,nx) :: lat_corner, lon_corner

integer :: ncid, dimid, varid, status        ! netcdf identifiers
integer :: dim_id_i, dim_id_j, dim_id_corner ! netcdf dimension identifiers

   proj%lat1     = 90.0 ! 1.665461E-01  ! 90.0 !
   proj%lon1     =  0.0 !-125.000000    ! 0.0 !
   proj%dx       = 23813.253175
   proj%dy       = 23813.253175
   proj%stdlon   = 280.0  !
   proj%truelat1 = 60.0
   proj%hemi     = 0.0
   proj%polei    = -999.9
   proj%polej    = -999.9
   proj%knowni   = -999.9
   proj%knownj   = -999.9
   proj%rsw      = -999.9
   proj%re_m     = 6371200.0
   proj%rebydx = proj%re_m / proj%dx

if (startloc == 'CENTER  ') then
  proj%knowni = real(nx)/2. !+ 0.5    ! + 0.5 needed to make symetric
  proj%knownj = real(ny)/2. !+ 0.5
else if (startloc == 'SWCORNER') then
  proj%knowni = 1.0
  proj%knownj = 1.0
else
  stop "no valid startloc"
end if
IF (proj%truelat1 .LT. 0.) THEN
  proj%hemi = -1.0 
ELSE
  proj%hemi = 1.0
ENDIF

if (proj%stdlon > 180.) proj%stdlon = proj%stdlon - 360.
if (proj%lon1 > 180.)     proj%lon1 = proj%lon1 - 360.

call set_ps(proj)

print *, proj%polei, proj%polej

i = 1
j = 1
call ijll_ps(i,j,proj,lat,lon)
print *, 'i=',i,'j=',j,'lat=',lat,'lon=',lon

i = 353
j = 353
call ijll_ps(i,j,proj,lat,lon)
print *, 'i=',i,'j=',j,'lat=',lat,'lon=',lon

i = 706
j = 706
call ijll_ps(i,j,proj,lat,lon)
print *, 'i=',i,'j=',j,'lat=',lat,'lon=',lon

i = 1
j = 706
call ijll_ps(i,j,proj,lat,lon)
print *, 'i=',i,'j=',j,'lat=',lat,'lon=',lon

do ix = 1, nx
do iy = 1, ny

  call ijll_ps(real(iy),real(ix),proj,lat_center(iy,ix),lon_center(iy,ix))
  call ijll_ps(iy-0.5,ix-0.5,proj,lat_corner(1,iy,ix),lon_corner(1,iy,ix))
  call ijll_ps(iy-0.5,ix+0.5,proj,lat_corner(2,iy,ix),lon_corner(2,iy,ix))
  call ijll_ps(iy+0.5,ix+0.5,proj,lat_corner(3,iy,ix),lon_corner(3,iy,ix))
  call ijll_ps(iy+0.5,ix-0.5,proj,lat_corner(4,iy,ix),lon_corner(4,iy,ix))

end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create the output filename and netcdf file (overwrite old)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  status = nf90_create("CMC_polar_stereo_grid.nc", NF90_NETCDF4, ncid)
    if (status /= nf90_noerr) call handle_err(status)

! Define dimensions in the file.

  status = nf90_def_dim(ncid, "idim"   , nx , dim_id_i)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "jdim"   , ny , dim_id_j)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "corner" ,  4 , dim_id_corner)
    if (status /= nf90_noerr) call handle_err(status)
  
! Define variables in the file.

  status = nf90_def_var(ncid, "latitude_center", NF90_FLOAT, (/dim_id_j, dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Latitude center of CMC grid")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999.0)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "longitude_center", NF90_FLOAT, (/dim_id_j, dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Longitude center of CMC grid")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999.0)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "latitude_corner", NF90_FLOAT, (/dim_id_corner, dim_id_j, dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Latitude corners of CMC grid (LL,LR,UR,UL)")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999.0)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_def_var(ncid, "longitude_corner", NF90_FLOAT, (/dim_id_corner, dim_id_j, dim_id_i/), varid)
    if (status /= nf90_noerr) call handle_err(status)

    status = nf90_put_att(ncid, varid, "long_name", "Longitude corners of CMC grid (LL,LR,UR,UL)")
      if (status /= nf90_noerr) call handle_err(status)
    status = nf90_put_att(ncid, varid, "missing_value", -9999.0)
      if (status /= nf90_noerr) call handle_err(status)

  status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "latitude_center", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lat_center)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "longitude_center", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lon_center)
    if (status /= nf90_noerr) call handle_err(status)
  
  status = nf90_inq_varid(ncid, "latitude_corner", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lat_corner)
    if (status /= nf90_noerr) call handle_err(status)

  status = nf90_inq_varid(ncid, "longitude_corner", varid)
    if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_var(ncid, varid , lon_corner)
    if (status /= nf90_noerr) call handle_err(status)
  
 status = nf90_close(ncid)

end program cmc_proj

  subroutine handle_err(status)
    use netcdf
    integer, intent ( in) :: status
 
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine handle_err

