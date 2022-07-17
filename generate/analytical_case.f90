subroutine cartesian_coordinates(x, y, z, lat, lon)
  double precision, intent(in) :: lat, lon
  double precision :: x, y, z

  ! Cartesian coordinates
  x = cos(lat)*cos(lon)
  y = cos(lat)*sin(lon)
  z = sin(lat)
!  print *, "xyz", x, y,z, "latlon", lat, lon

end subroutine

subroutine spherical_coordinates(lat1, lon1, x1, y1, z1)
  double precision, intent(in) :: x1, y1, z1
  double precision :: lat1, lon1, sinlon1
  double precision :: pi = 4.d0*atan(1.d0)
  double precision :: prec = 1e-13

  ! Compute the spherical coordinates angles in [-pi/2pi/2]
  lat1 = asin(z1) 

  !frac = y1/x1
  !lon1 = atan(frac)
  !if (x1/cos(lat1) < 0 .and. y1/cos(lat1) < 0) lon1 = pi + lon1 
  !if (x1/cos(lat1) < 0 .and. y1/cos(lat1) > 0) lon1 = pi - lon1 
  !lon1 = lon1 + 0.5*pi

  ! longitude in [0,pi]
  frac = x1/cos(lat1)

  ! For acos validity (should not happen)
  if (frac > 1.) frac = 1.
  if (frac < -1.) frac = -1.

  lon1 = acos(frac)
  sinlon1 = y1/cos(lat1)

  ! longitude in [0,2pi]
  if (sinlon1 < prec) then
    lon1 = 2.d0*pi - lon1
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Switching to the new coordinates, defined by euler angles (rotation axis and the 
!> perpendicular plane)
!> Input and output in [0,2pi]
!-------------------------------------------------------------------------------
subroutine new_coordinates(x1, y1, z1, x, y, z, alpha, beta)
  double precision, intent(in) :: x, y, z, alpha, beta
  double precision :: x1, y1, z1

  x1 = cos(beta)*cos(alpha)*x + cos(beta)*sin(alpha)*y - sin(beta)*z
  y1 = -sin(alpha)*x + cos(alpha)*y
  z1 = sin(beta)*cos(alpha)*x + sin(beta)*sin(alpha)*y + cos(beta)*z

end subroutine

subroutine spherical_to_old_cartesian(x, y, z, x1, y1, alpha, beta)
  double precision :: x, y, z
  double precision, intent(in) :: x1, y1, alpha, beta

  x = cos(beta)*cos(alpha)*x1 - sin(alpha)*y1
  y = cos(beta)*sin(alpha)*x1 + cos(alpha)*y1
  z = -sin(beta)*x1

end subroutine


!-------------------------------------------------------------------------------
!> Compute zonal or meridional winds at (lat, lon) for solid rotation around an
!> axis  defined by (alpha, beta).
!> Angle in input must be in radians.
!> No need to use dcos and the likes, the normal functions work in double
!> precision.
!-------------------------------------------------------------------------------
function find_wind_rotation(norm, lat, lon, alpha, beta, is_zonal) result (res)
  double precision, intent(in) :: lat, lon, norm
  double precision, intent(in) :: beta, alpha
  logical, intent(in) :: is_zonal
  double precision :: x, y, z
  double precision :: x1, y1, z1
  double precision :: lat1, lon1, sinlon1
  double precision :: u1, u1x, u1y
  double precision :: ux, uy, uz
  double precision :: frac,res

  ! Rotation around the z-axis, then around y'1 (euler notation)
  call new_spherical_coordinates(lat1, lon1, lat, lon, alpha, beta)

  ! In spherical coordinates, the speed are easy to compute : u is constant,
  ! v = 0. So in cartesian coordinates :
  u1 = norm*cos(lat1)
  u1x = -sin(lon1)*u1
  u1y = cos(lon1)*u1

  ! Switch back to first cartesian coordinates
  call spherical_to_old_cartesian(ux, uy, uz, u1x, u1y, alpha, beta)

  ! Finally, switch back to spherical coordinates (projection on (u,v))
  if (is_zonal) then
    res = -ux*sin(lon) + uy*cos(lon)
  else
    ! Winds are positive towards the south, so we must invert the sign
    res = ux*sin(lat)*cos(lon) + uy*sin(lon)*sin(lat) - uz*cos(lat)
  end if

end function

!> Angle in input must be in radians.
subroutine new_spherical_coordinates(lat1, lon1, lat, lon, alpha, beta)
  double precision, intent(in) :: lat, lon, beta, alpha
  double precision :: x, y, z
  double precision :: x1, y1, z1
  double precision :: lat1, lon1

  lat1 = lat
  lon1 = lon
  if (alpha == 0. .and. beta == 0) return

  ! Rotation around the z-axis, then around y'1 (euler notation)
  call cartesian_coordinates(x, y, z, lat, lon)
  call new_coordinates(x1, y1, z1, x, y, z, alpha, beta)
  call spherical_coordinates(lat1, lon1, x1, y1, z1)

end subroutine
