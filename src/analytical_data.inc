!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Analytical_data, used by Partitioning_class.
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Set analytically concentration and winds.
!-------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!> Winds are set analytically.
!-----------------------------------------------------------------------------
subroutine analytical_winds(test_case, part)
  type(Partition) :: part
  character(*) :: test_case

  if (test_case == "meridional") then
    call meridional_winds(part%grid)
  else if (test_case == "zonal") then
    call zonal_winds(part%grid)
  else if (test_case == "solid_rot") then
    call solid_rot_winds(part)
  else if (test_case == "testsuite_cv") then
    call cv_winds(0.d0, get_dt()*60,  part)
  else
    call print_error("Analytical wind not implemented",&
      "analytical_winds", fname_band)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Winds for zonal solid rotation
!-----------------------------------------------------------------------------
subroutine zonal_winds(grid)
  type(Band_Grid) :: grid
  integer :: i, j, k
  integer :: i_start, i_end
  integer :: j_start, j_end
  double precision :: U_0, coef
  double precision :: dlat, cur
  integer :: nb_lat2, nb_cells
  k = 1

  nb_lat2 = get_nb_lat2_Configuration()
  nb_cells = nb_cells_lat_sector(nb_lat2)
  ! Max is half of the smallest cells length
  U_0 = 60.d0/nb_cells;
  coef = pi/180.d0;
  call interior_lat_indices(i_start, i_end, grid)

  ! For each latitude line
  dlat = get_dlat(grid)
  do i = i_start, i_end
    cur = cell_center_lat(i, grid)
    call interior_lon_indices(j_start, j_end, grid, i)
    ! Interior cells boundaries
    do j = j_start, j_end+1
      grid%zonal_winds(k) = U_0*cos(cur*coef);
      k = k + 1
    end do
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Winds for meridional solid rotation
!-----------------------------------------------------------------------------
subroutine meridional_winds(grid)
  type(Band_Grid) :: grid
  double precision :: U_0

  U_0 = 1.d0/160d0
  grid%merid_winds = U_0
end subroutine

!> Time dependent winds according to Lauritzen et al.
!> For the model, we average winds over a time step
!> Input in seconds
subroutine cv_winds(t1, t2, part)
  type (Partition) :: part
  double precision, intent(in) :: t1, t2
  double precision :: T, R

  R = 6.3172e6 ! meters
  T = 12*24*3600 ! seconds

  call cv_winds_zonal(t1, t2, R, T, part%grid)
  call cv_winds_merid(t1, t2, R, T, part)

#ifdef NOWINDS_CORR
  if (part%id == 0) print *, "No winds corr"
  ! Do not correct winds for testing, but need to change the sing
  part%grid%merid_winds = -part%grid%merid_winds
#else
  call correct_winds(part)
#endif

end subroutine

!> Input in seconds
subroutine cv_winds_zonal(t1, t2, R, T, grid)
  type (Band_grid) :: grid
  double precision, intent(in) :: t1, t2
  double precision :: R, T
  double precision :: U_1, cur
  double precision :: lat, lon
  integer :: i, i_start, i_end
  integer :: j, j_start, j_end
  integer :: k

  k = 1
  grid%zonal_winds = UNDEFINED
  call interior_lat_indices(i_start, i_end, grid)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, grid, i)
    lat = cell_center_lat(i, grid)*pi/180.d0

    ! West boundary for the first cell
    lon = cell_west_lon(i, j_start, grid)*pi/180.d0
    !lon = lon - 2*pi/3
    grid%zonal_winds(k) = cv_zonal_wind(t1, t2, lat, lon, R, T)
    k = k + 1

    ! The rest is to the east
    do j = j_start, j_end
      lon = cell_east_lon(i, j, grid)*pi/180.d0
      !lon = lon - 2*pi/3
      grid%zonal_winds(k) = cv_zonal_wind(t1, t2, lat, lon, R, T)

      k = k + 1
    end do
  end do
end subroutine


!> Input in seconds
!> Input in radians, output in degree/s
function cv_zonal_wind(t1, t2, lat, lon, R, T) result(res)
  double precision, intent(in) :: t1, t2, R, T
  double precision :: lat, lon, coef1, coef2
  double precision :: cost, res, coef, time

  ! Winds are taken at half the timestep
  time = 0.5*(t2+t1)

  lon = lon - 2*pi*time/T
  coef1 = 10*R/T;
  coef2 = 2*pi*R/T;
  cost = cos(pi*time/T)
  res = coef1*sin(lon)**2 * sin(2*lat)*cost + coef2*cos(lat)

  ! Convert from m/s to degree/s
  res = res*180./(RADIUS*PI)

end function

!> Input in radians, output in degree/s
function cv_merid_wind(t1, t2, lat, lon, R, T) result(res)
  double precision, intent(in) :: t1, t2, R, T
  double precision :: lat, lon, coef1, coef2
  double precision :: cost, res, coef, time

  ! Winds are taken at half the timestep
  time = 0.5*(t2+t1)

  lon = lon - 2*pi*time/T
  coef1 = 10*R/T;
  ! Winds is taken at half the timestep
  cost = cos(pi*time/T)

  res = coef1* sin(2*lon)*cos(lat)*cost
  ! Convert from m/s to degree/s
  res = res*180./(RADIUS*PI)


end function

subroutine cv_winds_merid(t1, t2, R, T, part)
  type(Partition), target :: part
  type(Band_Grid), pointer :: grid
  double precision, intent(in) :: t1, t2, T, R
  double precision :: U_1, cur
  double precision :: lat, lon
  integer :: i, i_start, i_end
  integer :: j, j_start, j_end
  integer :: k, j_n1, j_n2
  integer :: l

  grid => part%grid
  call interior_lat_indices(i_start, i_end, grid)
  if (has_north_ghost_cells(grid)) i_start = i_start - 1
  if (.not. has_south_ghost_cells(grid)) i_end = i_end - 1

  k = 1
  grid%merid_winds = UNDEFINED
  do i = i_start, i_end
    !    call interior_lon_indices(j_start, j_end, grid, i)
    call lon_indices(j_start, j_end, part, i)
    lat = cell_south_lat(i, grid)*pi/180.d0

    do j = j_start, j_end

      call south_neighbour_cell_Partition(j_n1, j_n2, i, j, part, i+1)
      do l = j_n1, j_n2
        if (to_or_from_interior_cell(l, i, j, i+1, grid)) then
          lon = cell_interface_middle(i, j, i+1, l, grid)*pi/180.d0

          !lon = lon - 2*pi/3
          grid%merid_winds(k) = cv_merid_wind(t1, t2, lat, lon, R, T)
          !if (part%id == 3 .and. i == 2 .and. j < 160) then
          !  print *, grid%merid_winds(k)
          !end if
          k = k + 1
        end if
      end do
    end do
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Winds for solid rotation around an axis defined by its angles
!-----------------------------------------------------------------------------
subroutine solid_rot_winds(part)
  type(Partition) :: part

  double precision :: alpha, beta_0
  double precision :: U_0, U_1

  ! Degree/s
  U_0 = 0.0001d0*180.d0/pi
  !U_0 = 0.001*180./pi

  ! Rotation axis
  alpha = 0.d0
  beta_0 =0.5d0*pi

  call solid_rot_zonal_winds(U_0, alpha, beta_0, part%grid)
  call solid_rot_merid_winds(U_0, alpha, beta_0, part)
end subroutine

!-------------------------------------------------------------------------------
!> Set all zonal winds for solid rotation
!-------------------------------------------------------------------------------
subroutine solid_rot_zonal_winds(U_1, alpha, beta_0, grid)
  type (Band_grid) :: grid
  double precision :: alpha, beta_0
  double precision :: U_1, cur
  double precision :: lat, lon
  integer :: i, i_start, i_end
  integer :: j, j_start, j_end
  integer :: k

  call interior_lat_indices(i_start, i_end, grid)

  !print *, "U1 alpha beta0", U_1, alpha, beta_0
  k = 1
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, grid, i)
    lat = cell_center_lat(i, grid)*pi/180.d0

    ! West boundary for the first cell
    lon = cell_west_lon(i, j_start, grid)*pi/180.d0
    cur = find_wind_rotation(U_1, lat, lon, alpha, beta_0, .True., i, j)
    grid%zonal_winds(k) = cur
    k = k + 1

    ! The rest is to the east
    do j = j_start, j_end
      lon = cell_east_lon(i, j, grid)*pi/180.d0
      cur = find_wind_rotation(U_1, lat, lon, alpha, beta_0, .True., i,j)
      grid%zonal_winds(k) = cur

      k = k + 1
    end do
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Set all meridional winds for solid rotation
!-------------------------------------------------------------------------------
subroutine solid_rot_merid_winds(U_1, alpha, beta_0, part)
  type(Partition), target :: part
  type(Band_Grid), pointer :: grid
  double precision :: alpha, beta_0
  double precision :: U_1, cur
  double precision :: lat, lon
  integer :: i, i_start, i_end
  integer :: j, j_start, j_end
  integer :: k, j_n1, j_n2
  integer :: l

  grid => part%grid
  call interior_lat_indices(i_start, i_end, grid)
  if (has_north_ghost_cells(grid)) i_start = i_start - 1
  if (.not. has_south_ghost_cells(grid)) i_end = i_end - 1

  k = 1
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, grid, i)
    lat = cell_south_lat(i, grid)*pi/180.d0

    do j = j_start, j_end

      call south_neighbour_cell_Partition(j_n1, j_n2, i, j, part, i+1)
      do l = j_n1, j_n2
        lon = cell_interface_middle(i, j, i+1, l, grid)*pi/180.d0
        cur = find_wind_rotation(U_1, lat, lon, alpha, beta_0, .False., i,j)
        !        if (i == 99 .and. j == j_start) then
        !          print *, "j=",j,"merid winds", cur, lat, lon
        !        end if

        grid%merid_winds(k) = cur
        k = k + 1
      end do
    end do
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Compute zonal or meridional winds at (lat, lon) for solid rotation around an
!> axis.
!> Angle in input must be in radians.
!> No need to use dcos and the likes, the normal functions work in double
!> precision.
!-------------------------------------------------------------------------------
function find_wind_rotation(norm, lat, lon, alpha, beta_0, is_zonal, i, j) result (res)
  double precision, intent(in) :: lat, lon, norm
  double precision, intent(in) :: beta_0, alpha
  logical, intent(in) :: is_zonal
  double precision :: x, y, z
  double precision :: x1, y1, z1
  double precision :: lat1, lon1, sinlon1
  double precision :: u1, u1x, u1y
  double precision :: ux, uy, uz
  double precision :: frac,res
  ! For debug
  integer :: i, j

  ! Rotation around the z-axis, then around y'1 (euler notation)

  ! Cartesian coordinates
  x = cos(lat)*cos(lon)
  y = cos(lat)*sin(lon)
  z = sin(lat)

  ! Switching to the new coordinates, defined by the rotation axis and the 
  ! perpendicular plane
  x1 = cos(beta_0)*cos(alpha)*x + cos(beta_0)*sin(alpha)*y - sin(beta_0)*z
  y1 = -sin(alpha)*x + cos(alpha)*y
  z1 = sin(beta_0)*cos(alpha)*x + sin(beta_0)*sin(alpha)*y + cos(beta_0)*z

  ! Compute the spherical coordinates angles in [-pi/2pi/2]
  lat1 = asin(z1) 
  ! Special case : longitude in [0,2pi]
  frac = x1/cos(lat1)

  ! For acos validity
  if (frac > 1.) frac = 1.
  if (frac < -1.) frac = -1.

  lon1 = acos(frac)
  sinlon1 = y1/cos(lat1)

  if (sinlon1 < 0) then
    lon1 = 2.d0*pi - lon1
  end if

  ! In spherical coordinates, the speed are easy to compute : u is constant,
  ! v = 0. So in cartesian coordinates :
  u1 = norm*cos(lat1)
  u1x = -sin(lon1)*u1
  u1y = cos(lon1)*u1

  ! Switch back to first cartesian coordinates
  ux = cos(beta_0)*cos(alpha)*u1x - sin(alpha)*u1y
  uy = cos(beta_0)*sin(alpha)*u1x + cos(alpha)*u1y
  uz = -sin(beta_0)*u1x

  ! Finally, switch back to spherical coordinates (projection on (u,v))
  if (is_zonal) then
    res = -ux*sin(lon) + uy*cos(lon)
  else
    ! Winds are positive towards the south, so we must invert the sign
    res = ux*sin(lat)*cos(lon) + uy*sin(lon)*sin(lat) - uz*cos(lat)
  end if

end function


