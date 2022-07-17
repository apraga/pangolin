! Read netcdf file with unknown variables and dimensions
program read_data
  use netcdf
  use ScripInterface
  implicit none

  ! Fortran sucks so much we need to define constant string in the array declaration
  ! AND fill each string with spaces if needed
  character(len=14), dimension(12), parameter :: variables = (/ &
  "lat           ", &
  "lon           ", &
  "lev           ", &
  "A_array       ", &
  "B_array       ", &
  "SURFPRESSION  ", &
  "S0VENT_ZONAL  ", &
  "S0VENT_MERIDIE", &
  "S0VITESSE_VERT", &
  "S0TEMPERATURE ", &
  "S0POT_VORTICIT", &
  "S0HUMI.SPECIFI"/) 

  ! The integer is the position in the variables array
  enum, bind(C)
   enumerator :: ZONAL_WINDS = 7, MERID_WINDS = 8, VERT_WINDS = 9
   enumerator :: TEMPERATURE = 10
  end enum

  character(len=3), dimension(3), parameter :: dimensions = (/ &
    "lat", "lon", "lev"/)

  type PangoGrid
    real, allocatable :: lat(:)
    real, allocatable :: lon(:)
    real, allocatable :: lev(:)
    !real, allocatable :: zonal_winds(:,:,:)
    !real, allocatable :: merid_winds(:,:,:)
    real, allocatable :: temperature(:,:,:)
    ! Lon, lat, levels
    integer :: nb_lat
    integer :: nb_lon
    integer :: nb_lev = 1
    ! Total number
    integer :: nb_cells
  end type

  character(120) :: fname

  ! Change interpolation parameter here
  ! Grids names
  character(char_len) :: grid_f1    = 'coarse.nc'
  character(char_len) :: grid_f2    = 'fine.nc'
  ! Mapping files
  character(char_len) :: interp_f1  = 'interp1.nc'
  character(char_len) :: interp_f2  = 'interp2.nc'
  ! Interpolation type
  character(char_len) :: normalize_opt = 'frac'
  ! conservative, distwgt, bilinear
  character(char_len) :: map_opt       = 'conservative'

  type(PangoGrid) :: grid1, grid2
  

  fname = "/wkdir/pae1/emili/DATA_MCG/GLOB22/FM/ECMWF_FC_GLOB22L91/2008/NETCDF/"
  fname = trim(fname) // "FMGLOB22+2008090100.nc"

!  ! Generate analytical grids
  call new_analytical_grid(grid1,45,90)
  call write_file(grid1, "coarse_init.nc")
  call write_scrip_file(grid_f1, "coarse", grid1)

  call new_analytical_grid(grid2,90,180)
  call write_file(grid2, "fine_init.nc")
  call write_scrip_file(grid_f2, "fine", grid2)
!

  call create_mapping(grid_f1, grid_f2, interp_f1, interp_f2, map_opt, normalize_opt)
  call interpolate_data(grid1, grid2, interp_f1, normalize_opt)
  call write_file(grid2, "fine_interp.nc")

  call free_data(grid1)
  call free_data(grid2)
contains

! Create mapping between the grids with SCRIP
subroutine create_mapping(grid_f1, grid_f2, interp_f1, interp_f2,&
    map_opt, normalize_opt)
  character(char_len), intent(in) :: grid_f1, grid_f2
  character(char_len), intent(in) :: interp_f1, interp_f2
  character(char_len), intent(in) :: normalize_opt, map_opt

  print *, "Interpolating grids"
  call scrip_interpol(grid_f1, grid_f2, interp_f1, interp_f2, &
      map_opt, normalize_opt)
 
end subroutine

! Interpolate data from grid 1 to grid 2 using the coefs in interp1
! This means interpolate_data must have been called before !
! Arrayr must be unidimensional
subroutine interpolate_data(grid1, grid2, interp_f1, normalize_opt)
  type(PangoGrid) :: grid1, grid2
  character(char_len), intent(in) :: interp_f1
  character(char_len) :: normalize_opt
  double precision, allocatable :: temp1(:), temp2(:)

  ! Convert to 1D-array
  allocate(temp1(grid1%nb_cells))
  allocate(temp2(grid2%nb_cells))
  temp1 = reshape(grid1%temperature, (/grid1%nb_cells/))
  temp2 = 0.

  call interpolate_array(temp2, temp1, interp_f1, normalize_opt)

  ! Convert to multi-dimensional array
  grid2%temperature = reshape(temp2, (/grid2%nb_lon, grid2%nb_lat, grid2%nb_lev/))

  deallocate(temp1)
  deallocate(temp2)
end subroutine


! Create a regular lat-lon grid from analytical temperature
subroutine new_analytical_grid(grid, nb_lat, nb_lon)
  type(PangoGrid) :: grid
  integer, intent(in) :: nb_lat, nb_lon
  integer :: i, j
  real :: dlat, dlon, lat, lon
  real :: coef = PI/180.

  dlat = 180./nb_lat
  dlon = 360./nb_lon
  grid%nb_lat = nb_lat
  grid%nb_lon = nb_lon
  grid%nb_lev = 1
  grid%nb_cells = grid%nb_lat*grid%nb_lon*grid%nb_lev
  call allocate_data(grid)

  print *, "Generating grid ", nb_lat, nb_lon
  grid%lat =(/(90. -(i-0.5)*dlat, i = 1, nb_lat)/)
  grid%lon =(/((j-0.5)*dlon, j = 1, nb_lon)/)

  do i = 1, grid%nb_lat
    lat = grid%lat(i)
    do j = 1, grid%nb_lon
      lon = grid%lon(j)
      grid%temperature(j, i, 1) = cos(2*lat*coef)*cos(4*lon*coef)
    end do
  end do
end subroutine

! Create grid (with data information) from file
! Format : similar to the array variables
subroutine new_grid_from_file(fname, grid)
  character(*), intent(in) :: fname
  type(PangoGrid) :: grid
  integer :: ndim, nvars, nattr
  integer :: ncid, status, i, k
  character(NF90_MAX_NAME), allocatable :: names_dim(:), names_var(:)
  integer, allocatable :: nb(:)
  integer :: nb_elts(3)

  status = nf90_open(trim(fname), NF90_NOWRITE, ncid)
  call check_netcdf_err(status)

  status = nf90_inquire(ncid, ndim, nvars, nattr)
  call check_netcdf_err(status)

  allocate(names_dim(ndim))
  allocate(nb(ndim))
  allocate(names_var(nvars))

  ! Also check the dimension are ordered properly
  do i = 1, ndim
    status = nf90_inquire_dimension(ncid, i, names_dim(i), nb_elts(i))
    call check_netcdf_err(status)

    if (trim(names_dim(i)) /= dimensions(i)) then
      print *, "Wrong format : ", trim(names_var(i)), dimensions(i)
      call exit(1)
    end if
  end do
  grid%nb_lat= nb_elts(1)
  grid%nb_lon = nb_elts(2)
  grid%nb_lev = nb_elts(3)
  grid%nb_cells = grid%nb_lat*grid%nb_lon*grid%nb_lev

  print *, "latxlonxlev", nb_elts

  k = 1
  ! Also check the variables are ordered properly
  do i = 1, nvars
    status = nf90_inquire_variable(ncid, i, names_var(i))
    call check_netcdf_err(status)
    if (trim(names_var(i)) /= trim(variables(i))) then
      print *, "Wrong format : ", trim(names_var(i)), trim(variables(i))
      call exit(1)
    end if
  end do

  call allocate_data(grid)

  status = nf90_get_var(ncid, TEMPERATURE, grid%temperature)
  call check_netcdf_err(status)

  status = nf90_close(ncid)
  call check_netcdf_err(status)

  deallocate(names_dim)
  deallocate(names_var)
  deallocate(nb)
end subroutine

! Write data for scrip interpol
subroutine write_scrip_file(fname, title, grid)
  character(*), intent(in) :: fname, title
  type(PangoGrid) :: grid
  integer :: status
  integer :: ncid, ids(3), i, varids(6)
  integer :: dimids(3), sizes(6)
  character(50) :: strings(6)

  status = nf90_create(trim(fname), NF90_CLOBBER, ncid)
  call check_netcdf_err(status)

  ! Define metadata
  ! Must have the same size
  strings(1:3) = (/ &
    "grid_size   ", &
    "grid_corners", &
    "grid_rank   "/)
  sizes(1:3) = (/ grid%nb_cells, 4, 2/)
  do i = 1, 3
    status = nf90_def_dim(ncid, trim(strings(i)), sizes(i) , ids(i))
    call check_netcdf_err(status)
  end do

  ! Coordinates
  strings = (/ &
  "grid_dims      ", &
  "grid_center_lat", &
  "grid_center_lon", &
  "grid_imask     ", &
  "grid_corner_lat", &
  "grid_corner_lon" /)
  do i = 1, 6
    if (i == 1) then
      status = nf90_def_var(ncid, strings(i), NF90_INT, ids(3), varids(i))
    else if (i < 5) then
      status = nf90_def_var(ncid, strings(i), NF90_DOUBLE, ids(1), varids(i))
    else
      ! Beware, Fortran writes (grid corners, grid size)
      status = nf90_def_var(ncid, strings(i), NF90_DOUBLE, (/ids(2), ids(1)/), &
        varids(i))
    end if
    call check_netcdf_err(status)
  end do

  ! Units
  do i = 2,6
    if (i == 4) then
      status = nf90_put_att(ncid, varids(i), "units", "unitless")
    else
      status = nf90_put_att(ncid, varids(i), "units", "degrees")
    end if
    call check_netcdf_err(status)
  end do

  status = nf90_put_att(ncid, NF90_GLOBAL, "title", title)
  call check_netcdf_err(status)

  status = nf90_enddef(ncid)
  call check_netcdf_err(status)

  call write_scrip_coords(ncid, varids, grid)
  status = nf90_close(ncid)
  call check_netcdf_err(status)
end subroutine

subroutine write_scrip_coords(ncid, varids, grid)
  type(PangoGrid) :: grid
  integer, intent(in) :: ncid, varids(:)
  real, allocatable :: centers(:,:), corners(:,:,:)
  integer :: status
  real :: dlat, dlon
  integer :: i, j, k

  status = nf90_put_var(ncid, varids(1), (/grid%nb_lat, grid%nb_lon/))

  dlat = 180./grid%nb_lat
  dlon = 360./grid%nb_lon

  ! Store both lat and lon
  allocate(centers(2, grid%nb_cells))

 ! Center coordinates
  k = 1
  !do i = 1, grid%nb_lat
  do i = grid%nb_lat, 1, -1
    do j = 1, grid%nb_lon
      centers(1, k) = 90. -(i-0.5)*dlat ! lat
      centers(2, k) = (j-0.5)*dlon ! lon
      k = k + 1
    end do
  end do

  do i = 1, 2
    status = nf90_put_var(ncid, varids(i+1), centers(i,:))
    call check_netcdf_err(status)
  end do

  ! Mask size : 1 for plotting a cell
  centers(1,:) = 1
  status = nf90_put_var(ncid, varids(4), centers(1,:))
  call check_netcdf_err(status)

  ! Grid corners in counter-clockwise
  ! Beware, Fortran writes (grid corners, grid size)
  allocate(corners(2, 4, grid%nb_cells))
  k = 1
  !do i = 1, grid%nb_lat
  do i = grid%nb_lat, 1, -1
    do j = 1, grid%nb_lon
      corners(1, 1:2, k) = 90. -i*dlat
      corners(1, 3:4, k) = 90. -(i-1)*dlat ! lat

      corners(2, 4, k) = (j-1)*dlon ! lon
      corners(2, 2:3, k) = j*dlon 
      corners(2, 1, k ) = (j-1)*dlon
      k = k + 1
    end do
  end do

  do i = 1, 2
    status = nf90_put_var(ncid, varids(i+4), corners(i,:,:))
    call check_netcdf_err(status)
  end do

  deallocate(centers)
  deallocate(corners)

end subroutine


! Write grids centers and temperature
subroutine write_file(grid, fname)
  character(*), intent(in) :: fname
  type(PangoGrid) :: grid
  integer :: status, varids(4)
  integer :: ncid, ids(3), i
  integer :: dimids(3), nb_elts(3)
  character(char_len) :: names(3)

  print *, "Writing ", trim(fname)
  status = nf90_create(trim(fname), NF90_CLOBBER, ncid)
  call check_netcdf_err(status)

  ! Define metadata

  nb_elts = (/grid%nb_lat, grid%nb_lon, grid%nb_lev/)
  do i = 1, 3
    status = nf90_def_dim(ncid, dimensions(i), nb_elts(i), ids(i))
    call check_netcdf_err(status)
  end do

  ! Coordinates
  names = (/"lat", "lon", "lev"/)
  do i = 1, size(names)
    status = nf90_def_var(ncid, trim(names(i)), NF90_DOUBLE, ids(i), varids(i))
    call check_netcdf_err(status)
  end do

  ! Units
  status = nf90_put_att(ncid, varids(1), "units", "degrees_north")
  call check_netcdf_err(status)
  status = nf90_put_att(ncid, varids(2), "units", "degrees_east")
  call check_netcdf_err(status)
  !status = nf90_put_att(ncid, varids(3), "units", "level")
  !call check_netcdf_err(status)

  ! Data
  ! lon, lat, levels
  dimids = (/ids(2), ids(1), ids(3)/)
  status = nf90_def_var(ncid, "temperature", NF90_DOUBLE, dimids, varids(4))
  call check_netcdf_err(status)

  status = nf90_enddef(ncid)
  call check_netcdf_err(status)

  ! Define data
  status = nf90_put_var(ncid, varids(1), grid%lat)
  call check_netcdf_err(status)
  status = nf90_put_var(ncid, varids(2), grid%lon)
  call check_netcdf_err(status)
  ! dummy levels
  status = nf90_put_var(ncid, varids(3), grid%lev)
  call check_netcdf_err(status)

  status = nf90_put_var(ncid, varids(4), grid%temperature)
  call check_netcdf_err(status)

  status = nf90_close(ncid)
  call check_netcdf_err(status)
end subroutine


subroutine allocate_data(grid)
  type(PangoGrid) :: grid

  allocate(grid%lon(grid%nb_lon))
  allocate(grid%lat(grid%nb_lat))
  allocate(grid%lev(grid%nb_lev))
  allocate(grid%temperature(grid%nb_lon, grid%nb_lat, grid%nb_lev))
  !allocate(grid%zonal_winds(nlon, nlat, nlev))
  !allocate(grid%merid_winds(nlon, nlat, nlev))
end subroutine


subroutine free_data(grid)
  type(PangoGrid) :: grid

  if (allocated(grid%lon)) deallocate(grid%lon)
  if (allocated(grid%lat)) deallocate(grid%lat)
  !%if (allocated(grid%lev)) deallocate(grid%lev)

  if (allocated(grid%temperature)) deallocate(grid%temperature)
  !if (allocated(grid%zonal_winds)) deallocate(grid%zonal_winds)
  !if (allocated(grid%merid_winds)) deallocate(grid%merid_winds)
end subroutine

end program
