!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Band_Grid
!
!> @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Class for a generic grid defined in bands of latitude
!
!-------------------------------------------------------------------------------
module Band_Grid_class
use Data_check
use Global_grid_class
use Analytical_partitioning
use Parameters
use List_Tracers_class


implicit none


!> The grid is defined as bands. We store the positions of the boundaries
!> Latitude is defined in degree in [-90,90] and the longitude is in [0,360]
!> Also, the ghost cells are already contained in the grid (and not in a
!> separate array)
type Band_Grid
  !> Varition of latitude inside a cell at a given latitude (constant here). Positive on the
  !> northern hemisphere, negative on the southern.
  double precision,allocatable :: dlat(:)
  !> Variation of latitude inside a cell (constant at a given latitude)
  double precision,allocatable :: dlon(:)
  !> Stores the position of the first cell in the array at each latitude line. Allows for faster access
  !> for each cell. It is actually the sum of all previous nb_lon
  integer, allocatable :: first_cell(:)
  !> Latitude of the upper north line (in lat-lon plane)
  double precision :: lat_min = -1.d0
  !> Longitude of the leftmost cells (in lat-lon plane)
  double precision,allocatable :: lon_min(:)
  !> Number of east ghost cells
  integer, allocatable :: nb_ghosts_east(:)
  !> Number of west ghost cells
  integer, allocatable :: nb_ghosts_west(:)
  !> Number of latitudes
  integer :: nb_lat = -1
  !> Number of cells at a given latitude
  integer,allocatable :: nb_lon(:)

  !****************************************************************************!
  !            Analytical properties                      !
  !****************************************************************************!
  !> Gives the offset of the north neighbour for the first cell in the 
  !> current latitude line. This will give us the offset between the north
  !> neighbour and the current line. Beware, there are only nb_lat - 1
  integer, allocatable :: north_offset(:)
  !> Same as before, but for south neighbours.
  integer, allocatable :: south_offset(:)

  !> Tell us if there are ghost cells in a direction in base 2, north being the
  !> lowest indice, west being the fourth.
  integer :: ghost_cells = 0

  !> All data needed for multiple tracers (see List_Tracers_class)
  type(List_Tracers) :: tracers

  !****************************************************************************!
  !            MPI data                                                        !
  !****************************************************************************!
  ! mpi type (one for each instance)
  integer :: mpi_grid = mpi_datatype_null
  !> MPI type for ghost cells, that is, for data received from the neighbours.
  !> One type for each neighbour.
  integer, allocatable :: mpi_ghost_cells(:)
  !> MPI type for interior cells, that is, for data sent to neighbours.
  !> One type for each 
  integer, allocatable :: mpi_interior_cells(:)
  !> For zonal winds (only for analytical winds)
  integer :: mpi_zwinds_east = mpi_datatype_null
  integer :: mpi_zwinds_west = mpi_datatype_null

  !****************************************************************************!
  !              Field for advection                                           !
  !****************************************************************************!
  !> Air and tracer masses
  double precision, allocatable :: mass_air(:)
  !> Winds data  (with ghost cells) : current and previous/next for
  !> interpolation
  double precision, allocatable :: merid_winds(:)
  double precision, allocatable :: merid_winds_prev(:)
  double precision, allocatable :: merid_winds_next(:)
  double precision, allocatable :: zonal_winds(:)
  double precision, allocatable :: zonal_winds_prev(:)
  double precision, allocatable :: zonal_winds_next(:)

  !****************************************************************************!
  !              For (parallel) hdf5 I/O                                     !
  !****************************************************************************!
#ifdef HDF5
  !> Arrays of cell centers coordinates 
  double precision, allocatable :: center_lat(:)
  double precision, allocatable :: center_lon(:)
  !> Arrays of zonal winds coordinates 
  double precision, allocatable :: u_lat(:)
  double precision, allocatable :: u_lon(:)
  !> Arrays of merid winds coordinates 
  double precision, allocatable :: v_lat(:)
  double precision, allocatable :: v_lon(:)

  !> Needed for I/O using chunks of data
  !> It stores, the offset of the first cell/wind of the line in the global grid
  integer(hsize_t), allocatable :: ratio_offset(:)
  integer(hsize_t), allocatable :: zonal_offset(:)
  integer(hsize_t), allocatable :: merid_offset(:)
  !> Store the length of meridional winds on the line
  integer(hsize_t), allocatable :: merid_length(:)
#endif
   
end type

! For output_
character(*), parameter :: fname_band = "band_grid_class.F90"

contains 

!-----------------------------------------------------------------------------
!> Default constructor. Does not set number of ghosts
!-----------------------------------------------------------------------------
subroutine new_Band_Grid(this, n_lat, dlat, nb_lon, dlon, lat_min, lon_min,&
    ghost_cells)
  type(Band_Grid) :: this
  integer :: n_lat
  double precision :: dlat(:), dlon(:), lon_min(:), lat_min
  integer :: nb_lon(:)
  integer, optional :: ghost_cells

  this%nb_lat = n_lat
  call allocate_most_data(this, n_lat)

  this%dlon = dlon
  this%dlat = dlat
  this%lat_min = lat_min
  this%lon_min = lon_min
  this%nb_lon = nb_lon

  call set_first_cells_Band_grid(this)
  if (present(ghost_cells)) then
    this%ghost_cells = ghost_cells
  end if

  ! Initialize mpi type
  this%mpi_grid = mpi_datatype_null
end subroutine

!-----------------------------------------------------------------------------
!> Copy constructor.
!-----------------------------------------------------------------------------
subroutine new_Band_Grid_copy(this, other)
  type(Band_Grid) :: this, other
  integer :: n
  character(*), parameter :: func_name ="new_Band_grid_copy" 

  this%nb_lat = other%nb_lat
  this%ghost_cells = other%ghost_cells

  call check_allocated(this%dlat, "dlat", func_name, fname_band)
  call check_allocated(this%dlon, "dlon", func_name, fname_band)
  call check_allocated(this%first_cell, "first_cell", func_name, fname_band)
  call check_allocated(this%lon_min, "lon_min", func_name, fname_band)
  call check_allocated(this%nb_ghosts_east, "nb_ghosts_east", func_name, fname_band)
  call check_allocated(this%nb_ghosts_west, "nb_ghosts_west", func_name, fname_band)
  call check_allocated(this%nb_lon, "nb_lon", func_name, fname_band)
  call check_allocated(this%north_offset, "north_offset", func_name, fname_band)
  call check_allocated(this%south_offset, "south_offset", func_name, fname_band)

  n = this%nb_lat 
  call allocate_most_data(this, n)

  n = size(other%first_cell)
  allocate (this%first_cell(n))

  this%dlon = other%dlon
  this%dlat = other%dlat
  this%first_cell = other%first_cell
  this%lat_min = other%lat_min
  this%lon_min = other%lon_min
  this%nb_ghosts_west = other%nb_ghosts_west
  this%nb_ghosts_east = other%nb_ghosts_east
  this%nb_lon = other%nb_lon
  this%north_offset = other%north_offset
  this%south_offset = other%south_offset

  ! Initialize mpi type
  this%mpi_grid = mpi_datatype_null

end subroutine


!-----------------------------------------------------------------------------
!> Must compte size of meridional winds before (need partition info)
!-----------------------------------------------------------------------------
subroutine allocate_winds(size_merid, this)
  type (Band_Grid) :: this
  integer, intent(in) :: size_merid
  integer :: data_size
  character(*), parameter :: func_name ="allocate_winds" 


  call check_allocated(this%merid_winds, "merid winds", func_name, fname_band)
  call check_allocated(this%merid_winds_prev, "merid winds prev", func_name, &
    fname_band)
  call check_allocated(this%merid_winds_next, "merid winds next", func_name, &
    fname_band)
  call check_allocated(this%zonal_winds, "zonal winds", func_name, fname_band)
  call check_allocated(this%zonal_winds_prev, "zonal winds prev", func_name, &
    fname_band)
  call check_allocated(this%zonal_winds_next, "zonal winds next", func_name, &
    fname_band)
  
  data_size = size_winds_zonal(this)
  if (data_size < 0.) then
    call print_error("Negative zonal winds size", "allocate_winds", &
      fname_band)
  end if
  allocate (this%zonal_winds(data_size))
  allocate (this%zonal_winds_prev(data_size))
  allocate (this%zonal_winds_next(data_size))
  this%zonal_winds = undefined
  this%zonal_winds_prev = undefined
  this%zonal_winds_next = undefined

  if (size_merid < 0.) then
    call print_error("Negative meridional winds size", "allocate_winds", &
      fname_band)
  end if
  allocate (this%merid_winds(size_merid))
  allocate (this%merid_winds_prev(size_merid))
  allocate (this%merid_winds_next(size_merid))
  this%merid_winds = undefined
  this%merid_winds_prev = undefined
  this%merid_winds_next = undefined
end subroutine

!-----------------------------------------------------------------------------
!> Allocate data used in advection : slope, fluxes, mass
!-----------------------------------------------------------------------------
subroutine allocate_advection_data_Band_grid(size_merid, this)
  type (Band_Grid) :: this
  integer, intent(in) :: size_merid
  integer :: size_ratio, size_zonal
  character(*), parameter :: func_name ="allocate_advection_data_Band_grid" 

  size_ratio = sum(this%nb_lon)
  if (size_ratio < 0.) then
    call print_error("Negative ratio size", func_name, fname_band)
  end if

  call allocate_winds(size_merid, this) 

  size_zonal = size(this%zonal_winds)
  call new_List_Tracers(this%tracers, size_ratio, size_zonal, size_merid)

  call check_allocated(this%mass_air, "mass_air", func_name, fname_band)
  allocate (this%mass_air(size_ratio))

#ifdef HDF5
  call allocate_hdf5_data(this)
#endif
end subroutine

#ifdef HDF5
!-----------------------------------------------------------------------------
!> Allocate data used in advection : slope, fluxes, mass
!> Generation is done in partition_class
!-----------------------------------------------------------------------------
subroutine allocate_hdf5_data(this)
  type (Band_Grid) :: this
  integer :: nb_lat, nb_cells

  nb_lat = get_nb_lat_Band_grid(this)
  allocate(this%ratio_offset(nb_lat))
  allocate(this%zonal_offset(nb_lat))

  nb_cells = size(this%merid_length)
  allocate(this%merid_offset(nb_cells))
  !> merid_length is allocated in partition_class

  nb_cells = nb_cells_Band_grid(this)
  allocate(this%center_lat(nb_cells))
  allocate(this%center_lon(nb_cells))

  nb_cells = size(this%zonal_winds)
  allocate(this%u_lat(nb_cells))
  allocate(this%u_lon(nb_cells))

  nb_cells = size(this%merid_winds)
  allocate(this%v_lat(nb_cells))
  allocate(this%v_lon(nb_cells))

end subroutine

!> Cell center and zonal winds positions
subroutine generate_ratio_zonal_positions(this)
  type (Band_Grid) :: this
  double precision :: cur_lat
  integer :: i_start, i_end, i
  integer :: j_start, j_end, j
  integer :: nb_cells, k, k2

  nb_cells = nb_cells_Band_grid(this)

  k = 1
  k2 = 1
  call interior_lat_indices(i_start, i_end, this)
  do i = i_start, i_end
    cur_lat = cell_center_lat(i, this)
    call interior_lon_indices(j_start, j_end, this, i)
    do j = j_start, j_end
      this%center_lon(k) = cell_center_lon(i, j, this)
      this%center_lat(k) = cur_lat

      this%u_lat(k2) = cur_lat
      this%u_lon(k2) = cell_west_lon(i, j, this)
      k = k + 1
      k2 = k2 + 1
    end do

    this%u_lat(k2) = cur_lat
    this%u_lon(k2) = cell_west_lon(i, j_end+1, this)
    k2 = k2 + 1

  end do
end subroutine

#endif

#ifdef SYMMETRIC_SPLIT
subroutine cp_ratio_to_buffer(this, id)
  type(Band_grid) :: this
  integer, intent(in) :: id

  call cp_ratio_to_buffer_all(this%tracers, id)
end subroutine

subroutine cp_buffer_to_ratio(this, id)
  type(Band_grid) :: this
  integer, intent(in) :: id

  call cp_buffer_to_ratio_all(this%tracers, id)
end subroutine

subroutine average_with_buffer(this, id)
  type(Band_grid) :: this
  integer, intent(in) :: id

  call average_with_buffer_all(this%tracers, id)
end subroutine
#endif

!-----------------------------------------------------------------------------
!> Size of zonal winds
!-----------------------------------------------------------------------------
function size_winds_zonal(this) result(nb)
  type (Band_Grid) :: this
  integer :: nb
  integer :: i_start, i_end, i

  nb = 0
  call interior_lat_indices(i_start, i_end, this)
  do i = i_start, i_end
    nb = nb + nb_lon_Band_grid(i, this) + 1
  end do

end function

!!-----------------------------------------------------------------------------
!!> Size of meridional winds. We do not have north winds at the north pole, nor
!!> south winds at the south pole.
!!-----------------------------------------------------------------------------
!function size_winds_merid(this) result(nb)
!  type (Band_Grid) :: this
!  integer :: nb
!  integer :: i_start, i_end, i
!
!  nb = 0
!  call interior_lat_indices(i_start, i_end, this)
!  ! North and south extremities, if there are ghost cells
!  if (has_north_ghost_cells(this)) i_start = i_start - 1
!  if (.not. has_south_ghost_cells(this)) i_end = i_end - 1
!
!  do i = i_start, i_end
!    nb = nb + nb_interfaces_lat(i, this) 
!  end do
!
!end function

!-----------------------------------------------------------------------------
!> Sums the number of meridional interfaces up to line i (south latitude)
!-----------------------------------------------------------------------------
function sum_nb_interfaces_merid(i, this) result (nb)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  integer :: nb, k

  nb = 0
  if (i < 1 .or. i > this%nb_lat) return
  do k = 1, i
    nb = nb + nb_interfaces_lat(k, this) 
  end do

end function

!-----------------------------------------------------------------------------
!> Sums the number of zonal interfaces up to line i (south latitude)
!-----------------------------------------------------------------------------
function sum_nb_interfaces_zonal(i, this) result (nb)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  integer :: nb, k

  nb = 0
  if (i < 1 .or. i > this%nb_lat) return
  do k = 1, i
    nb = nb + nb_lon_Band_grid(k, this) + 1
  end do

end function


!-----------------------------------------------------------------------------
!> Number of cell interfaces between latitude line i and i+1. Only inside the current
!> grid.
!> Assume i and i+1 exists.
!-----------------------------------------------------------------------------
function nb_interfaces_lat(i, this) result (nb)
  type (Band_Grid) :: this
  integer :: nb
  integer, intent(in) :: i
  integer :: ratio

  if (is_latitude_just_before_equator_Band_grid(i, this)) then
    nb = nb_lon_Band_grid(i, this)
    return
  end if

  ! For a single partition, we go back to one sector
  ratio = 1
  if (has_single_partition()) ratio = 3

  ! + 1 for the lef interface to ghost cells
  nb = nb_lon_Band_grid(i, this)/ratio + nb_lon_Band_grid(i+1, this)/ratio + 1

  if (is_first_on_band_Band_grid(this)) nb = nb - 1
  if (is_last_on_band_Band_grid(this)) nb = nb - 1

  nb = nb*ratio
end function

!-----------------------------------------------------------------------------
!> Returns true if the latitude line i is just before the equator
!> Concurrent of is_latitude_just_before_equator
!> Uses float. TODO : optimize
!-----------------------------------------------------------------------------
function is_latitude_just_before_equator_Band_grid(i, this) result (equator)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  logical :: equator
  double precision :: lat, dlat

  lat = cell_center_lat(i, this)
  dlat = get_dlat(this)
  equator = lat > 0 .and. lat < dlat
end function


!-----------------------------------------------------------------------------
!> Concurrent version of is_first_on_band for a grid. Uses float.
!> TODO : optimise
!-----------------------------------------------------------------------------
function is_first_on_band_Band_grid(this) result(res)
  type (Band_Grid) :: this
  integer :: i_start, i_end
  integer :: j_start, j_end, k
  logical :: res
  double precision :: mid, dlon

  call interior_lat_indices(i_start, i_end, this)
  call interior_lon_indices(j_start, j_end, this, i_start)

  mid = cell_center_lon(i_start, j_start, this)
  dlon = this%dlon(i_start)
  res = .False.
  do k = 0, 2
    res = res .or. (mid > 120.d0*k .and. mid < 120.d0*k + dlon)
  end do

end function

!-----------------------------------------------------------------------------
!> Concurrent version of is_last_on_band for a grid. Uses float.
!> TODO : optimise
!-----------------------------------------------------------------------------
function is_last_on_band_Band_grid(this) result(res)
  type (Band_Grid) :: this
  integer :: i_start, i_end
  integer :: j_start, j_end, k
  logical :: res
  double precision :: mid, dlon

  call interior_lat_indices(i_start, i_end, this)
  call interior_lon_indices(j_start, j_end, this, i_start)

  mid = cell_center_lon(i_start, j_end, this)
  dlon = this%dlon(i_start)
  res = .False.
  do k = 1, 3
    res = res .or. (mid < 120.d0*k .and. mid > 120.d0*k - dlon)
  end do

end function


!-----------------------------------------------------------------------------
!> Initializes arrays for advection 
!-----------------------------------------------------------------------------
subroutine init_advection_data(this) 
  type(Band_Grid) :: this
  integer :: i, i_start, i_end
  integer :: j, j_start, j_end
  double precision :: m_i

  call interior_lat_indices(i_start, i_end, this)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, this, i)
    m_i = cell_area(i, j_start, this)
    do j = j_start, j_end
      call set_cell_air_mass(m_i, i, j, this)
    end do
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Number of cells except ghost cells
!-----------------------------------------------------------------------------
function nb_cells_Band_grid(this) result (nb)
  type(Band_Grid) :: this
  integer :: nb
  integer :: i_start, i_end, i

  nb = 0
  call interior_lat_indices(i_start, i_end, this)
  do i = i_start, i_end
    nb = nb + nb_lon_Band_grid(i, this)
  end do
end function

!-----------------------------------------------------------------------------
!> Set the position of first cells on each latitude line
!-----------------------------------------------------------------------------
subroutine set_first_cells_Band_grid(this)
  type(Band_Grid) :: this
  integer :: i

  allocate (this%first_cell(this%nb_lat))
  this%first_cell(1) = 1
  do i = 2, this%nb_lat
    this%first_cell(i) = sum(this%nb_lon(1:i-1)) + 1
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Allocate most of the data
!-----------------------------------------------------------------------------
subroutine allocate_most_data(this, nb_lat)
  type(Band_Grid) :: this
  integer, intent(in) :: nb_lat

  allocate (this%dlat(1))
  this%dlat = 0.
  allocate (this%nb_lon(nb_lat))
  this%nb_lon = 0
  allocate (this%dlon(nb_lat))
  this%dlon = 0.
  allocate (this%lon_min(nb_lat))
  this%lon_min = 0.

  ! True number of neighbours
  allocate (this%north_offset(nb_lat))
  this%north_offset = 0
  allocate (this%south_offset(nb_lat))
  this%south_offset = 0
  allocate (this%nb_ghosts_east(nb_lat))
  this%nb_ghosts_east = 0
  allocate (this%nb_ghosts_west(nb_lat))
  this%nb_ghosts_west = 0
end subroutine

!-----------------------------------------------------------------------------
!> Create a grid based on its index in a (lat-lon) matrix. Grid is in most case
!> square but the last band with a square number of partition can sometimes be 
!> rectangular hence the height and width. However, only the height varies,
!> while the width is constant. We then need the first latitude.
!> @param zon : partition zone
!> @param first_lat : starting latitude from the closest pole (beware of the 
!> minimal latitude)
!> @param nb_resized : number of resized partition on the current band (only 
!> one one half)
!> @param width : width of a square partition (before resizing)
!> @param i_pos : position in the partition matrix (latitude)
!> @param j_pos : position in the partition matrix (longitude)
!-----------------------------------------------------------------------------
subroutine new_Band_grid_from_pos(this, i_pos, j_pos, first_lat, height, &
    width, nb_resized, zone, ghost_cells)
  type(Band_Grid) :: this
  integer, intent(in) :: i_pos, j_pos, height
  integer, intent(in) :: width, first_lat
  integer, intent(in) :: zone, ghost_cells
  integer :: k, i_lat, nb_resized
  integer, allocatable :: list_firsts(:)
  double precision :: cur_dlon
  double precision :: dlat(1)
  integer :: nb_cells, nb_lat, nb
  integer :: first_lat_ghost

  dlat = get_dlat_Global_grid()

  ! Number of ghost cells in the 4 directions
  this%ghost_cells = ghost_cells

  nb_lat = height
  ! North-south ghost cells
  if (has_north_ghost_cells(this)) nb_lat = nb_lat + 1 
  if (has_south_ghost_cells(this)) nb_lat = nb_lat + 1
  this%nb_lat = nb_lat

  ! Ghost cells are in the array, so we allocate it here
  first_lat_ghost = first_lat
  ! On southern hemisphere, we revert the first lat
  if (.not. on_northern_hemisphere(zone)) then
    nb = get_nb_lat_Global_grid()
    ! Beware, we must not count the ghost cells now
    first_lat_ghost = nb - first_lat - height 
  end if

  ! Correct if there are north ghost cells 
  if (has_north_ghost_cells(this)) first_lat_ghost = first_lat_ghost - 1

  call allocate_most_data(this, nb_lat)

  ! Used for storing data used for neighbours positions
  allocate (list_firsts(nb_lat))

  this%dlat = dlat

  ! We always count from the north pole
  this%lat_min = 90 - first_lat_ghost*this%dlat(1)

  ! Southern hemisphere : start from nb_lat2 and decreas
  ! First set longitude without ghost cells
  do k = 1, this%nb_lat
    i_lat =  first_lat_ghost + k

    cur_dlon = get_dlon_Global_grid(i_lat)
    nb_cells = nb_cells_lat_sector(i_lat)

    call set_longitude_data(this, k, i_pos, j_pos, cur_dlon, nb_resized, &
      width, nb_cells, zone, list_firsts)
  end do

  ! Now add and store east/west neighbours
  call set_data_for_neighbours(this, list_firsts, first_lat_ghost, zone)
  ! Update grid spec for ghost
  call add_east_west_ghosts(this)


  deallocate (list_firsts)

  ! Initialize mpi type
  this%mpi_grid = mpi_datatype_null

end subroutine

!-------------------------------------------------------------------------------
!> Store the position of the first north and south neighbour line, and the
!> number of east and west ghost cells
!> @param list_firsts : list of the first cells on each latitude line
!> @param first_lat : first latitude (frome the north pole) of the partition
!-------------------------------------------------------------------------------
subroutine set_data_for_neighbours(this, list_firsts, first_lat, zone)
  type(Band_Grid) :: this
  integer, intent(in) :: list_firsts(:)
  integer, intent(in) :: first_lat, zone
  integer :: k, i_lat

  this%north_offset = 0
  this%south_offset = 0
  this%nb_ghosts_west = 0
  this%nb_ghosts_east = 0

  do k = 1, this%nb_lat
    i_lat = first_lat + k
    call set_nb_east_west_ghosts(i_lat, k, this, zone, list_firsts)
  end do

  do k = 1, this%nb_lat
    i_lat = first_lat + k
    call set_north_south_offsets(i_lat, k, this, zone, list_firsts)
  end do

end subroutine

!-------------------------------------------------------------------------------
!> Set the north and south offset
!> @param i_lat : latitude line from the north pole
!> @param k : line inside the grid
!-------------------------------------------------------------------------------
subroutine set_north_south_offsets(i_lat, k, this, zone, list_firsts)
  type (Band_Grid) :: this
  integer, intent(in) :: i_lat, k, zone
  integer, intent(in) :: list_firsts(:)
  logical :: not_first, not_last
  integer :: mid, neighb
  integer :: first
  integer :: tmp

  not_first = (.not. is_north_ghost(k, this) .and. k > 1)
  not_last = (.not. is_south_ghost(k, this) .and. k < this%nb_lat)

  ! Special case
  if (offsets_sector_extrem(not_first, not_last, k, this, list_firsts, zone)) &
    return

  mid = (nb_cells_lat_sector(i_lat) + 1)/2
  ! Don't forget west ghost cells
  first = list_firsts(k) - this%nb_ghosts_west(k)

  ! The offset is just the difference between the first cell on the neighbour
  ! line and the neihbour of the first cell on the current line
  ! Beware of the indices though for the offset
  if (not_first) then
    neighb = left_north_cell(first, mid, zone)
    tmp = neighb - list_firsts(k-1) + this%nb_ghosts_west(k-1)
    this%north_offset(k) = tmp
  end if

  if (not_last) then
    neighb = left_south_cell(first, mid, zone)
    tmp = neighb - list_firsts(k+1) + this%nb_ghosts_west(k+1)
    this%south_offset(k) = tmp
  end if
end subroutine

!-------------------------------------------------------------------------------
! Special case for north/south offset : sector extremity
!-------------------------------------------------------------------------------
function offsets_sector_extrem(not_first, not_last, k, this, list_firsts, zone) &
    result (is_done)
  type (Band_Grid) :: this
  integer, intent(in) :: k, zone
  integer, intent(in) :: list_firsts(:)
  logical, intent(in) :: not_first, not_last
  logical :: is_extr
  logical :: is_done, is_north
  logical :: to_north_ghost, to_south_ghost
  integer :: l
  ! North-south couples
  integer :: vals(4) = (/1, -2, -1, 0/)

  is_done = .False.
  is_extr = .True.
  do l = 2 , size(list_firsts)
    is_extr = is_extr .and. (list_firsts(1) == list_firsts(l))
  end do

  is_north = on_northern_hemisphere(zone)

  if (is_extr) then
    ! Default values
    if (is_north) then
      this%north_offset(k) = -1
      this%south_offset(k) = 0
    else
      this%south_offset(k) = -1
      this%north_offset(k) = 0
    end if

    to_north_ghost = not_first .and. (is_north_ghost(k-1, this))
    to_south_ghost = not_last .and. (is_south_ghost(k+1, this))

    ! From north ghost cells
    if (is_north .and.(is_north_ghost(k, this))) then
      this%south_offset(k) = vals(1)
    else if (.not. is_north .and. (is_south_ghost(k, this))) then
      this%north_offset(k) = vals(1)
    end if

    ! To north ghost cells
    if (is_north .and. to_north_ghost) then
      this%north_offset(k) = vals(2)
    else if (.not. is_north .and. to_south_ghost) then
      this%south_offset(k) = vals(2)
    end if

    ! To south ghost cells
    if (is_north .and. to_south_ghost) then
      this%south_offset(k) = vals(3)
    else if (.not. is_north .and. to_north_ghost) then
      this%north_offset(k) = vals(3)
    end if

    ! From south ghost cells
    if (is_north .and.(is_south_ghost(k, this))) then
      this%north_offset(k) = vals(4)
    else if (.not. is_north .and. (is_north_ghost(k, this))) then
      this%south_offset(k) = vals(4)
    end if

    is_done = .True.
  end if
end function

!-------------------------------------------------------------------------------
!> Set the number of east and west ghost cells at line k in the partition.
!> The number of ghost cells must be initialized to 0.
!> @parameter k : line inside the partition
!-------------------------------------------------------------------------------
subroutine set_nb_east_west_ghosts(i_lat, k, this, zone, list_firsts)
  type (Band_Grid) :: this
  integer, intent(in) :: i_lat, k, zone
  integer, intent(in) :: list_firsts(:)
  integer :: nb_cells, last

  ! No E/W ghost cells on north/south ghost line
  if (is_north_south_ghost(k, this)) return

  this%nb_ghosts_west(k) = nb_ghosts_side(.true., k, list_firsts, i_lat, zone, this)
  this%nb_ghosts_east(k) = nb_ghosts_side(.false., k, list_firsts, i_lat, zone, this)
  !if (k == 2) print *, "nb ghost west", this%nb_ghosts_west(k) 

  ! 1 cell for sector extremities
  if (list_firsts(k) == 1) then
    this%nb_ghosts_west(k) = 1
  end if
  nb_cells = nb_cells_lat_sector(i_lat)
  last = list_firsts(k) + this%nb_lon(k) - 1
  if (last == nb_cells) then
    this%nb_ghosts_east(k) = 1
  end if

  ! If no ghosts for north/south add 1 for east/west
  if (this%nb_ghosts_east(k) == 0) this%nb_ghosts_east(k) = 1
  if (this%nb_ghosts_west(k) == 0) this%nb_ghosts_west(k) = 1

end subroutine

!-------------------------------------------------------------------------------
!> Compute the number of ghosts on a side (either east or west) for north/south
!-------------------------------------------------------------------------------
function nb_ghosts_side(is_first, k, list_firsts, i_lat, zone, this) result (nb_ghosts)
  type (Band_Grid) :: this
  logical, intent(in):: is_first
  integer, intent(in) :: k, list_firsts(:)
  integer, intent(in) :: i_lat, zone
  integer :: extrema, nb_ghosts
  integer :: prev, next, extrema1, extrema2
  integer ::  cur
  integer :: i_start, i_end

  call interior_lat_indices(i_start, i_end, this)
  ! Special case : partition with 1 line
  if (i_start == i_end) then
    nb_ghosts = 1
    return
  end if

  cur = cur_cell_from_first(k, is_first, list_firsts, this)

  ! Compute the extrema by going to the previous (next) latitude and finding the
  ! south (north) neighbour
  extrema1 = int(UNDEFINED)
  if ( k < i_end) then
    next = cur_cell_from_first(k+1, is_first, list_firsts, this)
    extrema1 = north_extrema_cell(next, is_first, list_firsts, i_lat+1, zone)
  end if

  extrema2 = int(UNDEFINED)
  if ( k > i_start) then
    prev = cur_cell_from_first(k-1, is_first, list_firsts, this)
    extrema2 = south_extrema_cell(prev, is_first, list_firsts, i_lat-1, zone)
  end if

  if (extrema1 == UNDEFINED) then
    extrema = extrema2
  else if (extrema2 == UNDEFINED) then
    extrema = extrema1
  else
    if (is_first) then
      extrema = min(extrema1, extrema2)
    else
      extrema = max(extrema1, extrema2)
    end if
  end if

  if (is_first) then
    nb_ghosts = max(cur - extrema, 0)
  else
    nb_ghosts = max(extrema - cur, 0)
  end if

end function

!-------------------------------------------------------------------------------
!> Returns either the first or last cell
!-------------------------------------------------------------------------------
function cur_cell_from_first(k, is_first, list_firsts, this) result (cur)
  integer, intent(in) :: k, list_firsts(:)
  integer :: cur
  logical, intent(in) :: is_first
  type (Band_Grid) :: this

  if (is_first) then
    cur = list_firsts(k)
  else
    cur = list_firsts(k) + this%nb_lon(k) - 1
  end if
end function

!-------------------------------------------------------------------------------
!> Return either the leftmost or rightmost north neighbour
!> @param i_lat : current latitude 
!-------------------------------------------------------------------------------
function north_extrema_cell(cur, is_first, list_firsts, i_lat, zone) result (neighb)
  integer, intent(in) :: cur, i_lat, zone
  integer, intent(in) :: list_firsts(:)
  logical, intent(in) :: is_first
  integer :: neighb, mid, nb_cells
  integer :: nb_lat2

  ! Special case at the equator
  nb_lat2 = get_nb_lat2_Configuration()
  if (i_lat == nb_lat2 + 1) then
    neighb = cur
    return
  end if

  nb_cells = nb_cells_lat_sector(i_lat)
  mid = (nb_cells + 1)/2
  if (is_first) then
    neighb = left_north_cell(cur, mid, zone)
  else
    neighb = right_north_cell(cur, mid, nb_cells, zone)
  end if
end function

!-------------------------------------------------------------------------------
!> Return either the leftmost or rightmost north neighbour
!-------------------------------------------------------------------------------
function south_extrema_cell(cur, is_first, list_firsts, i_lat, zone) result (neighb)
  integer, intent(in) :: cur, i_lat, zone
  integer, intent(in) :: list_firsts(:)
  logical, intent(in) :: is_first
  integer :: neighb, mid, nb_cells
  integer :: nb_lat2

  ! Special case at the equator
  nb_lat2 = get_nb_lat2_Configuration()
  if (i_lat == nb_lat2) then
    neighb = cur
    return
  end if

  nb_cells = nb_cells_lat_sector(i_lat)
  mid = (nb_cells + 1)/2
  if (is_first) then
    neighb = left_south_cell(cur, mid, zone)
  else
    neighb = right_south_cell(cur, mid, nb_cells, zone)
  end if
end function


!-------------------------------------------------------------------------------
!> Correct nb_lon and lon_min for east/west ghost cells
!-------------------------------------------------------------------------------
subroutine add_east_west_ghosts(this)
  type (Band_Grid) :: this
  integer :: k

  do k = 1, this%nb_lat
    if (.not. is_north_south_ghost(k, this)) then
      this%nb_lon(k) = this%nb_lon(k) + this%nb_ghosts_west(k) + &
        this%nb_ghosts_east(k) 
      this%lon_min(k) = this%lon_min(k) - this%nb_ghosts_west(k)*this%dlon(k) 
    end if
  end do

end subroutine

!-------------------------------------------------------------------------------
!> Create a grid at position i on the leftover band. Also finishes to set some
!> @param height : partition height
!> @param width : partition width (except for the last one)
!> @param j : partition number on the band (no need to correct it)
!> @param n_left : number of partitions on the band
!> @param zone : zone indice
!-------------------------------------------------------------------------------
subroutine new_Band_grid_leftover(this, i_pos, j_pos, height, width, n_left, &
    zone, ghost_cells)
  type (Band_Grid) :: this
  integer, intent(in) :: ghost_cells, i_pos, j_pos
  integer, intent(in) :: width, zone
  integer :: height, nb_lat
  integer :: n_left, first_lat
  integer :: nb_lat2, k, nb_cells
  integer :: i_lat
  integer, allocatable :: list_firsts(:)

  double precision:: dlat(1)
  double precision  :: cur_dlon

  !write (*,*) "ghost cells for leftover", ghost_cells
  this%ghost_cells = ghost_cells

  nb_lat = height
  if (has_north_ghost_cells(this)) nb_lat = nb_lat + 1
  if (has_south_ghost_cells(this)) nb_lat = nb_lat + 1
  this%nb_lat = nb_lat

  call allocate_most_data(this, nb_lat)

  ! The current latitude (2 in the partition) is easy to compute as we are on
  ! the last band. We deduce the latitude line.
  nb_lat2 = nb_lat2_Global_grid()
  if (on_northern_hemisphere(zone)) then
    ! Nb lat -1 as there is one latitude line on the southern hemisphere
    first_lat = nb_lat2 - (nb_lat-1) + 1
  else
    ! On south hemisphere, we always have north ghost cells on the other side of
    ! the equator
    first_lat = nb_lat2
  end if

  ! First we define the normal grid, and then we add the ghost cells
  allocate (list_firsts(this%nb_lat))
  do k = 1, this%nb_lat
    i_lat =  first_lat + k - 1

    cur_dlon = get_dlon_Global_grid(i_lat)
    nb_cells = nb_cells_lat_sector(i_lat)

    call set_longitude_data_leftover(this, k, i_pos, j_pos, cur_dlon, &
      width, nb_cells, n_left, zone, list_firsts)
  end do

  !print *, "nb lon inside", this%nb_lon
  call set_data_for_neighbours(this, list_firsts, first_lat-1, zone)
  !print *, "nb ghosts", this%nb_ghosts_west, this%nb_ghosts_east 
  ! Update grid spec for ghost
  call add_east_west_ghosts(this)


  dlat = get_dlat_Global_grid()
  this%dlat = dlat
  ! Beware, on the northern hemisphere, south ghost cells are on the other side
  this%lat_min = (nb_lat -1)*dlat(1)
  ! On southern hemisphere, always starts from +dlat
  if (on_southern_hemisphere(zone)) this%lat_min = dlat(1)

  ! Initialize mpi type
  this%mpi_grid = mpi_datatype_null

  deallocate (list_firsts)
end subroutine

!-----------------------------------------------------------------------------
!> Convert global cell position to local
!> Float operations and comparison.
!> @parameter j_cur : cell in global grid (first sector)
!> @parameter sector : current sector (starts from 1)
!-----------------------------------------------------------------------------
function global_to_local_cell(j_cur, i, sector, this) result (j_pos)
  type (Band_Grid) :: this
  integer, intent(in) :: j_cur, i, sector
  logical :: res
  integer :: j_pos
  integer :: correct
  double precision :: center, lon_min, lon_max, dlon

  ! Find the center of the cell
  dlon = this%dlon(i)
  center = (j_cur - 0.5d0)*dlon

  ! Do not use west/east ghost cells at the global grid border
  lon_min = this%lon_min(i)
  correct = 0
  if (lon_min < 0) then
    lon_min = 0.
    correct = this%nb_ghosts_west(i)
  end if
  lon_max = lon_min + this%nb_lon(i)*dlon
  lon_max = min(lon_max, 360.)

  ! Longitude are translated to first sector
  lon_min = lon_min - (sector-1)*120.
  lon_max = lon_max - (sector-1)*120.
  !  if (j_cur == 1 .and. i == 2) then
  !    print *, "min max", lon_min, lon_max, "cur", center
  !    print *, "correct", correct
  !  end if
  ! Strict inegality as it is center (no need for precision)
  res = (center > lon_min .and. center < lon_max)
  j_pos = -1
  if (res) then
    j_pos = floor((center - lon_min)/dlon) + 1
    j_pos = j_pos + correct
  end if
end function

subroutine get_tracer_pointer(this, ptr, datatype, tracer)
  type(Band_Grid), target :: this
  integer, intent(in) :: datatype, tracer
  double precision, pointer :: ptr(:)

  call associate_tracer_pointer(this%tracers, ptr, datatype, tracer)
end subroutine

!-----------------------------------------------------------------------------
!> @param datatype : IS_ZWINDS, IS_MWINDS
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!-----------------------------------------------------------------------------
subroutine get_winds_pointer(this, ptr, datatype, array)
  type(Band_Grid), target :: this
  integer, intent(in) :: datatype, array
  double precision, pointer :: ptr(:)

  if (datatype == IS_ZWINDS) then
    if (array == IS_PREV) then
      ptr => this%zonal_winds_prev
    else if (array == IS_NEXT) then
      ptr => this%zonal_winds_next
    else
      ptr => this%zonal_winds
    end if

  else if (datatype == IS_MWINDS) then
    if (array == IS_PREV) then
      ptr => this%merid_winds_prev
    else if (array == IS_NEXT) then
      ptr => this%merid_winds_next
    else
      ptr => this%merid_winds
    end if


  else
    call print_error("Wrong datatype", "get_winds_pointer", fname_band)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Returns number of cells at latitude i (in the grid), without the ghost 
!> cells
!-----------------------------------------------------------------------------
function get_nb_lon(i, this) result(nb_lon)
  type (Band_Grid) :: this
  integer :: i, nb_lon
  integer :: nb_ghosts

  nb_lon = this%nb_lon(i)
  nb_ghosts = nb_ghost_cells_lat(i, this)
  nb_lon = nb_lon - nb_ghosts

end function


!-------------------------------------------------------------------------------
!> Tell us if there are ghost cells on the north boundary line
!-------------------------------------------------------------------------------
function has_north_ghost_cells(this) result(res)
  type (Band_Grid) :: this
  logical :: res

  res = mod(this%ghost_cells, 2) > 0
end function

!-------------------------------------------------------------------------------
!> Tell us if there are ghost cells on the east boundary line
!-------------------------------------------------------------------------------
function has_east_ghost_cells(this) result(res)
  type (Band_Grid) :: this
  logical :: res

  res = mod(this%ghost_cells / 2, 2) > 0
end function

!-------------------------------------------------------------------------------
!> Tell us if there are ghost cells on the south boundary line
!-------------------------------------------------------------------------------
function has_south_ghost_cells(this) result(res)
  type (Band_Grid) :: this
  logical :: res

  res = mod(this%ghost_cells / 4, 2) > 0
end function

!-------------------------------------------------------------------------------
!> Tell us if there are ghost cells on the west boundary line
!-------------------------------------------------------------------------------
function has_west_ghost_cells(this) result(res)
  type (Band_Grid) :: this
  logical :: res

  res = mod(this%ghost_cells / 8, 2) > 0
end function

!-------------------------------------------------------------------------------
!> Set longitude data (nb_lon, d_lon, lon_min). But does not contains east/west
!> ghost cells
!> Does not work for the leftover band. Another function does that : 
!> set_north_south_ghost in partition_class
!> Also, weast/east ghost cells are added later.
!> @param k : position in the array 
!> @param cur_dlon : current dlon
!> @param zone : partition zone
!> @param nb_resized : number of resized partition on the current band (only 
!> one one half)
!> @param width : width of a square partition (before resizing)
!> @param i_pos, j_pos : partition position in the array. j can be on any sector
!> @param nb_cells : number of cells on a sector for this latitude
!> @param list_firsts : store the first cell position (used later)
!-------------------------------------------------------------------------------
subroutine set_longitude_data(this, k, i_pos, j_pos, cur_dlon, nb_resized, width, &
    nb_cells, zone, list_firsts)
  type (Band_Grid) :: this
  integer, intent(in) :: k, width, nb_resized, nb_cells
  integer, intent(in) :: i_pos, j_pos, zone
  double precision, intent(in) :: cur_dlon
  integer :: list_firsts(:)
  integer :: j_corr, nb_parts
  integer :: first, length

  this%dlon(k) = cur_dlon

  nb_parts = nb_parts_on_band_sector(i_pos)
  j_corr = translate_jpos(j_pos, nb_parts)

  ! First we define the normal grid, and then we add the ghost cells
  ! Special case for north/south gohst cells
  if (is_north_south_ghost(k, this)) then
    call set_longitude_data_ghost(k, this, nb_resized, width, i_pos, j_corr, &
      nb_cells, zone, list_firsts)
  else
    ! For normal cells, we just call the function
    call find_first_length_cell(first, length, k, this, nb_resized, width, &
      i_pos, j_corr, nb_cells, zone)

    ! This does not contains east/west ghost cells
    this%nb_lon(k) = length
    ! Minus 1, as we want the west cell boundary
    this%lon_min(k) = (first - 1)*cur_dlon
    list_firsts(k) = first
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Compute the lon_min and nb_lon for ghost cells. The strategy is simple, we
!> compute the neighbours for the extremities
!> Only works on the first sector
!> @param k : latitude line in the grid
!> @param nb_cells : number of cells on a sector for this latitude
!> @param list_firsts : store the first cell position (used later)
!-------------------------------------------------------------------------------
subroutine set_longitude_data_ghost(k, this, nb_resized, width, i_pos, j_corr, &
    nb_cells, zone, list_firsts)
  type (Band_grid) :: this
  integer, intent(in) :: k, width, nb_resized, i_pos
  integer, intent(in) :: nb_cells, j_corr, zone
  integer :: list_firsts(:)
  integer :: k_next, nb_cells_next
  logical :: is_north
  integer :: first, length, last

  if (is_north_ghost(k, this)) then
    k_next = 2
    is_north = .True.
  else if (is_south_ghost(k, this)) then
    k_next = this%nb_lat - 1
    is_north = .False.
  end if
  nb_cells_next = nb_cells_next_lat_ghost(is_north, i_pos, zone, nb_cells)

  ! We find the extremities in global coordinates of normal cells
  call find_first_length_cell(first, length, k_next, this, nb_resized, width, &
    i_pos, j_corr, nb_cells_next, zone)
  last = first + length - 1

  ! And deduce the ghost cells
  call compute_longitude_data_ghost(k, first, last, this,  width, i_pos, &
    j_corr, nb_cells, zone, k_next, is_north, nb_cells_next, list_firsts)
end subroutine

!-------------------------------------------------------------------------------
!> Compute the lon_min and nb_lon for ghost cells on leftover band. 
!> Same strategy as "normal" bands.
!> Only works on the first sector
!> @param k : latitude line in the grid
!> @param nb_cells : number of cells on a sector for this latitude
!-------------------------------------------------------------------------------
subroutine set_longitude_data_ghost_leftover(k, this, width, i_pos, &
    j_corr, nb_cells, n_left, zone, list_firsts)
  type (Band_grid) :: this
  integer, intent(in) :: k, width, i_pos
  integer, intent(in) :: nb_cells, j_corr, zone
  integer, intent(in) :: n_left
  integer :: k_next, nb_cells_next
  logical :: is_north
  integer :: list_firsts(:)

  integer :: first, length, last

  if (is_north_ghost(k, this)) then
    k_next = 2
    is_north = .True.
  else if (is_south_ghost(k, this)) then
    k_next = this%nb_lat - 1
    is_north = .False.
  end if
  nb_cells_next = nb_cells_next_lat_ghost(is_north, i_pos, zone, nb_cells)

  ! We find the extremities in global coordinates
  call find_first_length_cell_leftover(first, length, this, width, &
    j_corr, nb_cells_next, n_left)
  last = first + length - 1

  call compute_longitude_data_ghost(k, first, last, this,  width, i_pos, &
    j_corr, nb_cells, zone, k_next, is_north, nb_cells_next, list_firsts)
end subroutine


!-------------------------------------------------------------------------------
!> Does not work for ghost cells ! (but used by its computation)
!> Only works on the first sector
!> @param j_pos : partition position on the band (first sector)
!-------------------------------------------------------------------------------
subroutine find_first_length_cell_leftover(first, length, this, width, j_pos,&
    nb_cells, n_left)
  type (Band_grid) :: this
  integer :: first, length
  integer, intent(in) :: width, j_pos, n_left
  integer, intent(in) :: nb_cells

  first = (j_pos-1)*width + 1

  ! Same width for all partitions, except the length
  if (j_pos /= n_left) then
    length = width
    ! length partition is a special case
  else
    length = nb_cells - (j_pos - 1)*width
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Does the true computing for nb_lon and lon_min (ghost cells).
!> Needs the neighbour first and last cells
!> @param first, last : neighbours inside the grid
!-------------------------------------------------------------------------------
subroutine compute_longitude_data_ghost(k, first, last, this,  width, i_pos, &
    j_corr, nb_cells, zone, k_next, is_north, nb_cells_next, list_firsts) 
  type (Band_Grid) :: this
  integer :: first, last, k
  integer :: list_firsts(:)
  integer, intent(in) :: nb_cells, zone
  integer :: mid, nb_cells_next
  integer :: k_next, width
  integer :: i_pos, j_corr
  integer :: ghost_first, ghost_last
  logical :: is_north

  ! Middle cell
  mid = (nb_cells_next + 1)/2

  if (is_north) then
    ! Special case for north ghost cells on the other side of the equator
    if (is_just_after_equator(i_pos) .and. k == 1) then
      ghost_first = first
      ghost_last = last
    else
      ghost_first = left_north_cell(first, mid, zone)
      ghost_last = right_north_cell(last, mid, nb_cells_next, zone)
    end if
  else
    ! Special case for south ghost cells on the other side of the equator
    if (is_just_before_equator(i_pos) .and. k == this%nb_lat) then
      ghost_first = first
      ghost_last = last
    else
      ghost_first = left_south_cell(first, mid, zone)
      ghost_last = right_south_cell(last, mid, nb_cells_next, zone)
    end if
  end if

  this%nb_lon(k) = ghost_last - ghost_first + 1
  this%lon_min(k) = (ghost_first-1)*this%dlon(k)
  list_firsts(k) = ghost_first

end subroutine

!-------------------------------------------------------------------------------
!> Set longitude data for leftover band
!-------------------------------------------------------------------------------
subroutine set_longitude_data_leftover(this, k, i_pos, j_pos, cur_dlon, &
    width, nb_cells, n_left, zone, list_firsts)
  type (Band_Grid) :: this
  integer, intent(in) :: k, width, nb_cells
  integer, intent(in) :: i_pos, j_pos, zone, n_left
  double precision, intent(in) :: cur_dlon
  integer :: j_corr, nb_parts
  integer :: first, length
  integer :: list_firsts(:)

  this%dlon(k) = cur_dlon

  nb_parts = nb_parts_on_band_sector(i_pos)
  j_corr = translate_jpos(j_pos, nb_parts)

  if (is_north_south_ghost(k, this)) then
    call set_longitude_data_ghost_leftover(k, this, width, i_pos, &
      j_corr, nb_cells, n_left, zone, list_firsts)
  else
    ! For normal cells, we just call the function
    call find_first_length_cell_leftover(first, length, this, width, &
      j_pos, nb_cells, n_left)
    this%nb_lon(k) = length
    this%lon_min(k) = (first-1)*cur_dlon
    list_firsts(k) = first
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Compute the first and last cell on latitude k. It does not take west/ghost
!> cells into account.
!> Does not work for leftover band
!> Does not work for ghost cells ! (but used by its computation)
!> Only works on the first sector
!> @param k : latitude on the grid
!> @param first, length : first cell on the latitude line, and the partition
!> length (first is NOT an offset)
!-------------------------------------------------------------------------------
subroutine find_first_length_cell(first, length, k, this, nb_resized, width, i_pos, &
    j_corr, nb_cells, zone)
  type(Band_grid) :: this
  integer, intent(in) :: nb_resized, width, j_corr, zone
  integer, intent(in) :: nb_cells, k, i_pos
  integer :: first, length
  integer :: j2, nb_parts, k2
  logical :: is_north
  integer :: min_left, offset

  ! There are nb_resized partitions of width+1 and j-nb_resized-1 of width
  offset = (nb_resized*(width + 1) + (j_corr - nb_resized - 1)*width)

  is_north = (on_northern_hemisphere(zone) .and. is_north_ghost(k, this))
  ! On south ghost cells for northern zoneere, the situation is the same as
  ! on the northern zoneere for north ghost cells
  is_north = is_north .or. (on_northern_hemisphere(zone) .and. &
    is_south_ghost(k, this)) 

  ! Middle
  if (middle_partition(i_pos, j_corr, zone)) then
    ! Northern hemisphere
    if (on_northern_hemisphere(zone)) then
      k2 = k
      if (has_north_ghost_cells(this)) k2 = k2 - 1
    else
      ! Inverted on southern hemisphere
      k2 = this%nb_lat - k + 1
      if (has_south_ghost_cells(this)) k2 = k2 - 1
    end if

    ! It is not an offset
    first = offset
    length = 2*k2-1
  else
    ! Left half
    j2 = j_corr
    min_left = (j_corr-1)*(width + 1) 

    ! Right half : we consider the symmetry of j. This changes the offset
    ! and min_left, as we start from nb_cells, instead of 1
    if (right_partition(i_pos, j_corr)) then
      ! The position is computed from the est boundary
      nb_parts = nb_parts_on_band_sector(i_pos)
      j2 = nb_parts - j_corr + 1
      ! As if we added 1 to j2
      offset = (nb_resized*(width + 1) + (j2 - nb_resized)*width)
      offset = nb_cells - offset
      min_left = j2*(width + 1)
      min_left = nb_cells - min_left

      ! For north ghost cells, the width does not correspond, so we must add 1
      if (is_north) offset = offset + 1
      if (is_north) min_left = min_left + 1
    end if

    ! If the grid has only resized partition to the left
    if (j2 < nb_resized + 1) then
      length = width + 1
      first = min_left
      ! Otherwise, there is a mix
    else
      length = width
      first = offset
    end if
  end if
  ! As it is not an offset 
  first = first + 1
end subroutine

!-------------------------------------------------------------------------------
!> Correct the partition extremities by adding east and west ghost
!>cells
!> Does not work for north/south ghost cells ! (but used by its computation)
!> Only works on the first sector
!-------------------------------------------------------------------------------
subroutine correct_with_east_west_ghost_cells(this)
  type (Band_grid) :: this
  integer :: nb, k
  integer :: i_start, i_end

  call interior_lat_indices(i_start, i_end, this)
  do k = i_start, i_end
    nb = this%nb_ghosts_west(k)
    this%lon_min(k) = this%lon_min(k) - this%dlon(k)*nb
    this%nb_lon(k) = this%nb_lon(k) + nb

    nb = this%nb_ghosts_east(k)
    this%nb_lon(k) = this%nb_lon(k) + nb
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Computes the total number of ghost cells laterally (east-west)
!-----------------------------------------------------------------------------
function nb_ghost_cells_lat(k, this) result(nb)
  type (Band_grid) :: this
  integer, intent(in) :: k
  integer :: nb

  nb = this%nb_ghosts_east(k)
  nb = nb + this%nb_ghosts_west(k)
end function

!-----------------------------------------------------------------------------
!> Returns true if the latitude line i is just after the equator
!> Use float !
!-----------------------------------------------------------------------------
function is_latitude_just_after_equator_float(i, this, zone) result (equator)
  type (Band_grid) :: this
  integer, intent(in) :: i, zone
  logical :: equator
  double precision :: lat

  lat = this%lat_min - (i - 0.5d0)*this%dlat(1)
  ! The middle is on the southern hemisphere, the prev is on the southern
  if (on_northern_hemisphere(zone)) then
    equator = (lat < 0) .and. (lat + this%dlat(1) > 0)
    ! Or the contrary
  else
    equator = (lat > 0) .and. (lat - this%dlat(1) < 0)
  end if


end function

!-------------------------------------------------------------------------------
!> Find the correction for minimal longitude. For non-north/south ghost cells,
!> the test is simple.
!> For north-south ghost cells, there are 2 cases.
!> @param k : position in the array
!> @param corr_east, corr_west : longitude correction
!> @param zone : partition zone
!> @param i, j : partition position
!-------------------------------------------------------------------------------
subroutine find_east_weast_corr(corr_east, corr_west, k, this, zone, i, j)
  type(Band_grid) :: this
  integer, intent(in) :: zone, k, i, j
  integer :: corr_east, corr_west

  ! Here we add the (eventual) west and east ghost cells
  corr_east = 0
  corr_west = 0
  ! Check if we are not on north/south ghost cells for boundary
  if (.not. is_north_south_ghost_for_boundary(k, this, i, j)) then
    ! East
    !if (has_east_ghost_cells(this)) corr_east = 1
    corr_east = 1
    ! West
    ! if (has_west_ghost_cells(this)) corr_west = 1
    corr_west = 1
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Compute the minimal longitude for south ghost cells on the leftover band
!> (south hemisphere). For that we compute the north neighbour of the west 
!> extremity.
!> @param : i_lat : latitude from the equator
!-------------------------------------------------------------------------------
subroutine find_lon_min_south_ghost_last(lon_min, grid, i, width, dlon, &
    i_lat, zone)
  type (Band_grid) :: grid
  integer, intent(in) :: i, width, i_lat
  integer, intent(in) :: zone
  integer :: mid
  integer :: cur, neighb
  integer :: nb_lat
  double precision :: lon_min(:)
  double precision, intent(in) :: dlon

  ! Find west extremity and its neighbour
  cur = (i - 1)*width + 1
  ! i_lat starts from the equator
  mid = get_nb_lat2_Global_grid() - i_lat + 1
  neighb = left_south_cell(cur, mid, zone)
  nb_lat = size(lon_min)
  lon_min(nb_lat) = (neighb - 1)*dlon
end subroutine


!-----------------------------------------------------------------------------
!> Allocate array for mpi types for borders 
!-----------------------------------------------------------------------------
subroutine allocate_mpi_borders(this, nb)
  type (Band_Grid) :: this
  integer :: nb

  allocate (this%mpi_ghost_cells(nb))
  allocate (this%mpi_interior_cells(nb))
  !write (*,*) "allocating ", nb, "mpi ghost cells"
  this%mpi_ghost_cells = mpi_datatype_null
  this%mpi_interior_cells = mpi_datatype_null
end subroutine

!-------------------------------------------------------------------------------
!> Set the MPI type for ratio or slope field as array have the same dimensions.
!> @param val : value to set
!-------------------------------------------------------------------------------
subroutine new_MPI_global_type(val, offsets, blocklens, this, nb)
  type (Band_Grid) :: this
  integer, intent(in) :: nb
  integer, allocatable :: blocklens(:), offsets(:)
  integer :: val

  call check_mpi_types(this%tracers, offsets, blocklens, 1)

  call new_indexed_mpi(nb, blocklens, offsets, val)

end subroutine

!-------------------------------------------------------------------------------
!> Create an MPI type for ghost cells and assing it to val
!-------------------------------------------------------------------------------
subroutine new_MPI_ghost_cells(val, offsets, blocklens, this, nb, direction, &
    start_prev, end_prev)
  type (Band_Grid) :: this
  integer, intent(in) :: nb
  integer, intent(in) :: direction
  integer, allocatable :: blocklens(:), offsets(:)
  logical :: cond
  integer :: val
  integer :: start_cur, end_cur
  integer :: start_prev, end_prev

  call check_mpi_types(this%tracers, offsets, blocklens, 0)

  ! Check the arrays does not overlap
  if (direction == NORTH .or. direction == SOUTH) then
    end_cur = offsets(1) + blocklens(1)
    start_cur = offsets(1) + 1

    ! Check at the junction
    cond = min(end_prev, end_cur) > max(start_prev, start_cur) - 1
    if (cond) then
      call print_error("Overlapping ghost cells", "new_MPI_ghost_cells", &
        fname_band)
    end if
    start_prev = offsets(1) + 1
    end_prev = start_prev + blocklens(1) - 1
  end if

  call new_MPI_global_type(val, offsets, blocklens, this, nb)

end subroutine

subroutine create_mpi_border_zwinds(this, direction)
  type (Band_grid) :: this
  integer, intent(in) :: direction
  integer, allocatable :: blocklens(:), offsets(:)
  integer :: i, i_start, i_end, prev, k, n
  integer :: val

  if (this%mpi_zwinds_west /= mpi_datatype_null .and. &
    this%mpi_zwinds_east /= mpi_datatype_null) return

  call interior_lat_indices(i_start, i_end, this)
  n = i_end - i_start + 1
  allocate(blocklens(n))
  allocate(offsets(n))

  k = 1
  prev = 0
  do i = i_start, i_end

    if (direction == WEST) then
      offsets(k) = prev
      prev = prev + nb_lon_Band_grid(i, this) + 1
    else
      prev = prev + nb_lon_Band_grid(i, this) 
      offsets(k) = prev
      prev = prev + 1
    end if

    k = k + 1
  end do

  blocklens = 1
  !print *, "offset", offsets(1:5)
  call new_MPI_global_type(val, offsets, blocklens, this, n)
  if (direction == WEST) then
    this%mpi_zwinds_west = val
  else
    this%mpi_zwinds_east = val
  end if
  deallocate(blocklens)
  deallocate(offsets)

end subroutine


!-----------------------------------------------------------------------------
!> Update/search starting indice of the data array for the ghost cells, with 
!> the following disposition. North/south are converted it from the local 
!> position (on a band) to a global (for the grid), while east/west are created
!>      N |           |
!>      -----------------
!>      W |           | E
!>        |           |
!>        |           |
!>      -----------------
!>      S |           |
!-----------------------------------------------------------------------------
subroutine start_indice_ghost(start, this, side) 
  type (Band_Grid) :: this
  integer :: start
  character(*) :: side

  ! North is untouched
  select case (side)
  case ("west")
    start = this%nb_lon(1) + 1
  case ("east")
    start = sum(this%nb_lon(1:2))
    ! South is updated
  case ("south")
    start = sum(this%nb_lon(1:this%nb_lat-1)) + start
  end select
end subroutine

!-----------------------------------------------------------------------------
!> Get starting indice of the data array for the border inside the partition
!-----------------------------------------------------------------------------
subroutine get_start_border_in(this, start, side)
  type (Band_Grid) :: this
  integer :: start
  character(*) :: side

  start = 1
  select case (side)
  case ("north")
    start = 1
  case ("east")
    start = this%nb_lon(1) + 1
  case ("south")
    start = this%nb_lon(1) + this%dlat(1) + 1
  case ("west")
    start = this%nb_lon(1) + this%nb_lat + this%nb_lon(this%nb_lat) + 1
  case default
    call print_error("Error : not a border", "get_start_border_in", &
      fname_band)
  end select
end subroutine

!-----------------------------------------------------------------------------
!> Creates and commit an new indexed mpi type for an array of double precision
!> @param mpi_type : storage in the class for the new MPI type
!-----------------------------------------------------------------------------
subroutine new_indexed_MPI(nb, blocklens, offsets, mpi_type)
  integer :: blocklens(:), offsets(:)
  integer :: nb, mpi_type, ierr

  call mpi_type_indexed(nb, blocklens, offsets, mpi_double_precision, mpi_type, ierr)
  call check_mpi_error(ierr, "create indexed type", "new_indexed_MPI", &
    fname_band)

  call mpi_type_commit(mpi_type, ierr)
  call check_mpi_error(ierr, "commit type", "new_indexed_MPI", &
    fname_band)

end subroutine

!-----------------------------------------------------------------------------
!> Find the center of the grid, that is the middle latitude and the middle
!> longitude at this latitude
!-----------------------------------------------------------------------------
subroutine get_center_Band_grid(this,mid_lat,mid_lon)
  type(Band_Grid) :: this
  integer :: i_min
  double precision :: mid_lat, mid_lon

  mid_lat = this%lat_min - 0.5d0*(this%nb_lat * this%dlat(1))
  ! Not exactly the middle
  i_min = this%nb_lat/2
  if (i_min == 0) i_min = 1
  mid_lon = this%lon_min(i_min) + 0.5d0*(this%nb_lon(i_min) * this%dlon(i_min))
end subroutine


!-----------------------------------------------------------------------------
!> MPI user-defined datatype : create the datatype from a grid
!> The grid must have all of its array allocated !
!-----------------------------------------------------------------------------
subroutine new_MPI_Band_grid(grid) 
  type(Band_Grid) :: grid
  integer :: ierr
  ! Length, displacement, type
  integer, parameter :: asize = 12
  integer :: mtype(asize), mlength(asize)
  integer(kind=mpi_address_kind) :: mlocation(asize)
  integer(kind=mpi_address_kind) :: start
  integer :: i, n_lat
  character(*), parameter :: func_name = "new_MPI_Band_grid"
  integer(kind=mpi_address_kind) :: loc

  !call print_Band_grid(grid)
  ! Uses addresses as we have allocatable arrays. 
  call mpi_get_address(grid, start, ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_band)

  call get_address(grid%dlat, mlocation(1), func_name, fname_band)
  call get_address(grid%dlon, mlocation(2), func_name, fname_band)

  call mpi_get_address(grid%dlon, loc, ierr)

  call get_address(grid%first_cell, mlocation(3), func_name, fname_band)
  call get_address(grid%ghost_cells, mlocation(4), func_name, fname_band)
  call get_address(grid%lat_min, mlocation(5), func_name, fname_band)
  call get_address(grid%lon_min, mlocation(6), func_name, fname_band)
  call get_address(grid%nb_lat, mlocation(7), func_name, fname_band)
  call get_address(grid%nb_lon, mlocation(8), func_name, fname_band)
  call get_address(grid%nb_ghosts_east, mlocation(9), func_name, fname_band)
  call get_address(grid%nb_ghosts_west, mlocation(10), func_name, fname_band)
  call get_address(grid%north_offset, mlocation(11), func_name, fname_band)
  call get_address(grid%south_offset, mlocation(12), func_name, fname_band)
  ! Concentration and winds are set by each process, no need to exchange

  do i = 1, size(mlength)
    mlocation(i) = mlocation(i) - start!mlocation(1)
  end do

  n_lat = size(grid%nb_lon)
  mlength(1) = size(grid%dlat)
  mlength(2) = size(grid%dlon)
  mlength(3) = size(grid%first_cell)
  mlength(4) = 1
  mlength(5) = 1
  mlength(6) = size(grid%lon_min)
  mlength(7) = 1
  mlength(8) = size(grid%nb_lon)
  mlength(9) = size(grid%nb_ghosts_east)
  mlength(10) = size(grid%nb_ghosts_west)
  mlength(11) = size(grid%north_offset)
  mlength(12) = size(grid%south_offset) 

  mtype(1) = mpi_double_precision
  mtype(2) = mpi_double_precision
  mtype(3) = mpi_integer
  mtype(4) = mpi_integer
  mtype(5) = mpi_double_precision
  mtype(6) = mpi_double_precision
  mtype(7) = mpi_integer
  mtype(8) = mpi_integer
  mtype(9) = mpi_integer
  mtype(10) = mpi_integer
  mtype(11) = mpi_integer
  mtype(12) = mpi_integer

  call mpi_type_create_struct(size(mlength), mlength, mlocation, mtype, grid%mpi_grid,&
    ierr)
  call check_mpi_error(ierr, "create struct type", "new_indexed_MPI", &
    fname_band)

  call mpi_type_commit(grid%mpi_grid, ierr)
  call check_mpi_error(ierr, "commit type", "new_indexed_MPI", &
    fname_band)

end subroutine

!-----------------------------------------------------------------------------
!> Print grid features (for debug)
!-----------------------------------------------------------------------------
subroutine print_Band_grid(this)
  type(Band_Grid) :: this
  integer :: i, i_max

  if (.not. allocated(this%dlat)) then
    write (*,*) "Grid not allocated"
  else
    write (*,*) "Nb lat :",this%nb_lat
    write (*,*) "lat min :",this%lat_min
    write (*,*) "dlat :",this%dlat(1)
    write (*,*) "ghost cells :",this%ghost_cells
    write (*,*) "(i,lon_min,nb lon,dlon, ghosts_w,ghosts_e, offset_n, offset_s)"
    i_max = min(5, this%nb_lat)
    do i=1,i_max!this%nb_lat
      write (*,'(i4,f10.5,i5,f10.5)', advance="no") i,this%lon_min(i),this%nb_lon(i),this%dlon(i)

      ! Conditional print
      call print_nb_ghosts(i, this)
      call print_offsets(i, this)

    end do
    if (i_max < this%nb_lat) then
      write (*,*) "..."
      i = this%nb_lat
      write (*,'(i4,f10.5,i5,f10.5)', advance="no") i,this%lon_min(i),this%nb_lon(i),this%dlon(i)
      call print_nb_ghosts(i, this)
      call print_offsets(i, this)
    end if
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Print conditionnaly the current number of east/west ghost cells at line i
!-----------------------------------------------------------------------------
subroutine print_nb_ghosts(i, this)
  type(Band_Grid) :: this
  integer :: i

  !if (.not. is_north_south_ghost(i, this)) then
  write (*,'(2i4)', advance="no") this%nb_ghosts_west(i), this%nb_ghosts_east(i)
  !else
  ! Blank space
  !write (*,'(8x)', advance="no")
  !end if
end subroutine

!-----------------------------------------------------------------------------
!> Print conditionnaly the current north/south offset at line i
!-----------------------------------------------------------------------------
subroutine print_offsets(i, this)
  type(Band_Grid) :: this
  integer :: i

  !if (i > 1) then
  write (*,'(i4)', advance="no") this%north_offset(i)
  !else
  !write (*,'(4x)', advance="no") 
  !end if

  !if (i < this%nb_lat) then
  write (*,'(i4)') this%south_offset(i)
  !else
  !write (*,'(4x)') 
  !end if
end subroutine

!-----------------------------------------------------------------------------
!> Write the cell ratio and coordinates of the cell center
!> dlat is supposed constant
!> Writing is done from the north
!-----------------------------------------------------------------------------
subroutine write_ratio_Band_grid(this, tracer, id_file)
  type(Band_Grid) :: this
  integer, intent(in) :: id_file, tracer
  integer :: i, j
  double precision :: cur_lat, cur_lon
  double precision :: dlat, dlon, cur
  character(10) :: fformat
  character(10) :: fformat2
  integer :: i_start, i_end
  integer :: j_start, j_end

  ! Does not print ghost cells
  call interior_lat_indices(i_start, i_end, this)
  cur_lat = get_lat_min(this)
  dlat = get_dlat(this)

  fformat = '(2'//DBLE_FORMAT//')'
  fformat2 = '('//DBLE_FORMAT//',a)'

  do i = i_start, i_end

    dlon = this%dlon(i)
    call interior_lon_indices(j_start, j_end, this, i)
    cur_lon = get_lon_min(i,this)

    ! Write a cell : boundary and data
    do j = j_start, j_end

      ! Write the cell data
      cur = get_cell_ratio(i, j, tracer, this)
      write (id_file, fformat2,advance="no") cur

      write (id_file,fformat) cur_lat - 0.5*dlat, cur_lon + 0.5*dlon
      cur_lon = cur_lon + dlon
    end do
    cur_lat = cur_lat - dlat
  end do

end subroutine

!-----------------------------------------------------------------------------
!> Write the partition id and the cell corners
!> @param id_part : partition id
!-----------------------------------------------------------------------------
subroutine write_id_Band_grid(this, id_part, id_file)
  type(Band_Grid) :: this
  integer, intent(in) :: id_part, id_file
  integer :: i, j
  double precision :: cur_lat, cur_lon
  double precision :: dlat, dlon
  character(10) :: fformat
  integer :: i_start, i_end
  integer :: j_start, j_end

  ! Does not print ghost cells
  call interior_lat_indices(i_start, i_end, this)
  cur_lat = get_lat_min(this)
  dlat = get_dlat(this)

  fformat = '(8'//DBLE_FORMAT//')'

  do i = i_start, i_end

    dlon = this%dlon(i)
    call interior_lon_indices(j_start, j_end, this, i)
    cur_lon = get_lon_min(i,this)

    ! Write a cell : boundary and data
    do j = j_start, j_end

      ! Write the partition id
      write (id_file, '(i10)',advance="no") id_part

      write (id_file,fformat) cur_lat, cur_lon, &
        cur_lat, cur_lon + dlon, &
        cur_lat - dlat, cur_lon + dlon, &
        cur_lat - dlat, cur_lon

      cur_lon = cur_lon + dlon

    end do
    cur_lat = cur_lat - dlat
  end do

end subroutine



!-----------------------------------------------------------------------------
!> Give the starting indices and minimal latitude for the interior cells (i.e 
!> without the ghost cells)
!-----------------------------------------------------------------------------
subroutine interior_lat_indices(i_start, i_end, this)
  type (Band_Grid) :: this
  integer :: i_start, i_end

  i_start = first_interior_lat(this) 
  i_end = last_interior_lat(this) 
end subroutine

!-----------------------------------------------------------------------------
!> First latitude indice for interior cells
!-----------------------------------------------------------------------------
function first_interior_lat(this) result(i_start)
  type (Band_Grid) :: this
  integer :: i_start

  ! Does not print ghost cells
  i_start = 1
  if (has_north_ghost_cells(this)) i_start = i_start + 1

end function

!-----------------------------------------------------------------------------
!> Last latitude indice for interior cells
!-----------------------------------------------------------------------------
function last_interior_lat(this) result(i_end)
  type (Band_Grid) :: this
  integer :: i_end, i_start

  i_start = first_interior_lat(this) 
  ! Cannot do it if only one band
  i_end = this%nb_lat
  if (i_end > i_start) then
    if (has_south_ghost_cells(this)) i_end = i_end - 1
  end if

end function


!-------------------------------------------------------------------------------
!> Compute interface length between cell (i,j) and (i_neighb, j_neighb) in
!> degree
!> It is a length on the sphere, not a difference of longitudes
!-------------------------------------------------------------------------------
function cell_interface_length(i, j, i_neighb, j_neighb, this) result(dlon)
  type (Band_grid) :: this
  integer, intent(in) :: i, j, i_neighb, j_neighb
  double precision :: lat, dlon

  dlon = min(cell_east_lon(i, j, this), cell_east_lon(i_neighb, j_neighb, this))
  dlon = dlon - max(cell_west_lon(i, j, this), cell_west_lon(i_neighb, j_neighb, this))

  if (i_neighb > i) then
    lat = cell_south_lat(i, this)
  else
    lat = cell_north_lat(i, this)
  end if

  dlon = dlon*cos(lat*pi/180.d0)
end function


!-------------------------------------------------------------------------------
!> Compute middle of the interface between cell (i,j) and (i_neighb, j_neighb)
!-------------------------------------------------------------------------------
function cell_interface_middle(i, j, i_neighb, j_neighb, this) result(middle)
  type (Band_grid) :: this
  integer, intent(in) :: i, j, i_neighb, j_neighb
  double precision :: lon, dlon, middle

  lon = max(cell_west_lon(i, j, this), cell_west_lon(i_neighb, j_neighb, this))
  dlon = min(cell_east_lon(i, j, this), cell_east_lon(i_neighb, j_neighb, this))
  dlon = dlon - lon

  middle = lon + 0.5d0*dlon

end function


!-----------------------------------------------------------------------------
!> Minimal north latitude (without ghost cells)
!-----------------------------------------------------------------------------
function get_lat_min(this) result (lat_min)
  type (Band_Grid) :: this
  double precision :: lat_min

  lat_min = this%lat_min
  if (has_north_ghost_cells(this)) lat_min = lat_min - this%dlat(1)
end function

!-----------------------------------------------------------------------------
!> Latitude (in degrees) of the cell center  (with ghost cells)
!-----------------------------------------------------------------------------
function cell_center_lat(i, this) result (lat)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  double precision :: lat, dlat

  dlat = get_dlat(this)
  lat = this%lat_min - (i-0.5d0)*dlat
end function

!-----------------------------------------------------------------------------
!> Latitude (in degrees) of the cell south border (with ghost cells)
!-----------------------------------------------------------------------------
function cell_south_lat(i, this) result (lat)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  double precision :: lat, dlat

  dlat = get_dlat(this)
  lat = this%lat_min - i*dlat
end function

!-----------------------------------------------------------------------------
!> Latitude (in degrees) of the cell north border (with ghost cells)
!-----------------------------------------------------------------------------
function cell_north_lat(i, this) result (lat)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  double precision :: lat, dlat

  dlat = get_dlat(this)
  lat = this%lat_min - (i-1)*dlat
end function

!-----------------------------------------------------------------------------
!> Longitude (in degrees) of the cell center  (with ghost cells)
!-----------------------------------------------------------------------------
function cell_center_lon(i, j, this) result (lon)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j
  double precision :: lon, dlon

  dlon = get_dlon(i, this)
  lon = this%lon_min(i) + (j-0.5d0)*dlon
end function


!-----------------------------------------------------------------------------
!> Give the starting indices  for the interior cells (i.e without the ghost
!> cells). If we are on a ghost line, return all cells.
!> @param k : latitude line
!-----------------------------------------------------------------------------
subroutine interior_lon_indices(j_start, j_end, this, k)
  type (Band_Grid) :: this
  integer :: j_start, j_end, k

  j_start = first_interior_lon(k, this)
  j_end = last_interior_lon(k, this)
end subroutine

function first_interior_lon(k, this) result(j_start)
  type (Band_Grid) :: this
  integer :: j_start, k

  j_start = this%nb_ghosts_west(k) + 1
end function

function last_interior_lon(k, this) result(j_end)
  type (Band_Grid) :: this
  integer :: j_end, k

  j_end = this%nb_lon(k) - this%nb_ghosts_east(k)
end function

!-----------------------------------------------------------------------------
!> Give the number of interior cells at line i inside the grid (without the 
!> ghost cells)
!-----------------------------------------------------------------------------
function nb_lon_Band_grid(i, this) result (nb_lon)
  type (Band_Grid) :: this
  integer :: nb_lon, i

  nb_lon = this%nb_lon(i)
  nb_lon = nb_lon - nb_ghost_cells_lat(i, this)
end function

!-----------------------------------------------------------------------------
!> Give the number of cells at line i inside the grid (with the 
!> ghost cells)
!-----------------------------------------------------------------------------
function true_nb_lon_Band_grid(this) result (nb_lon)
  type (Band_Grid) :: this
  integer :: nb_lon, i

  nb_lon = this%nb_lon(i)
end function

function get_nb_ghosts_west(i, this) result(nb)
  type (Band_Grid) :: this
  integer :: nb, i

  nb = this%nb_ghosts_west(i)
end function

function get_nb_ghosts_east(i, this) result(nb)
  type (Band_Grid) :: this
  integer :: nb, i

  nb = this%nb_ghosts_east(i)
end function

function get_nb_ghosts_north(this) result(nb)
  type (Band_Grid) :: this
  integer :: nb

  nb = 0
  if (has_north_ghost_cells(this)) then
    nb = this%nb_lon(1)
  end if
end function

function get_nb_ghosts_south(this) result(nb)
  type (Band_Grid) :: this
  integer :: nb

  nb = 0
  if (has_south_ghost_cells(this)) then
    nb = this%nb_lon(this%nb_lat)
  end if
end function

!-----------------------------------------------------------------------------
!> Give the number of interior latitude inside the grid (without the 
!> ghost cells)
!-----------------------------------------------------------------------------
function get_nb_lat_Band_grid(this) result (nb_lat)
  type (Band_Grid) :: this
  integer :: nb_lat

  nb_lat = this%nb_lat
  if (has_north_ghost_cells(this)) nb_lat = nb_lat - 1
  if (has_south_ghost_cells(this)) nb_lat = nb_lat - 1
end function

!-----------------------------------------------------------------------------
!> Give the first interior latitude inside the grid (without the 
!> ghost cells)
!-----------------------------------------------------------------------------
function first_lat_Band_grid(this) result (i_lat)
  type (Band_Grid) :: this
  integer :: i_lat

  i_lat = 1
  if (has_north_ghost_cells(this)) i_lat = i_lat + 1
end function

!-----------------------------------------------------------------------------
!> Give the last interior latitude inside the grid (without the 
!> ghost cells)
!-----------------------------------------------------------------------------
function last_lat_Band_grid(this) result (i_lat)
  type (Band_Grid) :: this
  integer :: i_lat

  i_lat = this%nb_lat
  if (has_south_ghost_cells(this)) i_lat = i_lat - 1
end function

!-----------------------------------------------------------------------------
!> Give the first interior longitude at latitude i inside the grid (without
!> the ghost cells)
!-----------------------------------------------------------------------------
function first_lon_Band_grid(i, this) result (i_lon)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  integer :: i_lon

  i_lon = 1
  i_lon = i_lon + this%nb_ghosts_west(i)
end function

!-----------------------------------------------------------------------------
!> Give the last interior longitude at latitude i inside the grid (without
!> the ghost cells)
!-----------------------------------------------------------------------------
function last_lon_Band_grid(i, this) result (i_lon)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  integer :: i_lon

  i_lon = this%nb_lon(i)
  i_lon = i_lon - this%nb_ghosts_east(i)
end function


!-----------------------------------------------------------------------------
! Returns the longitude number for the first cell of the i-th line in the grid
! Number must be > 0
!-----------------------------------------------------------------------------
function pos_lon_line_Band_grid(this, i) result(i_lon)
  type(Band_Grid) :: this
  integer :: i, i_lon

  i_lon = int(this%lon_min(i)/this%dlon(i)) + 1
end function

!-------------------------------------------------------------------------------
!> Returns the north neighbours for cells (i,j) inside the partition. Returns 
!> the position inside the grid. Used by north_neighbour_cell_Partition.
!> @param i : latitude line inside the partition
!> @param j : longitude on the latitude line inside the partition
!-------------------------------------------------------------------------------
subroutine north_neighbour_cell_Band_grid(j_start, j_end, i, j, this, i_neighb,&
    offset, i_pos, j_pos, zone)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, offset
  integer, intent(in) :: i_neighb, zone
  integer, intent(in) :: i_pos, j_pos
  integer :: offset2, j_start, j_end, mid
  logical :: done

  mid = middle_cell(i, this, zone)
  offset2 = offset
  done = .False.

  ! Default is [j, j+1]
  ! We use the offset to correct, according to the first cells.

  if (mid > 0) then
    ! If the middle is on the current partition
    ! At the middle, only 1 neighbour
    if (j == mid) then 
      j_start = j + offset2
      j_end = j_start
      done = .True.
      ! After the middle, we must modify the offset
    else if (j > mid) then
      offset2 = offset2 - 1
      j_start = j + offset2
      j_end = j + offset2 + 1
      done = .True.
    end if
  end if

  if (.not. done) then
    ! Default (middle is outside)
    j_start = j + offset2
    j_end = j + offset2 + 1
  end if

  ! With this strategy, some cells may be outside, so we take the first (or
  ! last) inside
  !  if (i == 41  .and. j == 1) then
  !    write (*,*) "befor", j_start, j_end, "mid", mid, "j", j, "offset", offset2
  !!    write (*,*) "i ineighb", i, i_neighb
  !  end if

  call correct_indices(i_neighb, j_start, this, i, i_pos, j_pos)
  call correct_indices(i_neighb, j_end, this, i, i_pos, j_pos)

  !  if (i == 1 .and. j == 1 .and. i_pos == 2 .and. j_pos == 3) then
  !    write (*,*) "after corr", j_start, j_end
  !  end if
end subroutine

!-------------------------------------------------------------------------------
!> Returns the south neighbours for cells (i,j) inside the partition. Returns also
!> the ghost cells.
!> Used by south_neighbour_cell_Partition
!> @param i_pos, j_pos : partition position in the matrix
!> @param i_neighb : neighbour latitude line (not always i+1, as we can be on
!> the south hemisphere)
!-------------------------------------------------------------------------------
subroutine south_neighbour_cell_Band_grid(j_start, j_end, i, j, this, &
    i_neighb, offset, i_pos, j_pos, zone)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, offset
  integer, intent(in) :: i_neighb, zone
  integer, intent(in) :: i_pos, j_pos
  integer :: offset2
  integer :: j_start, j_end
  integer :: mid
  logical :: done

  mid = middle_cell(i, this, zone)
  offset2 = offset
  done = .False.

  ! Default is [j, j+1]
  ! We use the offset to correct, according to the first cells.
  if (mid > 0) then
    ! If the middle is on the current partition
    ! At the middle, only 1 neighbour
    if (j == mid) then 
      j_start = j + offset2
      j_end = j_start + 2
      done = .True.
      ! After the middle, we must modify the offset
    else if (j > mid) then
      offset2 = offset2 + 1
      j_start = j + offset2
      j_end = j + offset2 + 1
      done = .True.
    end if
  end if

  if (.not. done) then
    ! Default (middle is outside)
    j_start = j + offset2
    j_end = j + offset2 + 1
  end if

  !    if (i == 1  .and. j == 89 .and. i_pos == 2 .and. j_pos == 1) then
  !      write (*,*) "befor", j_start, j_end, "mid", mid, "j", j, "offset", offset2
  !      write (*,*) "i ineighb", i, i_neighb
  !    end if

  ! With this strategy, some cells may be outside, so we take the first (or
  ! last) inside
  call correct_indices(i_neighb, j_start, this, i, i_pos, j_pos)
  call correct_indices(i_neighb, j_end, this, i, i_pos, j_pos)
  !  if (i == 1 .and. j == 1 .and. i_pos == 2 .and. j_pos == 3) then
  !    write (*,*) "after", j_start, j_end
  !  end if
end subroutine

!-------------------------------------------------------------------------------
!> If j is outside the grid, takes the closest indice
!-------------------------------------------------------------------------------
subroutine correct_indices(i_neighb, j, this, i, i_pos, j_pos)
  type (Band_Grid) :: this
  integer, intent(in) :: i_neighb
  integer, intent(in) :: i, i_pos, j_pos
  integer :: j

  !call correct_indices_ghost(j, i, this, i_neighb, i_pos, j_pos)
  call correct_indices_other(i_neighb, j, this, i_pos, j_pos)

end subroutine

!-----------------------------------------------------------------------------
!> Neighbour indices must be in 1 and the number of cells.
!-----------------------------------------------------------------------------
subroutine correct_indices_other(i_neighb, j, this, i_pos, j_pos)
  type (Band_Grid) :: this
  integer, intent(in) :: i_neighb
  integer, intent(in) :: i_pos, j_pos
  integer :: j

  ! New version, only one condition. A single partition case does not need ghost
  ! cells
  j = max(j, 1)
  ! For a single partition, we work only on the first sector
  if (has_single_partition()) then
    j = min(j, this%nb_lon(i_neighb)/3)
  else
    j = min(j, this%nb_lon(i_neighb))
  end if

  if (has_single_partition()) return

  ! We only have to check the sector boundaries
  if (is_second_cell_sector(i_neighb, j, i_pos, j_pos, this)) then
    j = max(j, 2)
  else if (is_beforelast_cell_sector(i_neighb, j, i_pos, j_pos, this)) then
    j = min(j, this%nb_lon(i_neighb)-1)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> If we are on north/south ghost cells and the partition is on a boundary, we
!> must adapt the indices (as there are no east/west ghost cells on north/south
!> ghost cells)
!-----------------------------------------------------------------------------
subroutine correct_indices_ghost(j, i, this, i_neighb, i_pos, j_pos)
  type (Band_Grid) :: this
  integer, intent(in) :: i_neighb
  integer, intent(in) :: i, i_pos, j_pos
  integer :: j
  logical :: cond

  ! From ghost cells, we need to add one (as there are less cells on the
  ! ghost line)
  if (is_north_south_ghost_for_boundary(i, this, i_pos, j_pos)) then
    ! North ghost cells to south
    cond = (i == 1 .and. i_neighb > i)
    ! south ghost cells to north
    cond = cond .or. (i == this%nb_lat .and. i_neighb < i)
    if (cond) j = j + 1
  end if

  ! To ghost cells, need to remove one
  if (is_north_south_ghost_for_boundary(i_neighb, this, i_pos, j_pos)) then
    cond = (i == 2 .and. i_neighb < i)
    ! South to south ghost cells 
    cond = cond .or. (i == this%nb_lat-1 .and. i_neighb > i)

    if (cond) j = j - 1
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Returns true if the cell is the westmost on the current sector
!>@param i, j : cell position in the partition
!>@param i_pos, j_pos : partition position
!-----------------------------------------------------------------------------
function is_westmost_on_sector(i, j, i_pos, j_pos, this) result(res)
  type (Band_grid) :: this
  integer, intent(in) :: i, j, i_pos, j_pos
  integer :: l, nb_parts
  logical ::res

  res = .False.
  ! Does not work for north and south ghost cells
  if (has_north_ghost_cells(this) .and. i == 1) return
  if (has_south_ghost_cells(this) .and. i == this%nb_lat) return

  ! The first true element is always the second (except for ghost cells, but we
  ! have eliminated this case)
  if (j == 2) then
    ! The partition must also be the first on the current band
    nb_parts = nb_parts_on_band_sector(i_pos)
    do l = 1, 3
      res = (j_pos == nb_parts*(l-1) + 1)
      if (res) return
    end do
  end if
end function

!-----------------------------------------------------------------------------
!> Returns false if we can adjust north or south ghost cells. We cannot for a
!> boundary partition (i.e the first on the sector).
!> @param i_pos, j_pos : partition position in the matrix array
!-----------------------------------------------------------------------------
function is_north_south_ghost_for_boundary(i, this, i_pos, j_pos) result &
    (res)
  type (Band_grid) :: this
  integer, intent(in) :: i
  integer, intent(in) :: i_pos, j_pos
  logical :: res

  res = .False.

  ! Only when left boundary we must increase of 1
  if (is_first_part_sector(i_pos, j_pos)) then
    res = is_north_south_ghost(i, this)
  end if
end function

!-----------------------------------------------------------------------------
!> Returns true if the current partition is the first on the sector
!> @param i_pos, j_pos : partition position in the matrix array
!-----------------------------------------------------------------------------
function is_first_part_sector(i_pos, j_pos) result (res)
  integer :: nb_parts
  integer, intent(in) :: i_pos, j_pos
  logical :: res

  nb_parts = nb_parts_on_band_sector(i_pos)
  res = (j_pos == 1 .or. j_pos == nb_parts+1 .or. j_pos == 2*nb_parts+1)

end function

!-----------------------------------------------------------------------------
!> Returns true if the current cell is the second on the sector. Returns false
!> for ghost cells.
!> @param i, j : cell position in the grid
!-----------------------------------------------------------------------------
function is_second_cell_sector(i, j, i_pos, j_pos, this) result(res) 
  type (Band_grid) :: this
  logical :: res
  integer, intent(in) :: i_pos, j_pos
  integer, intent(in) :: i, j

  res = .False.
  if (is_first_part_sector(i_pos, j_pos)) then
    if (.not. is_north_south_ghost(i, this)) then
      res = (j < 3)
    end if
  end if
  !  if (i_pos == 1 .and. j_pos == 1 .and. i == 2 .and. j == 2) then
  !    write (*,*) "milou is first", is_first_part_sector(i_pos, j_pos)
  !    write (*,*) "is north south", is_north_south_ghost(i, this)
  !    write (*,*) "res", res
  !  end if

end function 

!-----------------------------------------------------------------------------
!> Returns true if the current cell is before the last on the sector. Returns 
!>_false for ghost cells.
!> @param i, j : cell position in the grid
!-----------------------------------------------------------------------------
function is_beforelast_cell_sector(i, j, i_pos, j_pos, this) result(res) 
  type (Band_grid) :: this
  logical :: res
  integer, intent(in) :: i_pos, j_pos
  integer, intent(in) :: i, j
  integer :: nb_lon

  res = .False.
  if (is_last_part_sector(i_pos, j_pos)) then
    nb_lon = this%nb_lon(i)
    if (.not. is_north_south_ghost(i, this)) then
      res = (j > nb_lon - 2)
    end if
  end if
end function 

!-----------------------------------------------------------------------------
!> Returns true if the current partition is the last one on the sector
!> @param i_pos, j_pos : partition position in the matrix array
!-----------------------------------------------------------------------------
function is_last_part_sector(i_pos, j_pos) result (res)
  integer :: nb_parts
  integer, intent(in) :: i_pos, j_pos
  logical :: res

  nb_parts = nb_parts_on_band_sector(i_pos)
  res = (j_pos == nb_parts .or. j_pos == 2*nb_parts .or. j_pos == 3*nb_parts)

end function


function is_ghost_cell(i, j, this) result(res)
  type (Band_grid) :: this
  integer, intent(in) :: i, j
  integer :: nb_lon
  logical :: res

  res = is_north_south_ghost(i, this)
  res = res .or. (j >= 1 .and. j <= this%nb_ghosts_west(i))
  nb_lon = this%nb_lon(i)
  res = res .or. (j > nb_lon - this%nb_ghosts_east(i) .and. j <= nb_lon)

end function

!-------------------------------------------------------------------------------
!> Returns true if the cell (i,j) is a ghost cell interfacing with an interior
!> cell, or if (i,j) is an interior cell
!> Assume there is a neighbouring latitude.
!>@param cur : cell position at i_neighb
!-------------------------------------------------------------------------------
function to_or_from_interior_cell(cur, i, j, i_neighb, this) result(res)
  type (Band_grid) :: this
  integer, intent(in) :: cur, i, j, i_neighb
  integer :: j_next1, j_next2
  logical :: is_ghost, res

  if (has_single_partition()) then
    res =.True.
    return
  end if

  call interior_lon_indices(j_next1, j_next2, this, i_neighb)
  is_ghost = is_ghost_cell(i, j, this)
  res = is_ghost .and.  .not. is_north_south_ghost(i_neighb, this) .and. &
    (cur > j_next1-1 .and. cur < j_next2+1)
  ! or from interior cells
  res = res .or. (.not. is_ghost)
end function

function is_interior_cell(i, j, this) result(res)
  type (Band_grid) :: this
  integer, intent(in) :: i, j
  logical :: res
  integer :: j_start, j_end

  call interior_lon_indices(j_start, j_end, this, i)
  res = (j >= j_start .and. j <= j_end)

end function

!-------------------------------------------------------------------------------
!!> Returns true if the cell (i,j) is a ghost cell interfacing with an interior
!!> cell
!!> Assume there is a neighbouring latitude.
!!>@param cur : cell position at i_neighb
!!-------------------------------------------------------------------------------
!function from_ghost_to_interior(cur, i, j, i_neighb, this) result(res)
!  type (Band_grid) :: this
!  integer, intent(in) :: cur, i, j, i_neighb
!  integer :: j_next1, j_next2
!  logical :: is_ghost, res
!
!  call interior_lon_indices(j_next1, j_next2, this, i_neighb)
!  is_ghost = is_ghost_cell(i, j, this)
!  res = is_ghost .and.  .not. is_north_south_ghost(i_neighb, this) .and. &
!    (cur > j_next1-1 .and. cur < j_next2+1)
!end function
!

!-----------------------------------------------------------------------------
!> Returns true if we are on north or south ghost cells. 
!-----------------------------------------------------------------------------
function is_north_south_ghost(i, this) result (res)
  type (Band_grid) :: this
  integer, intent(in) :: i
  logical :: res

  ! Only when left boundary we must increase of 1
  res = is_north_ghost(i, this)
  res = res .or. is_south_ghost(i, this)
end function

!-----------------------------------------------------------------------------
!> Returns true if we are on north ghost cells. 
!-----------------------------------------------------------------------------
function is_north_ghost(i, this) result (res)
  type (Band_grid) :: this
  integer, intent(in) :: i
  logical :: res

  ! Only when left boundary we must increase of 1
  res = (has_north_ghost_cells(this) .and. i == 1)
end function

!-----------------------------------------------------------------------------
!> Returns true if we are on south ghost cells. 
!-----------------------------------------------------------------------------
function is_south_ghost(i, this) result (res)
  type (Band_grid) :: this
  integer, intent(in) :: i
  logical :: res

  ! Only when left boundary we must increase of 1
  res = (has_south_ghost_cells(this) .and. i == this%nb_lat)
end function

!-----------------------------------------------------------------------------
!> Returns true if the cell is the westmost on the current sector
!>@param i, j : cell position in the partition
!>@param i_pos, j_pos : partition position
!-----------------------------------------------------------------------------
function is_eastmost_on_sector(i, j, i_pos, j_pos, this) result(res)
  type (Band_grid) :: this
  integer, intent(in) :: i, j
  integer, intent(in) :: i_pos, j_pos
  logical ::res
  integer :: nb_parts, nb_lon, l

  res = .False.
  ! Does not work for north and south ghost cells
  if (has_north_ghost_cells(this) .and. i == 1) return
  if (has_south_ghost_cells(this) .and. i == this%nb_lat) return

  nb_parts = nb_parts_on_band_sector(i_pos)
  nb_lon = this%nb_lon(i)
  ! The last true element is always nb_lon-1 (except for ghost cells, but we
  ! have eliminated this case)
  if (j == nb_lon-1) then
    ! The partition must also be the first on the current band
    nb_parts = nb_parts_on_band_sector(i_pos)
    do l = 1, 3
      res = (j_pos == nb_parts*l)
      if (res) return
    end do
  end if

end function

!-----------------------------------------------------------------------------
!> Compute the global position of the cell on the line i with float
!approximation. No need for a precision, fortran does the rounding.
!-----------------------------------------------------------------------------
function global_cell_position_from_local(this, i, j) result(j_glob)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j
  double precision :: lon_west, ratio, dlon
  integer :: j_glob

  lon_west = cell_west_lon(i, j, this)
  dlon = this%dlon(i)
  ratio = lon_west / dlon

  j_glob = nint(ratio) + 1

  !  if (i == 1 .and. j == 1) then
  !    write (*,*) "j cur 2", j_glob, "lon west, dlon", lon_west, this%dlon(i)
  !    write (*,*) "lon min", this%lon_min(1), this%lon_min(2)
  !  end if

end function

!-----------------------------------------------------------------------------
!> Returns the indice of the middle cell on the line i for the current grid. If
!> the middle is not in the current partition, returns -1 if the middle is to
!> the east, -2 otherwise (to the west).
!> The middle cell is the one for the current sector.
!> @param i : latitude inside the grid
!> @param zone : partition zone
!-----------------------------------------------------------------------------
function middle_cell(i, grid, zone) result(mid)
  type (Band_Grid) :: grid
  integer, intent(in) :: i, zone
  integer :: mid, nb_lon
  integer :: sector
  double precision :: lon_min, lon_mid, dlon
  logical :: must_adjust

  nb_lon = grid%nb_lon(i)
  lon_min = grid%lon_min(i)
  dlon = grid%dlon(i)
  ! Half of a sector
  sector = sector_from_zone(zone)
  lon_mid = 60.d0 + 120.d0*(sector - 1)

  if (lon_min + nb_lon*dlon < lon_mid) then
    mid = -1
    return
  end if
  if (lon_min > lon_mid) then
    mid = -2
    return
  end if

  must_adjust = .False.
  ! For periodic boundary. We don't take into account these ghost cells
  if (lon_min < 0) then
    lon_min = 0
    must_adjust = .True.
  end if
  ! Must add 1 as the indice start from 1
  mid = floor((lon_mid - lon_min)/dlon) + 1
  if (must_adjust) mid = mid + 1
end function

!-----------------------------------------------------------------------------
! Check the indice is inside the current grid
!-----------------------------------------------------------------------------
function is_indice_inside(k, this, i) result(res)
  type(Band_Grid) :: this
  logical :: res
  integer, intent(in) :: k, i

  res = (k > 0 .and. k < this%nb_lon(i) + 1)
end function

!-----------------------------------------------------------------------------
! Returns a range [k_start, k_end] containing the north neighbours of the
! current cell (or partition) inside our reduced grid. This is computed
! analytically.
! Input and output are done only with the position inside the current latitude
! (or partition) line. Nb_cur and nb_next are the number of cells (or
! partitions) in the line.
! Beware, it only works inside a zone.
!-----------------------------------------------------------------------------
subroutine north_partition_reduced(j, nb_cur, nb_prev, k_start, k_end)
  integer, intent(in) :: j, nb_cur, nb_prev
  integer :: k_start, k_end

  k_start = floor(dble(nb_prev*(j-1))/nb_cur + 1) 
  k_end = ceiling(dble(nb_prev*j)/nb_cur) 

end subroutine

!-----------------------------------------------------------------------------
! Same as north
!-----------------------------------------------------------------------------
subroutine south_partition_reduced(j, nb_cur, nb_next, k_start, k_end)
  integer, intent(in) :: j, nb_cur, nb_next
  integer :: k_start, k_end

  k_start = floor(dble(nb_next*(j-1))/nb_cur + 1)
  k_end = ceiling(dble(nb_next*j)/nb_cur)
end subroutine

!-----------------------------------------------------------------------------
!> Find the cell whose south latitude is lat
!> Assume global latitude is constant
!-----------------------------------------------------------------------------
function find_global_lat(lat, this) result(lat_g)
  type(Band_Grid) :: this
  double precision, intent(in) :: lat
  integer :: lat_g

  lat_g = nint((90.d0 - lat)/get_dlat(this))
end function

!-----------------------------------------------------------------------------
!> Copies the data from the global grid (which contains numerical data). 
!> Global is an integer checking we parse all the cells of the grid. 
!> We need the previous latitude line (cannotget it analytically if dlat is 
!> not constant).
!> Data must be allocated.
!> global_ratio is a pointer to the global grid ratio data.
!-----------------------------------------------------------------------------
subroutine numerical_concentration(this)
  type(Band_Grid) :: this

  call print_error("No numerical concentration implemented", &
    "concentration_numerical", fname_band)
end subroutine

!-----------------------------------------------------------------------------
!> Minimal longitude without the ghost cells
!-----------------------------------------------------------------------------
function get_lon_min(k, this) result(lon_min)
  type (Band_Grid) :: this
  integer :: k, nb_ghosts
  double precision :: lon_min

  lon_min = this%lon_min(k)
  nb_ghosts = this%nb_ghosts_west(k)
  lon_min = lon_min + this%dlon(k)*nb_ghosts
end function

!-----------------------------------------------------------------------------
!> Returns longitude of the west boundary for the cell at position (i,j),
!> We suppose the caller know how to avoid the ghost cells
!-----------------------------------------------------------------------------
function cell_west_lon(i, j, this) result(pos)
  type (Band_Grid) :: this
  integer :: i, j
  double precision :: pos

  pos = this%lon_min(i) + this%dlon(i)*(j - 1)
end function

function cell_west_lon2(i, j, this) result(pos)
  type (Band_Grid) :: this
  integer :: i, j
  double precision :: pos

  pos = this%lon_min(i) + this%dlon(i)*(j - 1)
end function


!-----------------------------------------------------------------------------
!> Returns longitude of the east boundary for the cell at position (i,j)
!> We suppose the caller know how to avoid the ghost cells
!-----------------------------------------------------------------------------
function cell_east_lon(i, j, this) result(pos)
  type (Band_Grid) :: this
  integer :: i, j
  double precision :: pos

  pos = this%lon_min(i) + this%dlon(i)*j
end function

!-----------------------------------------------------------------------------
!> Get the variation of longitude in cells at latitude i. Constant at a given
!> latitude.
!-----------------------------------------------------------------------------
function get_dlon(i, this) result(dlon)
  type (Band_Grid) :: this
  integer :: i
  double precision :: dlon
  dlon = this%dlon(i)
end function

!-----------------------------------------------------------------------------
!> Get the variation of latitude in cells at latitude i. Constant for all 
!> cells for now.
!-----------------------------------------------------------------------------
function get_dlat(this) result(dlon)
  type (Band_Grid) :: this
  double precision :: dlon

  dlon = this%dlat(1)
end function

!-----------------------------------------------------------------------------
!> Set the cell tracer ratio at the position (i,j), the indices being strictly 
!> positive
!-----------------------------------------------------------------------------
subroutine set_cell_ratio(val, i, j, tracer, this) 
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, tracer
  double precision, intent(in) :: val
  integer :: pos

  pos = this%first_cell(i) + j - 1
  !  if (pos < 1 .or. pos > size(this%ratio)) then
  !    print *, "pos outside the grid", i, j
  !  end if
 
  call set_tracer_ratio(this%tracers, val, pos, tracer)

end subroutine

!-----------------------------------------------------------------------------
!> Set the cell slope at the position (i,j), the indices being strictly positive
!-----------------------------------------------------------------------------
subroutine set_cell_slope(val, i, j, tracer, this) 
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, tracer
  double precision, intent(in) :: val
  integer :: pos

  pos = this%first_cell(i) + j - 1
  call set_tracer_slope(this%tracers, val, pos, tracer)

end subroutine

!-----------------------------------------------------------------------------
!> Set the cell gradient at the position (i,j), the indices being strictly positive
!-----------------------------------------------------------------------------
subroutine set_cell_gradient(val, i, j, tracer, this) 
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, tracer
  double precision, intent(in) :: val
  integer :: pos

  pos = this%first_cell(i) + j - 1
  call set_tracer_gradient(this%tracers, val, pos, tracer)

end subroutine

!-----------------------------------------------------------------------------
!> Set a value for meridional winds at position k in the wind array
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!-----------------------------------------------------------------------------
subroutine set_merid_winds(val, k, array, this) 
  type(Band_Grid) :: this
  integer, intent(in) :: k, array
  double precision :: val

  if (array == IS_PREV) then
    this%merid_winds_prev(k) = val
  else if (array == IS_NEXT) then
    this%merid_winds_next(k) = val
  else
    this%merid_winds(k) = val
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Set a value for zonal winds at position k in the wind array
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!-----------------------------------------------------------------------------
subroutine set_zonal_winds(val, k, array, this) 
  type(Band_Grid) :: this
  integer, intent(in) :: k, array
  double precision :: val

  if (array == IS_PREV) then
    this%zonal_winds_prev(k) = val
  else if (array == IS_NEXT) then
    this%zonal_winds_next(k) = val
  else
    this%zonal_winds(k) = val
  end if
end subroutine



!-----------------------------------------------------------------------------
!> Set air mass  for cell (i,j), i, j > 0
!-----------------------------------------------------------------------------
subroutine set_cell_air_mass(val, i, j, this) 
  type (Band_Grid) :: this
  integer, intent(in) :: i, j
  double precision, intent(in) :: val
  integer :: pos

  pos = this%first_cell(i) + j - 1
  this%mass_air(pos) = val

end subroutine

!-----------------------------------------------------------------------------
!> Get a tracer  ratio for cell (i,j), i, j > 0
!-----------------------------------------------------------------------------
function get_cell_ratio(i, j, tracer, this) result(val)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, tracer
  double precision :: val
  integer :: pos

  pos = this%first_cell(i) + j - 1
  val = get_tracer_ratio(this%tracers, pos, tracer) 
end function

!-----------------------------------------------------------------------------
!> Get a tracer slope for cell (i,j), i, j > 0
!-----------------------------------------------------------------------------
function get_cell_slope(i, j, tracer, this) result(val)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, tracer
  double precision :: val
  integer :: pos

  pos = this%first_cell(i) + j - 1
  val = get_tracer_slope(this%tracers, pos, tracer) 
end function

!-----------------------------------------------------------------------------
!> Get a tracer gradient for cell (i,j), i, j > 0
!-----------------------------------------------------------------------------
function get_cell_gradient(i, j, tracer, this) result(val)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, tracer
  double precision :: val
  integer :: pos

  pos = this%first_cell(i) + j - 1
  val = get_tracer_gradient(this%tracers, pos, tracer) 
end function


!-----------------------------------------------------------------------------
!> Get air mass for cell (i,j), i, j > 0
!-----------------------------------------------------------------------------
function get_cell_air_mass(i, j, this) result(val)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j
  double precision :: val
  integer :: pos

  pos = this%first_cell(i) + j - 1
  val = this%mass_air(pos) 
end function

!-----------------------------------------------------------------------------
!> Area of the cell (i,j) : suppose dlat/2 << 1 for approximation
!> Returns an area in degrees
!-----------------------------------------------------------------------------
function cell_area(i, j, this) result(val)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j
  double precision :: val, coef
  double precision :: lat_n, lat_s

  coef = pi/180

  !lat = cell_center_lat(i, this)
  !val = this%dlon(i)*this%dlat(1)*cos(lat*coef)

  if (at_the_poles(i, this)) then
    lat_s = cell_south_lat(i, this)
    val = 2*PI*(1 - sin(lat_s*coef))/(3*coef)
  else
    lat_n = cell_north_lat(i, this)
    lat_s = cell_south_lat(i, this)
    val = this%dlon(i)*(sin(lat_n*coef) - sin(lat_s*coef))/coef
  end if
end function

function at_the_poles(i, this) result (res)
  type (Band_Grid) :: this
  integer, intent(in) :: i
  logical :: res

  res = .False.
  if (abs(cell_north_lat(i, this) - 90.) < DBLE_PREC) res = .True.
  if (abs(cell_south_lat(i, this) + 90.) < DBLE_PREC) res = .True.
end function

!-----------------------------------------------------------------------------
!> Check if the current cell is at the middle of the sector (float computation) 
!> @param : i, j : cell position in grid
!-----------------------------------------------------------------------------
function is_middle_float(i, j, this, zone) result(res)
  type (Band_Grid) :: this
  integer, intent(in) :: i, j, zone
  logical :: res
  double precision :: lon1, lon2, middle
  integer :: sector

  lon1 = cell_west_lon(i, j, this)
  lon2 = cell_east_lon(i, j, this)
  sector = sector_from_zone(zone)
  middle = 60.d0 + 120.d0*(sector-1)
  res = lon1 < middle .and. lon2 > middle
end function

!-----------------------------------------------------------------------------
!> Assuming the cell in inside a sector, return the global position
!-----------------------------------------------------------------------------
function local_to_global(i, j, this) result(pos)
  type (Band_Grid) :: this
  double precision :: center
  integer, intent(in)  :: i, j
  integer :: pos

  center = cell_center_lon(i, j, this)
  pos = floor(center/this%dlon(i))+1
end function

!-----------------------------------------------------------------------------
!> Assume constant dlat and round to nearest integer (the cast to int does not
!> work in some cases)
!-----------------------------------------------------------------------------
function first_global_lat(lat_min, grid) result(first_lat)
  type(Band_grid) :: grid
  double precision, intent(in) :: lat_min
  integer :: first_lat

  first_lat = nint((90.d0 - lat_min) / grid%dlat(1)) + 1
end function

!-----------------------------------------------------------------------------
!> Compute a tracer mass inside the grid (interior cells only)
!-----------------------------------------------------------------------------
function mass_tracer(tracer, this) result(mass)
  type (Band_Grid) :: this
  double precision :: ratio, mass, air_mass
  integer :: i, j, k, tracer
  integer :: i_start, i_end
  integer :: j_start, j_end

  mass = 0
  k = 0
  call interior_lat_indices(i_start, i_end, this)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, this, i)
    ! Concentration*cell area
    air_mass = cell_area(i, j_start, this)
    !if (i < 5) print *, "i=", air_mass
    do j = j_start, j_end
      ratio = get_cell_ratio(i, j, tracer, this)
      mass = mass + air_mass*ratio

      k = k + 1
    end do
  end do

end function

subroutine set_merid_flux(this, val, k_f, tracer)
  type(Band_Grid) :: this
  integer :: k_f, tracer
  double precision :: val

  call set_tracer_merid_flux(this%tracers, val, k_f, tracer)
end subroutine

subroutine set_zonal_flux(this, val, k_f, tracer)
  type(Band_Grid) :: this
  integer :: k_f, tracer
  double precision :: val

  call set_tracer_zonal_flux(this%tracers, val, k_f, tracer)
end subroutine

function get_merid_flux(this, k_f, tracer) result(val)
  type(Band_Grid) :: this
  integer :: k_f, tracer
  double precision :: val

  val = get_tracer_merid_flux(this%tracers, k_f, tracer)
end function

function get_zonal_flux(this, k_f, tracer) result(val)
  type(Band_Grid) :: this
  integer :: k_f, tracer
  double precision :: val

  val = get_tracer_zonal_flux(this%tracers, k_f, tracer)
end function

subroutine unset_gradients(this)
  type(Band_Grid) :: this

  call unset_gradients_all(this%tracers)
end subroutine

subroutine init_tracer_ratio(this)
  type(Band_Grid) :: this

  call init_tracer_ratio_all(this%tracers)
end subroutine


subroutine free_Band_grid(this)
  type(Band_Grid) :: this
  integer :: ierr

  call free_array(this%dlon, "dlon", fname_band)
  call free_array(this%lon_min, "lon_min", fname_band)
  call free_array(this%dlat, "dlat", fname_band)
  call free_array(this%first_cell, "first_cell", fname_band)
  call free_array(this%north_offset, "north_offset", fname_band)
  call free_array(this%south_offset, "south_offset", fname_band)
  call free_array(this%lon_min, "lon_min", fname_band)
  call free_array(this%mass_air, "mass_air", fname_band)
  call free_array(this%nb_lon, "nb_lon", fname_band)
  call free_array(this%nb_ghosts_east, "nb_ghosts_east", fname_band)
  call free_array(this%nb_ghosts_west, "nb_ghosts_west", fname_band)

  call free_array(this%merid_winds, "merid winds", fname_band)
  call free_array(this%merid_winds_prev, "merid winds prev", fname_band)
  call free_array(this%merid_winds_next, "merid winds next", fname_band)
  call free_array(this%zonal_winds, "zonal winds", fname_band)
  call free_array(this%zonal_winds_prev, "zonal winds prev", fname_band)
  call free_array(this%zonal_winds_next, "zonal winds next", fname_band)

#ifdef HDF5
  call free_array(this%center_lat, "center_lat", fname_band)
  call free_array(this%center_lon, "center_lon", fname_band)
  call free_array(this%u_lat, "u_lat", fname_band)
  call free_array(this%u_lon, "u_lon", fname_band)
  call free_array(this%ratio_offset, "ratio_offset", fname_band)
  call free_array(this%zonal_offset, "zonal_offset", fname_band)
  call free_array(this%merid_offset, "merid_offset", fname_band)
  call free_array(this%merid_length, "merid_length", fname_band)
#endif

  call free_List_Tracers(this%tracers)
  ! Free mpi arrays in the class
  call free_array_mpi_types(this%mpi_ghost_cells, "mpi_ghost_cells", fname_band)
  call free_array_mpi_types(this%mpi_interior_cells, "mpi_interior_cells",&
    fname_band)

  if (this%mpi_grid /= mpi_datatype_null) then
    call mpi_type_free(this%mpi_grid, ierr)
    call check_error(ierr,"trying to free mpi_grid",  "free_Band_grid", &
      fname_band)
  end if

end subroutine

end module
