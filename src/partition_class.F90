!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Partition
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> A partition is a grid with more information (id, neighbours...)
!
!-------------------------------------------------------------------------------

module Partition_class
use Band_grid_class
use Configuration_class
use Global_grid_class
use Analytical_partitioning
implicit none

public 
type Partition
  !> Id of the partitions (between 0 and _NB_PARTITIONS-1)
  integer :: id = -1
  !> Position in the partitioning matrix ( positive)
  integer :: i_pos = -1
  integer :: j_pos = -1
  !> If the partition in on the residue band (exists if the number of
  !> partitions is not a square)
  type(Band_grid) :: grid
  !> Contains the total variation of the height in the last band. Positive for
  !> square partitions, negative for the last band (we actually substracted 1 to
  !> the true value).
  !> If there is no last band, it is < 0 for the square partitions and contains
  !> the old height of the last band
  integer :: height_var = 0

  !> MPI type for partitioning
  integer :: mpi_partition

  !> As we sometimes resize the width of the partition, the extremities are 
  !> modified. This is stored here. In base 10, it tells us the number of
  !> inversion. First digit (<10) for north, second digit for south. 
  !> The digit is the number of inversion (1 for left, 2 for right, 3 for both)
  !> A north inversion means the partition is shifted toward the center of the
  !> sector. On the northern hemisphere, on the left, it means the partition has
  !> one more neighbour to the right (digit = 1 or 3) and one less neighbour to
  !> the left (digit = 2 or 3). On the left side, this is inverted.
  !> For south overlap, the neighbour is shift toward the extremities, so this
  !> is invereted.
  !> For the southern hemisphere, we invert only the digit for south with the
  !> one for norh, so the previous logic still holds
  integer :: overlap = 0

  ! Array of requests for all neighbours (send and receive in the same array)
  integer, allocatable :: requests(:)

  !> The zone in which the partition is in. Numbered from 1 to 6, 1-2-3 being 
  !> the northern hemisphere (west to east).
  integer :: zone = -1


end type


! For output
character(*), parameter :: fname_partition = "partition_class.F90"
contains 

#include "analytical_data.inc"
!-------------------------------------------------------------------------------
!> Default constructor. Id and grid must be set outside.
!-------------------------------------------------------------------------------
subroutine new_Partition(this, i_pos, j_pos, height_var, overlap, zone)
  type(Partition) :: this
  integer :: i_pos, j_pos
  !type(Band_grid) :: grid
  integer :: height_var, overlap, zone
  this%i_pos = i_pos
  this%j_pos = j_pos
  this%height_var = height_var
  this%overlap = overlap
  this%zone = zone

  ! Initialize mpi type
  this%mpi_partition = mpi_datatype_null

end subroutine

!-------------------------------------------------------------------------------
!> Copy constructor. Here we also copy the grid.
!-------------------------------------------------------------------------------
subroutine new_Partition_copy(this, other)
  type(Partition) :: this, other

  this%i_pos = other%i_pos
  this%j_pos = other%j_pos
  this%height_var = other%height_var
  this%overlap = other%overlap
  this%zone = other%zone
  call new_Band_grid_copy(this%grid, other%grid)

  ! Initialize mpi type
  this%mpi_partition = mpi_datatype_null
end subroutine


!-----------------------------------------------------------------------------
!> Set partition position and the grid
!> @param i : position in the partition matrix (latitude)
!> @param j : position in the partition matrix (longitude)
!> @param hemisph : 1 on northern hemisphere, -1 otherwise
!> @param first_lat : starting latitude from the closest pole (beware of the 
!> minimal latitude)
!> @param nb_resized : number of resized partition on the current band (only 
!> one one half)
!> @param width : width of a square partition (before resizing)
!-----------------------------------------------------------------------------
subroutine new_Partition_from_pos(this, i, j, first_lat, height, width, nb_resized,&
    zone, height_corr)
  type (Partition) :: this
  integer, intent(in) :: i, j, first_lat, height, zone
  integer, intent(in) :: width, nb_resized, height_corr
  integer :: ghost_cells

  call set_pos_Partition(this, i, j)
  this%zone = zone

  ! Number of directions for ghost cells
  ghost_cells = compute_ghost_cells_direction(this, i, j) 

  call new_Band_grid_from_pos(this%grid, i, j, first_lat, height, width, &
    nb_resized, zone, ghost_cells)

  ! Set misc data : height_var, and the ghost cells
  this%height_var = height_corr

  ! Finally, set first cells (could not do it before because of the ghost cells)
  call set_first_cells_Band_grid(this%grid)

  ! Initialize mpi type
  this%mpi_partition = mpi_datatype_null

end subroutine

!-------------------------------------------------------------------------------
!> Create a grid at position i on the leftover band. Also finishes to set some
!> values
!> @param height : partition height
!> @param width : partition width (except for the last one)
!> @param j_pos : partition number on the band (no need to correct it)
!> @param i_pos : partition band
!> @param n_left : number of partitions on the band
!> @param zone : zone indice
!-------------------------------------------------------------------------------
subroutine new_Partition_leftover(this, i_pos, j_pos, height, width, n_left, &
    zone)
  type(Partition) :: this
  integer :: height, n_left, j
  integer, intent(in) :: width, zone, i_pos, j_pos
  integer :: ghost_cells, corr_nb

  !write (*,*) "leftover width", width
  this%zone = zone
  ! Set the ghost cells
  call set_pos_Partition(this, i_pos, j_pos)
  ghost_cells = compute_ghost_cells_direction(this, i_pos, j_pos) 

  call new_Band_grid_leftover(this%grid, i_pos, j_pos, height, width, n_left,&
    zone, ghost_cells)

  ! Finally, set first cells (could not do it before because of the ghost cells)
  call set_first_cells_Band_grid(this%grid)

  ! Initialize mpi type
  this%mpi_partition = mpi_datatype_null

end subroutine

#ifdef HDF5
subroutine set_hdf5_data_Partition(this)
  type (Partition) :: this

  call set_hdf5_offset(this)
  call set_hdf5_offset_merid(this)
end subroutine

!-------------------------------------------------------------------------------
!> Set the offset for ratio and zonal winds (meridiional done later) for faster I/O
!> We have to read line by line as hdf5 does not allow for varying stride
!> merid_length is not set here, but when allocating the meridional winds
!> instead
!-------------------------------------------------------------------------------
subroutine set_hdf5_offset(this) 
  type (Partition), target :: this
  type (Band_Grid), pointer :: grid
  integer :: i, i_start, i_end
  integer :: k, tmp, i_glob, nb
  double precision :: lat_min
  integer :: sum_ratio_prev, first
  integer :: sum_zonal_prev

  grid => this%grid

  lat_min = get_lat_min(grid)
  i_glob = first_global_lat(lat_min, grid)
  call interior_lat_indices(i_start, i_end, grid)

  k = 1
  ! Store the sum of all previous latitudes
  sum_ratio_prev = sum_nb_cells(i_glob-1)
  sum_zonal_prev = sum_nb_zonal(i_glob-1)

  do i = i_start, i_end
    ! Store all the offsets at once
    first = first_interior_cell_position_global(i, this)-1
    grid%ratio_offset(k) = sum_ratio_prev + first
    grid%zonal_offset(k) = sum_zonal_prev + first

    nb = nb_cells_lat(i_glob)

    sum_ratio_prev = sum_ratio_prev + nb
    sum_zonal_prev = sum_zonal_prev + nb + 1

    i_glob = i_glob + 1
    k = k + 1
  end do

end subroutine

!-------------------------------------------------------------------------------
!> Set the offset for meridion winds for faster I/O
!-------------------------------------------------------------------------------
subroutine set_hdf5_offset_merid(this) 
  type (Partition), target :: this
  type (Band_Grid), pointer :: grid
  integer :: i, i_start, i_end
  integer :: k, tmp, i_glob, nb
  double precision :: lat_min
  integer :: first, sum_merid_prev

  grid => this%grid

  lat_min = get_lat_min(grid)
  i_glob = first_global_lat(lat_min, grid)
  call interior_lat_indices(i_start, i_end, grid)
  ! Read the interface with ghost cells
  if (has_north_ghost_cells(this%grid)) then
    i_start = i_start - 1
    i_glob = i_glob - 1
  end if
  if (.not. has_south_ghost_cells(this%grid)) i_end = i_end - 1

  k = 1
  ! Store the sum of all previous latitudes
  sum_merid_prev = sum_nb_merid(i_glob-1)

  !print *, "starts end", i_start, i_end
  do i = i_start, i_end
    !    if (this%id == 5) print *, "first cell", first, i
    first = first_merid_position_global(i, i_glob, this)

    grid%merid_offset(k) = sum_merid_prev + first-1
    !if (this%id == 2) then
      !print *, "offset",grid%merid_offset(k), "i", i, "first merid", first
    !end if

    nb = nb_merid_lat(i_glob)
    sum_merid_prev = sum_merid_prev + nb 

    i_glob = i_glob + 1
    k = k + 1
  end do
  !if (this%id == 5) then
  !  print *, "merid offset", grid%merid_offset
  !  print *, "merid length", grid%merid_length
  !end if

end subroutine

#endif

!-----------------------------------------------------------------------------
!> Idle processes are tagged with an id of -2
!-----------------------------------------------------------------------------
function is_idle_Partition(this) result (res)
  type (Partition) :: this
  logical :: res

  res = (this%id == -2)
end function


subroutine init_Partition(this,id)
  type(Partition) :: this
  integer :: id
  this%id = id

end subroutine

!-----------------------------------------------------------------------------
!> MPI user-defined datatype from already defined partition, and already defined
!> MPI type for band_grid. 
!-----------------------------------------------------------------------------
subroutine new_MPI_Partition(this)
  type(Partition) :: this
  integer :: ierr, mpi_band_grid
  ! Length, displacement, type
  integer :: mlength(7), mtype(7)
  ! Must use mpi_address_kind (MPI 2)
  integer(kind=mpi_address_kind) :: mlocation(7)
  integer(kind=mpi_address_kind) :: start
  integer :: i
  integer :: test
  character(*), parameter :: func_name = "new_MPI_Partition" 

  !call print_Partition(this)
  call mpi_get_address(this, start, ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_partition)

  call get_address(this%id, mlocation(1), func_name, fname_partition)
  call get_address(this%i_pos, mlocation(2), func_name, fname_partition)
  call get_address(this%j_pos, mlocation(3), func_name, fname_partition)

  call mpi_get_address(this%grid, mlocation(4), ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_partition)

  call get_address(this%height_var, mlocation(5), func_name, fname_partition)
  call get_address(this%overlap, mlocation(6), func_name, fname_partition)
  call get_address(this%zone, mlocation(7), func_name, fname_partition)

  ! Allocate arrays after mpi_get_address, as it is safer
  mlength = 1
  ! For an idle proc, we do not need a band grid, just send the integer
  if (is_idle_Partition(this)) then
    mtype = mpi_integer
  else
    call new_MPI_Band_grid(this%grid)
    mpi_band_grid = this%grid%mpi_grid
    mtype = (/mpi_integer, mpi_integer, mpi_integer, mpi_band_grid, &
      mpi_integer, mpi_integer, mpi_integer/)
  end if

  do i = 1, size(mlength)
    mlocation(i) = mlocation(i) - start
  end do

  call mpi_type_create_struct(size(mlength), mlength, mlocation, mtype, &
    this%mpi_partition, ierr)
  call check_mpi_error(ierr, "create struct", func_name, fname_partition)

  call mpi_type_commit(this%mpi_partition, ierr)
  call check_mpi_error(ierr, "commit type", func_name, fname_partition)
end subroutine


!-------------------------------------------------------------------------------
!> Compute minimal longitude for north and south ghost cells on a leftover band.
!> @param i : partition number of the leftover band.
!-------------------------------------------------------------------------------
subroutine lon_min_leftover_ghost_cells(lon_min, grid, i, width, dlon, nb_lat,&
    first_lat, zone)
  type (Band_grid) :: grid
  integer, intent(in) :: i, width,  nb_lat, first_lat
  integer, intent(in) :: zone
  integer :: nb_cells, mid, nb_lat2
  integer :: cur, neighb, i_lat
  double precision :: lon_min(:)
  double precision, intent(in) :: dlon(:)

  if (on_northern_hemisphere(zone)) then
    ! On the north hemisphere, there are always south ghost cells, but set later
    if (has_north_ghost_cells(grid)) then
      call find_lon_min_north_ghost_last(lon_min, grid, i, width, dlon(1), &
        first_lat)
    end if
    ! Other side of the equator
    lon_min(nb_lat) = lon_min(nb_lat - 1)

    ! Avoid periodic conditions on south ghost cells
    if (lon_min(nb_lat) < 0) lon_min(nb_lat) = 0
  else
    ! On the south hemisphere, same situation, but inverted.
    if (has_south_ghost_cells(grid)) then
      ! Always north ghost cells
      i_lat = nb_lat - 2
      !      write (*,*) "size lonmin", size(lon_min)
      call find_lon_min_south_ghost_last(lon_min, grid, i, width, dlon(nb_lat), &
        i_lat, zone)
    end if
    ! Other side of the equator
    lon_min(1) = lon_min(2)
    ! Avoid periodic conditions on south ghost cells
    if (lon_min(1) < 0) lon_min(1) = 0

  end if
end subroutine


!-------------------------------------------------------------------------------
!> Compute the minimal longitude for north ghost cells on the leftover band.
!> For that we compute the north neighbour of the west extremity
!-------------------------------------------------------------------------------
subroutine find_lon_min_north_ghost_last(lon_min, grid, i, width, dlon, &
    i_lat)
  type (Band_grid) :: grid
  integer, intent(in) :: i, width, i_lat
  integer :: nb_cells, mid
  integer :: cur, neighb
  double precision :: lon_min(:)
  double precision, intent(in) :: dlon

  ! Find west extremity and its neighbour
  cur = (i - 1)*width + 1
  ! mid is i_lat
  neighb = left_north_cell_northern(cur, i_lat)
  nb_cells = nb_cells_lat_sector(i_lat + 1)
  ! If cur is at the middle, neighb is too much to the left
  if (cur == (nb_cells + 1)/2) neighb = neighb + 1

  lon_min(1) = (neighb - 1)*dlon
end subroutine

!-------------------------------------------------------------------------------
!> True if the partition is on the boundary 
!-------------------------------------------------------------------------------
function is_partition_on_boundary(this) result(res)
  type (Partition) :: this
  logical :: res

  res = (this%j_pos == 1 .or. is_last_on_band(this))
end function


!-------------------------------------------------------------------------------
!> Store the directions in which there are ghost cells in base 2
!-------------------------------------------------------------------------------
function compute_ghost_cells_direction(this, i, j) result(ghost_cells)
  type (Partition) :: this
  integer :: ghost_cells, nb_lat
  integer :: nb_parts
  integer :: i, j

  nb_lat = nb_lat_Global_grid()
  ! Number of parts on the current band
  nb_parts = nb_parts_on_band_sector(i)
  !write (*,*) "nb parts on band", nb_parts, i

  ghost_cells = 0
  ! Always west ghost cells (periodicity)
  ghost_cells = ghost_cells + 2
  ! Same for east
  ghost_cells = ghost_cells + 8

  ! North
  if (.not. is_northmost(this)) ghost_cells = ghost_cells + 1
  ! South
  if (.not. is_southmost(this)) ghost_cells = ghost_cells + 4
end function

!-----------------------------------------------------------------------------
!> Returns true if the current latitude line is on the north pole
!-----------------------------------------------------------------------------
function is_north_pole(i, this) result (res)
  type (Partition) :: this
  integer, intent(in) :: i
  logical :: res

  ! No north ghost cells on first band
  res = (is_northmost(this) .and. i == 1)

end function

!-----------------------------------------------------------------------------
!> Returns true if the current latitude line is on the south pole
!-----------------------------------------------------------------------------
function is_south_pole(i, this) result (res)
  type (Partition) :: this
  integer, intent(in) :: i
  logical :: res

  ! No north ghost cells on first band
  res = (is_southmost(this) .and. i == this%grid%nb_lat)

end function


!-----------------------------------------------------------------------------
!> Returns true if the current partition is on the northmost band.
!> The partition zone must be set.
!-----------------------------------------------------------------------------
function is_northmost(this) result (res)
  type (Partition) :: this
  logical :: res
  res = .False.

  if (on_northern_hemisphere(this%zone)) res = (this%i_pos == 1)
end function

!-----------------------------------------------------------------------------
!> Returns true if the current partition is on the southmost band.
!> The partition zone must be set.
!-----------------------------------------------------------------------------
function is_southmost(this) result (res)
  type (Partition) :: this
  logical :: res
  integer :: nb_bands

  res = .False.
  nb_bands = get_total_nb_bands_Configuration_sector()
  ! Special case : only one partition or 3
  if (has_three_partitions() .or. has_single_partition()) res = .True.

  if (on_southern_hemisphere(this%zone)) then
    res = (this%i_pos == nb_bands)
  end if
end function

!-----------------------------------------------------------------------------
!> Returns true if the current band is on the northern hemisphere.
!-----------------------------------------------------------------------------
function band_on_northern_hemisphere(i_pos) result (res)
  integer, intent(in) :: i_pos
  logical :: res

  res = (i_pos < get_nb_bands_Configuration(1) + 1)
end function


!-----------------------------------------------------------------------------
!> MPI user-defined datatype for 2D advection : we create the interior cells 
!> (sent to the adjacent neighbours) and ghost cells (received from neihgbours)
!> The types are defined with the adequate offsets in order to be used 
!> directly with the start of the array (so mpi_indexed)
!> The types are actually in the grid, so we call the adequate constructor,
!> while most of the work is done here.
!-----------------------------------------------------------------------------
subroutine new_MPI_borders_Partition(this)
  type(Partition) :: this
  integer :: j, nb_neighb
  integer :: k, l
  integer :: i_neighb, j_start, j_end
  k = 1

  ! Don't forget to allocate as much mpi type as there are neighbours
  nb_neighb =  nb_neighbours_Partition(this)
  call allocate_mpi_borders(this%grid, nb_neighb)

  call north_partition_wrapper(this, i_neighb, j_start, j_end)
  call create_MPI_borders_Partition(this, k, i_neighb, j_start, j_end, NORTH)

  call east_partition_wrapper(this, i_neighb, j_start)
  call create_MPI_borders_Partition(this, k, i_neighb, j_start, j_start, EAST)

  call south_partition_wrapper(this, i_neighb, j_start, j_end)
  call create_MPI_borders_Partition(this, k, i_neighb, j_start, j_end, SOUTH)

  call west_partition_wrapper(this, i_neighb, j_start)
  call create_MPI_borders_Partition(this, k, i_neighb, j_start, j_start, WEST)

  call extend_MPI_borders_Partition(this)
end subroutine

!-------------------------------------------------------------------------------
!> Now that all MPI types are created, extend it for all tracer
!-------------------------------------------------------------------------------
subroutine extend_MPI_borders_Partition(this)
  type(Partition) :: this

  call new_MPI_interior_tracers(this%grid%tracers, this%grid%mpi_interior_cells)
  call new_MPI_ghost_tracers(this%grid%tracers, this%grid%mpi_ghost_cells)
end subroutine


!-----------------------------------------------------------------------------
!> Create the MPI type corresponding to a border of the current partition :
!> interior or ghost cells
!-----------------------------------------------------------------------------
subroutine create_MPI_borders_Partition(this, k, i_neighb, j_start, j_end, &
    direction)
  type(Partition) :: this
  integer :: i_neighb, j_start, j_end
  integer :: j, length, direction
  integer :: k, i_start
  integer :: length_prev, length_prev2
  integer :: start_prev, end_prev
  logical :: at_equator

  if (j_start > -1) then
    if (j_end < 0) j_end = j_start

    length_prev = 0
    length_prev2 = 0
    start_prev = 0
    end_prev = 0
    do j = j_start, j_end
      length = interface_length(this, i_neighb, j)
      !write (*,*) "inter length", length
      if (length < 0) then
        call print_error("Negative ghost interface length", &
          "create_MPI_borders_Partition", fname_partition)
      end if
      call create_mpi_ghost_cells(this, length, length_prev, k, &
        direction, this%id, start_prev, end_prev)
      ! Store the previous length if several north/souh neighbours
      length_prev = length_prev + length

      call interface_interior_cells(i_start, length, this, i_neighb, j)
      if (length < 0) then
        call print_error("Negative interior interface length", &
          "create_MPI_borders_Partition", fname_partition)
      end if

      call create_mpi_interior_cells(this, i_start, length, length_prev2, k, &
        direction, this%id)

      length_prev2 = length_prev2 + length
      k = k + 1
    end do
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Set the offset and block length for mpi interior cells.
!> Beware, the offset must start from 0.
!> Special cases for east and west as only 1 cell is used for East/West
!> transport but more can be used for North/South.
!-------------------------------------------------------------------------------
subroutine set_mpi_offset_length_interior_cells(offsets, blocklens, i_start, &
    length, length_prev, nb, this, direction)
  type (Partition), target :: this
  type (Band_grid), pointer :: grid
  integer, intent(in) :: direction, length, length_prev, nb
  integer, intent(in) :: i_start 
  integer :: blocklens(:), offsets(:)
  integer :: i, corr, i_corr
  integer :: nb_lat, length_prev_corr
  integer :: nb_int
  integer :: i1, i2
  corr = 0

  grid => this%grid
  ! North : first or second line (contiguous)
  if (direction == NORTH) then
    if (has_north_ghost_cells(grid)) corr = corr + grid%nb_lon(1) 
    offsets(1) = i_start + corr
    blocklens(1) = length

    ! East
  else if (direction == EAST) then
    ! If has northern ghost cells, we must correct the offset
    i_corr = 0
    if (has_north_ghost_cells(grid)) i_corr = 1
    call interior_lat_indices(i1, i2, grid)

    do i = 1, nb
      ! There can be more than one cells needed for the neighbor
      nb_int = nb_interior_cells_east(i1 + i-1, this)

      ! As the previous offset is on the previous line, we add nb_lon instead of
      ! nb_lon - 1
      offsets(i) = sum(grid%nb_lon(1:i+i_corr))- grid%nb_ghosts_east(i+i_corr)
      offsets(i) = offsets(i) - nb_int
      blocklens(i) = nb_int
    end do

    ! South : last line (contiguous)
  else if (direction == SOUTH) then
    blocklens(1) = length
    !! We sum all the cells up to the before last line
    nb_lat = last_lat_Band_grid(grid)
    offsets(1) = sum(grid%nb_lon(1:nb_lat-1))
    offsets(1) = offsets(1) + i_start

    ! West 
  else
    i_corr = 0
    if (has_north_ghost_cells(grid)) i_corr = 1

    call interior_lat_indices(i1, i2, grid)
    do i = 1, nb
      ! There can be more than one cells needed for the neighbor
      nb_int = nb_interior_cells_west(i1 + i-1, this)

      ! As the previous offset is on the previous line, we add nb_lon instead of
      ! nb_lon - 1
      offsets(i) = sum(grid%nb_lon(1:i+i_corr-1)) + grid%nb_ghosts_west(i+i_corr)
      blocklens(i) = nb_int
    end do

  end if
end subroutine

!-------------------------------------------------------------------------------
!> Create the mpi type for the ghost cells (exterior to the grid)
!> We only need the length of the interface
!> @param direction integer for the direction 0=north... 3=west
!> @param i_start : first cell on the current line
!> @param length : interface length
!> @param length_prev : sum of all previous interface length
!> @param start_prev : previous start of the sub array (for check)
!> @param end_prev : previous end of the sub array(for check)
!-------------------------------------------------------------------------------
subroutine create_mpi_ghost_cells(this, length, length_prev, k, &
    direction, id, start_prev, end_prev)
  type (Partition) :: this
  integer, intent(in) :: length,length_prev, direction, k
  ! For debug
  integer, intent(in) :: id
  integer :: i, sum_prev

  integer, allocatable :: blocklens(:), offsets(:)
  integer :: start_prev, end_prev
  integer :: nb, val


  ! North-south : data is contiguous (1 block)
  if (direction == NORTH .or. direction == SOUTH) then
    nb = 1
    ! East west : nb lat
  else
    nb = get_nb_lat_Partition(this)
  end if

  allocate (blocklens(nb))
  allocate (offsets(nb))

  call set_mpi_offset_length_ghost_cells(offsets, blocklens, length, &
    length_prev, nb, this, direction)
  !if (id == 16 .and. direction == NORTH) then
  !  if (id == 6 .and. direction == SOUTH) then
  !    print *, "ghost west offset for 6", offsets
  !  end if


  call new_MPI_ghost_cells(val, offsets, blocklens, this%grid, nb, direction, &
    start_prev, end_prev)
  this%grid%mpi_ghost_cells(k) = val

  deallocate (blocklens)
  deallocate (offsets)
end subroutine


!-------------------------------------------------------------------------------
!> Same as create_mpi_ghost_cells but for interior cells.
!> @param i_start : first cell on the current line
!> @param length : interface length
!> @param length_prev : sum of all previous interface length
!> @param at_equator : true if the current line is just before/after the equator,
!> and the ghost cells are on the other side
!-------------------------------------------------------------------------------
subroutine create_mpi_interior_cells(this, i_start, length, length_prev, k, &
    direction, id)
  type (Partition) :: this
  integer, intent(in) :: length, length_prev, direction, k
  ! For debug
  integer, intent(in) :: id
  integer :: i, sum_prev
  integer, intent(in) :: i_start
  integer :: nb, val
  integer, allocatable :: blocklens(:), offsets(:)

  ! North-south : data is contiguous (1 block)
  if (direction == NORTH .or. direction == SOUTH) then
    nb = 1
    ! East west : nb lat
  else
    nb = get_nb_lat_Partition(this)
  end if

  allocate (blocklens(nb))
  allocate (offsets(nb))

  call set_mpi_offset_length_interior_cells(offsets, blocklens, i_start, &
    length, length_prev,nb, this, direction)
  ! Create mpi for all fields, using the offset of the concentration array
  call new_MPI_global_type(val, offsets, blocklens, this%grid, nb)
  this%grid%mpi_interior_cells(k) = val

  !  if (id == 13 .and. direction == NORTH) then
  !    print *, "interior south offset for 6", offsets
  !  end if

  deallocate (blocklens)
  deallocate (offsets)
end subroutine

!-------------------------------------------------------------------------------
!> Returns all longitude indices, except for ghost cells at sector boundaries
!> Needs to be in partition_class, as we need to know the position on the band
!-------------------------------------------------------------------------------
subroutine lon_indices(j_start, j_end, this, k)
  type (Partition) :: this
  integer :: j_start, j_end, k

  j_start = 1
  j_end = this%grid%nb_lon(k)

  if (.not. has_single_partition()) then
    if (is_first_on_band(this)) j_start = first_interior_lon(k, this%grid)
    if (is_last_on_band(this)) j_end =  last_interior_lon(k, this%grid)
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Set the offset and block length for mpi ghost cells
!-------------------------------------------------------------------------------
subroutine set_mpi_offset_length_ghost_cells(offsets, blocklens, length, &
    length_prev, nb, this, direction)
  type (Partition), target :: this
  type (Band_grid), pointer :: grid

  integer, intent(in) :: direction, length, length_prev, nb
  integer :: blocklens(:), offsets(:)
  integer :: i, i_lat
  integer :: i_corr, corr
  integer :: i_start, i_end
  integer :: nb_ghosts

  grid => this%grid
  ! North : first line (contiguous)
  if (direction == NORTH) then
    offsets(1) = length_prev
    blocklens(1) = length
    !write (*,*) "offest north ghost", offsets(1)

    ! East
  else if (direction == EAST) then
    i_corr = 0
    if (has_north_ghost_cells(grid)) i_corr = 1

    do i = 1, nb
      ! As the previous offset is on the previous line, we add nb_lon instead of
      ! nb_lon - 1
      blocklens(i) = grid%nb_ghosts_east(i + i_corr)
      offsets(i) = sum(grid%nb_lon(1:i+i_corr)) - blocklens(i)
    end do

    ! South : last line (contiguous)
  else if (direction == SOUTH) then
    i_lat = grid%nb_lat
    blocklens(1) = length

    offsets(1) = sum(grid%nb_lon(1:i_lat-1)) + length_prev
    ! West 
  else
    offsets(1) = 0
    i_corr = 0
    ! Nb_lon with ghost cells
    if (has_north_ghost_cells(grid)) then
      offsets(1) = grid%nb_lon(1)
      i_corr = 1
    endif
    blocklens(1) = grid%nb_ghosts_west(1 + i_corr)

    ! An offset is deduced from the previous by adding the previous length line
    do i = 2, nb
      ! We only want one west ghost cells
      blocklens(i) = grid%nb_ghosts_west(i + i_corr)
      offsets(i) = sum(grid%nb_lon(1:i - 1 + i_corr))
    end do
  end if
  !  print *, "for ", id, "dir ", direction, "offset, length=", offsets, blocklens
end subroutine

!-----------------------------------------------------------------------------
!> Size of meridional winds in the partition. We do not have north winds at 
!> the north pole, nor south winds at the south pole.
!> We use this occasion to allocate merid_length and set it
!-----------------------------------------------------------------------------
function size_winds_merid(this) result(nb)
  type (Partition) :: this
  integer :: nb, k
  integer :: i_start, i_end, i

  nb = 0
  call interior_lat_indices(i_start, i_end, this%grid)
  ! North and south extremities, if there are ghost cells
  if (has_north_ghost_cells(this%grid)) i_start = i_start - 1
  if (.not. has_south_ghost_cells(this%grid)) i_end = i_end - 1

#ifdef HDF5
  allocate(this%grid%merid_length(i_end-i_start+1))
#endif
  !print *, "size merid length", size(this%grid%merid_length), i_start, i_end

  k = 1
  do i = i_start, i_end
#ifdef HDF5
    this%grid%merid_length(k) = nb_winds_merid_lat(i, this)
    nb = nb + this%grid%merid_length(k)
#else
    nb = nb + nb_winds_merid_lat(i, this)
#endif
    k = k +1
  end do

end function

!-------------------------------------------------------------------------------
!> Number of meridional winds at latitude i on the current partition
!> An interface is defined between 2 cells, with at least one of them being an 
!> interior cell
!-------------------------------------------------------------------------------
function nb_winds_merid_lat(i, this) result (nb)
  type (Partition) :: this
  integer, intent(in) :: i
  integer :: j_start, j_end, k
  integer :: j1, j2, j
  integer :: j_next1, j_next2
  integer :: nb

  nb = 0
  call interior_lon_indices(j_next1, j_next2, this%grid, i+1)

  !call interior_lon_indices(j_start, j_end, this%grid, i)
  call lon_indices(j_start, j_end, this, i)

  do j = j_start, j_end
    call south_neighbour_cell_Partition(j1, j2, i, j, this, i+1)

    ! Count only inside cells !
    do k = j1, j2
      if (to_or_from_interior_cell(k, i, j, i+1, this%grid)) nb = nb + 1
    end do
  end do


end function

!-------------------------------------------------------------------------------
!> Returns the first cell on latitude i, taking into accound the sector borders
!-------------------------------------------------------------------------------
function first_cell_position(i, this) result(istart)
  type (Partition) :: this
  integer, intent(in)  :: i
  integer :: pos, istart

  istart = 1
  if (has_single_partition()) then
    return
  end if

  ! Special case : sector boundary
  if (is_first_on_band(this) .and. &
    .not. is_north_south_ghost(i, this%grid)) istart = 2
end function

!-------------------------------------------------------------------------------
!> Returns the global position for the first cell on latitude i
!> Assumee constant dlon on the same latitude line
!-------------------------------------------------------------------------------
function first_cell_position_global(i, this) result(pos)
  type (Partition) :: this
  integer, intent(in)  :: i
  integer :: pos, istart

  if (has_single_partition()) then
    pos = 1
    return
  end if

  istart = 1
  ! Special case : sector boundary
  if (is_first_on_band(this) .and. &
    .not. is_north_south_ghost(i, this%grid)) istart = 2
  pos = local_to_global(i, istart, this%grid)
end function

!-------------------------------------------------------------------------------
!> Get the global position of the first meridion wind
!> @param i : latitude line in the partition
!> @param i_glob : latitude line in the global grid
!-------------------------------------------------------------------------------
function first_merid_position_global(i, i_glob, this) result(pos)
  type (Partition) :: this
  integer, intent(in) :: i, i_glob
  integer :: j, j1_loc, j2_loc, j1_glob, j2_glob
  integer :: l, pos, j_start
  integer :: local, global, first, last
!  integer :: local_contrib, global_contrib

  ! First ghost cell, except at the east boundary of the sector
  call lon_indices(first, last, this, i)

  pos = 0
  ! Count winds external to the partition with global topology
  ! Advance global position just before the first cell
  global = local_to_global(i, first, this%grid)
  do j = 1, global-1
    call south_neighbour_cell_Global(j1_glob, j2_glob, i_glob, j)
    if (j1_glob < 1 .or. j2_glob < 1) cycle
    pos = pos + j2_glob - j1_glob + 1
  end do

  ! Examine inside the partition (starting from ghost cells)
  j_start = first_interior_lon(i, this%grid)
  do j = first, j_start
    call south_neighbour_cell_Partition(j1_loc, j2_loc, i, j, this)
    call south_neighbour_cell_Global(j1_glob, j2_glob, i_glob, global)

    if (j1_glob < 1 .or. j2_glob < 1) then
      call print_error("no global south neighbour",&
      "first_merid_position_global", fname_partition)
    end if

    ! If local neighbours
    if (j1_loc > 0 .and. j2_loc > 0) then
      ! First we add global neighbours not counted by local topology
      pos = pos + (j2_glob - j1_glob + 1) - (j2_loc - j1_loc + 1)
      ! Then we examine local neighbours
      do l = j1_loc, j2_loc
        pos = pos + 1
        ! We stop at the first occurence
        if (to_or_from_interior_cell(l, i, j, i+1, this%grid)) return
        ! Otherwise, it is only a global neighbour
      end do
    else
      ! If no local neighbour, we still must update the global counter
      pos = pos + j2_glob - j1_glob + 1
    end if

    global = global + 1
  end do

end function

function first_interior_cell_position_global(i, this) result(pos)
  type (Partition) :: this
  integer, intent(in)  :: i
  integer :: pos, j_start

  if (has_single_partition()) then
    pos = 1
    return
  end if

  j_start = first_interior_lon(i, this%grid)
  pos = local_to_global(i, j_start, this%grid)
end function


!-------------------------------------------------------------------------------
!> Returns the global position for the last cell on latitude i
!-------------------------------------------------------------------------------
function last_cell_position_global(i, this) result(pos)
  type (Partition) :: this
  integer, intent(in)  :: i
  integer :: pos, iend
  double precision :: center

  if (has_single_partition()) then
    pos = this%grid%nb_lon(i)
    return
  end if

  iend = get_true_nb_lon_Partition(i, this)
  ! Special case : sector boundary for interior cells
  if (is_last_on_band(this) .and. .not. is_north_south_ghost(i, this%grid)) then
    iend = iend-1
  end if
  pos = local_to_global(i, iend, this%grid)

end function


!-----------------------------------------------------------------------------
!> Returns the number of interior cells inside sent to east neighbour
!-----------------------------------------------------------------------------
function nb_interior_cells_east(i, this) result (nb)
  type (Partition), target :: this
  type (Band_grid), pointer :: grid
  integer :: i, nb
  integer :: next_last, next_first
  integer :: last, first
  integer :: j_start, j_end, tmp
  integer :: nb_lon

  nb = 0
  grid => this%grid
  if (is_north_south_ghost(i, grid)) return

  nb = 1
  if (i < grid%nb_lat .and. .not. is_south_ghost(i+1, grid)) then
    !if (this%id == 37 .and. i == 3) then
    !  print *, "first loop", grid%nb_ghosts_east(i+1)
    !end if
    if (grid%nb_ghosts_east(i+1) < 2) goto 100

    ! Special case : we may have more than 1
    next_last = grid%nb_lon(i+1) - grid%nb_ghosts_east(i+1) + 1

    call north_neighbour_cell_Partition(j_start, tmp, i+1, next_last, this, i)
    nb_lon = grid%nb_lon(i) - grid%nb_ghosts_east(i)
    !if (this%id == 37 .and. i == 3) then
    !  print *, "(iter1) jstart j end ", j_start, nb_lon, "for", i
    !end if
    if (j_start <= nb_lon) then
      nb = max(nb, nb_lon - j_start + 1)
    end if
  end if

  100 if (i > 1 .and. .not. is_north_ghost(i-1, grid)) then
    if (grid%nb_ghosts_east(i-1) < 2) return

    ! Special case : we may have more than 1
    next_last = grid%nb_lon(i-1) - grid%nb_ghosts_east(i-1) + 1

    call south_neighbour_cell_Partition(j_start, tmp, i-1, next_last, this, i)
    nb_lon = grid%nb_lon(i) - grid%nb_ghosts_east(i)
    !if (this%id == 37 .and. i == 3) then
    !  print *, "(iter1) jstart j end ", j_start, nb_lon, "for", i
    !end if

    if (j_start <= nb_lon) then
      nb = max(nb, nb_lon - j_start + 1)
    end if
  end if

end function

!-----------------------------------------------------------------------------
!> Returns the number of interior cells inside sent to west neighbour
!-----------------------------------------------------------------------------
function nb_interior_cells_west(i, this) result (nb)
  type (Partition), target :: this
  type (Band_grid), pointer :: grid
  integer :: i, nb
  integer :: next_last, next_first
  integer :: last, first
  integer :: j_start, j_end, tmp
  integer :: nb_ghosts

  nb = 0
  grid => this%grid
  if (is_north_south_ghost(i, grid)) return

  nb = 1
  if (i < grid%nb_lat .and. .not. is_south_ghost(i+1, grid)) then
    if (grid%nb_ghosts_west(i+1) < 2) goto 100

    ! Special case : we may have more than 1
    next_first = grid%nb_ghosts_west(i+1)

    call north_neighbour_cell_Partition(tmp, j_end, i+1, next_first, this, i)
    nb_ghosts = grid%nb_ghosts_west(i)

    if (j_end > nb_ghosts) then
      nb = max(nb, j_end - nb_ghosts)
    end if
  end if

  100 if (i > 1 .and. .not. is_north_ghost(i-1, grid)) then
    if (grid%nb_ghosts_west(i-1) < 2) return

    ! Special case : we may have more than 1
    next_first = grid%nb_ghosts_west(i-1)

    call south_neighbour_cell_Partition(tmp, j_end, i-1, next_first, this, i)
    nb_ghosts = grid%nb_ghosts_west(i)

    if (j_end > nb_ghosts) then
      nb = max(nb, j_end - nb_ghosts)
    end if
  end if


end function



!-----------------------------------------------------------------------------
!> Accessors for mpi type. Interior or ghost cells are defined by side (0 and 1
!> respectively). Finding the correct indice inside the array is done by the 
!> caller
!-----------------------------------------------------------------------------
function get_mpi_border(this, side, neighb, direction) result (res)
  type (Partition) :: this
  integer :: direction
  integer :: neighb, res, side

  if (side == 0) then
#ifdef DEBUG
    if (neighb > size(this%grid%mpi_interior_cells)) then
      write (*,*) "cur", this%id, "indice", neighb, "side", side
      write (*,*) "direction", direction
      call print_error("Indice outside mpi interior cells array", &
        "get_mpi_border", fname_partition)
    end if
#endif

    res = this%grid%mpi_interior_cells(neighb)
    !write (*,*) "indice", neighb, "size", size(this%mpi_interior_cells)
  else
#ifdef DEBUG
    if (neighb > size(this%grid%mpi_ghost_cells)) then
      write (*,*) "cur", this%id, "indice", neighb, "side", side
      write (*,*) "direction", direction
      call print_error("Indice outside mpi ghost cells array", &
        "get_mpi_border", fname_partition)
    end if
#endif

    res = this%grid%mpi_ghost_cells(neighb)
    !write (*,*) "indice", neighb, "size", size(this%mpi_ghost_cells)
  end if
  !res = mpi_borders(pos)
end function

!-----------------------------------------------------------------------------
!> Creates and commit an new vector mpi type for an array of double precision
!-----------------------------------------------------------------------------
subroutine new_vector_MPI(blocklength, mpi_type)
  integer :: blocklength, mpi_type, ierr

  call mpi_type_contiguous(blocklength, mpi_double_precision, mpi_type, ierr)
  call check_mpi_error(ierr, "create mpi type", "new_vector_MPI", fname_partition)
  call mpi_type_commit(mpi_type, ierr)
  call check_mpi_error(ierr, "commit mpi", "new_vector_MPI", fname_partition)
end subroutine

!-----------------------------------------------------------------------------
!> Send dimensions for allocatable arrays (of the grid arrays from process 0 to process id)
!> @param id : id for receiving process
!> @param is_idle : true if we send to an idle process
!-----------------------------------------------------------------------------
subroutine send_dimensions_Partition(this, id, is_idle)
  type(Partition), target :: this
  type(Band_grid), pointer :: grid
  logical, intent(in) :: is_idle
  integer :: id, ierr, nb
  integer, allocatable :: dimensions(:)

  nb = 9
  allocate (dimensions(nb))

  ! Then the receiving process will know it is idle
  if (is_idle) then
    dimensions = -1
  else
    ! We need to instantiate the MPI type for each partition
    ! We only need to send 3 sizes
    grid => this%grid


    dimensions(1) = size(grid%dlat)
    dimensions(2) = size(grid%dlon)
    dimensions(3) = size(grid%first_cell)
    dimensions(4) = size(grid%north_offset)
    dimensions(5) = size(grid%south_offset)
    dimensions(6) = size(grid%lon_min)
    dimensions(7) = size(grid%nb_ghosts_east)
    dimensions(8) = size(grid%nb_ghosts_west)
    dimensions(9) = size(grid%nb_lon)
  end if

  ! Before sending the grid, we must send the dimensions of the arrays (as 
  ! they are allocatable). We sent it in a buffer
  call mpi_send(dimensions, nb, mpi_integer, id, 0, mpi_comm_world, ierr)
  call check_mpi_error(ierr, "send dimensions", "send_partition_data", &
    fname_partition)
end subroutine

!-----------------------------------------------------------------------------
!> Receive partition dimensions (partition + band_grid) from process id 
!> (should be 0  and allocate the corresponding data.
!-----------------------------------------------------------------------------
subroutine receive_dimensions_Partition(this, id)
  type(Partition), target :: this
  type(Band_grid), pointer :: cur
  integer :: ierr, id, nb
  integer :: status(mpi_status_size)
  integer :: mpi_band_grid
  integer, allocatable :: dimensions(:)
  logical :: is_idle

  ! Get data dimensions before receiving
  call mpi_probe(id, 0, mpi_comm_world, status, ierr);
  call check_mpi_error(ierr, "probing dimensions", "receive_dimensions_Partition", &
    fname_partition)
  call mpi_get_count(status, mpi_integer, nb, ierr);
  call check_mpi_error(ierr, "counting dimensions", "receive_dimensions_Partition", &
    fname_partition)

  if (nb /= 9) then
    call print_error("Dimensions should have 9 elements", &
      "receive_dimensions_Partition", fname_partition)
  end if
  allocate (dimensions(nb))

  ! Now we are ready to receive
  call mpi_recv(dimensions, nb, mpi_integer, id, 0, mpi_comm_world, status, &
    ierr)
  call check_mpi_error(ierr, "receive dimensions", "receive_partition_data", &
    fname_partition)

  is_idle = .False.
  if (dimensions(1) == -1) is_idle = .True.

  if (is_idle) then
    print *, "is idle"
    this%id = -2
  else
    ! Allocate the grid data
    allocate (this%grid%dlat(dimensions(1)))
    allocate (this%grid%dlon(dimensions(2)))
    allocate (this%grid%first_cell(dimensions(3)))
    allocate (this%grid%north_offset(dimensions(4)))
    allocate (this%grid%south_offset(dimensions(5)))
    allocate (this%grid%lon_min(dimensions(6)))
    allocate (this%grid%nb_ghosts_east(dimensions(7)))
    allocate (this%grid%nb_ghosts_west(dimensions(8)))
    allocate (this%grid%nb_lon(dimensions(9)))
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Send the current partition to the process dest
!> @param dest : target process
!-------------------------------------------------------------------------------
subroutine send_Partition(this, dest)
  type (Partition) :: this
  integer :: dest
  integer :: ierr, status(mpi_status_size)

  call mpi_send(this, 1, this%mpi_partition, dest, 0, mpi_comm_world, ierr)
  call check_mpi_error(ierr, "send partition", "send_Partition", &
    fname_partition)
end subroutine

!-------------------------------------------------------------------------------
!> Receive the partition from process source
!> @param source : sending process
!-------------------------------------------------------------------------------
subroutine receive_Partition(this, source)
  type (Partition) :: this
  integer :: source
  integer :: ierr, status(mpi_status_size)

  call mpi_recv(this, 1, this%mpi_partition, source, 0, mpi_comm_world, status, ierr)
  call check_mpi_error(ierr, "receive partition data", "receive_partition_data", &
    fname_partition)
end subroutine


!-----------------------------------------------------------------------------
!> Set the partition's position in the partitioning matrix
!-----------------------------------------------------------------------------
subroutine set_pos_Partition(this,i,j)
  type(Partition) :: this
  integer, intent(in) :: i,j
  this%i_pos = i
  this%j_pos = j
end subroutine

!-----------------------------------------------------------------------------
!> Create winds analytically
!-----------------------------------------------------------------------------
subroutine analytical_winds_Partition(this)
  type(Partition) :: this
  integer :: i_lat_prev
  character(LINE_WIDTH) :: test_case

  test_case = get_test_case()
  call analytical_winds(trim(test_case), this)

  !call check_ratio_is_set(this)
  call check_data(this)

end subroutine

subroutine check_data(this)
  type (Partition) :: this
  integer :: tracer

  call check_ratio_is_set(this)
  do tracer = 1, NB_TRACERS
    call check_ratio_undefined(this, tracer, IS_GRADIENT)
    call check_ratio_bounds(this, tracer)
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Check all interior cells are initialized
!-----------------------------------------------------------------------------
subroutine check_ratio_is_set(this)
  type (Partition) :: this
  integer :: tracer

  do tracer = 1, NB_TRACERS
    call check_ratio_undefined(this, tracer, IS_RATIO )
  end do

end subroutine



subroutine check_ratio_bounds(part, tracer)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: tracer
  integer :: i, j
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: ierr, rank, nb
  double precision :: cur
  character(8) :: nb_s, rank_s

  grid => part%grid
  call interior_lat_indices(i_start, i_end, grid)

  call mpi_comm_rank(mpi_comm_world, rank, ierr) 
  nb = 0
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, grid, i)

    do j = j_start, j_end
      cur = get_cell_ratio(i, j, tracer, grid)
      if (abs(cur) > 1+DBLE_PREC) nb = nb + 1
    end do
  end do

  if (nb > 0) then
    write (nb_s, '(i8)') nb
    write (rank_s, '(i8)') rank
    call print_warning("Cells with a ratio > or < 1 :"//nb_s//" rank="//rank_s)
  end if
end subroutine

subroutine check_ratio_undefined(part, tracer, array)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: tracer
  integer :: array, i, j
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: ierr, rank, nb
  double precision :: tmp, cur
  character(8) :: array_char = "ratio"
  character(8) :: nb_s, rank_s

  grid => part%grid
  if (array == IS_GRADIENT) array_char = "gradient"

  nb = 0
  call mpi_comm_rank(mpi_comm_world, rank, ierr) 
  call interior_lat_indices(i_start, i_end, grid)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, grid, i)
    do j = j_start, j_end
      if (array == IS_RATIO) then
        cur = get_cell_ratio(i, j, tracer, grid)
      else if (array == IS_GRADIENT) then
        cur = get_cell_gradient(i, j, tracer, grid)
      end if

      if (cur == UNDEFINED .and. rank == 0) then
        print *, "cell", i, j
      end if
      if (cur == UNDEFINED) nb = nb+1
    end do
  end do

  if (nb > 0) then
    write (nb_s, '(i8)') nb
    write (rank_s, '(i8)') rank

    call print_warning("Undefined "//array_char//" for "//nb_s//" cells, rank"//rank_s)
  end if
end subroutine


subroutine check_winds_div(this)
  type(Partition), target :: this
  type(Band_Grid), pointer :: grid
  integer :: i, i_start, i_end
  integer :: j, j_start, j_end
  integer :: k, k_n, k_p, l, nb_err
  double precision :: dlon, flux_merid, dlat
  double precision :: total

  grid => this%grid

  call interior_lat_indices(i_start, i_end, grid)
  ! Indices for the prev and next fluxes
  call starting_winds_indices(k_p, k_n, i_start, this)

  k = 2
  nb_err = 0

  dlat = get_dlat(grid)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, grid, i)

    call skip_first_merid_winds(k_n, k_p, i, this)
    do j = j_start, j_end

      flux_merid = sum_merid_fluxes_cell(i, j, k_p, k_n, this, grid)
      total = (grid%zonal_winds(k-1) - grid%zonal_winds(k))*dlat + flux_merid

      if (abs(total) > DBLE_PREC) then
        nb_err = nb_err + 1
      end if

      k = k + 1
    end do
    call skip_last_merid_winds(k_n, k_p, i, this)

    ! Switch to next line
    k = k + 1
  end do
  if (nb_err > 0) then
    print *, "partition", this%id, "nb errors", nb_err
    call print_warning("Winds div is not null")
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Starting winds indices for previous and next latitudes. Incremented as we go
!> Not the true indice, but updated later
!-------------------------------------------------------------------------------
subroutine starting_winds_indices(k_p, k_n, i_start, this)
  type (Partition) :: this
  integer :: k_p, k_n
  integer, intent(in) :: i_start
  k_p = 1
  k_n = 1

  if (i_start > 1 .and. has_north_ghost_cells(this%grid)) then
    k_n = nb_winds_merid_lat(i_start-1, this) + 1
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Increment meridional winds for the first interior cells
!-------------------------------------------------------------------------------
subroutine skip_first_merid_winds(k_n, k_p, i, this)
  type (Partition) :: this
  integer, intent(in) :: i
  integer :: k_n, k_p

  if (has_single_partition()) return
  call skip_next_merid_winds(k_n, i, i+1, EAST, this)
  call skip_next_merid_winds(k_p, i, i-1, EAST, this)
end subroutine

!-------------------------------------------------------------------------------
!> Increment meridional winds for skipping last cells
!-------------------------------------------------------------------------------
subroutine skip_last_merid_winds(k_n, k_p, i, this)
  type (Partition) :: this
  integer, intent(in) :: i
  integer :: k_n, k_p

  if (has_single_partition()) return
  call skip_next_merid_winds(k_n, i, i+1, WEST, this)
  call skip_next_merid_winds(k_p, i, i-1, WEST, this)
end subroutine

!-------------------------------------------------------------------------------
!> Skip next meridional winds, i.e winds from ghost to interior cells, by
!> incrementing winds indice.
!> @param dir : will skip either east or west section
!> @param i_next : previous or next latitude
!-------------------------------------------------------------------------------
subroutine skip_next_merid_winds(k, i, i_next, dir, this)
  type (Partition) :: this
  integer :: k, l, p
  integer :: j_first, j_last, j1, j2
  integer :: j_start, j_end
  integer, intent(in) :: i, i_next, dir

  call lon_indices(j_first, j_last, this, i)
  ! Interior only
  call interior_lon_indices(j_start, j_end, this%grid, i)

  if (dir == EAST) then
    do l = j_first, j_start - 1
      if (i_next < i) then
        call north_neighbour_cell_Partition(j1, j2, i, l, this, i_next)
      else
        call south_neighbour_cell_Partition(j1, j2, i, l, this, i_next)
      end if

      if (j1 < 1 .or. j2 < 1) cycle
      do p = j1, j2
        if (to_or_from_interior_cell(p, i, l, i_next, this%grid)) then
          k = k + 1
        end if
      end do
    end do

  else if (dir == WEST) then
    do l = j_end+1, j_last
      if (i_next < i) then
        call north_neighbour_cell_Partition(j1, j2, i, l, this, i_next)
      else
        call south_neighbour_cell_Partition(j1, j2, i, l, this, i_next)
      end if

      if (j1 < 1 .or. j2 < 1) cycle
      do p = j1, j2
        if (to_or_from_interior_cell(p, i, l, i_next, this%grid)) then
          k = k + 1
        end if
      end do
    end do
  else 
    call print_error("Wrong direction : east or west","skip_next_winds", &
      fname_partition)
  end if
end subroutine




!-------------------------------------------------------------------------------
!> Correct 2d winds in order to have a null divergence
!> Induces sequential order of operation
!-------------------------------------------------------------------------------
subroutine correct_winds(this)
  type(Partition) :: this
  integer :: neighb, ierr
  logical :: cond, status(mpi_status_size)

  ! Needed to correct winds in analytical version
  ! TODO avoid multiple call
  if (.not. has_single_partition()) then
    call create_mpi_border_zwinds(this%grid, WEST)
    call create_mpi_border_zwinds(this%grid, EAST)
  end if

  !print *, "Correcting winds"
  call correct_winds_merid(this)

  ! Sequential as we need the ghost cells
  if (.not. has_single_partition()) then
    ! First receive ghost cells
    if (this%j_pos > 1) then
      call west_neighbours(this, neighb)
      call mpi_recv(this%grid%zonal_winds, 1, this%grid%mpi_zwinds_west, &
        neighb, 0, mpi_comm_world, status, ierr)
      call check_error(ierr,"receiving ghost cells","correct_winds", &
        fname_partition)
    end if
  end if

  ! Then compute
  call correct_winds_zonal(this)

  ! And send
  if (.not. has_single_partition()) then
    if (.not. is_last_on_band_all(this)) then
      call east_neighbours(this, neighb)
      call mpi_ssend(this%grid%zonal_winds, 1, this%grid%mpi_zwinds_east, &
        neighb, 0, mpi_comm_world, ierr)
    end if
  end if

  ! Periodic boundaries
  if (.not. has_single_partition()) then
    if (is_last_on_band_all(this)) then
      call east_neighbours(this, neighb)
      call mpi_recv(this%grid%zonal_winds, 1, this%grid%mpi_zwinds_east, &
        neighb, 0, mpi_comm_world, status, ierr)

    else if (this%j_pos == 1) then
      call west_neighbours(this, neighb)
      call mpi_ssend(this%grid%zonal_winds, 1, this%grid%mpi_zwinds_west, &
        neighb, 0, mpi_comm_world, ierr)
    end if
  end if

  call check_winds_div(this)
end subroutine

!-------------------------------------------------------------------------------
!> Correct meridional winds for a divergence null. We want the sum of fluxes on
!> a latitude line to be null.
!> Input : positive merid winds are towards the north
!-------------------------------------------------------------------------------
subroutine correct_winds_merid(this)
  type(Partition), target :: this
  type(Band_Grid), pointer :: grid
  integer :: i, i_start, i_end
  integer :: nb_lat, cur_lat
  integer :: j, k, ierr
  double precision, allocatable :: div(:), cur_div(:)
  double precision :: dlat, dl, lat
  integer :: nb_procs

  grid => this%grid

  !print *, "before", grid%merid_winds(1:6)
  ! Change sign of winds : positive means towards south
  grid%merid_winds = -grid%merid_winds

  ! Sum of all flows on each latitude line, except the last
  nb_lat = get_nb_lat_Global_grid()
  allocate(div(nb_lat-1))
  allocate(cur_div(nb_lat-1))

  call interior_lat_indices(i_start, i_end, grid)
  if (is_southmost(this)) i_end = i_end - 1

  i_start=1
  ! Compute the divergence on all latitude lines
  call set_local_divergence_latitude(cur_div, i_start, i_end, this, grid)

  ! Update global vector
  div = 0
  call mpi_allreduce(cur_div, div , size(div), mpi_double_precision, mpi_sum,&
    mpi_comm_world, ierr)
  call check_error(ierr,"trying to sum all divergence","correct_winds_merid", &
    fname_partition)

  ! Divide by the line elements
  dlat = get_dlat(grid)
  do i = 1, nb_lat-1
    dl = 360.d0*cos((90.d0 -i*dlat)*pi/180.d0) 
    div(i) = div(i) / dl
  end do

  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)

  ! Correct winds
  k=1
  lat = cell_south_lat(i_start, grid)

  cur_lat = find_global_lat(lat, grid)

  do i = i_start, i_end
    do j = 1, nb_winds_merid_lat(i, this)

      grid%merid_winds(k) = grid%merid_winds(k)-div(cur_lat)
      k = k + 1
    end do
    cur_lat = cur_lat + 1
  end do

  deallocate(div)
  deallocate(cur_div)
end subroutine

!-------------------------------------------------------------------------------
!> Correct zonal winds such as the sum of fluxes inside a cell is null.
!> Requires meridional winds to be corrected first
!-------------------------------------------------------------------------------
subroutine correct_winds_zonal(this)
  type(Partition), target :: this
  type(Band_Grid), pointer :: grid
  integer :: i, i_start, i_end
  integer :: j, j_start, j_end
  integer :: k, k_n, k_p, l
  integer :: k_start
  double precision :: dlon, flux_merid, dlat
  double precision :: total, tmp
  integer :: nb_procs, ierr

  grid => this%grid

  call interior_lat_indices(i_start, i_end, grid)
  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)

  ! Indices for the prev and next fluxes
  call starting_winds_indices(k_p, k_n, i_start, this)
  k = 2
  k_start = 1

  dlat = get_dlat(grid)
  do i = i_start, i_end
    !do i = i_end, i_start, -1
    call interior_lon_indices(j_start, j_end, grid, i)

    call skip_first_merid_winds(k_n, k_p, i, this)

    do j = j_start, j_end

      flux_merid = sum_merid_fluxes_cell(i, j, k_p, k_n, this, grid)

      ! We must have zonal_west - zonal_east + merid_up - merid_down = 0
      ! Zonal flux = zonal winds * dlat
      grid%zonal_winds(k) = grid%zonal_winds(k-1) + flux_merid/dlat 

      k = k + 1
    end do

    call skip_last_merid_winds(k_n, k_p, i, this)

    ! For sequential version, periodic boundaries
    if (has_single_partition()) then
      grid%zonal_winds(k-1) = grid%zonal_winds(k_start)
    end if

    k_start = k
    k = k + 1
  end do

end subroutine

!-------------------------------------------------------------------------------
! Compute sum of all meridional fluxes for a cell
!-------------------------------------------------------------------------------
function sum_merid_fluxes_cell(i, j, k_p, k_n, this, grid) result (flux_merid)
  type (Partition) :: this
  type (Band_grid) :: grid
  integer, intent(in) :: i, j
  integer :: k_n, k_p, l
  integer :: j_p1, j_p2, j_n1, j_n2
  double precision :: dlon, flux_merid
  integer :: nb_procs, ierr
  double precision :: tmp

  call north_neighbour_cell_Partition(j_p1, j_p2, i, j, this, i-1)
  call south_neighbour_cell_Partition(j_n1, j_n2, i, j, this, i+1)
  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)

  flux_merid = 0
  if (j_p1 > 0 .and. j_p2 > 0) then
    do l = j_p1, j_p2
      if (to_or_from_interior_cell(l, i, j, i-1, grid)) then
        dlon = cell_interface_length(i, j, i-1, l, grid)
        flux_merid = flux_merid + grid%merid_winds(k_p)*dlon
      end if

      k_p = k_p + 1
    end do
  end if

  tmp = flux_merid

  if (j_n1 > 0 .and. j_n2 > 0) then
    do l = j_n1, j_n2
      if (to_or_from_interior_cell(l, i, j, i+1, grid)) then
        dlon = cell_interface_length(i, j, i+1, l, grid)
        flux_merid = flux_merid - grid%merid_winds(k_n)*dlon
      end if

      k_n = k_n + 1
    end do
  end if

end function


!-------------------------------------------------------------------------------
!> Compute part of the divergence on the latitude line 
!> @param cur_div : array storing the final divergence
!-------------------------------------------------------------------------------
subroutine set_local_divergence_latitude(cur_div, i_start, i_end, this, grid)
  type (Partition) :: this
  type (Band_grid) :: grid
  integer, intent(in) :: i_start, i_end
  integer :: i, neighb1, neighb2
  integer :: j, j_start, j_end
  integer :: l, k_p, k_n
  double precision :: dlon, lat
  double precision :: cur_div(:)
  integer :: cur_lat, nb_procs, ierr

  call starting_winds_indices(k_p, k_n, i_start, this)

  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  cur_div = 0.

  ! First global lat
  lat = cell_south_lat(i_start, grid)
  cur_lat = find_global_lat(lat, grid)

  do i = i_start, i_end
    if (cur_lat > size(cur_div)) then
      call print_error("latitude outside global array", &
      "set_local_divergence_latitude", fname_partition)
    end if

    call skip_first_merid_winds(k_n, k_p, i, this)
    call interior_lon_indices(j_start, j_end, grid, i)
    if (has_north_ghost_cells(grid) .and. i == 1) then
      j_start = 1
      j_end = grid%nb_lon(i)
    end if

    do j = j_start, j_end
      call south_neighbour_cell_Partition(neighb1, neighb2, i, j, this, i+1)
      do l = neighb1, neighb2
        if (to_or_from_interior_cell(l, i, j, i+1, grid)) then
          dlon = cell_interface_length(i, j, i+1, l, grid)
          cur_div(cur_lat) = cur_div(cur_lat) + grid%merid_winds(k_n)*dlon
          k_n = k_n + 1
        end if
      end do
    end do
    call skip_last_merid_winds(k_n, k_p, i, this)

    ! Hack : at the equator, divergence is computed twice, so we cancel (for the sum
    ! later) one of them 
    if (i == 1 .and. cur_lat == get_nb_lat2_Global_grid()) then
      cur_div(cur_lat) = 0
    end if
    cur_lat = cur_lat + 1
  end do

end subroutine


subroutine write_concentration(this, id_file, tracer)
  type(Partition) :: this
  integer, intent(in) :: id_file, tracer

  call write_ratio_Band_grid(this%grid, id_file, tracer)
end subroutine

subroutine write_partition_id(this, id_file)
  type(Partition) :: this
  integer, intent(in) :: id_file

  call write_id_Band_grid(this%grid, this%id, id_file)
end subroutine



!-----------------------------------------------------------------------------
!> Write the cell winds interpolated at the cell south boundary
!> Format : (u, v, lat, lon)
!> Done here as we need the neighbours
!-----------------------------------------------------------------------------
subroutine interpolate_winds_to_south(this, id_file)
  type(Partition), target :: this
  type(Band_grid), pointer :: grid
  integer, intent(in) :: id_file
  integer :: i, j, k
  character(10) :: fformat
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: k_p, k_n
  double precision :: u, v, lat, lon

  grid => this%grid

  fformat = '(4'//DBLE_FORMAT//')'

  call interior_lat_indices(i_start, i_end, grid)
  if (.not. has_south_ghost_cells(grid)) i_end = i_end - 1

  k_n = 1
  k = 1
  do i = i_start, i_end

    call interior_lon_indices(j_start, j_end, grid, i)
    lat = cell_south_lat(i, grid)

    do j = j_start, j_end
      lon = cell_center_lon(i, j, grid)
      u = interpolate_cell_zonal_winds_to_center(i, j, k, grid)
      v = interpolate_cell_merid_winds_to_south(i, j, k_n, this)

      write (id_file,fformat) u, v, lat, lon
    end do
    ! Switch to next line
    k = k + 1
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Write the cell winds interpolated at the cell center
!> Format : (u, v, lat, lon)
!> Done here as we need the neighbours
!-----------------------------------------------------------------------------
subroutine interpolate_winds_to_center(this, id_file)
  type(Partition), target :: this
  type(Band_grid), pointer :: grid
  integer, intent(in) :: id_file
  integer :: i, j, k
  character(10) :: fformat
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: k_p, k_n
  double precision :: u, v, lat, lon

  grid => this%grid

  fformat = '(4'//DBLE_FORMAT//')'
  ! Does not print ghost cells
  call interior_lat_indices(i_start, i_end, grid)

  call starting_winds_indices(k_p, k_n, i_start, this)
  k = 1
  do i = i_start, i_end

    call interior_lon_indices(j_start, j_end, grid, i)
    lat = cell_center_lat(i, grid)

    do j = j_start, j_end
      lon = cell_center_lon(i, j, grid)
      u = interpolate_cell_zonal_winds_to_center(i, j, k, grid)
      v = interpolate_cell_merid_winds_to_center(i, j, k_p, k_n, this)

      write (id_file,fformat) u, v, lat, lon
    end do
    ! Switch to next line
    k = k + 1
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Interpolate merid winds of cell at the south border
!> @param k : current position in the meridional winds array of south border
!-------------------------------------------------------------------------------
function interpolate_cell_merid_winds_to_south(i, j, k, this) result(v)
  type (Partition) :: this
  integer, intent(in) :: i, j
  integer :: k
  integer :: j_n1, j_n2
  integer :: nb_v, l
  double precision :: v
  double precision :: coef(3), length

  coef = 0.
  call south_neighbour_cell_Partition(j_n1, j_n2, i, j, this, i+1)
  v = 0.
  if (j_n1 > 0 .and. j_n2 > 0) then
    do l = j_n1, j_n2
      coef(l-j_n1+1) = cell_interface_length(i, j, i+1, l, this%grid)
    end do
    length = sum(coef)
    coef = coef/length

    do l = j_n1, j_n2
      !v = coef(l-j_n1+1) * this%grid%merid_winds(k)
      v = this%grid%merid_winds(k)
      k = k + 1
    end do
    v = v / (j_n2 - j_n1 + 1)
  end if

end function


!-------------------------------------------------------------------------------
!> Interpolate merid winds of cells borders into v
!> @param k_p, k_n : current position in the meridional winds array of north and
!> south border
!-------------------------------------------------------------------------------
function interpolate_cell_merid_winds_to_center(i, j, k_p, k_n, this) result(v)
  type (Partition) :: this
  integer, intent(in) :: i, j
  integer :: k_p, k_n
  integer :: j_p1, j_p2, j_n1, j_n2
  integer :: nb_v, l
  double precision :: v, v1, v2
  double precision :: coef(3), length

  call north_neighbour_cell_Partition(j_p1, j_p2, i, j, this, i-1)
  v = 0.
  coef = 0.
  v1 = 0.
  v2 = 0.
  if (j_p1 > 0 .and. j_p2 > 0) then
    ! Compute coefficients
    do l = j_p1, j_p2
      coef(l-j_p1+1) = cell_interface_length(i, j, i-1, l, this%grid)
    end do
    length = sum(coef)
    coef = coef/length

    ! Take the weighted sum
    do l = j_p1, j_p2
      v1 = coef(l-j_p1+1) * this%grid%merid_winds(k_p)
      k_p = k_p + 1
    end do

  end if

  coef = 0.
  call south_neighbour_cell_Partition(j_n1, j_n2, i, j, this, i+1)
  if (j_n1 > 0 .and. j_n2 > 0) then
    do l = j_n1, j_n2
      coef(l-j_n1+1) = cell_interface_length(i, j, i+1, l, this%grid)
    end do
    length = sum(coef)
    coef = coef/length
    !    if (i == 30) then
    !      print *, "j=", j, "coef=",coef
    !    end if

    do l = j_n1, j_n2
      v2 = coef(l-j_n1+1) * this%grid%merid_winds(k_n)
      k_n = k_n + 1
    end do
  end if

  v = 0.5d0*(v1 + v2)
end function

!-------------------------------------------------------------------------------
!> Interpolate zonal winds of cells borders into u
!> @param k : current position in the meridional winds array of west border
!-------------------------------------------------------------------------------
function interpolate_cell_zonal_winds_to_center(i, j, k, grid) result(u)
  type (Band_grid) :: grid
  integer, intent(in) :: i, j
  integer :: k
  double precision :: u

  u = (grid%zonal_winds(k) + grid%zonal_winds(k+1))/2.
  k = k + 1
end function



subroutine print_Partition(this)
  type(Partition) :: this

  write (*,'(a,i5)', advance="no") "Partition id=",this%id
  write (*,'(a,2i5,a)') " (i_pos,j_pos)=(",this%i_pos,this%j_pos," )"
  write (*,'(a,i5)') " Grid data size = ", size(this%grid%nb_lon)
  write (*,'(a,i5)') " overlap = ", this%overlap
  write (*,'(a,i5)') " zone = ", this%overlap
  call print_Band_grid(this%grid)

end subroutine


subroutine write_center_Partitition(this,id_file)
  type(Partition) :: this
  integer, intent(in) :: id_file
  double precision :: mid_lat, mid_lon

  call get_center_Band_grid(this%grid,mid_lat,mid_lon)
  write (id_file,'(f10.5,a,f10.5)',advance="no") mid_lat," ",mid_lon
end subroutine

!-----------------------------------------------------------------------------
!> Gives the west neighbour id >= 1 (-1 if none). 
!> This is based on an analytical partitioning.
!-----------------------------------------------------------------------------
subroutine west_neighbours(this, west)
  type(Partition) :: this
  integer :: west
  integer :: i_west, j_west
  west = -1
  call west_partition_wrapper(this, i_west, j_west)
  !if (this%id == 1) write (*,*) "west of", this%id,"at", i_west, j_west
  if (i_west > -1 .and. j_west > -1) then
    west = compute_id_from_pos(i_west, j_west)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Compute the id (>= 0) of the partition from its indices (> 0) in the matrix 
!> position
!>@param i, j : partition position in the matrix (in all sectors)
!-----------------------------------------------------------------------------
function compute_id_from_pos(i, j) result(id)
  integer, intent(in) :: i, j
  integer :: id, zone, j_corr
  integer :: nb_parts

  ! i and j > 0, while id >= 0
  nb_parts = nb_parts_on_band_sector(i)
  j_corr = translate_jpos(j, nb_parts) 
  ! The id is computed in the first sector 
  id = sum_partitions_sector(i-1) + j_corr - 1
  !write (*,*) "id in 1st" id 
  ! Must recompute the zone of the partition
  zone = compute_zone_from_pos(i, j, nb_parts) 

  nb_parts = get_nb_partitions_Configuration(1)
  ! And we must correct according to the zone
  if (on_northern_hemisphere(zone)) then
    id = id + (zone - 1)*nb_parts
    !if (i == 2 .and. j == 6) write (*,*) "fater", id
  else
    ! Must all the partitions on the northern hemisphere
    ! We only add the number of partitions on 2 sectors as the first is already
    ! counted by sum_partitions_sector
    id = id + 2*nb_parts
    nb_parts = get_nb_partitions_Configuration(4)
    id = id + (zone - 4)*nb_parts
  end if

end function

!-----------------------------------------------------------------------------
!> Recompute the zone from the partition position in the matrix
!>@param nb_parts : number of partitions on band i on a sector
!-----------------------------------------------------------------------------
function compute_zone_from_pos(i, j, nb_parts) result (zone)
  integer, intent(in) :: i, j, nb_parts
  integer :: zone

  ! Search the zone on the northern hemisphere
  zone = 1
  if (j > nb_parts) then
    zone = int((j-1) / nb_parts) + 1
  end if
  ! Correct if we are on the southern hemisphere
  if (i > get_nb_bands_Configuration(1)) then
    zone = zone + 3
  end if
end function

!-----------------------------------------------------------------------------
!> Computes the sector from the position in the matrix array. Used for neighbour
!> computation. 
!-----------------------------------------------------------------------------
function sector_from_pos(i, j) result(sector)
  integer, intent(in) :: i, j
  integer :: sector

  if (i < 1 .or. i > get_total_nb_bands_Configuration_sector()) then
    call print_error("Wrong partition position", "sector_from_pos",&
      fname_partition)
  end if

  sector = int((j-1)/nb_parts_on_band_sector(i)) + 1
  !write (*,*) "sector for", i,j,"is", sector
end function

!-----------------------------------------------------------------------------
!> Compute the west neighbour.
!> Contrary to north-south neighbours, we compute directly the global neighbour
!> @param this : current partition
!> @param i_east, j_east : position of the neighbour in the matrix
!-----------------------------------------------------------------------------
subroutine west_partition_wrapper(this, i_west, j_west)
  type(Partition) :: this
  integer :: i_west, j_west
  integer :: nb_parts

  i_west = this%i_pos
  j_west = -1

  if (this%j_pos > 1) then
    j_west = this%j_pos - 1
  else
    if (has_single_partition()) return

    nb_parts = nb_parts_on_band(i_west)
    ! Periodic boundary conditions
    j_west = nb_parts
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Gives the east neighbour id >= 0 (-1 if none). 
!> This is based on an analytical partitioning.
!-----------------------------------------------------------------------------
subroutine east_neighbours(this, east)
  type(Partition) :: this
  integer :: east, i_east, j_east
  east = -1

  call east_partition_wrapper(this, i_east, j_east)
  if (i_east > -1 .and. j_east > -1) then
    east = compute_id_from_pos(i_east, j_east)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Gives the east neighbour indices inside the position matrix
!> Contrary to north-south neighbours, we compute directly the global neighbour
!> @param this : current partition
!> @param i_east, j_east : position of the neighbour in the matrix
!----------------------------------------------------------------------------
subroutine east_partition_wrapper(this, i_east, j_east)
  type(Partition) :: this
  integer :: i_east, j_east
  integer :: nb_parts
  integer :: nb_bands

  j_east = -1
  i_east = this%i_pos

  nb_parts = 3*nb_parts_on_band_sector(i_east)
  !if (this%id == 5) write (*,*) "nb parts for 5", nb_parts
  if (this%j_pos < nb_parts) then
    j_east = this%j_pos + 1
  else
    if (has_single_partition()) return

    ! Periodic boundary conditions
    j_east = 1
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Gives the north neighbours ids >= 0 (-1 if none). 
!> This is based on an analytical partitioning.
!> The neighbours are between k_start and k_end
!-----------------------------------------------------------------------------
subroutine north_neighbours(this, id_start, id_end)
  type(Partition) :: this
  integer :: i_north, j_start, j_end
  !integer :: sum_prev
  integer :: id_start, id_end
  id_start = -1
  id_end = -1

  call north_partition_wrapper(this, i_north, j_start, j_end)
  !if (this%id == 5) write (*,*) "neighb north", j_start, j_end
  !if (this%id == 16) write (*,*) "north=", j_start, j_end
  if (i_north > -1 .and. j_start > -1) then
    id_start = compute_id_from_pos(i_north, j_start)
  end if
  if (i_north > -1 .and. j_end > -1) then
    id_end = compute_id_from_pos(i_north, j_end)
  end if
  !if (this%id == 29) write (*,*) "id after", id_start, id_end

end subroutine

!-----------------------------------------------------------------------------
!> Gives the indice in the mpi border array
!> @param direction : an integer from 0 to 3 for north to west
!> @param neighb : the neighbour, between 1 and the number of neighbours in this 
!> direction
!> The call must ensure there is a neighbour.
!-----------------------------------------------------------------------------
function border_indice(this, direction, neighb) result (i)
  type(Partition) :: this
  integer :: neighb, i, direction

  select case (direction)
  ! North
case (0)
  i = neighb
  ! East
case (1)
  i = nb_north_neighbours(this) + 1
  ! South
case (2)
  i = nb_north_neighbours(this) + neighb
  !if (this%id == 2) write (*,*) "nb neighb", nb_neighbours_Partition(this), "neigh", neighb
  ! Ignore the case of a single partition 
  if (.not. has_single_partition()) i = i + 1
  !if (this%id == 2) write (*,*) "res", i, "dir", direction
  ! West
case (3)
  i = nb_neighbours_Partition(this)
case default
  call print_error("Wrong border indice", "border_indice", fname_partition)
end select
end function

!-------------------------------------------------------------------------------
!> Returns the neighbour according to a integer (0 = north, 3 = west)
!-------------------------------------------------------------------------------
subroutine neighbour_from_direction(neighb, neighb2, this, direction)
  type (Partition) :: this
  integer :: neighb, neighb2, direction

  select case(direction)
  case (0)
    call north_neighbours(this, neighb, neighb2)
  case (1)
    call east_neighbours(this, neighb)
    neighb2 = neighb
  case (2)
    call south_neighbours(this, neighb, neighb2)
  case (3)
    call west_neighbours(this, neighb)
    neighb2 = neighb
  case default
    call print_error("Not a correct direction", "neighbour_from_direction", &
      fname_partition)
  end select
end subroutine 

!-----------------------------------------------------------------------------
!> Wrapper for north_partition_sector. Once we computed the neighbour, we correct
!> the position as the analytical partitioning only gives the neighbours on 
!> the first sector
!-----------------------------------------------------------------------------
subroutine north_partition_wrapper(this, i_north, j_start, j_end)
  type(Partition) :: this
  integer :: i_north, j_start, j_end
  integer :: sector, nb_parts

  call north_partition_sector(this, i_north, j_start, j_end)
  if (j_start > 0 .and. j_end > 0) then
    sector = sector_from_zone(this%zone)
    nb_parts = nb_parts_on_band_sector(i_north)
    j_start = j_start + (sector - 1)*nb_parts
    j_end = j_end + (sector - 1)*nb_parts
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Find the north neighbour with the indices inside the position matrix
!> The neighbours will range from (i_north, j_start) to (i_north, j_end)
!> Returns the neighbours on the first sector, so you have to correct the
!> results
!-----------------------------------------------------------------------------
subroutine north_partition_sector(this, i_north, j_start, j_end)
  type(Partition) :: this
  integer :: i_north, j_start, j_end
  integer :: i, j, j_corr
  integer :: nb_parts, nb_parts_prev

  i_north = -1
  j_start = -1
  j_end = -1

  ! No north neighbour. Must be done first.
  if (is_northmost(this)) return

  i = this%i_pos
  nb_parts = nb_parts_on_band_sector(i)
  j_corr = translate_jpos(this%j_pos, nb_parts)

  i_north = i - 1

  ! Last band on north hemisphere is a special case
  if (on_northern_hemisphere(this%zone) .and. is_leftover_band(i)) then 
    call north_south_partition_from_leftover(this, j_start, j_end, 1)
    !if (this%id == 73) write (*,*) "plouf", j_start, j_end, this%i_pos
    return
  end if

  ! If the previous band is on the other side of the equator. Must happen before
  ! the following test (as both can returns true).
  if (is_just_after_equator(this%i_pos)) then
    call north_south_partition_equator(this, j_start, j_end, i-1)
    return
  end if

  if (on_southern_hemisphere(this%zone) .and. is_leftover_band(i-1)) then 
    call north_south_partition_last_square(this, j_start, j_end, -1)
    return
  end if

  !if (this%i_pos == 3) write (*,*) "blarf", j_start, j_end
  !if (this%id == 7) write (*,*) "blarf", j_start, j_end

  ! Number of partitions on the band 
  nb_parts_prev = nb_square_parts_on_band(i-1)
  ! Only local so must update

  ! Southern hemisphere : north relation becomes south (actually, it is the
  ! same)
  if (on_northern_hemisphere(this%zone)) then 
    call north_partition_reduced(j_corr, nb_parts, nb_parts_prev, j_start, j_end)
  else
    call south_partition_reduced(j_corr, nb_parts, nb_parts_prev, j_start, j_end)
    !if (this%id == 65) write (*,*) "there after", j_start, j_end, this%overlap, i,j
  end if

  call correct_neighbours_resized(j_start, j_end, this, 1, j_corr)
end subroutine

!-----------------------------------------------------------------------------
!> If square partitions have been resized, we must correct the neighbours
!> The integer tell us which digit to use, therefore the direction (1 for 
!> north, 2 for south). See the definition of the class for more details.
!> @param j_corr : position on the band (first sector only)
!-----------------------------------------------------------------------------
subroutine correct_neighbours_resized(j_start, j_end, this, pos, j_corr)
  type (Partition) :: this
  integer :: j_start, j_end
  integer, intent(in) :: pos, j_corr
  integer :: nb_parts, mid, nb
  integer :: dir_left, dir_right

  dir_left = 1
  ! For south, the translation is in the other direction
  if (pos == 2) dir_left = -1

  nb_parts = nb_parts_on_band_sector(this%i_pos)
  mid = (nb_parts + 1)/2
  ! On the other half, the translation is in the other direction
  if (j_corr > mid) dir_left = -dir_left

  dir_right = dir_left
  if (j_corr == mid) dir_right = -dir_left

  if (on_southern_hemisphere(this%zone)) then
    dir_left = -dir_left
    dir_right = -dir_right
  end if

  ! Southern hemisphere, the direction of south and north are inverted
  nb = this%overlap / (10**(pos-1))
  nb = mod(nb, 10)
  ! Remove/add a left neighbour
  if (nb == 1 .or. nb == 3) j_start = j_start + dir_left
  ! Add/remove a right neighbour
  if (nb == 2 .or. nb == 3) j_end = j_end + dir_right

end subroutine

!-----------------------------------------------------------------------------
!> Based on the analytical partitioning, we can find the interface length for
!> the ghost cells of the neighbour of a partition. It computes the starting 
!>_indice inside the current grid and the length.
!> To do it analytically avoids to store the list of neighbours.
!> We need the matrix position of the neighbour.
!> We assume there is a neighbour.
!-----------------------------------------------------------------------------
function interface_length(this,i_neighb, j_neighb) result(i_length) 
  type (Partition) :: this
  integer :: i_start, j_corr
  integer :: i_neighb, j_neighb
  integer :: i, i_length
  integer :: i_start_loc
  integer :: nb_parts

  i = this%i_pos
  i_start = -1

  ! The neighbour column can be outside the first sector, so we translate it
  ! (easier to do it here than in the following functions)
  nb_parts = nb_parts_on_band_sector(i_neighb)
  j_corr = translate_jpos(j_neighb, nb_parts)

  ! Get the length of the interface and a local starting indice
  ! North neighbour
  if (i_neighb < i) then
    call interface_north_partition(i_start, i_start_loc, i_length, this, &
      j_corr) 
    ! No need to update north (local is already global)
    ! South neighbour
  else if (i_neighb > i) then
    call interface_south_partition(i_start, i_start_loc, i_length, this, &
      j_corr) 
    ! West or east (they are the same)
  else
    i_length = get_nb_lat_Partition(this)
  end if
end function

!-----------------------------------------------------------------------------
!> Same as interface_length but for interior cells
!> Returns the position of the first cells and the interface length
!>@param i_neighb, j_neighb : position in the partition matrix
!> j_neighb will be translated in the first sector
!-----------------------------------------------------------------------------
subroutine interface_interior_cells(i_start, i_length, this, i_neighb, &
    j_neighb ) 
  type (Partition) :: this
  integer :: i_neighb, j_neighb
  integer :: i, i_length, i_start
  integer :: nb_parts, j_corr

  i = this%i_pos

  ! North neighbour
  if (i_neighb < i) then
    call interface_interior_north_partition(i_start, i_length, this, j_neighb)
    ! South neighbour
  else if (i_neighb > i) then
    !i_length = interface_south_length_interior_cells(this, j_neighb)
    call interface_interior_south_partition(i_start, i_length, this, j_neighb)
    ! West or east
  else
    i_length = get_nb_lat_Partition(this)
    i_start = 1
  end if
end subroutine

!
!-----------------------------------------------------------------------------
!> Compute the width of a partition on the last band and the previous width (if
!> last band)
!> All variables refer to the leftover band
!-----------------------------------------------------------------------------
subroutine width_leftover(width_prev, width, j_pos, nb_cells, &
    nb_parts_band) 
  integer :: width, width_prev, j_pos
  integer :: nb_cells, nb_parts_band

  width_prev = nb_cells / nb_parts_band
  width = width_prev
  ! For the last partition, width is different
  if (j_pos > nb_parts_band - 1) then
    width = nb_cells  - (j_pos - 1)*width
  end if
end subroutine

!-----------------------------------------------------------------------------
!> From the leftover band, we retrieve the width of the "square" partitions
!> before an eventual resizing
!-----------------------------------------------------------------------------
function width_from_leftover(i_pos) result(width)
  integer :: width
  integer :: nb_lat, i_pos

  nb_lat = get_nb_lat2_Global_grid()
  ! We know there is a leftover band, so easy to find the height (we know we
  ! must substract 1

  width = (nb_lat-1) / (i_pos - 1)
end function

!-----------------------------------------------------------------------------
!> Find the interface length between a partition and its northern neighbour (for
!> the neighbour).
!> Also returns the first indice of the ghost cells inside the current partition grid.
!> Beware, the second indice is only local (for the corresponding band), while
!> the first correspond to the latitude band
!> Depends entirely of the analytical partitioning
!> Same as south, except that the discussion about the leftover band is on the
!> width, not the left/right neighbour cells.
!> @param j_neighb : neighbour column. Can be on any sector, we correct it
!> anyway.
!-----------------------------------------------------------------------------
subroutine interface_north_partition(i_start, i_start_loc, i_length, this, &
    j_neighb)
  type (Partition) :: this
  integer :: j_neighb, i_neighb, height_var
  integer :: i,  left, right
  integer :: nb_cells, j_neighb_corr
  integer :: mid, nb_resized
  integer :: neighb_left, neighb_right, width
  integer :: i_start, i_length, i_start_loc
  integer :: i_lat, height
  integer :: nb_cells_prev, nb_lat2
  logical :: is_top
  integer :: j1, j2, nb_parts

  double precision :: dlon
  i = this%i_pos

  i_lat = get_first_lat_Partition(this)
  dlon =get_dlon(i_lat, this%grid)
  nb_cells = nb_cells_from_dlon_sector(dlon)

  j1 = -2
  j2 = 1
  mid = (nb_cells + 1)/2
  height_var = this%height_var
  ! First, we get the current partition extremities (position on the previous latitude line)
  ! Compute the extremities
  ! Tells us which latitude to choose on a leftover band and the middle
  ! partition respectively
  is_top = .False.

  if (on_northern_hemisphere(this%zone)) then
    is_top = .True.
  end if

  call partition_extremities(left, right, width, height_var, this, &
    nb_cells, is_top)

  ! Then we must find the neighbour for these extremities 
  ! At the equator, the cells are in 1-to-1 correspondance, so no need to
  ! compute the neighbours
  nb_cells_prev = nb_cells_north_last(i, this%zone, nb_cells)
  if (.not. is_latitude_just_after_equator(2, this)) then
    left = left_north_cell(left, mid, this%zone) 
    right = right_north_cell(right, mid, nb_cells_prev, this%zone) 
  end if

  ! Finally, we need the extremities for the neighbours (only the leftmost and
  ! rightmost)
  ! Number resized on the previous partition
  i_neighb = i - 1
  nb_parts = nb_parts_on_band_sector(i_neighb)
  j_neighb_corr = translate_jpos(j_neighb, nb_parts)
  call neighbour_extremities(neighb_left, neighb_right,i_neighb, j_neighb_corr,&
    width, this, nb_cells, height_var)

  ! The length is then given by a min
  i_start = max(left, neighb_left)
  i_length = min(right, neighb_right) - i_start  + 1
  !if (this%id == j1 .and. j_neighb == j2) then
  !  write (*,*) "left right", left, neighb_left, right, neighb_right
  !  write (*,*) "ilength", i_length, "for", this%id, "neigh", j_neighb
  !end if
  i_start_loc = local_start_interface(left, neighb_left)

end subroutine

!-----------------------------------------------------------------------------
!> Find the interior interface length between a partition and its northern neighbour (for
!> the neighbour). 
!> Also returns the position of the first cell inside the current partition grid.
!> @param j_neighb : neighbour column. Can be on any sector, we correct it
!> anyway.
!-----------------------------------------------------------------------------
subroutine interface_interior_north_partition(i_start, i_length, this, j_neighb)
  type (Partition) :: this
  integer :: j_neighb, i_neighb, height_var
  integer :: i,  left, right
  integer :: nb_cells, j_neighb_corr
  integer :: mid
  integer :: neighb_left, neighb_right, width
  integer :: i_start, i_length, i_start_loc
  integer :: i_lat, height, nb_cells_next
  logical :: is_top
  integer :: j1, j2, nb_parts

  double precision :: dlon
  i = this%i_pos

  i_lat = get_first_lat_Partition(this)
  dlon = get_dlon(i_lat, this%grid)
  nb_cells = nb_cells_from_dlon_sector(dlon)

  mid = (nb_cells + 1)/2
  height_var = this%height_var
  ! First, we get the current partition extremities (position on the previous latitude line)
  ! Compute the extremities
  ! Tells us which latitude to choose on a leftover band and the middle
  ! partition respectively
  is_top = .False.

  if (on_northern_hemisphere(this%zone)) then
    is_top = .True.
  end if

  ! Get partition extremities
  call partition_extremities(left, right, width, height_var, this, &
    nb_cells, is_top)

  ! Then the neighbour extremities (only the leftmost and rightmost)
  ! Number resized on the previous partition
  i_neighb = i - 1
  nb_parts = nb_parts_on_band_sector(i_neighb)
  j_neighb_corr = translate_jpos(j_neighb, nb_parts)
  call neighbour_extremities(neighb_left, neighb_right,i_neighb, j_neighb_corr,&
    width, this, nb_cells, height_var)

  ! Finally, we take the north neighbour of these
  ! At the equator, the cells are in 1-to-1 correspondance, so no need to
  ! compute the neighbours
  nb_cells_next = nb_cells_north_last(i, this%zone, nb_cells)
  mid = (nb_cells_next + 1)/2
  if (.not. is_latitude_just_after_equator(2, this)) then
    neighb_left = left_south_cell(neighb_left, mid, this%zone) 
    neighb_right = right_south_cell(neighb_right, mid, nb_cells_next, this%zone) 
  end if

  ! The length is then given by a min
  i_start = max(left, neighb_left)
  i_length = min(right, neighb_right) - i_start  + 1
  i_start = cast_local_start_interface(left, neighb_left, i_lat, this)

end subroutine

!-----------------------------------------------------------------------------
!> Find the interior interface length between a partition and its southern 
!> neighbour
!> Also returns the position of the first cell inside the current partition grid.
!> @param j_neighb : neighbour column. Can be on any sector, we correct it
!> anyway.
!-----------------------------------------------------------------------------
subroutine interface_interior_south_partition(i_start, i_length, this, j_neighb)
  type (Partition) :: this
  integer :: j_neighb, i_neighb, height_var
  integer :: i,  left, right
  integer :: nb_cells, j_neighb_corr
  integer :: mid
  integer :: neighb_left, neighb_right, width
  integer :: i_start, i_length, i_start_loc
  integer :: i_lat, height, nb_cells_next
  logical :: is_top
  integer :: j1, j2, nb_parts

  double precision :: dlon
  i = this%i_pos

  i_lat = get_last_lat_Partition(this)
  dlon = get_dlon(i_lat, this%grid)
  nb_cells = nb_cells_from_dlon_sector(dlon)

  mid = (nb_cells + 1)/2
  height_var = this%height_var
  ! First, we get the current partition extremities (position on the previous latitude line)
  ! Compute the extremities
  ! Tells us which latitude to choose on a leftover band and the middle
  ! partition respectively
  is_top = .False.

  if (on_southern_hemisphere(this%zone)) then
    is_top = .True.
  end if

  ! Get partition extremities
  call partition_extremities(left, right, width, height_var, this, &
    nb_cells, is_top)

  ! Then the neighbour extremities (only the leftmost and rightmost)
  ! Number resized on the previous partition
  i_neighb = i + 1
  nb_parts = nb_parts_on_band_sector(i_neighb)
  j_neighb_corr = translate_jpos(j_neighb, nb_parts)
  call neighbour_extremities(neighb_left, neighb_right,i_neighb, j_neighb_corr,&
    width, this, nb_cells, height_var)

  ! Finally, we take the north neighbour of these
  ! At the equator, the cells are in 1-to-1 correspondance, so no need to
  ! compute the neighbours
  nb_cells_next = nb_cells_south_last(i, this%zone, nb_cells)
  mid = (nb_cells_next + 1)/2
  if (.not. is_latitude_just_before_equator(i_lat, this)) then
    neighb_left = left_north_cell(neighb_left, mid, this%zone) 
    neighb_right = right_north_cell(neighb_right, mid, nb_cells_next, this%zone) 
  end if

  !  if (this%id == 2) then
  !    write(*,*) "left right", left, right, neighb_left, neighb_right
  !    write(*,*) "nb cells", nb_cells, nb_cells_next
  !  end if
  ! The length is then given by a min
  i_start = max(left, neighb_left)
  i_length = min(right, neighb_right) - i_start  + 1
  i_start = cast_local_start_interface(left, neighb_left, i_lat, this)

end subroutine



!-----------------------------------------------------------------------------
!> Width of a square partition before the resize. We use the current partition
!> but it works also for the neighbours if we are not on the leftover band
!-----------------------------------------------------------------------------
function width_before_resize(this) result(width)
  type (Partition) :: this
  integer :: width, nb_bands, height_var
  ! Avoid ghost cells
  width = get_nb_lat_Partition(this)
  !rite (*,*) "nblat", this%grid%nb_lat
  ! The neighbour width is the same as the current (before the resize)
  !write (*,*) "height var", this%height_var

  if (has_modified_height(this%i_pos, this%height_var, this%zone)) then
    width = width - 1
  end if
end function

!-----------------------------------------------------------------------------
!> Number of partitions resized in width on the current band from height_var and
!> the latitude position
!> Partitions are resized in height if abs(height_var) >= 1, and resized in
!> width also if abs(height_var) >= 2
!> @param zone : zone indice
!-----------------------------------------------------------------------------
function nb_resized_width(i_pos, zone) result(nb_resized)
  integer, intent(in) :: i_pos, zone
  integer :: nb_resized 
  integer :: nb_square_bands, tmp
  integer :: nb_bands_resized 
  integer :: nb_bands1, nb_squares
  integer :: nb_tot

  ! As the number of resized bands is stored, no need to compute it again
  nb_bands_resized = get_nb_bands_resized_Configuration(zone)
  ! Half number of partitions on the last band
  nb_resized = -1
  if (nb_bands_resized > 0) then
    nb_bands1 = get_nb_bands_Configuration(1)
    nb_squares = get_nb_bands_square_Configuration(zone)
    ! We compute the distance to the last square band (leftover is ignored)
    if (on_northern_hemisphere(zone)) then
      if (i_pos < nb_squares+1) nb_resized = nb_squares - i_pos 
    else
      nb_tot = get_total_nb_bands_Configuration_sector()
      if (i_pos > nb_tot - nb_squares) nb_resized = i_pos - (nb_tot - nb_squares+1)
    end if

    ! Must be >=0 and <= nb_bands_resized-2
    if (nb_resized > -1 .and. nb_resized < nb_bands_resized -1) then
      ! The number decreases from last band, maximum is nb_bands_resized -1
      nb_resized = nb_bands_resized - 1 - nb_resized
    else
      nb_resized = 0
    end if
  end if
  if (nb_resized < 0) nb_resized = 0

end function

!-----------------------------------------------------------------------------
!> Find the interface length between a partition and its southern neighbour.
!> Also returns the first indice of the ghost cells inside the current 
!> partition grid and on the current latitude line .
!> Beware, the indice is only local (for the corresponding band)
!> Depends entirely of the analytical partitioning
!> Same as north, except that the discussion about the leftover band is on the
!> left/right neighbour cells, not the width.
!> @param j_neighb : neighbour column. Can be in any sector
!-----------------------------------------------------------------------------
subroutine interface_south_partition(i_start, i_start_loc, i_length, this, &
    j_neighb) 
  type (Partition) :: this
  integer :: j_neighb, i_neighb, nb_resized
  integer :: i, left, right, nb_cells
  integer :: mid, nb_parts_band
  integer :: neighb_left, neighb_right, width
  integer :: i_start, i_length, i_start_loc
  integer :: height_var, i_lat
  integer :: height, nb_cells_next
  double precision :: dlon
  logical :: is_top
  integer :: j1, j2
  integer :: nb_parts, j_neighb_corr

  i = this%i_pos

  j1 = -2
  j2 = 1
  ! Number of cells on the last band (for both hemispheres) without the ghost
  ! cells
  i_lat = get_last_lat_Partition(this) 
  dlon = get_dlon(i_lat, this%grid)
  nb_cells = nb_cells_from_dlon_sector(dlon)
  mid = (nb_cells + 1)/2

  ! First, we get the current partition extremities (position on the previous latitude line)
  i_neighb = i + 1

  ! Tells us which latitude to choose on a leftover band and the middle
  ! partition respectively
  is_top = .False.

  if (on_southern_hemisphere(this%zone)) then
    is_top = .True.
  end if

  ! Compute the extremities
  height_var = this%height_var
  call partition_extremities(left, right, width, height_var, this, &
    nb_cells, is_top)
  if (this%id == j1 .and. j_neighb == j2) then
    write (*,*) "extremities", left, right
  end if

  ! Then we must find the neighbour of the extremities, unchanged at the
  ! equator.
  if (.not. is_latitude_just_before_equator(i_lat, this)) then
    left = left_south_cell(left, mid, this%zone) 
    right = right_south_cell(right, mid, nb_cells, this%zone) 
  end if
  if (this%id == j1 .and. j_neighb == j2) write (*,*) "extremities neighb", left, right

  ! Finally, we compute the extremites for the neighbours
  nb_parts = nb_parts_on_band_sector(i_neighb)
  j_neighb_corr = translate_jpos(j_neighb, nb_parts)

  call neighbour_extremities(neighb_left, neighb_right,i_neighb, j_neighb_corr,&
    width, this, nb_cells, height_var)

  if (this%id == j1 .and. j_neighb == j2) then
    write (*,*) "neighb, corr", j_neighb, j_neighb_corr
    write (*,*) "neighb extremities", neighb_left, neighb_right
    write (*,*) "ijneighb width nbcells", i_neighb, j_neighb, width, nb_cells
  end if
  !write (*,*) "south neighb for extremities", left, right
  ! The length is then given by a difference
  i_start = max(left, neighb_left)
  i_length = min(right, neighb_right) - i_start + 1
  ! Local indice (inside the grid)
  i_start_loc = local_start_interface(left, neighb_left)
end subroutine

!-------------------------------------------------------------------------------
!> Compute the extremities for a neighbour partition. 
!> @param neighb_left : left extremity
!> @param neighb_right : left extremity
!> @param this : current partition
!> @param i_neighb : band on which the neighbour partition is
!> @param j_neighb : position on this band
!> @param nb_cells : number of cells on the current line (not the neighbour !)
!-------------------------------------------------------------------------------
subroutine neighbour_extremities(neighb_left, neighb_right, i_neighb, &
    j_neighb, width, this, nb_cells, height_var)
  type (Partition) :: this
  integer :: neighb_left, neighb_right
  integer, intent(in) :: i_neighb, j_neighb
  integer, intent(in) ::  nb_cells, height_var
  integer :: i, width, nb_cells_next
  logical :: is_top, is_north, cond
  integer :: hemisph, neighb_hemisph
  integer :: height, nb_parts_band
  integer :: nb_resized, nb_lat2
  integer :: next_zone, next_hemisph
  integer :: i_lat, neighb_zone

  i = this%i_pos
  nb_parts_band = nb_parts_on_band_sector(i_neighb)
  nb_lat2 = get_nb_lat2_Global_grid()
  !write (*,*) "nb resized neighb", nb_resized

  if (i_neighb > i) then
    nb_cells_next = nb_cells_south_last(i, this%zone, nb_cells)
  else
    nb_cells_next = nb_cells_north_last(i, this%zone, nb_cells)
  end if

  ! Find neighbour zone and hemisphere
  next_hemisph = -1
  hemisph = 1
  ! The zone on the southern hemisphere (on first sector)
  next_zone = 4
  is_north = .True.
  if (on_southern_hemisphere(this%zone)) then
    next_hemisph = 1
    next_zone = 1
    hemisph = -1
    is_north = .False.
  end if
  ! Default : both are on the same hemisphere

  ! Zone for the neighbour is by default on the the same hemisphere
  neighb_zone = this%zone
  neighb_hemisph = hemisph
  cond = (is_just_before_equator(i) .and. i_neighb > i)
  ! It changes only before the equator and for a south neighbour, or after
  ! and for a north neighbour
  cond = cond .or. ( is_just_after_equator(i) .and. i_neighb < i)
  if (cond) then
    neighb_zone = next_zone
    neighb_hemisph = next_hemisph
  end if


  ! We want the top of the triangular partition only if both are on northern
  ! hemisphere and for a south relation
  ! or both on southern and for a north relation
  cond = (is_north .and. (.not.  is_just_before_equator(i)) .and. i_neighb > i) 
  cond = cond .or. (.not. is_north .and. (.not.  is_just_after_equator(i))&
    .and. i_neighb < i) 
  is_top = .False.
  if (cond) is_top = .True.

  if (is_leftover_band(i_neighb)) then
    ! Recompute leftover height
    i_lat = leftover_neighbour_extrema_lat( i_neighb, neighb_hemisph)
    !if (this%id == 2) write (*,*) "ilat", i_lat
    call left_right_cell_partition_leftover(neighb_left, neighb_right, &
      j_neighb, nb_parts_band, nb_cells_next, i_lat)
  else
    nb_resized = nb_resized_width(i_neighb, neighb_zone)
    ! If we are on a leftover band, must retrieve the width of a square
    ! partition
    if (is_leftover_band(this%i_pos)) then
      width = width_square_partition(neighb_zone) 
    end if


    call left_right_cell(neighb_left, neighb_right, i_neighb, j_neighb, &
      nb_resized, neighb_zone,  height_var, nb_parts_band, nb_cells_next, width, is_top)
    !if (this%id == 5) write (*,*) "inside neighb extremities", neighb_left, neighb_right
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Compute the neighbour latitude extrema. Used for computing the leftover width
!> so we need the latitude the furthest from the equator
!> @param neighb_hemisph : hemisphere for the neighbour
!-------------------------------------------------------------------------------
function leftover_neighbour_extrema_lat(i_neighb, neighb_hemisph) result(i_lat)
  integer,intent(in) :: i_neighb
  integer,intent(in) :: neighb_hemisph
  integer :: i_lat, nb_lat2
  logical :: north_lat

  nb_lat2 = get_nb_lat2_Global_grid()
  ! Default is a south relation (we want the north lat of the neighbour)
  !  north_lat = .True.

  ! Neighbour is on northern hemisphere
  if (band_on_northern_hemisphere(i_neighb)) then
    ! South latitude is the equator
    i_lat = recompute_leftover_height(neighb_hemisph) 
    i_lat = nb_lat2 - i_lat + 1
    ! Neighbour is on southern hemisphere
  else
    ! North latitude is the equator
    i_lat = recompute_leftover_height(neighb_hemisph) 
    i_lat = nb_lat2 + i_lat
  end if

end function

!-------------------------------------------------------------------------------
!>  Compute the partition extremities. Also returns width and height var
!> @param left[inout] : left extremity (cell number)
!> @param right[inout] : right extremity (cell number)
!> @param width[inout] : width of the partition
!> @param height_var[inout] : height_var of the partition
!> @param is_top : select the top or bottom for middle partition
!-------------------------------------------------------------------------------
subroutine partition_extremities(left, right, width, height_var, this, &
    nb_cells, is_top)
  type (Partition) :: this
  integer :: left, right
  integer :: width, height_var
  integer, intent(in) :: nb_cells
  logical, intent(in) :: is_top
  integer :: i, j, nb_lat2
  integer :: i_lat, j_corr
  integer :: nb_parts_band, nb_resized

  i = this%i_pos
  ! Number of partitions on the band
  nb_parts_band = nb_parts_on_band_sector(i)

  j = this%j_pos
  ! Translate to the first sector
  j_corr = translate_jpos(j, nb_parts_band)

  ! There can be a leftover on the southern hemisphere
  if (is_leftover_band(i)) then
    ! We need the latitude of the extrema (near the pole) to compute the width
    nb_lat2 = get_nb_lat2_Global_grid()

    if (on_northern_hemisphere(this%zone)) then
      i_lat = nb_lat2 + get_nb_lat_Partition(this)
    else
      i_lat = nb_lat2 - get_nb_lat_Partition(this) + 1
    end if

    call left_right_cell_partition_leftover(left, right, j_corr, nb_parts_band,&
      nb_cells, i_lat)
    ! Compute the width of square partitions from the last band
    width = width_from_leftover(this%i_pos)
  else
    ! Retrieve the partition width (before an eventual resize)
    width = width_before_resize(this)
    height_var = this%height_var
    nb_resized = nb_resized_width(i, this%zone)

    call left_right_cell(left, right, i, j_corr, nb_resized, this%zone, &
      height_var, nb_parts_band, nb_cells, width, is_top)
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Small function for computing the local starting indice of a neighbour
!-------------------------------------------------------------------------------
function local_start_interface(left, neighb_left) result(i_start)
  integer :: i_start
  integer :: left, neighb_left
  ! And the indice is either the extrema, or the difference (as we want the
  ! position on the current band grid)
  if (left > neighb_left) then
    i_start = left
  else
    i_start = neighb_left - left + 1
  end if
end function

!-------------------------------------------------------------------------------
!> Cast global indice on the sector to local starting indice
!-------------------------------------------------------------------------------
function cast_local_start_interface(left, neighb_left, i, this) result(i_start)
  type (Partition) :: this
  integer, intent(in) :: i
  integer :: i_start, nb_ghosts
  integer :: left, neighb_left

  nb_ghosts = get_nb_ghosts_west(i, this%grid)
  ! And the indice is either the extrema, or the difference (as we want the
  ! position on the current band grid)
  if (left > neighb_left) then
    i_start = nb_ghosts
  else
    i_start = neighb_left - left + nb_ghosts
  end if
end function


!-----------------------------------------------------------------------------
!> Length of the northern ghost cells interface for partition (i,j). Works on
!> northern hemisphere.
!> Idea : the interface is of the same length as the first line with eventually
!> 1 cell at the right and 1 at the left
!> @param i : band number of the partition
!> @param j : position of the partition on the band (must be on the first
!sector)
!> @param width : number of cells
!> @param dlon : variation of longitude
!> @param this : current partition
!-----------------------------------------------------------------------------
function nb_ghost_cells_north(i, j, width, this, dlon) result (nb)
  type (Partition) :: this
  integer, intent(in) :: i, j, width
  integer :: nb, nb_cells, nb_parts_band
  integer :: left, right, mid
  double precision, intent(in) :: dlon
  integer :: nb_lat2, height

  nb = width
  ! If we are on the leftover band, it is more tricky. We need to get the
  ! extremities and their neighbours
  if (is_leftover_band(i)) then
    ! Extremities closest to the pole
    height = 0
    if (on_northern_hemisphere(this%zone)) height = get_nb_lat_Partition(this)

    nb_cells = nb_cells_from_dlon_sector(dlon)
    mid = (nb_cells + 1)/2
    nb_parts_band = nb_parts_on_band_sector(i)
    call left_right_cell_partition_leftover(left, right, j, nb_parts_band, &
      nb_cells, height)
    ! Neighbours for extremities
    left = left_north_cell_northern(left, mid) 
    right = right_north_cell_northern(right, mid, nb_cells - 2) 
    nb = right - left + 1

    ! Other wise, we can do it analytically
  else
    ! Only at the center we have 1 neighbour
    if (i == j) then
      nb = 1
    else
      nb = nb
    end if
  end if

end function

!-----------------------------------------------------------------------------
!> Length of the southern ghost cells interface for partition (i,j)
!> Same idea as north, but easier as we do not deal with the leftover band
!> @param i, j : partition position in the matrix
!> @param zone : partition zone
!-----------------------------------------------------------------------------
function nb_ghost_cells_south(i, j, width, zone) result (nb)
  integer, intent(in) :: i, j, width
  integer, intent(in) :: zone
  integer :: nb

  nb = width
  ! Only at the center we have 2 neighbours on each side
  if (middle_partition(i, j, zone)) then
    nb = nb + 2
    ! Other wise, only one neighbour to add
  else 
    nb = nb + 1
  end if
end function

!-----------------------------------------------------------------------------
!> Same as interface_north_length_interior_cells_leftover. We cannot factorize it 
!> though.
!-----------------------------------------------------------------------------
function interface_south_length_interior_cells_leftover(this, j_neighb) result (i_length)
  type (Partition) :: this
  integer :: i_length, j_neighb
  integer :: nb_cells, mid, i_start_loc
  integer :: left, right, nb_parts_band
  integer :: nb_cells_next
  double precision :: dlon
  integer :: i_start, i_end, neighb_left, neighb_right
  integer :: width_prev, last_lat
  integer :: i_lat, nb_lat, nb_lat2
  integer :: height_var, width
  logical :: is_top

  last_lat = get_last_lat_Partition(this)
  dlon = get_dlon(last_lat, this%grid)
  ! Number of cells on the last latitude
  nb_cells = nb_cells_from_dlon_sector(dlon)
  ! Mid on the next
  ! Find the number of cells on the south latitude line
  nb_cells_next = nb_cells_south_last(this%i_pos, this%zone, nb_cells)
  mid = (nb_cells_next + 1)/2

  ! Get the interface extremitieswith the neighbour
  call interface_south_partition(i_start, i_start_loc, i_length, this, j_neighb)
  i_end = i_start + i_length - 1

  !if (this%id == 0) then
  !  write (*,*) "interface", i_start, i_start_loc, i_length
  !end if

  ! Compute the neighbours for these extremities, except at the equator
  if (is_latitude_just_before_equator(last_lat, this)) then
    neighb_left = i_start
    neighb_right = i_end
  else
    neighb_left = left_north_cell(i_start, mid, this%zone) 
    ! Beware of the indices (mid, nbcells)
    neighb_right = right_north_cell(i_end, mid, nb_cells, this%zone) 
  end if

  ! Extremities for the current partition
  ! Number of cells on the last band (for botho hemispheres)
  nb_lat = get_last_lat_Partition(this) 
  dlon = get_dlon(nb_lat, this%grid)
  nb_cells = nb_cells_from_dlon_sector(dlon)

  ! Tells us which latitude to choose on the middle partition 
  is_top = .False.
  if (on_southern_hemisphere(this%zone)) then
    is_top = .True.
  end if
  ! Compute the extremities
  call partition_extremities(left, right, width, height_var, this, &
    nb_cells, is_top)

  !  if (this%id == 0) then
  !    write (*,*) "jneighb", j_neighb
  !    write (*,*) "nb_cells next", nb_cells_next
  !    write (*,*) "south interior leftover", left, right, neighb_left, neighb_right
  !  end if
  ! We do the difference
  i_length = min(right, neighb_right) - max(left, neighb_left) + 1

end function

!-----------------------------------------------------------------------------
!> Compute the width of a square partition (before resizing)
!-----------------------------------------------------------------------------
function width_square_partition(zone) result(width)
  integer, intent(in) :: zone
  integer :: width, nb_bands_square
  integer :: nb_bands, nb_lat2
  type (Configuration) :: config

  nb_bands_square = get_nb_bands_square_Configuration(zone)
  config = get_Configuration()
  nb_lat2 = get_nb_lat2_Global_grid()
  ! If leftover band, we reduce the number of latitudes
  if (nb_bands_square < get_nb_bands_Configuration(zone)) then
    nb_lat2 = nb_lat2 - 1
  end if
  ! Find the height before resizing (also the width)
  width = nb_lat2 / nb_bands_square
end function


!-----------------------------------------------------------------------------
!> Gives the neighbours position from the leftover band. On the northern
!> hemisphere, the neighbour is at the north, otherwise at the south. 
!> The indice correspond to the position on the previous band and are > 0
!> Glbal position is done later
!> @param hemisph : positive on the northern hemisphere, otherwise negative.
!-----------------------------------------------------------------------------
subroutine north_south_partition_from_leftover(this, j_start, j_end, hemisph)
  type(Partition) :: this
  integer, intent(in) :: hemisph
  integer :: i, j, mid, j_start, j_end
  integer :: nb_parts_neighb, nb_parts_left
  integer :: nb_lat, nb_cells_neighb, neighb_left
  integer :: left, right, neighb_right
  integer :: mid_left, mid_right, nb_cells
  integer :: width_neighb, width, neighb_lat
  integer :: nb_parts, nb_lon, first_lat
  integer :: j_corr, i_neighb

  i = this%i_pos
  j = this%j_pos

  ! Translate j_pos to the first sector
  nb_parts = nb_parts_on_band_sector(i)
  j_corr = translate_jpos(j, nb_parts)

  ! Find number of partitions on cur and neighbour band
  i_neighb = i - 1
  ! Southern hemisphere : the neighbour is at the south
  if (hemisph < 0) i_neighb = i + 1
  call leftover_nb_parts_cur_neighb(nb_parts_neighb, nb_parts_left, i_neighb, i)

  ! Last lat for neighbour, then we go toward the pole
  neighb_lat = leftover_extrema_lat(this, hemisph)
  if (hemisph > 0) then
    neighb_lat = neighb_lat - 1
  else
    neighb_lat = neighb_lat + 1
  end if

  ! Find width, mid, nb cells
  call neighbour_info_from_leftover(width_neighb, mid_left, mid_right, mid,&
    nb_cells, nb_cells_neighb, neighb_lat, this, hemisph)
  width = nb_cells / nb_parts_left

  ! Inside the partition (numbering starts from 1)
  left = (j_corr - 1)*width + 1
  neighb_left = left_north_cell_northern(left, mid)
  j_start = find_partition_from_cell_from_last(this, neighb_left, mid_left, &
    mid_right, width_neighb, nb_parts_neighb, nb_cells_neighb)

  ! Same for right
  if (hemisph > 0) then
    first_lat = get_first_lat_Partition(this)
  else
    ! Southern hemisphere : the neighbour is at the south
    first_lat = get_last_lat_Partition(this)
  end if

  nb_lon = get_nb_lon_Partition(first_lat, this)
  right = left + nb_lon - 1
  ! North cell works on both hemispheres.
  neighb_right = right_north_cell_northern(right, mid, nb_cells_neighb)

  j_end = find_partition_from_cell_from_last(this, neighb_right, mid_left, mid_right, &
    width_neighb, nb_parts_neighb, nb_cells_neighb)
  !write (*,*) j_start,j_end
  if (this%id == -5) then
    !    write (*,*) "width", width
    write (*,*) "left right", left, right
    write (*,*) "neighb left right", neighb_left, neighb_right
    write (*,*) "jstart jend", j_start, j_end
    !    write (*,*) "nb_cells_neigb", nb_cells_neighb
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Computes necessery information for finding north or south partition from
!> leftover
!-------------------------------------------------------------------------------
subroutine neighbour_info_from_leftover(width_neighb, mid_left, mid_right, &
    mid, nb_cells, nb_cells_neighb, neighb_lat, this, hemisph)
  type (Partition) :: this
  integer, intent(in) :: neighb_lat, hemisph
  integer, intent(out) :: width_neighb, mid_left, mid
  integer, intent(out) :: mid_right, nb_cells, nb_cells_neighb
  integer :: height

  ! Get some info
  height = partition_unmodified_height(this)
  ! Retrieve the normal height 
  width_neighb = height

  ! First compute the leftmost northern neighbour cell
  if (hemisph > 0) then
    ! Northern hemisphere : the neighbour is north, so we go towards the equator
    nb_cells = nb_cells_lat_sector(neighb_lat + 1)
  else
    nb_cells = nb_cells_lat_sector(neighb_lat - 1)
  end if
  mid = (nb_cells + 1)/2

  ! Then find the partition 
  ! Retrieve the total number of cells on the neighbour band
  nb_cells_neighb = nb_cells_lat_sector(neighb_lat)
  ! Middle cell as the center partition is a special case
  mid_left = (nb_cells_neighb + 1)/2
  ! Minus half the width (exactly the extremity)
  mid_left = mid_left - (height - 1)
  ! Right extremity
  mid_right = mid_left + 2*height - 2
end subroutine

!-------------------------------------------------------------------------------
!> Find the partition containing the cell (all cases).
!> The partition position is given in the latitude band.
!> Works only from a last band and not for a middle part !
!> @param mid_left : left extremity of the middle partition (neighbour band)
!> @param mid_right : right extremity of the middle partition (neighbour band)
!> @param width_neighb : width on the neighbour band
!> @param nb_parts_neighb : number of partitions on neighbour band
!> @param nb_cells : number of cells on neighbour latitude line
!> @param neighb : the cell for which we need to know what partition it's in
!-------------------------------------------------------------------------------
function find_partition_from_cell_from_last(this, neighb, mid_left, mid_right, &
    width_neighb, nb_parts_neighb, nb_cells) result(id_j)
  type (Partition) :: this
  integer :: mid_left, mid_right, neighb
  integer :: width_neighb, nb_parts_neighb
  integer :: id_j, nb_cells
  ! Left part
  if (neighb < mid_left + 1) then
    id_j = find_extr_partition_from_cell(this, neighb, width_neighb)
    ! Middle part 
  else if (neighb < mid_right + 1) then
    ! Mid of the neighbious line
    id_j = (nb_parts_neighb + 1)/2
    ! Right part is like the left, with a variable change
  else
    neighb = nb_cells - neighb + 1
    id_j = find_extr_partition_from_cell(this, neighb, width_neighb)
    ! Don't forget to switch back
    id_j = nb_parts_neighb - id_j + 1
  end if
end function

!-----------------------------------------------------------------------------
!> On northern hemisphere, computes the last latitude line for the north
!> neighbour of the leftover band.
!> On southern hemisphere, computes the first latitude line for the south
!> neighbour of the leftover band.
!-----------------------------------------------------------------------------
function leftover_extrema_lat(this, hemisph) result (neighb_lat)
  type (Partition) :: this
  integer, intent(in) :: hemisph
  integer :: neighb_lat, nb_lat2

  nb_lat2 = get_nb_lat2_Global_grid()
  neighb_lat = get_nb_lat_Partition(this)
  if (hemisph > 0) then
    neighb_lat = nb_lat2 - neighb_lat + 1
  else
    neighb_lat = nb_lat2 + neighb_lat! + 1
  end if
end function

!-------------------------------------------------------------------------------
!> Find the number of parts on the current leftover band and on the neighbour
!> band (which must containes square partitions.
!-------------------------------------------------------------------------------
subroutine leftover_nb_parts_cur_neighb(nb_neighb, nb_cur, i_neighb, i_cur)
  integer, intent(in) :: i_neighb, i_cur
  integer, intent(out) :: nb_neighb, nb_cur

  nb_neighb = nb_square_parts_on_band(i_neighb)
  nb_cur = nb_parts_on_band_sector(i_cur)
end subroutine

!-----------------------------------------------------------------------------
!> Find the last band square partition containing the cell. This works for a 
!> partition on the other hemisphere.
!> This is an alternative to find_partition_from_cell_from_last as we do not have
!> information on the partition.
!> @param j_cell : cell number
!> @param i : band indice
!> @param width : width of square partition before resize
!-----------------------------------------------------------------------------
function find_last_band_square_partition_from_cell(j_cell, i, width, zone) &
    result(j_neighb)
  type (Partition) :: this
  integer, intent(in) :: width, i, zone
  integer :: j_cell, j2
  integer :: nb_corr, j_neighb
  integer :: height, resized_right
  integer :: mid_right, mid_left
  integer :: nb_lat2, nb_resized
  integer :: nb_cells, nb_parts
  integer :: i_lat, nb_bands_resized
  integer :: a1, a2

  nb_bands_resized = get_nb_bands_resized_Configuration(zone)
  ! Half number of partitions on the last band
  nb_resized = 0
  if (nb_bands_resized > 0) nb_resized = nb_bands_resized - 1
  ! Left extremity of resized partitions
  resized_right = nb_resized*(width + 1)
  height = width
  nb_parts  = nb_parts_on_band_sector(i)
  ! Resizing is only of 1
  if (nb_resized > 0) height = width + 1

  nb_lat2 = get_nb_lat2_Global_grid()
  ! First partition lat
  if ( on_northern_hemisphere(zone)) then
    i_lat = nb_lat2 - height + 1
  else
    i_lat = nb_lat2 + height
  end if
  ! Number of cells near the pole (serves as a left frontier for the middle
  ! partition)
  nb_cells = nb_cells_lat_sector(i_lat)
  ! Middle partitions extremities
  ! Left does not change when the latitude line varies
  mid_left = (nb_cells + 1)/2
  ! We have 2*height - 1 cells at the middle
  mid_right = mid_left + 2*height - 2

  !write (*,*) "ilat nb cells", i_lat, nb_cells
  a1 = -3
  a2 = -45
  !if (j_cell == 1 .or. j_cell == 45) write (*,*) "midleft midright", mid_left, mid_right

  j2 = j_cell
  ! If on the right side, change the variable
  nb_cells = nb_cells_lat_sector(nb_lat2)
  if (j_cell > mid_right) j2 = nb_cells - j_cell + 1

  ! Left size, inside resized partitions
  if (j2 < resized_right + 1) then
    j_neighb = (j2 - 1)/ (width + 1) + 1
    ! Left size, inside normal partitions
  else if (j2 < mid_left) then
    j_neighb = (j2 - 1 - resized_right) / width + nb_resized + 1
    ! Middle partition
  else if (j2 < mid_right + 1) then
    j_neighb  = (nb_parts + 1)/2
  end if
  if (i == a1 .and. j_cell == a2) then
    write (*,*) "j_cell", j_cell, "j2", j2
    !write (*,*) "j_neighb", j_neighb
    write (*,*) "mid_left", mid_left, "ilat", i_lat
    !write (*,*) "resize righ", resized_right
    !write (*,*) "nb resized", nb_resized
  end if

  ! If on the right side, adjust partition number
  if (j_cell > mid_right) j_neighb = nb_parts - j_neighb + 1
end function

!-----------------------------------------------------------------------------
!> Find the partition containing the cell (works for only for extremities). 
!> Conditions : 
!> * the neighbour is in the same hemisphere.
!> * works only from a last band and not for a middle part !
!> * works only for square partitions.
!-----------------------------------------------------------------------------
function find_extr_partition_from_cell(this, j_cell, width) result(j_neighb)
  type (Partition) :: this
  integer :: j_cell, j2
  integer :: nb_corr, j_neighb
  integer :: width, width_corr
  width_corr = width + 1

  ! Compute the number of corrected partition (with larger width) on a side
  ! Retrive it for the modification of last band height, as it is done only for
  ! a modification larger than 1
  ! Then number of partitions modified on last line (on one half)
  ! On last band, heightvar is -nb_corr - 2 (as partition width is modified only
  ! if height_var > 1)
  nb_corr = 0
  if (this%height_var < -2) nb_corr =  - this%height_var - 2

  ! If the neighbour is actually inside the stack of modified 
  if (j_cell < nb_corr * width_corr + 1) then
    ! Numbering starts from 1 so we must substract 1 to j_cell
    j_neighb = int((j_cell-1)/ width_corr) + 1
    ! Otherwise must it is in normal width
  else
    j2 = j_cell - nb_corr * width_corr
    j_neighb = int((j2-1)/ width) + 1
    ! Don't forget to add the offset
    j_neighb = j_neighb + nb_corr
  end if

end function 

!-----------------------------------------------------------------------------
!> Search the height of an unmodified square partition for north neighbours on last band
!-----------------------------------------------------------------------------
function partition_unmodified_height(this) result(height)
  type (Partition) :: this
  integer :: height
  integer :: nb_bands, nb_lat, nb_lat_square

  nb_bands = get_nb_bands_square_Configuration(this%zone)
  nb_lat = get_nb_lat2_Global_grid()

  nb_lat_square = get_nb_lat_Partition(this)
  nb_lat_square = nb_lat - nb_lat_square
  height = nb_lat_square / nb_bands
end function

!-----------------------------------------------------------------------------
!> Gives the north or south neighbours position >= 1 (-1 if none) for a band 
!> on one side of the equator, the neighbouring band being on the other side.
!> Works also for leftover partitions.
!> @param i_neighb : band containing the neighbour partition
!-----------------------------------------------------------------------------
subroutine north_south_partition_equator(this, j_start, j_end, i_neighb)
  type(Partition), target :: this
  integer, intent(in) :: i_neighb
  integer :: j_start, j_end
  integer :: left, right
  integer :: nb_lat2, nb_cells, nb_parts
  integer :: width, width_prev
  integer :: zone, nb_lat
  integer :: i_lat
  integer :: height
  logical :: is_top

  nb_lat2 = get_nb_lat2_Global_grid()
  nb_parts = nb_parts_on_band_sector(this%i_pos)
  nb_cells = nb_cells_lat_sector(nb_lat2)

  ! First, we get the current partition extremities (position on the previous
  ! latitude line)
  ! We must also retrieve the partition width (before an eventual resize)

  ! Compute the extremities
  ! Tells us which latitude to choose on the middle partition respectively
  ! Always bottom here
  is_top = .False.

  ! Compute the partition extremities
  call partition_extremities(left, right, width, this%height_var, this,&
    nb_cells, is_top)

  ! The neighbours are the same on each side of the equator
  ! So we only need to find which partition it's in
  if (is_leftover_band(i_neighb)) then
    ! Find width of the other band
    if (on_northern_hemisphere(this%zone)) then
      width = compute_width_leftover(i_neighb, -1)
    else
      width = compute_width_leftover(i_neighb, 1)
    end if
    nb_parts = nb_parts_on_band_sector(i_neighb)
    ! On the other side of the equator, the position on the latitude line is the
    ! same
    ! Find to what partitions they belong
    !if (this%id == 5) write (*,*) "left right", left, right
    j_start = find_partition_from_cell_leftover(left, nb_parts, width)
    j_end = find_partition_from_cell_leftover(right, nb_parts, width)
    !if (this%id == 5) write (*,*) "jend", j_end
  else
    ! Find width of the other band
    zone = 4
    if (.not. on_northern_hemisphere(this%zone)) zone = 1
    width = width_square_partition(zone)
    !if (this%id == 2) write (*,*) "width", width

    j_start = find_last_band_square_partition_from_cell(left, i_neighb, width, zone)
    j_end = find_last_band_square_partition_from_cell(right, i_neighb, width, zone)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Re-compute analytically the width of a partition on a leftover band
!> @param i_neighb : band indice
!> @param hemisph : hemisphere
!-----------------------------------------------------------------------------
function compute_width_leftover(i_neighb, hemisph) result (width)
  integer, intent(in) :: i_neighb
  integer :: width, nb_lat2, height
  integer :: nb_cells, nb_parts, i_lat
  integer :: nb_bands_resized
  integer :: hemisph, nb_bands_square
  integer :: side

  nb_lat2 = get_nb_lat2_Global_grid()
  height = recompute_leftover_height(hemisph)
  ! Find latitude indice
  if (hemisph > 0) then
    i_lat = nb_lat2 - height + 1
  else
    i_lat = nb_lat2 + height
  end if
  nb_cells = nb_cells_lat_sector(i_lat)
  nb_parts = nb_parts_on_band_sector(i_neighb)
  width = nb_cells / nb_parts
  !if (i_neighb == 3) then
  !      write (*,*) "!! nb_cells there", nb_cells
  !      write (*,*) "!! height", height, "ilat", i_lat
  !      write (*,*) "!! width", width
  !    end if

end function

!-----------------------------------------------------------------------------
!> Recomputes analytically the height of the leftover
!-----------------------------------------------------------------------------
function recompute_leftover_height(hemisph) result (height)
  integer, intent(in) :: hemisph
  integer :: height, nb_lat2
  integer :: side, nb_bands_square
  integer :: nb_bands_resized

  nb_lat2 = get_nb_lat2_Global_grid()
  if (hemisph > 0) then
    nb_bands_square = get_nb_bands_square_Configuration(1)
    ! Height of the leftover
    side = (nb_lat2 - 1) /nb_bands_square
    height = nb_lat2 - nb_bands_square*side

    nb_bands_resized = get_nb_bands_resized_Configuration(1)
    if (nb_bands_resized > 0) height = height - nb_bands_resized
  else
    nb_bands_square = get_nb_bands_square_Configuration(4)
    side = (nb_lat2 - 1) /nb_bands_square
    height = nb_lat2 - nb_bands_square*side

    nb_bands_resized = get_nb_bands_resized_Configuration(4)
    if (nb_bands_resized > 0) height = height - nb_bands_resized
  end if
end function

!-----------------------------------------------------------------------------
!> Gives the south neighbours position >= 1 (-1 if none) for the last band
!> containing square partitions on northern hemisphere. 
!> Works also for finding north partitions for last square band on southern
!> hemisphere.
!-----------------------------------------------------------------------------
subroutine north_south_partition_last_square(this, j_start, j_end, hemisph)
  type(Partition), target :: this
  integer, intent(in) :: hemisph
  integer :: i, j, j_start, j_end
  integer :: nb_parts, nb_parts_left, nb_partitions
  integer :: nb_cells, nb_cells_next, width
  double precision :: dlon
  integer :: left, right, i_lat
  integer :: j1, j2, mid, width_prev

  i = this%i_pos
  nb_parts = nb_parts_on_band_sector(i)
  j = translate_jpos(this%j_pos, nb_parts)


  ! Find number of cells on the next latitude line (use float operations)
  if (hemisph > 0) then
    i_lat = get_last_lat_Partition(this)
  else
    i_lat = get_first_lat_Partition(this)
  end if
  dlon = get_dlon(i_lat, this%grid) 
  nb_cells = nb_cells_from_dlon_sector(dlon)

  ! Either way, we goes toward the equator
  nb_cells_next = nb_cells  + 2
  width_prev = get_nb_lat_Partition(this)
  !write (*,*) "height var", this%height_var
  ! If the partition has been extended, it is only by 1
  if (this%height_var /= 0) width_prev = width_prev - 1

  ! Width of last band partitions
  if (hemisph > 0) then
    nb_parts_left = nb_parts_on_band_sector(i + 1)
  else
    nb_parts_left = nb_parts_on_band_sector(i - 1)
  end if
  width = nb_cells_next / nb_parts_left

  call left_right_cell_partition(left, right, this, nb_parts, nb_cells,&
    width_prev, .False.)

  ! Find the extrema neighbours cells
  mid = (nb_cells + 1)/2
  j1 = left_south_cell_northern(left, mid)
  j2 = right_south_cell_northern(right, mid)

  ! Find to what partitions they belong
  j_start = find_partition_from_cell_leftover(j1, nb_parts_left, width)
  j_end = find_partition_from_cell_leftover(j2, nb_parts_left, width)
end subroutine

!-----------------------------------------------------------------------------
!> Wrapper for left_right_cell
!-----------------------------------------------------------------------------
subroutine left_right_cell_partition(left, right, this, nb_parts, &
    nb_cells, width, is_top)
  type (Partition) :: this
  integer :: left, right, nb_resized
  integer :: width, nb_parts, nb_cells
  logical :: is_top

  nb_resized = nb_resized_width(this%i_pos, this%zone)
  !if (this%id == 107)  write (*,*) "nbresize", nb_resized
  !write (*,*) "id", this%id,"ij", this%i_pos, this%j_pos
  call left_right_cell(left, right, this%i_pos, this%j_pos, nb_resized, &
    this%zone, this%height_var, nb_parts, nb_cells, width, is_top)
end subroutine

!-----------------------------------------------------------------------------
!> Finds analytically the position of the left and rightmost cell on the band 
!> (before the leftover) for the current partition. For rectangular partitions,
!> left and right are the same for north and south. At the middle, we use a
!> logical to tell us which to choose.. Position is given on the
!> latitude line. Done analytically as it avoids to use floats. 
!> Does not need a partition (only 3 integer).
!> We need the width on the current band
!> Does not take ghost cells into account.
!> @param width : width of a square partition
!> @param is_top : true if we want the top of a triangular partition
!> @param nb_parts : nb of partitions on the current band (only in one sector)
!-----------------------------------------------------------------------------
subroutine left_right_cell(left, right, i, j, nb_resized, zone, height_var, nb_parts, &
    nb_cells, width, is_top)
  integer, intent(in) :: i, j, width, nb_resized, zone
  integer :: left, right, right_tmp
  integer :: j2, j_mid, j_corr
  integer :: height_var
  integer :: nb_parts, nb_cells
  ! Different situation at the middle if we consider the top or bottom
  logical :: is_top
  integer :: a1, a2

  a1 = -2
  a2 = 3

  if (nb_parts < 0) then
    call print_error("nb parts < 0", "left_right_cell", fname_partition)
  end if

  ! Translate j_pos to the first sector
  j_corr = translate_jpos(j, nb_parts)
  j2 = j_corr

  ! Middle indice
  j_mid = (nb_parts + 1)/2
  if (j_corr > j_mid) j2 = nb_parts - j_corr  + 1

  ! If inside the resized partitions (stacked at the left)
  if (j2 < nb_resized + 1) then
    left = (j2 - 1)*(width + 1) + 1
    ! Here the width is width + 1
    right = left + width
  else
    left = nb_resized*(width + 1) + (j2 - nb_resized - 1)*width + 1
    right = left + width - 1
  end if

  ! Special case at the middle. Do not forget to substract one
  ! Also, the height may have been modified (add 1 in this case)
  if (j_corr == j_mid) then
    ! At the top, only one cell
    if (is_top) then
      right = left
    else
      ! At the bottom, it can be resized
      ! Width is 2*width - 1 at the center (triangular structure)
      right = left + 2*width - 1 - 1
      if (has_modified_height(i, height_var, zone)) then
        ! If there is one more latitude band, there are 2 more cells at the
        ! middle
        right = left + 2*width 
      end if
    end if
  end if

  ! If on the right part
  if (j_corr > j_mid) then
    ! We must switch left and right on the other half
    right_tmp = right
    right = nb_cells - left  + 1
    left = nb_cells - right_tmp + 1
  end if
  if (i == a1 .and. j == a2) then
    write (*,*) "left right at the end", left, right, "nbcells", nb_cells
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Returns true if the partitions has been resized in height. Works only for
!> "squared" partitions
!> @param i_pos : current band number
!> @param height_var : variation of partition height
!> @param zone : zone indice
!-----------------------------------------------------------------------------
function has_modified_height(i_pos, height_var, zone) result(res)
  integer, intent(in) :: i_pos, height_var, zone
  integer :: nb_square_bands, i_pos2
  integer :: nb_bands_resized
  integer :: nb_bands1, nb_bands_tot
  integer :: start
  logical :: res

  res = .False.
  nb_square_bands = get_nb_bands_square_Configuration(zone)
  ! Number of square bands resized actually
  nb_bands_resized = get_nb_bands_resized_Configuration(zone)
  nb_bands1 = get_nb_bands_Configuration(1)
  if (nb_bands_resized < 1) return

  ! North hemisphere
  if (i_pos < nb_bands1 + 1) then
    res = (i_pos > (nb_square_bands - nb_bands_resized))
  else
    ! South hemisphere
    nb_bands_tot = get_total_nb_bands_Configuration_sector()
    start = (nb_bands_tot - nb_square_bands)
    res = (i_pos > start .and. i_pos < start + nb_bands_resized + 1)
  end if
  !    if (i_pos == 11) then
  !      write (*,*) "has modified height", i_pos, ":", res
  !      write (*,*) "nb bands1, nb square, nbresized", nb_bands1, nb_square_bands, nb_bands_resized
  !      write (*,*) "zone", zone
  !    end if
end function

!-----------------------------------------------------------------------------
!> Find in which partition the cell number j_cell is in.
!> Works only for a cell on the leftover partition band
!-----------------------------------------------------------------------------
function find_partition_from_cell_leftover(j_cell, nb_parts_left, width) result (part)
  integer :: part
  integer, intent(in) :: j_cell, nb_parts_left, width

  if (j_cell > (nb_parts_left - 1)*width) then
    part = nb_parts_left
  else
    ! We want part such as
    ! (part - 1)*w + 1 <= j_cell < part*w + 1
    ! <=> part <= (j_cell - 1)/w + 1 < part + 1
    part = int(dble(j_cell-1)/width) + 1
  end if
end function

!-----------------------------------------------------------------------------
!> Finds analytically the position of the left and rightmost cell on the last 
!> band. Equivalent of left_right_cell_partition
!> @param nb_cells : number of cells at the neighbour latitude (or last on south
!> hemisphere)
!> @param i_lat : latitude at the extrema
!-----------------------------------------------------------------------------
subroutine left_right_cell_partition_leftover(left, right, j, nb_parts, &
    nb_cells, i_lat)
  integer :: left, right, width, width_prev
  integer, intent(in) :: nb_parts, nb_cells, i_lat 
  integer :: nb_cells_extr
  integer :: j, nb_lat2

  nb_cells_extr = nb_cells_lat_sector(i_lat)
  call width_leftover(width_prev, width, j, nb_cells_extr, nb_parts)
  ! Last latitude is always on the equator

  !write (*,*) "leftover j", j
  left = (j-1)*width_prev + 1

  if (j < nb_parts) then
    right = left + width - 1
    ! Last partition has a different width (the rest)
  else
    right = nb_cells
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Numerical way to retrieve the number of cells at latitude line number i for a
!> partition (float operations).
!> Number of cells on a sector
!-----------------------------------------------------------------------------
function nb_cells_from_dlon_sector(dlon) result(nb_cells)
  double precision, intent(in) :: dlon
  integer :: nb_cells

  nb_cells = int(120.d0/dlon)
  ! Check
  if (nb_cells*dlon < 120.d0 - DBLE_PREC) then
    nb_cells = nb_cells + 1
  else if (nb_cells*dlon > 120.d0 + DBLE_PREC) then
    nb_cells = nb_cells - 1
  end if
end function

!-----------------------------------------------------------------------------
!> Tells us if the partition is on the right boundary on a group at the leftover
!> band
!-----------------------------------------------------------------------------
function is_right_boundary(j,r,q) result(res)
  integer, intent(in) :: j,r,q
  integer :: j2
  logical :: res

  res = .False.
  if (j <= r*(q+1)) then
    res = (mod(j,q+1) == 0)
  else
    j2 = j - r*(q+1)
    res = (mod(j2,q) == 0)
  endif
end function

!-----------------------------------------------------------------------------
!> Tells us if the partition is on the left boundary on a group at the leftover
!> band
!-----------------------------------------------------------------------------
function is_left_boundary(j,r,q) result(res)
  integer, intent(in) :: j,r,q
  integer :: j2
  logical :: res

  res = .False.
  if (j <= r*(q+1)) then
    res = (mod(j,q+1) == 1)
  else
    j2 = j - r*(q+1)
    if (q > 1) then 
      res = (mod(j2,q) == 1)
    else
      res = .True.
    end if
  endif
end function


!-----------------------------------------------------------------------------
! Returns the number of partitions on last band for each side of the central
! partition
!-----------------------------------------------------------------------------
subroutine get_coef_leftover_last_band(coef,rest,is_even,nb_squares,n_left)
  integer, intent(in) :: nb_squares,n_left
  logical, intent(in) :: is_even
  integer :: coef, rest

  if (is_even) then
    coef = (nb_squares-1)/(n_left-2)
    rest = mod(nb_squares-1,n_left-2)
  else
    coef = (nb_squares-1)/(n_left-1)
    rest = mod(nb_squares-1,n_left-1)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Gives the south neighbours ids >= 1 (-1 if none). 
!> This is based on an analytical partitioning.
!-----------------------------------------------------------------------------
subroutine south_neighbours(this, id_start, id_end)
  type(Partition) :: this
  integer :: id_start, id_end
  integer :: i_south, j_start, j_end
  integer :: tmp

  id_start = -1
  id_end = -1
  call south_partition_wrapper(this , i_south, j_start, j_end)
  if (i_south > -1 .and. j_start > -1) then
    id_start = compute_id_from_pos(i_south, j_start)
  end if
  if (i_south > -1 .and. j_end > -1) then
    id_end = compute_id_from_pos(i_south, j_end)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Wrapper for south_partition_sector. Once we computed the neighbour, we correct
!> the position as the analytical partitioning only gives the neighbours on 
!> the first sector
!-----------------------------------------------------------------------------
subroutine south_partition_wrapper(this, i_south, j_start, j_end)
  type(Partition) :: this
  integer :: i_south, j_start, j_end
  integer :: sector, nb_parts

  call south_partition_sector(this, i_south, j_start, j_end)
  if (j_start > 0 .and. j_end > 0) then
    sector = sector_from_zone(this%zone)
    nb_parts = nb_parts_on_band_sector(i_south)
    j_start = j_start + (sector - 1)*nb_parts
    j_end = j_end + (sector - 1)*nb_parts
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Gives the south neighbours indices (> 0) in the matrix position
!> This is based on an analytical partitioning.
!> Returns the neighbours on the first sector, so you have to correct the
!> results
!-----------------------------------------------------------------------------
subroutine south_partition_sector(this, i_south, j_start, j_end)
  type(Partition) :: this
  integer :: i_south, j_start, j_end
  integer :: i, j, j_corr
  integer :: nb_parts, nb_parts_next

  i = this%i_pos
  !write (*,*) "starting south"
  i_south = -1
  j_start = -1
  j_end = -1

  ! No neighbours after the south pole
  if (is_southmost(this)) return

  i_south = i + 1
  ! If the next band is on the other side of the equator
  if (is_just_before_equator(this%i_pos)) then
    call north_south_partition_equator(this, j_start, j_end, i+1)
    !if (this%id == 2) write (*,*) "res", j_start, j_end
    return
  end if

  ! If the next band is the leftover band on northern hemisphere
  if (on_northern_hemisphere(this%zone) .and. is_leftover_band(i+1)) then
    !  if (this%id == 5) write (*,*) "blouf", "this i", this%i_pos, "nex", i+1
    call north_south_partition_last_square(this, j_start, j_end, 1)
    return
  end if

  ! If the current is leftover on the southern hemisphere
  if (on_southern_hemisphere(this%zone) .and. is_leftover_band(i)) then
    call north_south_partition_from_leftover(this, j_start, j_end, -1)
    !if (this%id == 43) write (*,*) "kikoo", j_start, j_end
    return
  end if

  !if (this%id == 4) write (*,*) "otherwise"
  ! Otherwise, square partitions
  ! Number of partitions on the band 
  nb_parts = nb_square_parts_on_band(i)
  nb_parts_next = nb_square_parts_on_band(i+1)
  j = this%j_pos
  j_corr = translate_jpos(j, nb_parts)

  ! Southern hemisphere : south relation becomes north (actually, it is the
  ! same)
  if (on_northern_hemisphere(this%zone)) then
    call south_partition_reduced(j_corr, nb_parts, nb_parts_next, j_start, j_end)
  else
    call north_partition_reduced(j_corr, nb_parts, nb_parts_next, j_start, j_end)
  end if

  call correct_neighbours_resized(j_start, j_end, this, 2, j_corr)
end subroutine

!-----------------------------------------------------------------------------
!> Revert north and south for overlap
!-----------------------------------------------------------------------------
subroutine revert_overlap(this)
  type (Partition) :: this
  integer :: o1, o2

  ! Reverting means the number are exchanged
  !if (this%id == 79) write (*,*) "overlap before", this%overlap
  o2 = this%overlap/10
  o1 = this%overlap - o2*10

  this%overlap = 10*o1 + o2
  !  if (this%id == 70) write (*,*) "overlap after", this%overlap
end subroutine

!-----------------------------------------------------------------------------
!> Check if the next partition is on the leftover band in our partitioning. Different
!> from the last band.
!> @param i_pos : band indice
!> @param zone : zone indice
!-----------------------------------------------------------------------------
function next_is_leftover_band(i_pos, zone) result(res)
  integer, intent(in) :: i_pos, zone
  logical :: res

  res = is_leftover_band(i_pos + 1)
end function 

!-----------------------------------------------------------------------------
!> Computes the number of neighbours for the current partition
!-----------------------------------------------------------------------------
function nb_neighbours_Partition(this) result(nb)
  type(Partition) :: this
  integer :: nb
  integer :: cur, cur2
  nb = 0

  nb = nb + nb_north_neighbours(this)
  nb = nb + nb_east_neighbours(this)
  nb = nb + nb_south_neighbours(this)
  nb = nb + nb_west_neighbours(this)
end function

!-----------------------------------------------------------------------------
!> Computes the number of north neighbours
!-----------------------------------------------------------------------------
function nb_north_neighbours(this) result(nb)
  type(Partition) :: this
  integer :: nb
  integer :: cur, cur2
  nb = 0

  if (is_northmost(this)) return

  call north_neighbours(this, cur, cur2)
  !if (this%i_pos == 4) write (*,*) "north neigh", cur, cur2
  if (cur > -1 .and. cur2 > -1) then
    nb = nb + cur2 - cur + 1
  else if (cur > -1 .or. cur2 > -1) then
    nb = nb + 1
  end if
end function

!-----------------------------------------------------------------------------
!> Computes the number of south neighbours
!-----------------------------------------------------------------------------
function nb_east_neighbours(this) result(nb)
  type(Partition) :: this
  integer :: nb
  integer :: cur
  nb = 0

  call east_neighbours(this, cur)
  if (cur > -1) nb = nb + 1
end function

!-----------------------------------------------------------------------------
!> Computes the number of south neighbours
!-----------------------------------------------------------------------------
function nb_west_neighbours(this) result(nb)
  type(Partition) :: this
  integer :: nb
  integer :: cur
  nb = 0

  call west_neighbours(this, cur)
  if (cur > -1) nb = nb + 1
end function

!-----------------------------------------------------------------------------
!> Computes the number of south neighbours
!-----------------------------------------------------------------------------
function nb_south_neighbours(this) result(nb)
  type(Partition) :: this
  integer :: nb
  integer :: cur, cur2
  nb = 0

  if (is_southmost(this)) return

  call south_neighbours(this, cur, cur2)
  if (cur > -1 .and. cur2 > -1) then
    nb = nb + cur2 - cur + 1
  else if (cur > -1 .or. cur2 > -1) then
    nb = nb + 1
  end if
end function


!-----------------------------------------------------------------------------
!> This "slope" tells us if the first cell on each line goes toward the west      
!> (posivite) or the west (negative). Used for computing cells neighbours in
!> tests.
!> We cannot use the first longitude due to float approximation.
!> Slope is always 1, except on last partition on the leftover band.
!-----------------------------------------------------------------------------
function pseudo_slope_Partition(this) result(slope)
  type (Partition) :: this
  integer :: mid, slope
  integer :: i_pos, j_pos
  integer :: nb_parts

  i_pos = this%i_pos
  nb_parts = nb_parts_on_band_sector(i_pos)
  j_pos = translate_jpos(this%j_pos, nb_parts)

  ! Slope is positive on the last band (always) or on the left side of the
  ! "squared" partitions
  mid = (nb_parts_on_band_sector(i_pos) + 1)/2
  !write (*,*) "pseudo slope for", this%id, this%i_pos
  if (is_leftover_band(i_pos)) then
    slope = 1
    ! The slope is larger on the last partition, if it is to the right of the
    ! middle, so if there are at least 3 partitions on the last band
    if (is_last_on_band(this) .and. j_pos > 2) slope = 2
  else if (j_pos < mid + 1) then
    slope = 1
  else
    slope = -1
  end if
  if (j_pos == 1) slope = 0
end function

!-----------------------------------------------------------------------------
!> Cell Neighbour on the other side of the equator
!-----------------------------------------------------------------------------
subroutine neighbour_equator(j_start, j_end, i, j2, i_neighb, sector, single, this)
  type (Partition) :: this
  integer :: j_start, j_end
  integer, intent(in) :: i, i_neighb, j2, sector
  logical, intent(in) :: single
  integer :: nb_ghosts

  ! No west/east ghost cells on the north ghost cells
  ! We correct by substracting (to ghost, the default) or adding (from ghost) 
  nb_ghosts = this%grid%nb_ghosts_west(i)
  if (i == 1 .or. i == this%grid%nb_lat) then
    nb_ghosts = -this%grid%nb_ghosts_west(i_neighb)
  end if

  j_start = j2 - nb_ghosts
  j_end = j2 - nb_ghosts
  ! Do not take ghost cells into account
  if (single .or. has_three_partitions()) then
    j_start = j2 
    j_end = j2
  end if

  if (single) call correct_cell(j_start, sector, i, this, 1)
  if (single) call correct_cell(j_end, sector, i, this, 1)
end subroutine

!-----------------------------------------------------------------------------
!> Returns the north neighbours for cells (i,j) inside the partition. Returns only
!> the indice on the previous line. Also works (and include) the ghost cells.
!> The caller must check there is a neighbour.
!> @param i : latitude line inside the partition
!> @param j : cell position on the latitude line inside the partition
!-----------------------------------------------------------------------------
subroutine north_neighbour_cell_Partition(j_start, j_end, i, j, this, &
    i_neighb)
  type (Partition) :: this
  integer, intent(in) :: i, j
  integer, intent(in), optional ::  i_neighb
  integer :: j_start, j_end
  integer :: i_neighb2, offset
  integer :: mid, slope
  integer :: nb_parts, j_pos
  integer :: sector, j2
  logical :: single

  j_start = -1
  j_end = -1
  if (i == 1) return

  ! For a single partition, we must go back to one sector
  j2 = j
  single = has_single_partition()
  if (single) call correct_cell(j2, sector, i, this, 0)

  ! At the equator, 1-to-1 correspondance
  if (is_latitude_just_after_equator(i, this)) then
    call neighbour_equator(j_start, j_end, i, j2, i-1, sector, single, this)
    return
  end if

  nb_parts = nb_parts_on_band(this%i_pos)                                                  
  j_pos = translate_jpos(this%j_pos, nb_parts) 
  offset = this%grid%north_offset(i)

  ! After the equator (don't forget the single partition case) : revert cases
  ! We call directly the band_grid function to avoid loops
  if (on_southern_hemisphere(this%zone) .or. i > get_nb_lat2_Global_grid()) then
    ! With the correct neighbour
    call south_neighbour_cell_Band_grid(j_start, j_end, i, j2, this%grid, &
      i-1, offset, this%i_pos, j_pos, this%zone)
  else
    i_neighb2 = i - 1
    if (present(i_neighb)) i_neighb2 = i_neighb

    call north_neighbour_cell_Band_grid(j_start, j_end, i, j2, this%grid, &
      i_neighb2, offset, this%i_pos, j_pos, this%zone)
  end if

  if (single) call correct_cell(j_start, sector, i-1, this, 1)
  if (single) call correct_cell(j_end, sector, i-1, this, 1)

end subroutine

!-----------------------------------------------------------------------------
!> Correct cell position for single partition case. We want to go back to the
!> first sector for the computation and then adjust it.
!-----------------------------------------------------------------------------
subroutine correct_cell(j, sector, i, this, increase)
  type (Partition) :: this
  integer, intent(in) :: i, increase
  integer :: sector, j, nb_cells

  nb_cells = this%grid%nb_lon(i)/3
  if (increase == 0) then
    sector = 0
    do while (j > nb_cells) 
      j = j - nb_cells
      sector = sector + 1
    end do
  else
    j = j + nb_cells*sector
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Returns the south neighbours for cells (i,j) inside the partition. Returns only
!> the indice on the previous line. Also includes the ghost cells.
!> The caller must check there is a neighbour.
!> @param i_neighb : optional parameter if we want to specify another
!> neighbour latitude (used in south hemisphere)
!-----------------------------------------------------------------------------
subroutine south_neighbour_cell_Partition(j_start, j_end, i, j, this,&
    i_neighb)
  type (Partition) :: this
  integer, intent(in) :: i, j
  integer, intent(in), optional :: i_neighb
  integer :: j_start, j_end
  integer :: i_neighb2, offset
  integer :: nb_parts, j_pos
  integer :: j2, sector
  logical :: single

  j_start = -1
  j_end = -1
  if (i == this%grid%nb_lat) return

  ! For a single partition, we must go back to one sector
  j2 = j
  single = has_single_partition()
  if (single) call correct_cell(j2, sector, i, this, 0)

  ! At the equator, 1-to-1 correspondance
  if (is_latitude_just_before_equator(i, this)) then
    call neighbour_equator(j_start, j_end, i, j2, i+1, sector, single, this)
    !if (i == 90 .and. j == 3) print *, "inside", j_start, j_end
    return
  end if

  nb_parts = nb_parts_on_band(this%i_pos)                                                  
  j_pos = translate_jpos(this%j_pos, nb_parts) 
  offset = this%grid%south_offset(i)

  if (on_southern_hemisphere(this%zone) .or. i > get_nb_lat2_Global_grid()) then
    call north_neighbour_cell_Band_grid(j_start, j_end, i, j2, this%grid, &
      i+1, offset, this%i_pos, j_pos, this%zone)

  else
    i_neighb2 = i + 1
    if (present(i_neighb)) i_neighb2 = i_neighb

    call south_neighbour_cell_Band_grid(j_start, j_end, i, j2, this%grid, &
      i_neighb2, offset, this%i_pos, j_pos, this%zone)
  end if

  !  if (this%id == 1 .and. i == 1 .and. j == 89) then
  !if (i == 81 .and. j == 1) then
  !  print *, "inside", j_start, j_end, "for", i, j
  !end if


  if (single) call correct_cell(j_start, sector, i+1 , this, 1)
  if (single) call correct_cell(j_end, sector, i+1, this, 1)

end subroutine

!-----------------------------------------------------------------------------
!> Returns true if the latitude line i is just before the equator
!-----------------------------------------------------------------------------
function is_latitude_just_before_equator(i, this) result (equator)
  type (Partition) :: this
  integer, intent(in) :: i
  logical :: equator
  integer :: i_lat, nb_lat2

  nb_lat2 = get_nb_lat2_Global_grid()
  i_lat = get_last_lat_Partition(this) 

  ! Nb_lat is always the line before south ghost cells at the equator 
  equator = (is_just_before_equator(this%i_pos) .and. i == i_lat) 
  ! Or if the partition is after and it is the true first line
  equator = equator .or. (is_just_after_equator(this%i_pos) &
    .and. is_north_ghost(i, this%grid)) 
  ! Or for a single partition
  equator = equator .or. (has_single_partition() .and. i == nb_lat2)
  ! Or for 3 partitions
  equator = equator .or. (has_three_partitions() .and. i == nb_lat2)
end function


!-----------------------------------------------------------------------------
!> Returns true if the latitude line i is just after the equator
!> @param i : latitude line inside the partition
!-----------------------------------------------------------------------------
function is_latitude_just_after_equator(i, this) result (equator)
  type (Partition) :: this
  integer, intent(in) :: i
  logical :: equator
  integer :: nb_lat, nb_lat2
  integer :: i_lat

  nb_lat2 = get_nb_lat2_Global_grid()
  i_lat = get_true_nb_lat_Partition(this)
  ! If the partition is after and it is the first line

  ! Always ghost cells at the equator
  equator = (is_just_after_equator(this%i_pos) .and. i == 2)
  ! Or if the partition is before and it is the true last line (ghost)
  equator = equator .or. (is_just_before_equator(this%i_pos) &
    .and. is_south_ghost(i, this%grid)) 
  ! Or for a single partition
  equator = equator .or. (has_single_partition() .and. i == nb_lat2 + 1) 
  ! Or for three partitions
  equator = equator .or. (has_three_partitions() .and. i == nb_lat2 + 1) 
end function

!-----------------------------------------------------------------------------
!> Returns the east neighbour for cell (i,j) inside the partition. Returns only
!> the indice on the current line. Also includes the ghost cells.
!> Does not work for cells outside the grid.
!> The caller must check there is a neighbour.
!-----------------------------------------------------------------------------
subroutine east_neighbour_cell_Partition(j_neighb, j)
  integer, intent(in) :: j
  integer :: j_neighb

  j_neighb = j + 1
end subroutine

!-----------------------------------------------------------------------------
!> Returns the west neighbour for cell (i,j) inside the partition. Returns only
!> the indice on the current line. Also includes the ghost cells.
!> Does not work for cells outside the grid.
!> The caller must check there is a neighbour.
!-----------------------------------------------------------------------------
subroutine west_neighbour_cell_Partition(j_neighb, j)
  integer, intent(in) :: j
  integer :: j_neighb

  j_neighb = j - 1
end subroutine

!-----------------------------------------------------------------------------
!> Number of neighbours for the current cell inside the grid. Does not 
!> include the ghost cells.
!> @param i : latitude line inside the partition
!> @param j : longitude on the latitude line inside the partition
!-----------------------------------------------------------------------------
function nb_neighbours_cell_Partition(this, i, j) result(nb)
  type(Partition), target :: this
  type(Band_grid), pointer :: grid
  integer, intent(in) :: i, j
  integer :: j_start, j_end
  integer :: nb

  nb = 0
  grid => this%grid
  if (j > 1) then 
    call west_neighbour_cell_Partition(j_start, j)
    if (is_indice_inside(j_start, grid, i)) nb = nb + 1
  end if
  if (j < grid%nb_lon(i)) then 
    call east_neighbour_cell_Partition(j_start, j)
    if (is_indice_inside(j_start, grid, i)) nb = nb + 1
  end if
  ! North
  nb = nb + nb_neighbours_north_cell_Partition(i, j, this)
  ! South
  nb = nb + nb_neighbours_south_cell_Partition(i, j, this)

end function

!-----------------------------------------------------------------------------
!> Number of north neighbours for the current cell inside the grid. Does not 
!> include the ghost cells.
!> @param i, j : cell coordinates inside the partition
!-----------------------------------------------------------------------------
function nb_neighbours_north_cell_Partition(i, j, this) result(nb)
  type(Partition), target :: this
  integer, intent(in) :: i, j
  integer :: j_start, j_end, nb

  nb = 0
  if (i > 1) then 
    call north_neighbour_cell_Partition(j_start, j_end, i, j, this)
    ! Indices are already checked
    nb = j_end - j_start + 1
  end if
end function

!-----------------------------------------------------------------------------
!> Number of south neighbours for the current cell inside the grid. Does not 
!> include the ghost cells.
!> @param i, j : cell coordinates inside the partition
!-----------------------------------------------------------------------------
function nb_neighbours_south_cell_Partition(i, j, this) result(nb)
  type(Partition), target :: this
  integer, intent(in) :: i, j
  integer :: j_start, j_end, nb

  nb = 0
  if (i < this%grid%nb_lat) then 
    call south_neighbour_cell_Partition(j_start, j_end, i, j, this)
    ! Indices are already checked
    nb = j_end - j_start + 1
  end if
end function

!-----------------------------------------------------------------------------
!> Returns number of cells at latitude i (in the grid), without the ghost 
!>_cells
!-----------------------------------------------------------------------------
function get_nb_lon_Partition(i, this) result(nb_lon)
  type (Partition) :: this
  integer :: i, nb_lon
  integer :: nb_ghosts

  nb_lon = get_nb_lon(i, this%grid)
end function

!-----------------------------------------------------------------------------
!> Returns number of latitudes for the current grid, without the ghost cells
!-----------------------------------------------------------------------------
function get_nb_lat_Partition(this) result(nb)
  type (Partition) :: this
  integer :: nb

  nb = get_nb_lat_Band_grid(this%grid)
end function

!-----------------------------------------------------------------------------
!> Returns number of latitudes for the current grid, with the ghost cells
!-----------------------------------------------------------------------------
function get_true_nb_lat_Partition(this) result(nb)
  type (Partition) :: this
  integer :: nb

  nb = this%grid%nb_lat
end function

!-----------------------------------------------------------------------------
!> Returns number of longitudes for the current grid, with the ghost cells
!-----------------------------------------------------------------------------
function get_true_nb_lon_Partition(i, this) result(nb)
  type (Partition) :: this
  integer, intent(in) :: i
  integer :: nb

  nb = this%grid%nb_lon(i)
end function

!-----------------------------------------------------------------------------
!> Returns the first latitude for the current grid, without the ghost cells
!-----------------------------------------------------------------------------
function get_first_lat_Partition(this) result(i_lat)
  type (Partition) :: this
  integer :: i_lat

  i_lat = first_lat_Band_grid(this%grid)
end function

!-----------------------------------------------------------------------------
!> Returns the last latitude for the current grid, without the ghost cells
!-----------------------------------------------------------------------------
function get_last_lat_Partition(this) result(i_lat)
  type (Partition) :: this
  integer :: i_lat

  i_lat = last_lat_Band_grid(this%grid)
end function

!-----------------------------------------------------------------------------
!> Returns the last latitude for the current grid, with the ghost cells
!-----------------------------------------------------------------------------
function get_true_last_lat_Partition(this) result(i_lat)
  type (Partition) :: this
  integer :: i_lat

  i_lat = this%grid%nb_lat
end function

!-----------------------------------------------------------------------------
!> Returns the first longitude for the current grid, without the ghost cells
!-----------------------------------------------------------------------------
function get_first_lon_Partition(i, this) result(i_lon)
  type (Partition) :: this
  integer, intent(in):: i
  integer :: i_lon

  i_lon = first_lon_Band_grid(i, this%grid)
end function

!-----------------------------------------------------------------------------
!> Returns the last longitude for the current grid, without the ghost cells
!-----------------------------------------------------------------------------
function get_last_lon_Partition(i, this) result(i_lon)
  type (Partition) :: this
  integer, intent(in):: i
  integer :: i_lon

  i_lon = last_lon_Band_grid(i, this%grid)
end function

!-----------------------------------------------------------------------------
!> Returns true if the partition is the last on the current band (on all 
!> sectors)
!> Works also for 3 partitions.
!-----------------------------------------------------------------------------
function is_last_on_band_all(this) result(res)
  type (Partition) :: this
  integer :: nb_parts
  logical :: res

  nb_parts = nb_parts_on_band_sector(this%i_pos)
  ! One partition on the band means we must go to the extremities
  res = (this%j_pos == 3*nb_parts)
end function

!-----------------------------------------------------------------------------
!> Returns true if the partition is the first on the current band (on a sector)
!> Works also for 3 partitions.
!-----------------------------------------------------------------------------
function is_first_on_band(this) result(res)
  type (Partition) :: this
  integer :: nb_parts
  logical :: res

  nb_parts = nb_parts_on_band_sector(this%i_pos)
  res = (this%j_pos == 1)
  res = res .or. (this%j_pos == nb_parts+1)
  res = res .or. (this%j_pos == 2*nb_parts+1)
end function


!-----------------------------------------------------------------------------
!> Returns true if the partition is the last on the current band (on a sector)
!> Works also for 3 partitions.
!-----------------------------------------------------------------------------
function is_last_on_band(this) result(res)
  type (Partition) :: this
  integer :: nb_parts
  logical :: res

  nb_parts = nb_parts_on_band_sector(this%i_pos)
  res = (this%j_pos == nb_parts)
  res = res .or. (this%j_pos == 2*nb_parts)
  res = res .or. (this%j_pos == 3*nb_parts)
end function

!-----------------------------------------------------------------------------
!> Check all merid winds are initialized
!-----------------------------------------------------------------------------
subroutine check_merid_winds_are_set(this)
  type(Partition) :: this
  double precision :: cur
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: i, j, k 

  if (this%id == 0) print *, "skipping merid winds check (partition_class)"
  return

  k = 1
  call interior_lat_indices(i_start, i_end, this%grid)
  if (has_north_ghost_cells(this%grid)) i_start = i_start - 1
  if (.not. has_south_ghost_cells(this%grid)) i_end = i_end - 1

  do i = i_start, i_end
    call lon_indices(j_start, j_end, this, i)
    do j = 1, nb_winds_merid_lat(i, this)
      cur = this%grid%merid_winds(k) 

      if (cur == undefined) then
        print *, "ijk", i, j, k, "id", this%id
        call print_error("Not all merid winds were initialized in the grid", &
          "set_data_Band_grid", fname_partition)
      end if
      k = k + 1
    end do
  end do
end subroutine



!-----------------------------------------------------------------------------
!> Set a value for meridional winds at position k in the wind array
!-----------------------------------------------------------------------------
subroutine set_merid_winds_Partition(val, k, array, this) 
  type(Partition) :: this
  integer, intent(in) :: k, array
  double precision :: val

  call set_merid_winds(val, k, array, this%grid) 
end subroutine

!-----------------------------------------------------------------------------
!> Returns the total number of cells (without ghost cells)
!-----------------------------------------------------------------------------
function nb_cells_Partition(this) result(nb)
  type(Partition) :: this
  integer :: nb, i
  integer :: i_start, i_end

  nb = 0
  i_start = get_first_lat_Partition(this)
  i_end = get_last_lat_Partition(this)
  do i = i_start, i_end
    nb = nb + get_nb_lon_Partition(i, this)
  end do

end function

subroutine allocate_advection_data(this)
  type(Partition) :: this
  integer :: size_merid

  size_merid = size_winds_merid(this)
  call allocate_advection_data_Band_grid(size_merid, this%grid)

#ifdef HDF5
  ! Define the winds and cell positions
  call generate_ratio_zonal_positions(this%grid)
  call generate_merid_positions(this)
#endif

end subroutine

#ifdef HDF5
!-------------------------------------------------------------------------------
!> We need some partition information so it cannot be done in band_grid
!-------------------------------------------------------------------------------
subroutine generate_merid_positions(this)
  type (Partition), target :: this
  type (Band_grid), pointer :: grid
  double precision :: cur_lat
  integer :: i_start, i_end, i
  integer :: j_start, j_end, j
  integer :: j_n1, j_n2, l
  integer :: nb_cells, k, k2

  k = 1
  grid => this%grid
  ! Write north interface (easier that way)
  call interior_lat_indices(i_start, i_end, grid)
  if (has_north_ghost_cells(grid)) i_start = i_start - 1
  if (.not. has_south_ghost_cells(grid)) i_end = i_end - 1

  do i = i_start, i_end
    ! Takes ghost cells too
    call lon_indices(j_start, j_end, this, i)
    cur_lat = cell_south_lat(i, grid)
    do j = j_start, j_end
      call south_neighbour_cell_Partition(j_n1, j_n2, i, j, this, i+1)
      do l = j_n1, j_n2
        if (to_or_from_interior_cell(l, i, j, i+1, this%grid)) then
          grid%v_lat(k) = cur_lat
          grid%v_lon(k) = cell_interface_middle(i, j, i+1, l, grid)
          k = k + 1
        end if
      end do
    end do
  end do

end subroutine
#endif

!-----------------------------------------------------------------------------
!> Allocate array of requests once for all neighbours (even though we may use
!> less than all)
!-----------------------------------------------------------------------------
subroutine init_requests(this)
  type (Partition) :: this
  integer :: nb_neighb

  ! Initi arrays for MPI communications
  nb_neighb = nb_neighbours_Partition(this) 

  ! We send only one array at the time
  if (nb_neighb > 0) then
    allocate (this%requests(2*nb_neighb))
    this%requests = mpi_request_null
  end if
end subroutine


subroutine free_Partition(this)
  type(Partition) :: this
  integer :: ierr

  call free_Band_grid(this%grid)

  if (this%mpi_partition /= mpi_datatype_null) then
    call mpi_type_free(this%mpi_partition, ierr)
    call check_error(ierr,"trying to free mpi_partition","free_Partition", &
      fname_partition)
  end if

  if (allocated(this%requests)) deallocate (this%requests)

end subroutine

!-------------------------------------------------------------------------------
!> Start asynchronous MPI exchange for either zonal or meridional advection.
!> We post sending and receiving requests. For meridional advection, we must
!> exchange with the 4 directions. 
!> @param datatype : cell ratio, slope or gradient
!> @param merid : true for meridional advection
!> @param zonal : true for zonal advection
!-------------------------------------------------------------------------------
subroutine start_mpi_exchange(this, datatype, zonal, merid)
  type (Partition) :: this
  logical, intent(in) :: zonal, merid
  integer, intent(in) :: datatype
  integer :: l

  this%requests = MPI_REQUEST_NULL

  l = 1
  if (zonal) then
    call start_receive(this, 3, this%requests, l, datatype)
    call start_receive(this, 1, this%requests, l, datatype)
  end if

  if (merid) then
    call start_receive(this, 0, this%requests, l, datatype)
    call start_receive(this, 2, this%requests, l, datatype)
  end if

  if (zonal) then
    call start_send(this, 1, this%requests, l, datatype)
    call start_send(this, 3, this%requests, l, datatype)
  end if

  if (merid) then
    call start_send(this, 2, this%requests, l, datatype)
    call start_send(this, 0, this%requests, l, datatype)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Receive interior cells from neighbour and store it in ghost cells.
!> Non-blocking
!> @param datatype : cell ratio, slope or gradient
!> @param direction : 0 for north, 3 for west
!-----------------------------------------------------------------------------
subroutine start_receive(this, direction, request, l, datatype)
  type(Partition), target :: this
  integer :: neighb, neighb2
  integer :: request(:), ierr
  integer :: mpi_border, direction
  integer :: status(mpi_status_size)
  integer :: k, l, p
  integer :: i_pos, tsize
  integer, intent(in) :: datatype

  ! If there is a neighbour
  call neighbour_from_direction(neighb, neighb2, this, direction)
  !  write (*,*) "neighbours of ",this%id, neighb, neighb2
  if (neighb < 0) neighb = neighb2
  if (neighb2 < 0) neighb2 = neighb
  if (neighb > -1 .or. neighb2 > -1) then

    ! Indice in the MPI border for the neighbour (the first and only)
    p = 1
    do k = neighb, neighb2

      i_pos = border_indice(this, direction, p)
      call start_receive_data(this%grid%tracers, i_pos, datatype, k, request, l)

      p = p + 1
    end do
  end if
end subroutine


!-----------------------------------------------------------------------------
!> Send interior cells to the ghost cells of the neighbour.
!> Non-blocking
!> @param datatype : cell ratio, slope or gradient
!> @param cells : true for sending cells, false for slope
!-----------------------------------------------------------------------------
subroutine start_send(this, direction, request, l, datatype)
  type(Partition), target :: this
  integer, intent(in) :: datatype
  integer :: request(:)
  integer :: neighb, neighb2, ierr
  integer :: mpi_border, direction
  integer :: k, l, p
  integer :: i_pos

  ! If there is a neighbour
  call neighbour_from_direction(neighb, neighb2, this, direction)
  if (neighb < 0) neighb = neighb2
  if (neighb2 < 0) neighb2 = neighb
  if (neighb > -1 .or. neighb2 > -1) then
    p = 1

    do k = neighb, neighb2
      ! Indice in the MPI border for the neighbour l
      i_pos = border_indice(this, direction, p)
      !print *, "sending ", datatype, "to ", k, "from", this%id

      call start_send_data(this%grid%tracers, i_pos, datatype, k, request, l)
      p = p + 1
    end do
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Wait for all requests. We don't use mpi_wait as it will block onto the
!> northmost process (which we want to avoid).
!> Instead, we loop on all and periodically checks each one of them
!-----------------------------------------------------------------------------
subroutine wait_for_all(requests)
  integer :: requests(:), nb
  !logical, allocatable :: completed(:)
  integer :: ierr
  integer, allocatable :: status(:,:)
  !integer :: nb_completed, i, nb

  nb = size(requests)
  allocate (status(mpi_status_size,nb))

  call mpi_waitall(nb, requests, status, ierr)
  deallocate (status)

end subroutine

end module
