!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Partitioning
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> List of partitions for MPI
!> Contains also the partitioning, which operates on the global grid (defind
!> with some integers, and not with an array of cells
!
!-------------------------------------------------------------------------------

module Partitioning_class
use Profiling
use Configuration_class
use Global_grid_class
use Partition_class
use Band_grid_class
use Analytical_partitioning
#ifdef HDF5
use hdf5
#endif

implicit none

public 
type Partitioning
  !> Number of partitions on each zone, numbered from the west to the east, and
  !> north to south
  integer, allocatable :: nb_parts(:)
  !> List of partitions
  type(Partition), allocatable :: parts(:)
end type

! For output
character(*), parameter :: fname_list = "partitioning_class.F90"


contains

! Must use this include for preprocessing io.inc also
#include "io.inc"

!-----------------------------------------------------------------------------
!> Constructor
!-----------------------------------------------------------------------------
subroutine new_Partitioning(this)
  type(Partitioning) :: this
  integer :: ierr, rank
  integer :: nb_procs

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call check_mpi_error(ierr, "get processus rank", "new_Partitioning",&
    fname_list)

  ! Create global grid and send it to other process
  if (rank == 0) call create_global_grid()
  call broadcast_Global_grid()

  ! Partition the grid
  if (rank == 0) call partitioning_root(this)

  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  call check_mpi_error(ierr, "get number processes", "new_Partitioning",&
    fname_list)

  if (nb_procs > 1) then
    ! Update the number of partitions as it can have changed
    call update_Configuration()

    call broadcast_partitioning(this, rank)

    ! Don't forget to keep only the adequate number of partitions in the master
    ! process
    if (rank == 0) call free_master_Partitioning(this)
  end if

  ! Now we can create the data for each process
  call partition_data(this)
  !call print_Partitioning(this)

  ! Define MPI type for advection
  if (nb_procs > 1) call new_MPI_borders_Partitioning(this)

end subroutine

subroutine create_global_grid()
  integer :: nb_lat2, tmp
  integer(kind=k12) :: t_start

  nb_lat2 = get_nb_lat2_Configuration()
  if (is_data_numerical()) then
    t_start = get_t_start()
    tmp = nb_lat2_from_file(1)
    if (nb_lat2 /= tmp) then
      print *, "nb in file :", tmp, "nb in config", nb_lat2
      call print_error("Wrong nb_lat2 : file has different nb_lat2",&
        "create_global_grid", fname_list)
    end if

  end if
  call new_Global_grid(nb_lat2)
end subroutine

!-----------------------------------------------------------------------------
!> Create partitioning on root process
!-----------------------------------------------------------------------------
subroutine partitioning_root(this)
  type(Partitioning) :: this
  type (Configuration) :: config
  integer :: total, tmp

  config = get_Configuration()
  ! First element contains all the partitions
  total = get_nb_partitions_Configuration(1)

  ! 3 identical sectors on each hemispheres
  ! Numbered from west to east, and north to south

  if (total == 1) then
    allocate (this%nb_parts(1))
    this%nb_parts(1) = 1
  else if (total == 3) then
    allocate (this%nb_parts(3))
    this%nb_parts(:) = 1
  else
    allocate (this%nb_parts(6))
    ! South hemisphere has the first half
    tmp = total / 3
    this%nb_parts(4:6) = tmp / 2
    ! North hemisphere has the rest
    this%nb_parts(1:3) = tmp - (tmp/ 2)
  end if

  call set_nb_partitions_Configuration(this%nb_parts)

  call allocate_Partitioning(this)

  ! Special case : only one partition for all sectors
  if (total == 1) then
    call single_partition(this)
    ! Special case : south sectours are empty
  else if (total == 3) then
    call three_partitions(this)
  else
    call multiple_partitions(this)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Allocate the partitioning
!-----------------------------------------------------------------------------
subroutine allocate_Partitioning(this)
  type (Partitioning) :: this
  integer :: i, nb_parts

  nb_parts = sum(this%nb_parts)

  allocate(this%parts(nb_parts))
  ! Numbered by sector, with normal ordering in a sector
  do i = 1, nb_parts
    call init_Partition(this%parts(i), i-1)
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Create MPI type for borders cells
!-----------------------------------------------------------------------------
subroutine new_MPI_borders_Partitioning(this)
  type(Partitioning) :: this
  integer :: i

  do i = 1, sum(this%nb_parts)
    call new_MPI_borders_Partition(this%parts(i))
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Creates a single partition covering all the globe. Everything is done here
!> (easier).
!-----------------------------------------------------------------------------
subroutine single_partition(this)
  type(Partitioning), target :: this
  type(Partition), pointer :: part
  integer :: nb_lat, i
  double precision :: lat_min, dlat(1)
  integer, allocatable :: nb_lon(:)
  double precision, allocatable :: lon_min(:), dlon(:)
  integer :: nb_bands(6), nb_lat2

  call print_mesg("Single partition.")
  nb_lat = get_nb_lat_Global_grid()
  dlat = get_dlat_Global_grid()
  lat_min = 90.

  nb_bands = (/1, 0, 0, 0, 0, 0/)
  call set_nb_bands_Configuration(nb_bands)
  part => this%parts(1)
  call new_Partition(part, 1, 1, 0, 0, 1)

  allocate (nb_lon(nb_lat))
  allocate (lon_min(nb_lat))
  allocate (dlon(nb_lat))
  lon_min = 0.
  do i = 1, nb_lat
    nb_lon(i) = nb_cells_lat(i)
    dlon(i) = 360.d0/nb_lon(i)
  end do

  ! We use the default structure constructor
  call new_Band_Grid(part%grid, nb_lat, dlat, nb_lon, dlon, lat_min, &
    lon_min)
  ! Set by hand the rest
  part%grid%nb_ghosts_east = 0
  part%grid%nb_ghosts_west = 0

  nb_lat2 = get_nb_lat2_Global_grid()
  part%grid%north_offset(1:nb_lat2) = -1
  part%grid%north_offset(nb_lat2+1:nb_lat) = 0

  part%grid%south_offset(1:nb_lat2) = 0
  part%grid%south_offset(nb_lat2+1:nb_lat) = -1


end subroutine

!-----------------------------------------------------------------------------
!> Creates 3 partitions covering all the globe. It is a special case as there
!> are no sectors on the southern hemisphere. Everything is done here (easier).
!-----------------------------------------------------------------------------
subroutine three_partitions(this)
  type(Partitioning), target :: this
  type(Partition), pointer :: part, cur
  integer :: nb_lat, i
  integer :: nb_lat2
  double precision :: lat_min, tmp, dlat(1)
  integer, allocatable :: nb_lon(:)
  double precision, allocatable :: lon_min(:), dlon(:)
  integer :: nb_bands(6)

  nb_lat = get_nb_lat_Global_grid()
  dlat = get_dlat_Global_grid()
  lat_min = 90.

  nb_bands = (/1, 1, 1, 0, 0, 0/)
  call print_mesg("Three partitions")
  call set_nb_bands_Configuration(nb_bands)
  call set_nb_bands_square_Configuration(nb_bands)
  part => this%parts(1)
  call new_Partition(part, 1, 1, 0, 0, 1)

  allocate (nb_lon(nb_lat))
  allocate (lon_min(nb_lat))
  allocate (dlon(nb_lat))
  lon_min = 0.
  do i = 1, nb_lat
    tmp = dble(nb_cells_lat_sector(i))
    dlon(i) = 120.d0/tmp
    ! There are always east and west ghost cells
    nb_lon(i) = int(tmp) + 2
    lon_min(i) = -dlon(i)
  end do

  ! We use the default structure constructor
  ! Don't forget to set the ghost cells (east-west neighbours)
  call new_Band_Grid(part%grid, nb_lat, dlat, nb_lon, dlon, lat_min, &
    lon_min, 10)

  ! Set by hand the rest
  part%grid%nb_ghosts_east = 1
  part%grid%nb_ghosts_west = 1

  nb_lat2 = get_nb_lat2_Global_grid()
  part%grid%north_offset(1:nb_lat2) = -1
  part%grid%north_offset(nb_lat2+1:nb_lat) = 0

  part%grid%south_offset(1:nb_lat2) = 0
  part%grid%south_offset(nb_lat2+1:nb_lat) = -1

  ! Copy this partition to the other two sectors
  do i = 2, 3
    cur => this%parts(i)
    call new_Partition_copy(cur, part)
    ! Just lon min changes
    cur%grid%lon_min = part%grid%lon_min + 120.*(i-1)
    !print '(a,f15.7,a,f15.7)', "init", part%grid%lon_min(41), "copy", cur%grid%lon_min(41)
    ! And the zone and j_pos
    cur%zone =  i
    cur%j_pos = i
  end do
  !call print_Partitioning(this)
end subroutine


!-----------------------------------------------------------------------------
!> Creates partitions for all sectors
!-----------------------------------------------------------------------------
subroutine multiple_partitions(this)
  type(Partitioning), target :: this
  integer :: i

  ! Create partitions on the first sector : north and south hemisphere
  ! North
  call partition_grid(this, 1)
  ! South
  call partition_grid(this, 4)

  ! Do not forget to update the nb of parts in the list of partitions (in case
  ! it has changed)
  do i = 2, 3
    this%nb_parts(i) = this%nb_parts(1)
    this%nb_parts(3 + i) = this%nb_parts(4)
  end do

  ! Then we just copy this partitioning to other sectors. We only need to change
  ! minimal longitude
  do i = 2, 3
    call copy_to_other_sector(this, i, 1)
    call copy_to_other_sector(this, 3 + i, 4)
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Create partitioning in other sector by simple copy
!> @param zone : current zone
!> @param prev : previous zone
!-----------------------------------------------------------------------------
subroutine copy_to_other_sector(this, zone, prev)
  type(Partitioning), target :: this
  type(Partition), pointer :: part, part_prev
  integer ::  zone, prev, i
  integer :: i_start, i_end, i_prev
  integer :: nb_parts, sector, nb_bands
  integer :: offset

  ! Offset for lon_min
  if (prev == 1) then
    offset = (zone-1)*120
    i_prev = 1
  else
    offset = (zone-4)*120
    i_prev = sum(this%nb_parts(1:3)) + 1
  endif
  i_start = sum(this%nb_parts(1:zone-1)) + 1
  i_end = i_start + this%nb_parts(zone) - 1

  ! We only want the sector
  sector = sector_from_zone(zone)

  do i = i_start, i_end
    part => this%parts(i)
    part_prev => this%parts(i_prev)
    ! Create new partition by simple copiy
    call new_Partition_copy(part, part_prev)

    ! Only lon_min and j_pos change
    part%grid%lon_min = part%grid%lon_min + offset
    nb_parts = nb_parts_on_band_sector(part_prev%i_pos)*(sector - 1)
    ! j_pos must be smooth between the sectors
    part%j_pos = part_prev%j_pos + nb_parts

    ! And of course the zone
    part%zone = zone
    i_prev = i_prev + 1
  end do

  ! Don't forget to copy the number of bands
  nb_bands = get_nb_bands_Configuration(prev)
  call set_nb_bands_Configuration_zone(nb_bands, zone)
  nb_bands = get_nb_bands_square_Configuration(prev)
  call set_nb_bands_square_Configuration_zone(nb_bands, zone)
  nb_bands = get_nb_bands_resized_Configuration(prev)
  call set_nb_bands_resized_Configuration_zone(nb_bands, zone)

end subroutine


!-----------------------------------------------------------------------------
!> Compute the number of cells in all partition
!-----------------------------------------------------------------------------
function nb_cells_Partitioning(this) result (nb)
  type (Partitioning) :: this
  integer :: nb, i

  nb = 0
  do i = 1, get_nb_parts_Partitioning(this)
    nb = nb + nb_cells_Band_grid(this%parts(i)%grid)
  end do

end function


!-----------------------------------------------------------------------------
!> Attach a partition to each process
!> The process 0 contains all the partitions. Other process contains only 
!> one partition here
!-----------------------------------------------------------------------------
subroutine broadcast_partitioning(this, rank)
  type(Partitioning) :: this
  integer :: rank, ierr

  ! Wait the master has set the data before sending
  call mpi_barrier(mpi_comm_world, ierr)
  call check_mpi_error(ierr, "mpi barrier", "broadcast_partitioning", &
    fname_config)

  ! Non master processes only have one partition each
  if (rank > 0) then
    allocate (this%nb_parts(1))
    this%nb_parts = 1
    allocate(this%parts(1))
  end if

  ! First, exchange array dimensions
  call broadcast_dimensions(this, rank)

  ! Then define MPI type and exchange data 
  call broadcast_data(this, rank)

end subroutine

!-----------------------------------------------------------------------------
!> Send or receive dimensions in order to allocate them (used for creating the
!> grid).
!> @parameter rank : rank of the current process
!-----------------------------------------------------------------------------
subroutine broadcast_dimensions(this, rank)
  type(Partitioning) :: this
  type(Partition) :: dummy
  integer :: rank
  integer :: nb_parts, i
  integer :: nb_procs, ierr

  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  call check_mpi_error(ierr, "get number processes", "broadcast_partitioning", &
    fname_config)

  ! Send (if we are on the master) or receive (other processes)
  if (rank == 0) then
    nb_parts = get_nb_parts_Partitioning(this)
    ! We send active partitions
    print *, "nb parts", nb_parts, "nb procs", nb_procs
    do i = 2, nb_parts 
      call send_dimensions_Partition(this%parts(i), i-1, .False.)
    end do
    ! Idle procs (if we remove the last band for example)
    do i = nb_parts + 1, nb_procs
      call send_dimensions_Partition(dummy, i-1, .True.)
    end do
  else
    ! And we will allocate the grid there
    call receive_dimensions_Partition(this%parts(1), 0)
  end if
end subroutine


!-----------------------------------------------------------------------------
!> Send or receive dimensions in order to allocate them (used for creating the
!> grid).
!-----------------------------------------------------------------------------
subroutine broadcast_data(this, rank)
  type(Partitioning), target :: this
  integer :: rank
  integer :: nb_parts, i

  ! First create the needed MPI type for all partitions
  nb_parts = get_nb_parts_Partitioning(this)
  do i = 1, nb_parts
    ! Idle procs have their own MPI types
    call new_MPI_Partition(this%parts(i))
  end do

  if (rank == 0) then
    nb_parts = get_nb_parts_Partitioning(this)
    ! We send active partitions
    do i = 2, nb_parts
      ! Send the partition. Beware, for gfortran, we cannot use local variables
      ! but variables outside the subroutine scope
      call send_Partition(this%parts(i), i-1)
    end do
  else
    if (.not. is_idle_Partition(this%parts(1))) then
      call receive_Partition(this%parts(1), 0)
    end if
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Set winds and concentration for each process
!-----------------------------------------------------------------------------
subroutine partition_data(this)
  type(Partitioning), target :: this
  type(Partition), pointer :: cur
  integer :: k
  integer :: rank, ierr
  integer(kind=k12) :: t_start

  call mpi_barrier(mpi_comm_world, ierr)
  call check_mpi_error(ierr, "mpi barrier", "partition_data", &
    fname_config)
  call mpi_comm_rank(mpi_comm_world, rank, ierr) 

  if (is_idle_Partition(this%parts(1))) return

  call allocate_data_Partitioning(this)
#ifdef HDF5
  call set_hdf5_data_Partitioning(this)
#endif

  if (NO_RUN) then
    call print_warning("Not reading data")
    return
  end if

  if (rank == 0) call print_mesg("Numerical ratio (always)")
  t_start = get_t_start()
  call read_ratio(1, this)
  ! Read winds later, once we have dt

  if (rank == 0) then
    if (is_data_analytical()) then 
      call print_mesg("Analytical winds")
    else 
      call print_mesg("Numerical winds")
    end if
  end if

  ! Reading winds as an estimation, corrected later during advection
  do k = 1, get_nb_parts_Partitioning(this)
    cur => this%parts(k)
    if (is_data_analytical()) then
      call cv_winds(0.d0, 0.d0, cur)
    else
      call read_winds(1, 1, cur)
    end if
    !    call analytical_winds_Partition(cur)
  end do

  call check_data_is_set_Partitioning(this)
end subroutine

!-------------------------------------------------------------------------------
!> Check interior cells
!-------------------------------------------------------------------------------
subroutine check_data_is_set_Partitioning(this)
  type(Partitioning) :: this
  integer :: k

  do k = 1, get_nb_parts_Partitioning(this)
    call check_ratio_is_set(this%parts(k))
    call check_merid_winds_are_set(this%parts(k))
  end do
end subroutine

subroutine allocate_data_Partitioning(this)
  type(Partitioning), target :: this
  integer :: k

  do k = 1, get_nb_parts_Partitioning(this)
    call allocate_advection_data(this%parts(k))
  end do
end subroutine

#ifdef HDF5
subroutine set_hdf5_data_Partitioning(this)
  type(Partitioning), target :: this
  integer :: k

  do k = 1, get_nb_parts_Partitioning(this)
    call set_hdf5_data_Partition(this%parts(k))
  end do
end subroutine
#endif


!-----------------------------------------------------------------------------
!> Sum all the cells on the line nb i for all partitions left of the 
!> current partition
!> id is the global partition indice
!-----------------------------------------------------------------------------
function nb_cells_prev_lon_partition(this,i, id, j_pos) result (nb_prev_lon)
  type(Partitioning) :: this
  integer, intent(in) :: i, id, j_pos
  integer :: j, k_prev, nb_prev_lon

  nb_prev_lon = 0
  do j = 1, j_pos
    ! Indice of the previous partition
    k_prev = id - (j_pos - j)
    nb_prev_lon = nb_prev_lon + this%parts(k_prev)%grid%nb_lon(i)
  end do
end function

!-----------------------------------------------------------------------------
!> Partitioning procedure for the grid
!> The grid is filled with the largest "square" partitions (i.e, with the
!> same number of points
!> @param zone : zone indice
!-----------------------------------------------------------------------------
subroutine partition_grid(this, zone)
  type(Partitioning) :: this
  integer, intent(in) :: zone
  integer :: side, nb_bands, n_left
  integer :: height_corr, rest
  integer :: nb_waiting, nb_lat2
  integer :: nb_bands_square, nb_parts
  ! 0 is the default
  height_corr = 0

  nb_parts = this%nb_parts(zone)
  ! We set and get the closest square 
  call compute_nb_bands_Configuration(zone, nb_parts)
  nb_bands_square = get_nb_bands_square_Configuration(zone) 
  nb_bands = get_nb_bands_Configuration(zone) 

  ! Then we have the side of a square partition
  ! Except for the last band
  call compute_side_square_parts(side, n_left, rest, nb_bands_square, nb_parts)

  ! First we fill the last band as its height can vary
  if (n_left > 0) then
    call leftover_band(this, nb_bands_square, side, n_left, zone, height_corr)
  else
    ! Reset the number of bands then
    nb_bands = nb_bands_square
    call set_nb_bands_Configuration_zone(nb_bands, zone)

    ! Special case : the number of partitions is a square but its square root does
    ! not divides the number of latitudes
    nb_lat2 = get_nb_lat2_Global_grid()
    nb_waiting = nb_lat2 - nb_bands_square*side 
    if (nb_waiting > 0) then
      write (*,*) "Square number but missing ",nb_waiting, "latitudes... We add them."
      height_corr = -nb_waiting
    end if
  end if

  ! Then set the squares (which can be rectangles as the height rest is splitted
  ! into several latitudes)
  call split_to_semi_squares(height_corr, this, nb_bands_square, side, &
    zone, n_left)

  !nb_bands = get_nb_partitions_Configuration(zone) 
  !write (*,*) "parts on zone", zone, nb_bands
  !write (*,*) "modified nb bands , square", nb_bands, nb_bands_square, "zone", zone
end subroutine


!-----------------------------------------------------------------------------
!> Split global grid to smaller grid with the same number of elements
!> All are "squared", except at the middle, where it is triangular
!> Also, we may extend the partitions by 1 latitude if height_corr is > 0 for 
!> a more balanced load the partition on the last band.
!> Finally, they may not be squared if the number of
!> bands does not divides exactly the number of latitudes. This is a pseudo last
!> band as the number of partitions is double precisionly a square number
!> @param k_start : first position in the partition array
!> @param nb_bands : number of bands with square partitions
!> @param side : height for square partitions
!> @param zone : zone indice
!> @param n_left : number of partitions on an eventual leftover band
!-----------------------------------------------------------------------------
subroutine split_to_semi_squares(height_corr, this, nb_bands, side, zone,&
    n_left)
  type(Partitioning), target :: this
  integer,  intent(in) :: nb_bands, side, zone
  integer ::  height_corr
  integer :: i_lat, n_left
  integer :: k, dk, i_corr
  integer :: i_pos, di_pos

  !> With a grid gaining 2 cells at each latitude, we have the following result :
  !> If we split the grid in bands of equal height k, there are (2i+1)
  !> rectangular grid with k^2 cells in each, i being the band indice (from 0)
  ! So nb_squares is 2i-1
  ! Find the starting indice
  if (on_northern_hemisphere(zone)) then 
    ! All nb parts may not have been updated yet, so we use only the first
    ! sector
    k = (zone-1)*this%nb_parts(1) + 1
    !k = sum(this%nb_parts(1:zone-1)) + 1
    ! i_corr must be outside the current hemisphere by default
    i_corr = get_nb_bands_Configuration(zone) + 1
  else
    ! Start from the south pole
    ! Same here
    k = 3*this%nb_parts(1) + this%nb_parts(4)
    i_corr = get_total_nb_bands_Configuration() + 1
  end if
  !write (*,*) "k here", k, "zone", zone

  ! Creating the partition is done the same way on both hemispheres, we just
  ! change the position
  i_lat = 0
  i_pos = 1
  di_pos = 1
  dk = 1

  if (.not. on_northern_hemisphere(zone)) then
    ! We suppose the 3 sectors identical
    i_pos = get_nb_bands_Configuration(1) + get_nb_bands_Configuration(4) 
    dk = -1
    di_pos = -1
  end if

  !write (*,*) "heightcorr", height_corr, "icorr", i_corr
  if (height_corr /= 0 ) then
    i_corr = nb_bands - abs(height_corr) + 1
    call set_nb_bands_resized_Configuration_zone(abs(height_corr), zone)
  end if
  write (*,*) "Correcting height of ", height_corr, abs(height_corr)

  ! Splitting is done there
  call true_semi_square_splitting(this, k, dk, i_lat, i_pos, side, &
    nb_bands, height_corr, i_corr, di_pos, zone)

end subroutine

!-------------------------------------------------------------------------------
!> Split a part of the global grid into pseud square partitions. We do not worry 
!> about hemispheres, as the indices are computed for that by the caller.
!-------------------------------------------------------------------------------
subroutine true_semi_square_splitting(this, k, dk, first_lat, i_pos, side, &
    nb_bands, height_corr, i_corr, di_pos, zone)
  type(Partitioning), target :: this
  integer, intent(in) :: nb_bands, i_corr, dk
  integer, intent(in) :: di_pos, side
  integer, intent(in) :: height_corr, zone
  integer :: height, width
  logical :: next_first_resize
  integer :: i, nb_resized
  integer :: nb_squares, first_lat, i_pos, j
  integer :: k, j_start, j_end
  integer :: i_pos2

  nb_resized = 0
  height = side
  width = side
  ! Number of squares on the current band
  nb_squares = 1

  ! i_pos2 always start from 1 (for overlap)
  i_pos2 = 1
  do i = 1, nb_bands
    ! Tell us if the next line will be the first to have resized partitions
    next_first_resize = (i == (i_corr))
    ! If we must extends the partitions, only extends height for 1
    if (i > i_corr - 1) then
      ! nb_resized is the number of partitions also resized in width (one only
      ! half). It only happens when at least 2 partitions bands are modified
      nb_resized = 0
      if (i > i_corr) then
        nb_resized = (i - i_corr - 1) + 1
      end if
      height = side + 1
    end if

    ! On the southern hemisphere, we fill the array starting from the end
    j_start = 1
    j_end = nb_squares
    if (dk < 0) then
      j_start = nb_squares
      j_end = 1
    end if

    do j = j_start, j_end, dk
      ! Set the partition position and the grid
      call new_Partition_from_pos(this%parts(k), i_pos, j, first_lat, height, &
        width, nb_resized, zone, height_corr)

      ! And the width modification (for neighbours)
      call set_partition_overlap(this%parts(k), i_pos2, j , nb_resized, &
        nb_squares, nb_bands, next_first_resize)
      ! On the southern hemisphere, north and south are reverted
      if (.not. on_northern_hemisphere(zone)) then
        call revert_overlap(this%parts(k))
      end if

      k = k + dk
    end do

    nb_squares = nb_squares + 2
    first_lat = first_lat + height
    i_pos = i_pos + di_pos
    i_pos2 = i_pos2 + 1
  end do
  nb_squares = nb_squares - 2
end subroutine

!-----------------------------------------------------------------------------
!> Find the north south neighbours is easy. For a partition, it is the partition
!> on the line before. However, we have also a single cell either to the left or
!to the right. If the partitions are resized, this single cell may be on the
!> other side. This is what this function sets : an overlap which tells us if the
!> single cell is inverted.
!> For a partition j on a band, we set if its (north-south) overlap is 
!> inverted. 1 for the left overlap and 2 for the right (both being 3)
!> There can be 2 inversions if we are after the number of resized
!> @param i : band number (always on northern hemisphere)
!-----------------------------------------------------------------------------
subroutine set_partition_overlap(cur, i, j, nb_resized, nb_squares, nb_bands, &
    next_first_resize)
  type (Partition) :: cur
  integer, intent(in) :: i, j, nb_resized, nb_squares, nb_bands
  logical, intent(in) :: next_first_resize
  integer :: overlap, mid, j2
  integer :: left, right
  overlap = 0

  ! Not for the last band of course
  if (i > nb_bands + 1) return

  ! North and south start with at least one partition resized in width
  if (nb_resized < 1) then
    ! If the line is just before the first, everything is interverted for the
    ! south
    if (next_first_resize) then
      overlap = 30
      ! Only one overlap is inverted at the extremities
      if (j == 1) overlap = 20
      if (j == nb_squares) overlap = 10
    end if
    cur%overlap = overlap
    !if (overlap > 0) write (*,*) "overlap for", i, j, overlap
    return
  end if

  mid = (nb_squares + 1)/2
  ! North overlap start with the first line which has been modified and up to
  ! the center partition (idem for the other half)
  if(j /= mid) then
    j2 = j
    left = 1
    right = 2

    ! If we are on the other side, just change the variable, principe is the
    ! same by symetry
    if (j > mid) then
      j2 = nb_squares - j + 1
      ! Reversed of the other half
      left = 2
      right = 1
    end if

    ! Overlap with the partition starting from the "diagonal"
    if (j2 > nb_resized - 1) then
      overlap = right
      ! Add the left after the diagonal
      if (j2 > nb_resized) overlap = overlap + left

      ! Special cases : if we are just before/after the middle, we do not change
      ! the right/left boundary
      if (j == mid - 1) overlap = 1
      if (j == mid + 1) overlap = 2

      ! Only one overlap is inverted at the extremities
      if (j == 1) overlap = 2
      if (j == nb_squares) overlap = 1
      !  if (i == 4) write (*,*) "init overlap", overlap
    end if
  end if

  ! South overlap start with the first line but ends at the band before the last
  if (i < nb_bands) then
    ! Also at the middle though
    j2 = j
    left = 10
    right = 20
    if (j > mid) then
      ! On the other side, roles are reversed
      j2 = nb_squares - j + 1
      left = 20
      right = 10
    end if

    ! Overlap with the partition after the "diagonal"
    if (j2 > nb_resized ) then
      overlap = overlap + right
      ! Add only left overlap after the diagonal + 1
      if (j2 > nb_resized + 1) overlap = overlap + left

      ! Only one overlap is inverted at the extremities
      if (j == 1) overlap = overlap - left
      if (j == nb_squares) overlap = overlap - right
    end if

  end if
  cur%overlap = overlap
  !if (overlap > 0) write (*,*) "overlap for", i, j, overlap
end subroutine

!-----------------------------------------------------------------------------
!> Compute the last band height. If <= 0, then there is no last band, which 
!> contains the rest of the division in theorically equal squares
!> Also allocates the partitioning and create the last band
!> @param side : height of a square partition
!> @param zone : zone indice
!-----------------------------------------------------------------------------
subroutine leftover_band(this, nb_bands, side, n_left, zone, height_corr)
  ! We must declare this as a target for using a pointer on it
  type(Partitioning), target :: this
  integer, intent(in) :: nb_bands, side, n_left, zone
  integer :: height_corr, nb_lat
  integer :: width, width_left
  integer :: first_lat, height
  integer :: nb_bands_square, i_neighb
  integer :: nb_lat2

  nb_lat2 = nb_lat2_Global_grid()
  if (on_northern_hemisphere(zone)) then
    first_lat = nb_bands*side
    height = nb_lat2 - first_lat
    ! First lat is for ghost cells
    i_neighb = first_lat + 1
  else
    ! Must add one on the south hemisphere
    nb_lat = nb_lat_Global_grid()
    first_lat = nb_lat2 + 1
    height = nb_lat - nb_bands*side - first_lat + 1
    !write (*,*) "height", height
    i_neighb = nb_lat2 + height
  end if

  call split_last(width, width_left, i_neighb, n_left)
  !write (*,*) "width left", width, "ineighb", i_neighb
  !write (*,*) "nleft", n_left

  call adapt_height_leftover_band(height_corr, height, 0.99d0, side, width, &
    width_left, first_lat, n_left, zone, this)
  !write (*,*) "last height before and after", height, height - height_corr
  height = height - height_corr

  ! If we remove the last band
  if (height == 0) then
    this%nb_parts(zone) = this%nb_parts(zone) - n_left
  end if

  ! Now we can allocate
  if (height > 0) then
    call create_leftover_band(this, nb_bands, height, height_corr, width, &
      n_left, nb_lat, zone)
  else
    ! No last band so we save it
    height_corr = -height_corr
    nb_bands_square = get_nb_bands_square_Configuration(zone) 
    call set_nb_bands_Configuration_zone(nb_bands_square, zone)
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Create last band partitions (setting are now done)
!-------------------------------------------------------------------------------
subroutine create_leftover_band(this, nb_bands, height, height_corr, width, &
    n_left, nb_lat, zone)
  type(Partitioning) :: this
  integer :: height, nb_lat, height_corr
  integer :: j, k, n_left, zone
  integer :: nb_bands, width
  integer :: cur_band

  ! Now we create the last band partitions with the new height and eventual
  ! south ghost cells (north are always here)

  ! Find the starting indice
  if (on_northern_hemisphere(zone)) then
    ! Up from the equator
    k = sum(this%nb_parts(1:zone)) - n_left + 1
    cur_band = nb_bands + 1
  else 
    ! After the equator
    k = sum(this%nb_parts(1:3)) + 1!n_left
    cur_band = get_total_nb_bands_Configuration_sector() - nb_bands
  end if

  ! All partitions have the same number of cells, except the last one
  do j=1, n_left 
    !write (*,*) "k in leftover", k
    call new_Partition_leftover(this%parts(k), cur_band, j, height, width, &
      n_left, zone)

    this%parts(k)%height_var = -height_corr - 1
    k = k + 1
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Modify the height of the last band according to the load ratio to equilibrate 
!>_the load
!-------------------------------------------------------------------------------
subroutine adapt_height_leftover_band(height_corr, height, load_ratio, side, &
    width, width_left, first_lat, n_left, zone, this)
  type(Partitioning) :: this
  integer :: height
  integer :: first_lat
  integer :: side, width, width_left, n_left
  integer, intent(in) :: zone
  double precision, intent(in) :: load_ratio
  double precision :: limit, ratio
  integer :: load_last, load_last_rest
  integer :: load_square, height_corr
  integer :: i, di, nb_lat2
  integer :: nb_parts
  load_last_rest = 0

  ! nb for square partitions
  load_square = side*side

  call max_load_last(load_last, width, width_left, height)
  !write (*,*) "max load on last", load_last

  height_corr = 0
  ! Only continue if the last band is overloaded
  if (load_last < load_square + 1) then
    write (*,*) "The last band is not the limiting factor, OK"
    return
  endif

  write (*,*) "The last band is overloaded, trying to get load ratio under ", load_ratio
  nb_lat2 = nb_lat2_Global_grid()
  ratio = dble(load_last)/load_square
  !  first_lat = nb_lat_Global_grid() - height + 1
  ! Not the "true" first lat, but the one which will be needed to compute the
  ! width (in southern hemisphere, it is the last)
  if (on_northern_hemisphere(zone)) then
    first_lat = nb_lat2 - height + 1
    di = 1
  else
    first_lat = nb_lat2 + height
    di = -1
  end if

  limit = load_ratio - 0.00001
  if (ratio  > limit) then
    write (*,'(a,f13.4,a, f13.4)') "Ratio too far away from input (", ratio,")"
    ! We try to remove one line after the other, as we must split the new line
    ! again...
    do i = 1, height -1
      call split_last(width, width_left, first_lat + i*di, n_left)
      ! Get new load
      call max_load_last(load_last, width, width_left, height -i)
      ! Square partitions are now rectangular and just extends
      ! But at the center, it is still square
      load_square = (side + i)**2
      ratio = dble(load_last)/load_square
      if (ratio < limit) then
        height_corr = i
        exit
      endif
    end do

    if (height_corr > 0) then
      write (*,*) "We can satisfy the load limit. Load found :", ratio
      write (*,*) "Will remove", height_corr, "latitude lines on last band"
    else
      write (*,*) "Cannot satisfy load margin, last band should be removed"
      write (*,*)  "but this functionality is not yet implemented"
      !write (*,*) "Cannot satisfy load margin, removing last band entirely."
      !height_corr = height

      !! Remove from global configuration also
      !nb_parts = get_nb_partitions_Configuration(zone)
      !nb_parts = nb_parts - n_left
      !call set_nb_partitions_Configuration_zone(nb_parts, zone)
    end if
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Small function for splitting a line in the last band and get the quotient and
!> rest
!-----------------------------------------------------------------------------
subroutine split_last(width, width_left, i, n_left)
  integer :: width, width_left, n_left
  integer :: nb_cells
  integer, intent(in) :: i

  nb_cells = nb_cells_lat_sector(i)
  width = nb_cells / n_left
  width_left = nb_cells - (n_left - 1)*width
end subroutine


!-----------------------------------------------------------------------------
!> Compute the maximal load (i.e number of cells) for partitions on
!> the last band. 
!-----------------------------------------------------------------------------
subroutine max_load_last(load, width, width_left, height)
  integer, intent(in) :: width, width_left, height
  integer :: load
  integer :: nb_last_square, nb_last_rest
  nb_last_rest = 0

  ! The larger has the larger width
  nb_last_square = height*width
  ! If triangular, we use grid property (+2 at each new latitude)
  if (width_left > 0) nb_last_rest = height*width_left + 2*(height - 1)

  !write (*,*) "width, widthleft", width, width_left
  !write (*,*) "load last band", nb_last_square, nb_last_rest
  load = max(nb_last_square, nb_last_rest)
end subroutine

!-----------------------------------------------------------------------------
!> Find the number of longitude for the i-th partition on the last band
!> The aim is to continue the partitioning of the previous band
!> There are at most nb_prev+1 partitions (otherwise a new band is created)
!-----------------------------------------------------------------------------
subroutine get_nb_lon_leftover_band(i, lat_prev, nb_lon, nb_left, side, nb_prev)
  integer, intent(in) :: i, lat_prev, side, nb_left, nb_prev
  integer, intent(inout) :: nb_lon(:)
  integer :: j, j_start, j_end, mid, nb_lon_mid_prev
  integer :: q, r
  logical :: special_mid 

  ! Nb of cells in the middle partition at previous band
  nb_lon_mid_prev = nb_cells_lat_sector(lat_prev) - (nb_prev-1)*side
  q = nb_prev/nb_left
  r = modulo(nb_prev, nb_left)

  mid = (nb_prev+1)/2
  special_mid = .False.
  ! If nb_prev + 1 partitions, the central one is split in two, with the left
  ! partition being the smallest
  if (nb_left > nb_prev) then
    special_mid = .True.
    q = 1
    r = 0
  end if

  ! The rest of the division is set starting from the right
  if (i < r + 1) then
    j_start = (i-1)*(q+1) + 1
    j_end = j_start + q 
  else
    j_start = r*(q+1) + (i-r-1)*q + 1
    j_end = j_start + q - 1
  end if
  do j = 1, size(nb_lon)
    nb_lon(j) = create_chunk(j_start, j_end, j, nb_lon_mid_prev, side, mid, &
      special_mid)
  end do
end subroutine

!-----------------------------------------------------------------------------
!> This function set the width of a partition on the last band. The idea is to
!> use a width multiple of the previous partition, so either from the square
!> partition or the central one
!-----------------------------------------------------------------------------
function create_chunk(j_start, j_end, j, nb_lon_mid_prev, side, mid, &
    special_mid) result(width)
  integer, intent(in) :: j_start, j_end, j, nb_lon_mid_prev
  integer, intent(in) :: side, mid
  integer :: k, width, tmp
  logical, intent(in) :: special_mid

  width = 0
  do k=j_start, j_end
    if (k == mid) then
      tmp = nb_lon_mid_prev + 2*j
      if (special_mid) then
        width = width + tmp/2
      else
        width = width + tmp
      end if
      elseif (k == mid+1 .and. special_mid) then
      tmp = nb_lon_mid_prev + 2*j
      width = width + (tmp - (tmp/2))
    else
      width = width + side
    end if
  end do
  return
end function

!-----------------------------------------------------------------------------
!> print partitioning to stdout and grid (debug)
!-----------------------------------------------------------------------------
subroutine print_Partitioning(this)
  type(Partitioning) :: this
  integer :: i, rank, ierr

  call mpi_comm_rank(mpi_comm_world, rank, ierr) 
  call check_mpi_error(ierr, "process rank", "print_Partitioning", &
    fname_config)

  do i=1, sum(this%nb_parts)
    call print_Partition(this%parts(i))
  end do
end subroutine


!-------------------------------------------------------------------------------
!> Compute the side for square partitions (most of the partitioning).
!> Return also the number of leftover partitions and the number of latitude left
!-------------------------------------------------------------------------------
subroutine compute_side_square_parts(side, n_left, rest, nb_bands, nb_parts)
  integer, intent(in) :: nb_bands, nb_parts
  integer :: n_left, side
  integer :: rest
  integer :: nb_lat

  rest = 0
  nb_lat = get_nb_lat2_Global_grid()
  ! If the number of partitions is a square, partitioning is exact
  ! Otherwise, there is a band with the rest
  n_left = nb_parts - nb_bands ** 2
  if (n_left <= 0) then
    side = int(nb_lat/nb_bands)
    rest = nb_lat - nb_bands*side
  else
    side = int((nb_lat-1)/nb_bands)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> The master process contains all partitions, we remove them except for the
!> first
!-----------------------------------------------------------------------------
subroutine free_master_Partitioning(this)
  type(Partitioning) :: this
  type(Partitioning) :: tmp
  integer :: k, nb_parts
  integer :: nb_procs, ierr

  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  call check_mpi_error(ierr, "get number processes", "free_master_Partitioning", &
    fname_list)

  nb_parts = get_nb_parts_Partitioning(this)

  ! OK for 1 process though.
  if (nb_procs == 1) return

  ! We suppose there is enough processes
  if (nb_procs < nb_parts)  then
    call print_error("Not enough processes for the partitions",&
      "free_master_Partitioning", fname_list)
  end if

  ! Free all partitions by the master, except the first
  do k = 2, nb_parts
    call free_Partition(this%parts(k))
  end do

  ! Resize list
  allocate (tmp%parts(1))
  tmp%parts(1) = this%parts(1) 
  deallocate (this%parts)
  allocate (this%parts(1))
  this%parts(1) = tmp%parts(1)
  deallocate (tmp%parts)

  deallocate (this%nb_parts)
  allocate (this%nb_parts(1))
  this%nb_parts = 1
end subroutine

!-----------------------------------------------------------------------------
!> Returns all number of partitions
!-----------------------------------------------------------------------------
function get_nb_parts_Partitioning(this) result (nb)
  type(Partitioning) :: this
  integer :: nb

  nb = sum(this%nb_parts)
end function

!-----------------------------------------------------------------------------
!> Total tracer mass (only for master process)
!-----------------------------------------------------------------------------
subroutine mass_tracer_Partitioning(mass, this) 
  type(Partitioning) :: this
  double precision, allocatable :: mass(:), mass_loc(:)
  integer :: i, ierr, k

  allocate(mass_loc(NB_TRACERS))
  mass_loc = 0
  do i = 1, get_nb_parts_Partitioning(this)
    do k = 1, nb_tracers
      mass_loc(k) = mass_loc(k) + mass_tracer(k, this%parts(i)%grid)
    end do
  end do

  mass = 0
  call mpi_reduce(mass_loc, mass, nb_tracers, mpi_double_precision, &
    mpi_sum, 0, mpi_comm_world, ierr)

end subroutine

end module
