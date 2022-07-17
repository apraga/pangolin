!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Tests_advection_2d
!
!> @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Tests the advection in the 4 directions : east, west, north, south
!
!-------------------------------------------------------------------------------

module Tests_advection_2d
use Pangolin_run
use Partitioning_class
use Advection_2d

! Show the details of borders ghosts
!logical :: show_all = .true. !.false.
logical :: show_all = .false.
integer :: TRACER = 1
contains

!-------------------------------------------------------------------------------
! Start all the tests
!-------------------------------------------------------------------------------
subroutine run_advection_2d_tests(nb_fails)
  integer :: rank, nb_fails
  integer :: nb_parts, ierr

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  !write (*,*) "rank", rank

  nb_parts = get_nb_parts_Partitioning(list_parts)
  call init_advection(list_parts, nb_parts)

  call west_to_east_advection_2d(nb_fails, rank)
  call east_to_west_advection_2d(nb_fails, rank)
  call to_poles_advection_2d(nb_fails, rank)
end subroutine

!-------------------------------------------------------------------------------
!> Advection to and from the closest pole. We set a concentration of 1 at the 
!> equator and propagate it. Needs 2 pass to fill every cell at the leftover band
!> The result will be a grid full of 1.
!------------------------------------------------------------------------------
subroutine to_poles_advection_2d(nb_fails, rank)
  type(Partition), pointer :: part
  integer :: nb, i, nb_procs
  integer :: rank, ierr, status(mpi_status_size)
  integer :: nb_fails, nb_fails_i
  logical, allocatable :: results(:)

  !call init_array_result(results, rank)

  nb = get_nb_parts_Partitioning(list_parts)

  ! Set and init data for comparaison
  call set_data_equator(nb)

  if (rank == 0) then
    write (*,'(a)', advance="no") "Checking equator to poles advection..."
  end if
  nb_fails_i = nb_fails
  do i = 1, nb
    part => list_parts%parts(i)

    if (on_northern_hemisphere(part%zone)) then
      call test_meridional_advection_to_north(part)
    else
      call test_meridional_advection_to_south(part)
    end if

  end do

  ! Need a second pass
  do i = 1, nb
    part => list_parts%parts(i)
    if (on_northern_hemisphere(part%zone)) then
      call test_meridional_advection_to_south(part)
    else
      call test_meridional_advection_to_north(part)
    end if
    call check_ratio(nb_fails, part, 2)
  end do

  call check_if_failed(nb_fails - nb_fails_i, rank)
end subroutine

!-------------------------------------------------------------------------------
!> Advection to the east. We set a concentration of 1 to the westmost cell and
!> propagate it to the west, cell by cell. The result will be a grid full of 1.
!-------------------------------------------------------------------------------
subroutine west_to_east_advection_2d(nb_fails, rank)
  type(Partition), pointer :: part
  integer :: nb, i, nb_procs
  integer :: rank, ierr, status(mpi_status_size)
  integer :: nb_fails, nb_fails_i
  logical, allocatable :: results(:)

  !call init_array_result(results, rank)

  nb = get_nb_parts_Partitioning(list_parts)

  ! Set data for comparing
  call set_data_west(nb)

  if (rank == 0) then
    write (*,'(a)', advance="no") "Checking west to east advection..."
  end if
  nb_fails_i = nb_fails
  do i = 1, nb
    part => list_parts%parts(i)
    call test_zonal_advection_west_to_east(part)

    call check_ratio(nb_fails, part, 1)
  end do

  call check_if_failed(nb_fails - nb_fails_i, rank)
end subroutine

!-------------------------------------------------------------------------------
!> Advection to the west. We set a concentration of 1 to the eastmost cell and
!> propagate it to the east, cell by cell. The result will be a grid full of 1.
!-------------------------------------------------------------------------------
subroutine east_to_west_advection_2d(nb_fails, rank)
  type(Partition), pointer :: part
  integer :: nb, i
  integer :: rank, ierr, status(mpi_status_size)
  integer :: nb_fails, nb_fails_i
  logical, allocatable :: results(:)

  !call init_array_result(results, rank)

  nb = get_nb_parts_Partitioning(list_parts)

  ! Set data for comparing
  call set_data_east(nb)

  !if (rank > 0) return

  if (rank == 0) then
    write (*,'(a)', advance="no") "Checking east to west advection..."
  end if
  do i = 1, nb
    part => list_parts%parts(i)
    call test_zonal_advection_east_to_west(part)

    call check_ratio(nb_fails, part, 3)
  end do

  call check_if_failed(nb_fails - nb_fails_i, rank)
end subroutine


!-------------------------------------------------------------------------------
! Init the the result array for master process
!-------------------------------------------------------------------------------
subroutine init_array_result(results, rank)
  integer :: rank, nb_procs, ierr
  logical, allocatable :: results(:)

  ! Array of results 
  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  allocate (results(nb_procs))
end subroutine

!-------------------------------------------------------------------------------
! Sums the number of fails (for a single test) and send it to master process,
! which will display the eventul error message
!-------------------------------------------------------------------------------
subroutine check_if_failed(nb_failed_loc, rank)
  integer :: i, rank, ierr
  integer, intent(in) :: nb_failed_loc
  logical :: failed
  integer :: results

  ! Sums all number of fails
  call mpi_reduce(nb_failed_loc, results, 1, mpi_integer, mpi_sum, 0, &
    mpi_comm_world, ierr)
  if (rank == 0) then
    if (results > 0) then
      write (*,*) results, "fails"
    else
      write (*,*) "ok"
    end if
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Copies data from north to south
!-----------------------------------------------------------------------------
subroutine propagate_to_south_inside(part)
  type(Partition), pointer :: part
  integer :: i_start, i_end, j_start, j_end
  integer :: i, j, k, nb_lat
  double precision :: val
  integer :: k_start, k_end

  nb_lat = get_true_nb_lat_Partition(part)
  if (has_south_ghost_cells(part%grid)) nb_lat = nb_lat-1
  ! Do not update south ghost cells
  do i = 1, nb_lat - 1
    ! For the first line, read all cells (ghost cells or first lat)
    if (is_north_ghost(i, part%grid)) then
      j_start = 1
      j_end = get_true_nb_lon_Partition(i, part)
    else
      call interior_lon_indices(j_start, j_end, part%grid, i)
    end if
    do j = j_start, j_end
      val = get_cell_ratio(i, j, TRACER, part%grid)
      ! Get south neighbour
      call south_neighbour_cell_Partition(k_start, k_end, i, j, part)
      ! Copies data to neighbours
      do k = k_start, k_end
        call set_cell_ratio(val, i + 1, k, TRACER, part%grid)
      end do

    end do
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Copies data from south to north
!-----------------------------------------------------------------------------
subroutine propagate_to_north_inside(part)
  type(Partition), pointer :: part
  integer :: i_start, i_end, j_start, j_end
  integer :: i, j, k, nb_lat
  double precision :: val
  integer :: k_start, k_end

  call interior_lat_indices(i_start, i_end, part%grid)
  nb_lat = get_true_nb_lat_Partition(part)
  ! Do not update north ghost cells
  do i = nb_lat, i_start + 1, -1
    ! For the ghost cells, read all cells
    if (is_south_ghost(i, part%grid)) then
      j_start = 1
      j_end = get_true_nb_lon_Partition(i, part)
    else
      call interior_lon_indices(j_start, j_end, part%grid, i)
    end if
    do j = j_start, j_end
      ! Get value
      val = get_cell_ratio(i, j, TRACER, part%grid)
      ! Get south neighbour
      call north_neighbour_cell_Partition(k_start, k_end, i, j, part)
      ! Copies data to neighbours
      do k = k_start, k_end
        call set_cell_ratio(val, i - 1, k, TRACER, part%grid)
      end do
    end do
  end do
end subroutine



!-----------------------------------------------------------------------------
!> Copies data from west to east
!-----------------------------------------------------------------------------
subroutine propagate_to_east_inside(part)
  type(Partition), pointer :: part
  integer :: i_start, i_end, j_start, j_end
  integer :: i, j, nb_lon
  integer :: extrema, nb_ghosts
  double precision :: val

  call interior_lat_indices(i_start, i_end, part%grid)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, part%grid, i)
    nb_lon = part%grid%nb_lon(i)
    ! Copy from west ghost cells
    if (j_start > 1) then
      nb_ghosts = part%grid%nb_ghosts_west(i) 
      extrema = 1
      ! Special case : we have more than 1 west ghost cells, but we only need
      ! one
      if (nb_ghosts > 1) extrema = 1 + nb_ghosts - 1
      val = get_cell_ratio(i, extrema, TRACER, part%grid)
 
      call set_cell_ratio(val,  i, j_start, TRACER, part%grid)
    end if

    ! Update the interior
    do j = j_start + 1, j_end 
      val = get_cell_ratio(i, j - 1, TRACER, part%grid)
      call set_cell_ratio(val,  i, j, TRACER, part%grid)
    end do
  end do
  !write (*,*) "after", part%grid%ratio
end subroutine

!-----------------------------------------------------------------------------
!> Copies data from east to west
!-----------------------------------------------------------------------------
subroutine propagate_to_west_inside(part)
  type(Partition), pointer :: part
  integer :: i_start, i_end, j_start, j_end
  integer :: i, j, nb_lon
  integer :: nb_ghosts, extrema
  double precision :: val

  call interior_lat_indices(i_start, i_end, part%grid)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, part%grid, i)
    nb_lon = part%grid%nb_lon(i)
    ! Copy from west ghost cells
    if (j_end < nb_lon) then
      nb_ghosts = part%grid%nb_ghosts_east(i) 
      extrema = nb_lon
      ! Special case : we have more than 1 west ghost cells, but we only need
      ! one
      if (nb_ghosts > 1) extrema = nb_lon - nb_ghosts + 1
      val = get_cell_ratio(i, extrema, TRACER, part%grid)
      call set_cell_ratio(val,  i, j_end, TRACER, part%grid)
    end if

    ! Update the interior
    do j = j_end - 1, j_start, -1
      !write (*,*) "j", j
      val = get_cell_ratio(i, j + 1, TRACER, part%grid)
      call set_cell_ratio(val,  i, j, TRACER, part%grid)
    end do
  end do
  !if (part%i_pos > 1) write (*,*) "after", part%grid%ratio
end subroutine


!-----------------------------------------------------------------------------
!> Check all the data for a partition. Ratio must be set to 1 at the interior.
!> Ghost cells must be set at 0, except in the direction opposite to the
!> advection (given in input).
!-----------------------------------------------------------------------------
subroutine check_ratio(nb_fails, part, direction)
  type(Partition), pointer :: part
  integer, intent(inout) :: nb_fails
  double precision :: val
  integer :: i_start, i_end, j_start, j_end
  integer :: i, j, nb_lon
  integer :: nb_errors, nb_errors2
  integer :: direction
  logical :: has_neighbour

  call interior_lat_indices(i_start, i_end, part%grid)
  nb_errors = 0
  nb_errors2 = 0

  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, part%grid, i)

    ! All cells inside must be set to 1
    do j = j_start, j_end
      val = get_cell_ratio(i, j, TRACER, part%grid )

      has_neighbour = .True.

      ! Value is updated only if there is a neighbour
      if (val /= 1.d0) then
#ifdef VERBOSE
         write (*,*) "wrong interior cell data at", i,j, "for", part%id
#endif
        nb_errors = nb_errors + 1
      end if
    end do
  end do

  !nb_errors2 = check_ghost_cells_ratio(direction, part, i_start, i_end)

  if (nb_errors > 0 .or. nb_errors2 > 0) then
#ifdef VERBOSE
    write (*,*) "Partition", part%id, ":"
    write (*,*) "Nb errors at the interior", nb_errors
    write (*,*) "Nb errors at the ghost cells", nb_errors2
#endif
    nb_fails = nb_fails + nb_errors + nb_errors2
  end if
end subroutine

!-------------------------------------------------------------------------------
! Check the ghost cells are either to 0 or 1
!-------------------------------------------------------------------------------
function check_ghost_cells_ratio(direction, part, i_start, i_end) &
    result (nb_errors)
  type (Partition) :: part
  integer, intent(in) :: direction, i_start, i_end
  double precision :: val
  integer :: i, j_start, j_end, nb_lat
  integer :: j, nb_lon, nb_errors
  integer :: val_north, val_east, val_south, val_west

  val_north = 0
  val_east = 0
  val_south = 0
  val_west = 0
  ! Ghost cells must be 1 in the opposition propagation direction
  if (direction == 0) then
    val_south = 1
    if (.not. is_last_on_band(part)) val_east = 1
    if (.not. is_first_on_band(part)) val_west = 1
    ! East
  else if (direction == 1) then
    val_west = 1
    ! South
  else if (direction == 2) then
    val_north = 1
    if (.not. is_last_on_band(part)) val_east = 1
    if (.not. is_first_on_band(part)) val_west = 1
    ! West
  else if (direction == 3) then
    val_east = 1
  end if

  nb_errors = 0
  nb_lat = part%grid%nb_lat

  ! North ghost cells
  do i = 1, i_start - 1
    nb_lon = part%grid%nb_lon(i)
    do j = 1, nb_lon
      val = get_cell_ratio(i, j, TRACER, part%grid)
      if (val /= val_north) then
        write (*,*) "error in north ghost for ", part%id, "has", val, "instead of", val_north
        nb_errors = nb_errors + 1
      end if
    end do
  end do

  ! South ghost cells
  do i = i_end + 1, nb_lat
    nb_lon = part%grid%nb_lon(i)
    do j = 1, nb_lon
      val = get_cell_ratio(i, j, TRACER, part%grid)
      if (val /= val_south) then
        write (*,*) "error in south ghost for ", part%id, "at", i, j
        nb_errors = nb_errors + 1
      end if
    end do
  end do

  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, part%grid, i)
    nb_lon = part%grid%nb_lon(i)

    ! East ghost cells
    do j = j_end + 1, nb_lon
      val = get_cell_ratio(i, j, TRACER, part%grid)
      if (val /= val_east) then
#ifdef DEBUG
        write (*,*) "error in east ghost for ", part%id, " at ", i, j 
        print *, "val", val, "instead of ", val_east
#endif
        nb_errors = nb_errors + 1
      end if
    end do

    ! West ghost cells
    do j = 1, j_start - 1
      val = get_cell_ratio(i, j, TRACER, part%grid)
      if (val /= val_west) then
        nb_errors = nb_errors + 1
#ifdef DEBUG
        write (*,*) "error in west ghost for ", part%id, " at ", i, j 
        !write (*,*) "error ghost west at", i,j," for ", part%id
#endif
      end if
    end do
  end do

end function

!-----------------------------------------------------------------------------
!> Set the first cell on the western side of the grid to 1, the rest is 
!> supposed to be 0
!-----------------------------------------------------------------------------
subroutine set_data_west(nb)
  type(Partition), pointer :: part
  integer :: i_start, i_end, j_start, j_end
  integer :: k, nb, i
  double precision :: val

  do k = 1, nb
    part => list_parts%parts(k)
    ! Just to be sure
    call init_tracer_ratio(part%grid)

    if (part%j_pos == 1) then
      call interior_lat_indices(i_start, i_end, part%grid)

      do i = i_start, i_end
        call set_cell_ratio(1.d0,  i, 1, TRACER, part%grid)
      end do
    end if
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Set the last latitude line on the grid to 1, the rest to 0.
!-----------------------------------------------------------------------------
subroutine set_data_equator(nb)
  type(Partition), pointer :: part
  integer :: j, j_start, j_end
  integer :: k, nb, nb_lon, nb_bands
  integer :: i_lat
  double precision :: val

  do k = 1, nb
    part => list_parts%parts(k)

    ! Just to be sure
    call init_tracer_ratio(part%grid)

    i_lat = -1
    if (is_just_before_equator(part%i_pos)) then
      i_lat = part%grid%nb_lat
    end if

    if (is_just_after_equator(part%i_pos)) then
      i_lat = 1
    end if

    if (i_lat > -1) then
      ! set the ghost cells
      !write (*,*) "setting lat", i_lat, "lon", j_start, j_end, "for part", part%id
      do j = 1, part%grid%nb_lon(i_lat)
        call set_cell_ratio(1.d0,  i_lat, j, TRACER, part%grid)
      end do
    end if
  end do
end subroutine


!-----------------------------------------------------------------------------
!> Set the last cell on the eastern side of the grid to 1, the rest is 
!> supposed to be 0
!-----------------------------------------------------------------------------
subroutine set_data_east(nb)
  type(Partition), pointer :: part
  integer :: i_start, i_end, j_start, j_end
  integer :: k, nb, nb_lon, i
  double precision :: val

  do k = 1, nb
    part => list_parts%parts(k)

    ! Just to be sure
    call init_tracer_ratio(part%grid)
    if (is_last_on_band_all(part)) then
      call interior_lat_indices(i_start, i_end, part%grid)
      do i = i_start, i_end
        ! Set ghost cells
        nb_lon = part%grid%nb_lon(i)
        call set_cell_ratio(1.d0,  i, nb_lon, TRACER, part%grid)
      end do
    end if
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Zonal advection for tests from east to west
!-----------------------------------------------------------------------------
subroutine test_zonal_advection_east_to_west(part)
  type (Partition), pointer :: part
  integer :: direction, l
  integer :: ierr
  integer :: status(mpi_status_size)

  ! Receive from east
  ! Last partition on band must not receive (not periodic boundaries) 
  if (.not. is_last_on_band_all(part)) then
    l = 1
    call start_receive(part, 1, part%requests, l, IS_RATIO)
    call wait_for_all(part%requests)
  end if

  ! propagate it inside first if not receiving
  call propagate_to_west_inside(part)

  ! Do not send to west for first partition
  if (part%j_pos > 1) then
    l = 1
    call start_send(part, 3, part%requests, l, IS_RATIO)
    call wait_for_all(part%requests)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Zonal advection for tests from west to east
!-----------------------------------------------------------------------------
subroutine test_zonal_advection_west_to_east(part)
  type (Partition), pointer :: part
  integer :: direction, nb_parts
  integer :: ierr, l
  integer :: status(mpi_status_size)

  ! Receive from east
  ! First partition on band must not receive (not periodic boundaries)
  if (part%j_pos > 1) then
    l = 1
    call start_receive(part, 3, part%requests, l, IS_RATIO)
    call wait_for_all(part%requests)
  end if

  ! propagate it inside first if not receiving
  call propagate_to_east_inside(part)

  ! Do not send to east for last partition
  nb_parts = nb_parts_on_band_sector(part%i_pos)
  if (.not. is_last_on_band_all(part)) then
    l = 1
    call start_send(part, 1, part%requests, l, IS_RATIO)
    call wait_for_all(part%requests)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Meridional advection for tests from the Equator to south
!-----------------------------------------------------------------------------
subroutine test_meridional_advection_to_south(part)
  type (Partition), pointer :: part
  integer :: direction, l
  integer :: ierr
  integer :: status(mpi_status_size)

  ! First partition on band must not receive (not periodic boundaries)
  if (.not. is_just_after_equator(part%i_pos)) then
    l = 1
    call start_receive(part, 0, part%requests, l, IS_RATIO)
    call wait_for_all(part%requests)
  end if

  ! propagate it inside first if not receiving
  call propagate_to_south_inside(part)
!
  ! Do not send to south for last partition
  if (.not. is_southmost(part)) then
    l = 1
    call start_send(part, 2, part%requests, l, IS_RATIO)
    call wait_for_all(part%requests)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Meridional advection for tests from the Equator to north
!-----------------------------------------------------------------------------
subroutine test_meridional_advection_to_north(part)
  type (Partition), pointer :: part
  integer :: direction
  integer :: ierr, l
  integer :: status(mpi_status_size)

  if (.not. is_just_before_equator(part%i_pos)) then
    l = 1
    call start_receive(part, 2, part%requests, l, IS_RATIO)
    call wait_for_all(part%requests)
  end if

  ! propagate it inside first if not receiving
  call propagate_to_north_inside(part)

  ! Do not send to east for last partition
  if (.not. is_northmost(part)) then
    l = 1
    call start_send(part, 0, part%requests, l, IS_RATIO)
    call wait_for_all(part%requests)
  end if

end subroutine

end module
