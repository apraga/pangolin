!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Tests_borders
!
!> @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Performs all cells tests.
!
!-------------------------------------------------------------------------------

module Tests_cells
use Tests_output
use Pangolin_run
use Partitioning_class

contains

!-------------------------------------------------------------------------------
! Start all the tests
!-------------------------------------------------------------------------------
function run_cells_tests() result (nb_fails)
  integer :: nb_fails

  nb_fails = no_cell_isolated()
  nb_fails = compare_cells_neighbours() + nb_fails

end function

!-------------------------------------------------------------------------------
! Check all cells for all partitions have at least one neighbour
!-------------------------------------------------------------------------------
function no_cell_isolated() result(nb_fails)
  integer :: nb_parts, nb
  integer :: nb_fails
  integer :: i, j, k
  type(Partition), pointer :: cur

  nb_parts = get_nb_parts_Partitioning(list_parts) 
  nb_fails = 0
  write (*, '(a)', advance="no") "Checking no cell is isolated..."
  do k = 1, nb_parts
    cur => list_parts%parts(k)
    
    do i = 1, cur%grid%nb_lat
      do j = 1, cur%grid%nb_lon(i)
        nb = nb_neighbours_cell_Partition(cur, i, j)
        if (nb == 0) then
!          write (*,*) 
!          write (*,*) "Failed : one cell is alone at", i,j
          nb_fails = nb_fails + 1
        end if
      end do
    end do

  end do
  call check_nb_fails(nb_fails, nb_parts)
end function

!-------------------------------------------------------------------------------
!> Test the functions computing the neighbours. We compare the number of
!> neighbours to a numer computed numerically
!-------------------------------------------------------------------------------
function compare_cells_neighbours() result(nb_fails)
  type (Partition), pointer :: part
  type (Band_grid), pointer :: cur
  integer :: nb_fails
  integer :: j_start, j_end
  integer :: j_start2, j_end2
  integer :: i, j, k, nb_parts

  nb_fails = 0
  write (*, '(a)', advance="no") "Checking number of north-south neighbours for each cell..."
  nb_parts = get_nb_parts_Partitioning(list_parts) 
  do k = 1, nb_parts
    part => list_parts%parts(k)
    cur => part%grid

    do i = 1, cur%nb_lat
      call interior_lon_indices(j_start, j_end, cur, i)
      !j_start = 1
      !j_end = cur%nb_lon(i)

      do j = j_start, j_end
        ! North neighbours
        if (i > 1) then
          call north_neighbour_cell_Partition(j_start, j_end, i, j, part)
          call neighbours_north_anal(j_start2, j_end2, part, i, j)
          call compare_indices(j_start, j_start2, j_end, j_end2, i, j,&
            part, 0, nb_fails)
        end if

       ! South
       if (i < cur%nb_lat) then
         call south_neighbour_cell_Partition(j_start, j_end, i, j, part)
         call neighbours_south_anal(j_start2, j_end2, part, i, j)
         call compare_indices(j_start, j_start2, j_end, j_end2, i, j,&
            part, 2, nb_fails)
        end if
      end do
    end do
  end do

  call check_nb_fails(nb_fails, nb_parts)
end function

!-------------------------------------------------------------------------------
!> Compare numerical neighbours to analytical version
!-------------------------------------------------------------------------------
subroutine compare_indices(j_start, j_start2, j_end, j_end2, i, j, part, &
    direction, nb_fails)
    type (Partition) :: part
  integer, intent(in) :: j_start, j_start2, j_end, j_end2
  integer, intent(in) :: i, j
  integer :: nb_fails
  logical :: failed
  integer :: direction
  character(5) :: dir

  if (direction == 0) dir = "north"
  if (direction /= 0) dir = "south"

  if (j_start /= j_start2 .or. j_end /= j_end2) then
    ! If all 4 indices are negative (or 0) it is still OK, as we are
    ! outside anyway
    if (j_start > 0 .and. j_end > 0 .and. j_start2 > 0 .and. j_end2 > 0) then
if (i == 1 .and. j == 2) then
  print *, "id", part%id
  write (*,*) j_start, j_start2, j_end, j_end2, "dir", direction
end if
!print *, i, j
      nb_fails = nb_fails + 1
    end if
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Neighbour is on the other side of the equator
!-------------------------------------------------------------------------------
subroutine neighbours_equator(j_start, j_end, i, j, i_neighb, part)
  type (Partition), target :: part
  integer, intent(in) :: i, j, i_neighb
  integer :: j_start, j_end
  integer :: nb_ghosts

  j_start = j
  j_end = j

  ! For 6 or more partitions, there is a separation at the equator
  if (get_total_nb_partitions_Configuration() > 3) then
    nb_ghosts = part%grid%nb_ghosts_west(i)
    ! From ghosts, we must add instead of substracting
    if (i == 1 .or. i == part%grid%nb_lat) then
      nb_ghosts = -part%grid%nb_ghosts_west(i_neighb)
    end if

    j_start = j_start  - nb_ghosts
    j_end = j_end - nb_ghosts
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Analytical north neighbours
!> Compute the number of neighbours analytically for the cell (i,j) in the 
!> current grid. We cannot do it numerically due to float approximation.
!-------------------------------------------------------------------------------
subroutine neighbours_north_anal(j_start, j_end, part, i, j) 
  type (Partition), target :: part
  type (Band_grid), pointer :: grid
  integer, intent(in) :: i, j
  integer :: nb, j_start, j_end
  integer :: nb_cells, mid
  integer :: j1, j2, i_neighb
  integer :: nb_lat2, sector, j_glob
  double precision :: prec
  logical :: equator, is_north
  logical :: is_south, south_hemisph

  south_hemisph = on_southern_hemisphere(part%zone)
  ! Special case, the current latitude line is just after the equator.

  if (is_latitude_just_after_equator(i, part)) then
    call neighbours_equator(j_start, j_end, i, j, i-1,part)
    return
  end if
 
  prec = 1e-5
  grid => part%grid
  ! We switch to the global position of the cell in the latitude line (float
  ! computation)
  j_glob = global_cell_position_from_local(grid, i, j)


  !if (i == 1) print *, "j=",j,"glob before", j_glob
  ! Find the neighbours (much easier)

  nb_cells = nb_cells_from_dlon_sector(grid%dlon(i))
  ! We translate into the first sector
  sector = 1
  do while (j_glob > nb_cells)
    sector = sector + 1
    j_glob = j_glob - nb_cells
  end do
  !if (i == 1) print *, "j=",j,"glob after", j_glob

  mid = (nb_cells + 1)/2
  i_neighb = i - 1

  ! For a single partition, this case is not taken into account by
  ! left_north_cell
  nb_lat2 = get_nb_lat2_Global_grid()
  if (i > nb_lat2) then
    j_start = left_south_cell_northern(j_glob, mid)
    j_end = right_south_cell_northern(j_glob, mid)
  else
    j_start = left_north_cell(j_glob, mid, part%zone)
    j_end = right_north_cell(j_glob, mid, nb_cells, part%zone)
  end if

!  if (i == 1 .and. j == 2 .and. part%id == 1) then
!    print *, "global", j_start, j_end
!  end if

  call revert_to_local(j_start, j_end, i_neighb, sector, nb_cells, i, j, &
    part, grid)
   
end subroutine

!-------------------------------------------------------------------------------
!> Once the cell neighbours are computed globally, we switch back to local scope
!> Cells positions are on the first sector
!-------------------------------------------------------------------------------
subroutine revert_to_local(j_start, j_end, i_neighb, sector, nb_cells, i, j, &
    part, grid)
  type (Band_grid) :: grid
  type (Partition) :: part
  integer :: j_start, j_end, nb_cells
  integer, intent(in) :: i_neighb, sector
  integer, intent(in) :: i, j

  ! Here, we revert to true global on the first sector
  j_start = global_to_local_cell(j_start, i_neighb, sector, grid) 
  j_end = global_to_local_cell(j_end, i_neighb, sector, grid) 
  !if (i == 1 .and. j == 2 .and. part%id == 1) then
  !  print *, "local start end", j_start, j_end
  !end if

  ! If both are outside, return
  if (j_start < 1 .and. j_end < 1) return

  ! If the cells are outside the partition, set it to the min or max
  if (j_start < 1) j_start = 1
  if (j_end < 1) j_end = grid%nb_lon(i_neighb)

  ! Eventually correct (boundary conditions)
  call correct_indices_other(i_neighb, j_start, part%grid, part%i_pos, part%j_pos)
  call correct_indices_other(i_neighb, j_end, part%grid, part%i_pos, part%j_pos)

end subroutine

!-------------------------------------------------------------------------------
!> Analytical south neighbours
!> Compute the number of neighbours analytically for the cell (i,j) in the 
!> current grid. We cannot do it numerically due to float approximation.
!-------------------------------------------------------------------------------
subroutine neighbours_south_anal(j_start, j_end, part, i, j) 
  type (Partition), target :: part
  type (Band_grid), pointer :: grid
  integer, intent(in) :: i, j
  integer :: j_start, j_end
  integer :: nb_cells, mid
  integer :: i_neighb, nb_lat2 
  integer :: sector, j_glob
  double precision :: prec
  logical :: equator, is_north
  logical :: is_south, south_hemisph
  integer :: idref, iref, jref

  ! For debug
  idref = -11
  iref = 2
  jref = 2
  south_hemisph = on_southern_hemisphere(part%zone)
  ! Special case, the current latitude line is just before the equator.
  if (is_latitude_just_before_equator(i, part)) then
    call neighbours_equator(j_start, j_end, i, j, i+1, part)
    return
  end if
 
  prec = 1e-5
  grid => part%grid
  ! We switch to the global position of the cell in the latitude line (float
  ! computation). 
  j_glob = global_cell_position_from_local(grid, i, j)

  nb_cells = nb_cells_from_dlon_sector(grid%dlon(i))
 
  ! However we stay in the first sector
  sector = 1
  do while (j_glob > nb_cells)
    sector = sector + 1
    j_glob = j_glob - nb_cells
  end do

  nb_lat2 = get_nb_lat2_Global_grid()
 
  ! Find the neighbours (much easier)
  mid = (nb_cells + 1)/2
  i_neighb = i + 1
 
  ! For a single partition, this case is not taken into account by
  ! left_north_cell
  if (i > nb_lat2) then
    j_start = left_north_cell_northern(j_glob, mid)
    j_end = right_north_cell_northern(j_glob, mid, nb_cells)
  else
    j_start = left_south_cell(j_glob, mid, part%zone)
    j_end = right_south_cell(j_glob, mid, nb_cells, part%zone)
  end if

  ! Revert to correct sector
!  if (i == 1 .and. j == 2 .and. part%id == 1) then
!    print *, "global ", j_start, j_end
!  end if
  call revert_to_local(j_start, j_end, i_neighb, sector, nb_cells, i, j, &
    part, grid)
!   if (i == 1 .and. j == 2) then
!    print *, "global start end", j_start, j_end
!  end if


end subroutine


!-------------------------------------------------------------------------------
!> Compute the number of neighbours analytically for the cell (i,j) in the 
!> current grid. We cannot do it numerically due to float approximation.
!-------------------------------------------------------------------------------
function nb_neighbours_north_south_anal2(this, i, j, i_neighb) result (nb)
  type (Partition) :: this
  integer, intent(in) :: i, j, i_neighb
  integer ::  nb, slope, mid
  double precision :: lon_west, lon_east, prec

  mid = middle_cell(i, this%grid, this%zone)
  !if (i == 3 .and. j == 46) then
    !write (*,*) "nb anal", nb, i, i_neighb
  !end if

  ! Single precision
  prec = 1e-5
  ! Property of the grid : 2 neighbours except at the center, on the borders on
  ! the grid (for north), and on the borders of the sector (for both, but not at
  ! the same time)
  if (j == mid) then
    ! North neighbour
    if (i_neighb == i -1) then
      nb = 1
    ! South neighbour
    else if (i_neighb == i +1) then
      nb = 3
      ! If the middle is on the extremity of the partition, only 1
      if (j == 1) nb = 1
    else 
      call print_error("Wrong direction", "nb_neighbours_north_south_anal2", "tests_cells.f90")
    end if
  else
    ! 2 by default
    nb = 2
    slope = pseudo_slope_Partition(this)

    ! North special cases
    if (i_neighb == i -1) then
      ! With a "slope" toward the west, there can only be 1 north neighbour at the
      ! west extremity of the partition
      if (j == 1 .and. slope > 0) then
        nb = 1
        return
      end if
      
      ! At the east extremity, there should be 2, but on a triangular
      ! partition, only 1
      if (j == this%grid%nb_lon(i) .and. this%j_pos == this%i_pos) then
        nb = 1
        return
      end if
 
      ! On the border of the sector, only one north neighbour
      lon_west = cell_west_lon(i, j, this%grid)
      lon_east = cell_east_lon(i, j, this%grid)
      if (abs(lon_west) < prec) nb = 1
      if (abs(lon_east - 120.) < prec) nb = 1
    end if

    ! South special case : right extremity of the partition
    ! With a "slope" toward the east, there can only be 1 south neighbour at the
    ! east extremity of the partition
    if (i_neighb == i + 1) then
      if (j == this%grid%nb_lon(i) .and. slope < 0 ) nb = 1
    end if
       
  end if

end function


end module
