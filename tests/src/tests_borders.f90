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
!
!-------------------------------------------------------------------------------

module Tests_borders
use Pangolin_run
use Partitioning_class
use Band_grid_class
use Tests_output

! Show the details of borders ghosts
!logical :: show_all = .true. !.false.
logical :: show_all = .false.
contains

!-------------------------------------------------------------------------------
! Start all the tests
!-------------------------------------------------------------------------------
function run_borders_tests() result(nb_fails)
  integer :: nb_fails

  nb_fails = total_nb_ghost_cells()
  nb_fails = size_ghost_cells() + nb_fails
  nb_fails = extremities_ghost_cells() + nb_fails
end function

!-------------------------------------------------------------------------------
! Compute the total number of ghost cells for each partition's neighbours and 
! check nothing has been forgotten
!-------------------------------------------------------------------------------
function total_nb_ghost_cells() result(nb_fails)
  type(Partition), pointer :: cur
  integer :: nb_parts, i
  integer :: nb_fails


  nb_fails = 0
  write (*, '(a)', advance="no") "Checking total number of ghost cells ..."
  nb_parts = get_nb_parts_Partitioning(list_parts)
  do i = 1,nb_parts
    cur => list_parts%parts(i)
    ! North
    call compare_nb_ghosts(cur, 1, nb_fails)
    ! South
    call compare_nb_ghosts(cur, 0, nb_fails)
    ! No need for east-west
  end do
  call check_nb_fails(nb_fails, nb_parts)

end function

!-------------------------------------------------------------------------------
! Compare the total number of ghost cells in a direction given by a numerical 
! value
!-------------------------------------------------------------------------------
subroutine compare_nb_ghosts(cur, direction, nb_fails)
  type(Partition), pointer:: cur
  integer :: i_neighb, j_neighb, j_neighb2
  integer :: nb_fails
  integer :: nb_ref, direction
  logical :: has_neighb

  ! Here we find all the neighbour's position   
  if (direction == 1) then
    call north_partition_wrapper(cur, i_neighb, j_neighb, j_neighb2)
    has_neighb = i_neighb > 0 .and. (j_neighb > 0 .or. j_neighb2 > 0) 
    !write (*,*) "north neighb", j_neighb, j_neighb2

  else if (direction == 0) then
    call south_partition_wrapper(cur, i_neighb, j_neighb, j_neighb2)
    has_neighb = i_neighb > 0 .and. (j_neighb > 0 .or. j_neighb2 > 0) 
    !write (*,*) "south neighb", j_neighb, j_neighb2

  else 
    write (*,*) "Wrong direction"
    return
  end if

  ! Here, we get the number of cells in the adequate direction (numerically),
  ! which serves as a reference
  if (has_neighb) then
    !write (*,*) "direction", direction
    nb_ref = nb_ghosts_num_north_south(cur, direction)

    ! And we compare it to the function value
    call compare_nb(cur, nb_ref, i_neighb, j_neighb, j_neighb2,  direction, &
      nb_fails)
  else
    if (show_all) write (*,*) "no ", direction, " neighbour."
  end if
end subroutine

!-------------------------------------------------------------------------------
! Compare the total number of ghost cells to the analytical version
!-------------------------------------------------------------------------------
subroutine compare_nb(cur, nb_ref, i_neighb, j_neighb, j_neighb2, &
    direction, nb_fails)
  type(Partition), pointer:: cur
  integer :: nb_ghosts, nb_neighbours
  integer :: i_neighb, i_length
  integer :: j_neighb, j_neighb2
  integer :: nb_ref, j, direction
  integer :: nb_fails
  logical :: failed

  nb_ghosts = 0
  ! Sum the ghost cells
  !if (cur%id == 16) write (*,*) "j ", j_neighb, j_neighb2
  do j = j_neighb, j_neighb2
    !call interface_indices(cur, i_start, i_length, i_neighb, j)
    i_length = interface_length(cur, i_neighb, j)
    nb_ghosts = nb_ghosts + i_length
    if (show_all) write (*,*) "nb_ghosts loc", i_length
  end do

  ! Compare it to a numerical value
  if (nb_ghosts /= nb_ref) then
    failed = .True.
    if (direction == 1) write (*,*) "Failed : wrong north ghost cells."
    if (direction == 0) write (*,*) "Failed : wrong south ghost cells."
    nb_fails = nb_fails + 1
    write (*,*) "for", cur%id, "at", cur%i_pos, cur%j_pos
    write (*,*) "code gives", nb_ghosts, "instead of ", nb_ref
  end if
end subroutine

!-------------------------------------------------------------------------------
! Compute the total number of northern or southern neighbours by searching all
! the cells in the previous latitude line
! According to direction, we choose the north neighbour (1) or the south (0)
!-------------------------------------------------------------------------------
function nb_ghosts_num_north_south(cur, direction) result(nb)
  type (Partition), pointer :: cur
  integer :: direction
  integer :: nb_bands_max, part_line
  logical :: check, nb
  integer :: k_start, k_end, lat_cur
  double precision :: lon_min, lon_max
  integer :: j_start, j_end

  nb = 0

  ! Check we can have a neighbour 
  ! North
  if (direction == 1) then 
    check = (cur%i_pos > 1)
  ! South
  else
    nb_bands_max = get_total_nb_bands_Configuration() + 1
    check = (cur%i_pos < nb_bands_max)
  end if
  if (.not. check) return

  ! Indices of the current partition
  call interior_lat_indices(k_start, k_end, cur%grid)
  ! North latitude line
  if (direction == 1) then
    lat_cur = k_start 
    part_line = cur%i_pos - 1
    ! South
  else
    lat_cur = k_end
    part_line = cur%i_pos + 1
  end if
  ! Get the extrema of the current partition
  call interior_lon_indices(j_start, j_end, cur%grid, lat_cur)
  lon_min = cell_west_lon(lat_cur, j_start, cur%grid)
  lon_max = cell_east_lon(lat_cur, j_end, cur%grid)

  ! We examine all the partitions on the previous line in the partitioning
  nb = nb_neighbours_north_south(cur, lon_min, lon_max, part_line, direction)
   
end function

!-------------------------------------------------------------------------------
! Find the number of cells of the neighbour which intersect the current
! partition (only for north-south)
!-------------------------------------------------------------------------------
function nb_neighbours_north_south(cur, lon_min, lon_max, part_line, &
    direction) result(nb)
  integer :: j, nb_parts, lat_next
  integer, intent(in) :: part_line, direction
  double precision, intent(in) :: lon_min, lon_max
  type (Partition), pointer :: cur, neighb
  integer :: nb, id, k
  logical :: is_inside
  double precision :: lon_west, lon_east
  integer :: nb_tmp, k_start, k_end
  integer :: k_start_n, k_end_n
  double precision :: prec

  ! Precision
  prec = 1e-4
  nb_parts = nb_parts_on_band(part_line)
  !write (*,*) "nb parts", nb_parts

  nb = 0
  do j = 1, nb_parts
    ! part_line is the line for the neighbour
    id = compute_id_from_pos(part_line, j)
    neighb => list_parts%parts(id+1)

    ! Find the neighbour latitude line
    call interior_lat_indices(k_start_n, k_end_n, neighb%grid)
    ! North
    if (direction == 1) then
      lat_next = k_end_n
    else
      lat_next = k_start_n
    end if

    nb_tmp = 0
    call interior_lon_indices(k_start, k_end, neighb%grid, lat_next)

    ! Ignore ghost cells
    do k = k_start, k_end
      lon_west = cell_west_lon(lat_next, k, neighb%grid)
      lon_east = cell_east_lon(lat_next, k, neighb%grid)
      ! Condition : the current cell must intersect the partition interface
      ! i.e west boundary must be < east border of the partition and
      ! east boundary must be > west border of the partition and
      is_inside = lon_west < lon_max - prec .and. lon_east > lon_min + prec
      if (is_inside) then
        nb = nb + 1
        nb_tmp = nb_tmp + 1
      end if
      ! Exit if we are over the last
      if (lon_west > lon_max + prec) then
        if (show_all) write (*,*) "nb_tmp", nb_tmp
        return
      end if
    end do
    if (show_all) write (*,*) "nb_tmp", nb_tmp

  end do
end function

!-------------------------------------------------------------------------------
!> Compare the size of the ghost cells to the interior cells for a neighbours.
!> They should be equals.
!-------------------------------------------------------------------------------
function size_ghost_cells() result(nb_fails)
  type(Partition), pointer :: cur
  integer :: i
  integer :: nb_fails, nb_parts
  logical :: failed
 
  write (*, '(a)', advance="no") "Checking the number of ghost cells for"//&
    "each neighbour ..."
  nb_fails = 0

  nb_parts = get_nb_parts_Partitioning(list_parts) 
  do i = 1, nb_parts
    cur => list_parts%parts(i)

    call size_ghost_cells_north(cur, nb_fails)
    call size_ghost_cells_south(cur, nb_fails)
  end do

  call check_nb_fails(nb_fails, nb_parts)
end function

!------------------------------------------------------------------------------
!> Compare the size of the north ghost cells to the south interior cells 
!> for a neighbours.
!------------------------------------------------------------------------------
subroutine size_ghost_cells_north(cur, nb_fails)
  type (Partition), pointer :: cur, neighb
  integer :: k, k_start, k_end
  integer :: length, length2
  integer :: i_start, i_start_loc
  integer :: nb_fails

  call north_neighbours(cur, k_start, k_end)
  if (k_start > -1) then
    ! For each neighbour, we compute the north ghost cells interface length
    ! and compare it to the south interior cells of the neighbours
    do k = k_start, k_end
      neighb => list_parts%parts(k+1)

      call interface_north_partition(i_start, i_start_loc, length, cur, &
        neighb%j_pos)
      !length2 = interface_south_length_interior_cells(neighb, cur%j_pos)
      call interface_interior_south_partition(i_start, length2, neighb, &
        cur%j_pos)

      if (length /= length2) then
      !  write (*,*) "Failed : number of north ghost cells different from the ",&
      !    "interior :", length, "!=", length2
      !  write (*,*) "cur ", cur%id, "neighb", k
        nb_fails = nb_fails + 1
      end if
    end do
  end if
end subroutine

!------------------------------------------------------------------------------
!> Same as size_ghost_cells_north but for south
!------------------------------------------------------------------------------
subroutine size_ghost_cells_south(cur, nb_fails)
  type (Partition), pointer :: cur, neighb
  integer :: k_start, k_end, k
  integer :: length, length2
  integer :: i_start, i_start_loc
  integer :: nb_fails

  call south_neighbours(cur, k_start, k_end)
  if (k_start > -1) then
    do k = k_start, k_end
      neighb => list_parts%parts(k+1)

      call interface_south_partition(i_start, i_start_loc, length, cur, &
        neighb%j_pos)
      !length2 = interface_north_length_interior_cells(neighb, cur%j_pos)
      call interface_interior_north_partition(i_start, length2, neighb, &
        cur%j_pos)
      !write (*,*) "length", length, length2
      if (length /= length2) then
!        write (*,*) "Failed : number of south ghost cells different from the ",&
!          "interior :", length, "!=", length2
!        write (*,*) "cur ", cur%id, "neighb", k
        nb_fails = nb_fails + 1
      end if
    end do
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Compute the north and south extremities for the nourth/south ghost cells. 
!> Except for the first and last cell on the line, the west extremities must be 
!> the min and the east, the max
!> Must remove ghost cells 
!-------------------------------------------------------------------------------
function extremities_ghost_cells() result(nb_fails)
  type(Partition), pointer :: cur
  integer :: nb_parts, i
  integer :: nb_fails

  nb_fails  = 0
  write (*, '(a)', advance="no") "Checking N/S ghost cells extremities..."
  nb_parts = get_nb_parts_Partitioning(list_parts) 
  do i = 1,nb_parts
    cur => list_parts%parts(i)

    call check_extr_ghost_cells(cur, nb_fails)
  end do

  call check_nb_fails(nb_fails, nb_parts)
end function
 
!-------------------------------------------------------------------------------
!> Here we compare truly the extremities
!-------------------------------------------------------------------------------
subroutine check_extr_ghost_cells(cur, nb_fails)
  type(Partition), pointer :: cur
  integer :: nb_fails, nb_lat
  logical :: cond, at_equator

  if (has_north_ghost_cells(cur%grid)) then
    at_equator = is_just_after_equator(cur%i_pos)
    call compare_extr_ghost(cur, 1, 2, nb_fails, at_equator)
  end if

  if (has_south_ghost_cells(cur%grid)) then
    nb_lat = cur%grid%nb_lat 
    at_equator = is_just_before_equator(cur%i_pos)
    call compare_extr_ghost(cur, nb_lat, nb_lat-1, nb_fails, at_equator) 
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Compare the extremities for ghost cells at line i_ghost (inside the grid) in
!> respect to i_lat
!> @param at_equator : true if the ghost cells is on the other side of the
!> equator
!-------------------------------------------------------------------------------
subroutine compare_extr_ghost(cur, i_ghost, i_lat, nb_fails, at_equator) 
  type(Partition), pointer :: cur
  integer, intent(in) :: i_ghost, i_lat
  logical, intent(in) :: at_equator
  integer :: i_pos, j_pos, j_corr
  integer :: nb_fails, nb_lon
  integer :: j_start, j_end, nb_parts
  double precision :: lon_ghost, lon
  logical :: fail, special, tmp
  double precision :: prec

  ! Precision 
  prec = 0.01

  i_pos = cur%i_pos
  j_pos = cur%j_pos
  nb_parts = nb_parts_on_band_sector(i_pos)
  j_corr = translate_jpos(j_pos, nb_parts)

  ! First cell
  lon_ghost = cell_west_lon(i_ghost, 1, cur%grid)

  ! Gives the first and last interior cells
  call interior_lon_indices(j_start, j_end, cur%grid, i_lat)
  ! First cell on the next lat
  lon = cell_west_lon(i_lat, j_start, cur%grid)

  ! Fail when ghost is before normal cells
  ! 3 cases : the first partition on the sector, or if the second cell is at the
  ! middle, or at the equator

  special = at_equator
  !special = special .or. is_first_part_sector(i_pos, j_pos)
  !special = special .or. (is_middle_float(i_lat, 2, cur%grid, cur%zone) &
  !.and. .not. middle_partition(i_pos, j_pos, cur%zone))
  if (special) then
    fail = lon > lon_ghost + prec
  ! Fail when ghost is after normal cells
  else
    fail = lon < lon_ghost - prec
  end if
  if (fail) then
!    write (*,*) "failed extremity (first cell) at lat", i_ghost
!    write (*,*) "for part at", cur%i_pos, cur%j_pos
  tmp =  middle_partition(i_pos, j_pos, cur%zone)
    nb_fails = nb_fails + 1
  end if

  ! Last cell
  nb_lon = get_true_nb_lon_Partition(i_ghost, cur)
  lon_ghost = cell_east_lon(i_ghost, nb_lon, cur%grid)
  nb_lon = get_true_nb_lon_Partition(i_lat, cur)
  ! Last cell not ghost
  lon = cell_east_lon(i_lat, j_end, cur%grid)

  ! Special casese : at the equator or on the last sector, or the cell is at the
  ! middle
  special = .False.
!  special = at_equator
!  special = special .or. is_last_part_sector(i_pos, j_pos)
!  special = special .or. (is_middle_float(i_lat, nb_lon-1, cur%grid, cur%zone) &
!    .and. .not. middle_partition(i_pos, j_pos, cur%zone))
  
  ! Fail when ghost is after normal cells
  if (special) then
    fail = lon < lon_ghost - prec
  ! Fail when ghost is before normal cells
  else
    fail = lon > lon_ghost + prec
  end if
  if (fail) then
    !write (*,*) "failed extremity (last cell) at lat", i_ghost
    !write (*,*) "for part at", cur%i_pos, cur%j_pos
    !!write (*,*) "at eq", at_equator,"test", is_just_before_equator(cur%i_pos)
    !write (*,*) "cur ghost", lon, lon_ghost
    !write (*,*) "is mid", is_middle_float(i_lat, nb_lon-1, cur%grid, cur%zone)
    !tmp =  middle_partition(i_pos, j_pos, cur%zone)
    !write (*,*) "is not mid", tmp
 
    nb_fails = nb_fails + 1
  end if


end subroutine

end module
