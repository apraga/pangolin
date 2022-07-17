!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Analytical partitioning.
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Contains all functions and subroutines used by the analytical
!> partitioning.
!
!-------------------------------------------------------------------------------
module Analytical_partitioning
use Configuration_class
use Global_grid_class
implicit none

contains 

!-----------------------------------------------------------------------------
!> Return the number of partitions in the band i, whatever the band.
!> Works for all sectors
!> @param zone : zone indice
!> @param i : band number
!-----------------------------------------------------------------------------
function nb_parts_on_band(i) result(nb_parts)
  integer, intent(in) :: i
  integer :: nb_parts

  nb_parts = 3*nb_parts_on_band_sector(i)
end function

!-----------------------------------------------------------------------------
!> Return the number of partitions in the band i, whatever the band.
!> Works for both hemispheres, but on the first sector only !
!> @param zone : zone indice
!> @param i : band number
!-----------------------------------------------------------------------------
function nb_parts_on_band_sector(i) result(nb_parts)
  integer :: nb_parts, nb_partitions
  integer, intent(in) :: i
  integer :: i2
  integer :: nb_bands, nb_bands_square2
  integer :: nb_bands1, nb_bands_square1

  ! No parts on negative latitude bands
#ifdef DEBUG
  if (i < 1) then
    call print_error("Negative band number", "nb_parts_on_band_sector", &
      "analytical_partitioning.F90")
  end if
#endif

  nb_partitions = get_nb_partitions_Configuration(1)
  nb_bands1 = get_nb_bands_Configuration(1)

  ! On the northern hemisphere
  if (i < nb_bands1 + 1) then
    nb_bands_square1 = get_nb_bands_square_Configuration(1)
    ! If it is the leftover band
    if (i > nb_bands_square1) then
      nb_parts = nb_partitions - sum_partitions_sector(i-1)
    else
      nb_parts = nb_square_parts_on_band(i)
    end if
    ! On the southern hemisphere
  else
    nb_bands_square2 = get_nb_bands_square_Configuration(4)
    nb_bands = get_total_nb_bands_Configuration_sector()
    nb_partitions = get_nb_partitions_Configuration(4)
    ! Indice from the south pole
    i2 = nb_bands - i + 1
    ! If it is the leftover band
    if (i2 > nb_bands_square2) then
      nb_parts = nb_partitions - sum_partitions_sector(i2 - 1)
    else
      nb_parts = nb_square_parts_on_band(i)
    end if
  end if

  if (nb_parts < 0) then
    call print_error("local nb parts < 0", "nb_parts_on_band_sector", &
      "analytical_partitioning.f90")
  end if
end function 

!-----------------------------------------------------------------------------
!> Check if our partitioning has a "last band", meaning a leftover of partititions
!-----------------------------------------------------------------------------
function has_leftover_band(zone) result(res)
  logical :: res
  integer, intent(in) :: zone
  integer :: nb_bands, nb_bands_square

  nb_bands_square = get_nb_bands_square_Configuration(zone)
  nb_bands = get_nb_bands_Configuration(zone)
  res = (nb_bands > nb_bands_square)

end function

!-----------------------------------------------------------------------------
!> Returns the sum of all partitions up to the ith line on the adequate sector. 
!> For that, we need to compute the current sector
!> @param i : band number
!-----------------------------------------------------------------------------
function sum_partitions(i) result(n)
  integer, intent(in) :: i
  integer :: sector,  n

  n = 3*sum_partitions_sector(i)
end function

!-----------------------------------------------------------------------------
!> Returns the sum of all partitions up to the ith line on the first sector. 
!> On the south hemisphere, it sums from the south pole.
!> Work also for the leftover band.
!> @param i : band number
!-----------------------------------------------------------------------------
function sum_partitions_sector(i) result(n)
  integer, intent(in) :: i
  integer :: n, i2, nb_bands1
  integer :: nb_square1, nb_square2

  nb_square1 = get_nb_bands_square_Configuration(1)
  nb_bands1 = get_nb_bands_Configuration(1)
  ! Inside the square partitioning (northern hemisphere)
  if (i < nb_square1 + 1) then
    n = i**2
    ! On the leftover band (northern hemisphere)
  else if (i == nb_bands1) then
    n = get_nb_partitions_Configuration(1)
    ! On the leftover band (southern hemisphere)
  else 
    ! On the leftover band (southern hemisphere)
    if (has_leftover_band(4) .and. (i == nb_bands1 + 1)) then
      n = get_total_nb_partitions_sector_Configuration()
      nb_square2 = get_nb_bands_square_Configuration(4)
      n = n - nb_square2**2
      ! Inside the square partitioning (southern hemisphere)
    else
      ! Previous latitude band
      i2 = get_total_nb_bands_Configuration_sector() - i
      n = get_total_nb_partitions_sector_Configuration()
      ! Sum is a square, but starting from south pole
      n = n - i2**2
    end if
  end if
end function
!-----------------------------------------------------------------------------
!> Check if the partition is on the leftover band in our partitioning. Different
!> from the last band.
!-----------------------------------------------------------------------------
function is_leftover_band(i_pos) result(res)
  integer, intent(in) :: i_pos
  logical :: res
  integer :: nb_bands1

  ! Find ourselves the hemisphere
  res = .False.
  nb_bands1 = get_nb_bands_Configuration(1)
  if (i_pos < nb_bands1 + 1) then
    nb_bands1 = get_nb_bands_square_Configuration(1)
    if (has_leftover_band(1)) res = (i_pos == nb_bands1 + 1)
  else
    if (has_leftover_band(4)) res = (i_pos == (nb_bands1 + 1))
  end if
end function 


!-----------------------------------------------------------------------------
!> Small function for finding the zone from the band latitude.
!-----------------------------------------------------------------------------
function zone_from_ipos(i) result (zone)
  integer, intent(in) :: i
  integer :: zone

  zone = 1
  if (i > get_nb_bands_Configuration(1)) zone = 2
end function

!-----------------------------------------------------------------------------
!> Returns true if the partition is on the middle of the current band (for
!> square partitions).
!> @param i : band number
!> @param j : partition number on the band (for any sector)
!> @param zone : partition zone
!-----------------------------------------------------------------------------
function middle_partition(i, j, zone) result (res)
  integer, intent(in) :: i, j, zone
  logical :: res
  integer :: nb_bands, i2
  integer :: sector, nb_parts

  ! Number of bands on northern hemisphere
  nb_bands = get_nb_bands_Configuration(1)
  i2 = i
  if (i > nb_bands) then
    ! Must change the variable on south hemisphere
    nb_bands = get_total_nb_bands_Configuration_sector()
    i2 = nb_bands - i + 1
  end if

  ! We must adapt i2 to the sector
  nb_parts = nb_parts_on_band_sector(i2)
  sector = sector_from_zone(zone)
  i2 = i2 + nb_parts*(sector-1)

  res = (j == i2)
end function

!-----------------------------------------------------------------------------
!> Translate the column position of the partition into the first sector if needs
!> to be
!> @parameter j :: columun position
!> @parameter nb_parts :: number of partition on the current band (first sector)
!-----------------------------------------------------------------------------
function translate_jpos(j, nb_parts) result(j_corr)
  integer, intent(in) :: j, nb_parts
  integer :: j_corr

  j_corr = j
  ! We translate j to the first sector (as it is our reference)
  do while (j_corr > nb_parts)
    j_corr = j_corr - nb_parts
  end do
end function 

!-------------------------------------------------------------------------------
!> Find number of cells on the south latitude line for the current 
!> partition.
!> @param nb_cells : number of cells on the last latitude (current partition)
!-------------------------------------------------------------------------------
function nb_cells_south_last(i_pos, zone, nb_cells) result (nb_cells_next)
  integer, intent(in) :: i_pos, zone, nb_cells
  integer :: nb_cells_next

  ! At the equator, same number of cells
  if (is_just_before_equator(i_pos)) then
    nb_cells_next = nb_cells
  else
    if (on_northern_hemisphere(zone)) then
      nb_cells_next = nb_cells + 2
    else
      nb_cells_next = nb_cells - 2
    end if
  end if
end function

!-------------------------------------------------------------------------------
!> Find number of cells on latitude north of the last latitude for the current 
!> partition. 
!> @param nb_cells : number of cells on the last latitude (current partition)
!-------------------------------------------------------------------------------
function nb_cells_north_last(i_pos, zone, nb_cells) result (nb_cells_next)
  integer, intent(in) :: i_pos, zone,nb_cells
  integer :: nb_cells_next

  ! At the equator, same number of cells
  if (is_just_after_equator(i_pos)) then
    nb_cells_next = nb_cells
  else
    if (on_northern_hemisphere(zone)) then
      nb_cells_next = nb_cells - 2
    else
      nb_cells_next = nb_cells + 2
    end if
  end if
end function

!-------------------------------------------------------------------------------
!> Find number of cells on the next north/south latitude from ghost cells
!> @param is_first : true if the latitude line is the first (otherwise, the
!last)
!> @param nb_cells : number of cells on the current latitude (current partition)
!-------------------------------------------------------------------------------
function nb_cells_next_lat_ghost(is_first, i_pos, zone, nb_cells) result (nb_cells_next)
  integer, intent(in) :: i_pos, zone,nb_cells
  logical, intent(in) :: is_first
  logical :: at_equator
  integer :: nb_cells_next

  ! Check if the ghost cells are just on the other side of the equator
  at_equator = is_first .and. is_just_after_equator(i_pos)
  at_equator = at_equator .or. (.not. is_first .and. &
    is_just_before_equator(i_pos))
  ! At the equator, same number of cells
  if (at_equator) then
    nb_cells_next = nb_cells
  else
    if (on_northern_hemisphere(zone)) then
      ! Default is towards the south
      nb_cells_next = nb_cells + 2
      if (.not. is_first) nb_cells_next = nb_cells - 2
    else
      nb_cells_next = nb_cells - 2
      if (.not. is_first) nb_cells_next = nb_cells + 2
    end if
  end if
end function


!-----------------------------------------------------------------------------
!> Returns if the partition is on the band just before the equator.
!-----------------------------------------------------------------------------
function is_just_before_equator(i_pos) result(res)
  integer, intent(in) :: i_pos
  logical :: res
  integer :: nb_bands

  nb_bands = get_nb_bands_Configuration(1)
  ! Must have more than 1 partition
  res = (i_pos == nb_bands .and. .not. has_single_partition())
  !if (this%id == 1) write (*,*) "nb bands 1", nb_bands
end function 

!-----------------------------------------------------------------------------
!> Returns if the partition is on the band just after the equator.
!-----------------------------------------------------------------------------
function is_just_after_equator(i_pos) result(res)
  logical :: res
  integer, intent(in) :: i_pos
  integer :: nb_bands, nb

  nb_bands = get_nb_bands_Configuration(1)
  nb = get_total_nb_bands_Configuration_sector()
  res = (i_pos == nb_bands + 1 .and. nb > 1)
end function 

!-----------------------------------------------------------------------------
!> Returns true if the partition is on the right of the current band.
!> @param i : band number
!> @param j : partition number on the band
!-----------------------------------------------------------------------------
function right_partition(i, j) result (res)
  integer, intent(in) :: i, j
  logical :: res
  integer :: nb_bands, i2

  ! Number of bands on northern hemisphere
  nb_bands = get_nb_bands_Configuration(1)
  i2 = i
  if (i > nb_bands) then
    ! Must change the variable on south hemisphere
    nb_bands = get_total_nb_bands_Configuration_sector()
    i2 = nb_bands - i + 1
  end if
  res = (j > i2)
end function

!-----------------------------------------------------------------------------
!> Return the number of partitions in the band i on a sector. Works only for the
!> partitions in [1,N^2].
!> Works also for the southern hemisphere.
!-----------------------------------------------------------------------------
function nb_square_parts_on_band(i) result(nb_parts)
  integer :: nb_parts, nb_bands1
  integer, intent(in) :: i
  integer :: i2

  nb_bands1 = get_nb_bands_Configuration(1)
  i2 = i
  ! On southern hemisphere, change the variable
  if (i > nb_bands1) then
    i2 = get_total_nb_bands_Configuration_sector() - i + 1
  end if
  nb_parts = 2*i2 - 1
end function 

!-----------------------------------------------------------------------------
!> Find the sector from the zone
!-----------------------------------------------------------------------------
function sector_from_zone(zone) result(offset)
  integer, intent(in) :: zone
  integer :: offset

  ! Returns 0 for sector 1 and 4
  if (zone < 4) then
    offset = zone
  else
    offset = zone - 3
  end if
end function

!-----------------------------------------------------------------------------
!> Returns true if the current zone is on the northern hemisphere.
!-----------------------------------------------------------------------------
function on_northern_hemisphere(zone) result (res)
  integer, intent(in) :: zone
  logical :: res

  res = (zone < 4)
end function

!-----------------------------------------------------------------------------
!> Returns true if the current zone is on the south hemisphere.
!-----------------------------------------------------------------------------
function on_southern_hemisphere(zone) result (res)
  integer, intent(in) :: zone
  logical :: res

  res = (.not. on_northern_hemisphere(zone))
end function

!-----------------------------------------------------------------------------
!> Check if the partition is on the last band in our partitioning (last band
!> here means either the last band in the square partitioning or the last band
!> with a rest
!> @param zone : current partition zone
!> @param i_pos : partition band
!-----------------------------------------------------------------------------
function is_on_last_band(i_pos, zone) result(res)
  integer, intent(in) :: i_pos, zone
  logical :: res
  integer :: nb_bands, nb_bands_square
  integer :: i_pos2

  nb_bands_square = get_nb_bands_square_Configuration(zone)
  nb_bands = get_nb_bands_Configuration(zone)

  i_pos2 = i_pos
  ! On southern hemisphere, we start from the south pole
  if (on_southern_hemisphere(zone)) then
    i_pos2 = get_total_nb_bands_Configuration_sector() - i_pos + 1
  end if

  if (has_leftover_band(zone)) then
    res = (i_pos2 == nb_bands)
  else
    res = (i_pos2 == nb_bands_square)
  end if
end function 

!-------------------------------------------------------------------------------
!> Wrapper for left_south_cell. On the south hemisphere, south and north
!> relations are inverted.
!> !! Warning !! it determines the hemisphere based on the zone (does not work for a
!> single partition, so better call the underlying function)
!> !! Warning !! does not work at the equator
!> @param mid : the middle cell on the current line
!> @param cur : current cell (in global coordinates)
!-------------------------------------------------------------------------------
function left_south_cell(cur, mid, zone) result (j)
  integer, intent(in) :: cur, mid, zone
  integer :: j

  if (on_northern_hemisphere(zone)) then
    j = left_south_cell_northern(cur, mid)
  else
    j = left_north_cell_northern(cur, mid)
  end if
end function

!-------------------------------------------------------------------------------
!> Wrapper for left_north_cell. On the south hemisphere, south and north
!> relations are inverted.
!> !! Warning !! it determines the hemisphere based on the zone (does not work for a
!> single partition, so better call the underlying function)
!> !! Warning !! does not work at the equator
!> @param mid : the middle cell on the current part
!-------------------------------------------------------------------------------
function left_north_cell(cur, mid, zone) result (j)
  integer, intent(in) :: cur, mid, zone
  integer :: j

  if (on_northern_hemisphere(zone)) then
    j = left_north_cell_northern(cur, mid)
  else
    j = left_south_cell_northern(cur, mid)
  end if
end function

!-------------------------------------------------------------------------------
!> Wrapper for right_south_cell_northern. On the south hemisphere, south and north
!> relations are inverted.
!> !! Warning !! it determines the hemisphere based on the zone (does not work for a
!> single partition, so better call the underlying function)
!> !! Warning !! does not work at the equator
!> @param mid : the middle on the current latitude line (single sector)
!> @param last : the number of cells on the current latitude line
!-------------------------------------------------------------------------------
function right_south_cell(cur, mid, last, zone) result (j)
  integer, intent(in) :: cur, mid, last, zone
  integer :: j

  if (on_northern_hemisphere(zone)) then
    j = right_south_cell_northern(cur, mid)
  else
    j = right_north_cell_northern(cur, mid, last)
  end if
end function

!-------------------------------------------------------------------------------
!> Wrapper for right_north_cell_northern. On the south hemisphere, south and north
!> relations are inverted.
!> !! Warning !! it determines the hemisphere based on the zone (does not work for a
!> single partition, so better call the underlying function)
!> !! Warning !! does not work at the equator
!> @param mid : the middle on the current part
!> @param last : the number of cells on the current latitude band
!-------------------------------------------------------------------------------
function right_north_cell(cur, mid, last, zone) result (j)
  integer, intent(in) :: cur, mid, last, zone
  integer :: j

  if (on_northern_hemisphere(zone)) then
    j = right_north_cell_northern(cur, mid, last)
  else
    j = right_south_cell_northern(cur, mid)
  end if
end function


!-----------------------------------------------------------------------------
!> Fast computation of the rightmost south neighbour of the cell number cur
!> Warning : only works on northern hemisphere !
!> Returns only the position on the next lat
!-----------------------------------------------------------------------------
function right_south_cell_northern(cur, mid) result (j)
  integer, intent(in) :: cur, mid
  integer :: j
  if (cur < mid ) then
    j = cur + 1
  else 
    j = cur + 2
  end if

end function

!-----------------------------------------------------------------------------
!> Fast computation of the leftmost south neighbour of the cell number cur
!> Returns only the position on the next lat
!> Warning : only works on northern hemisphere !
!> @param mid : the middle of the current latitude line.
!-----------------------------------------------------------------------------
function left_south_cell_northern(cur, mid) result (j)
  integer, intent(in) :: cur, mid
  integer :: j
  if (cur < mid + 1) then
    j = cur
  else 
    j = cur + 1
  end if
end function

!-----------------------------------------------------------------------------
!> Fast computation of the leftmost north neighbour of the cell number cur
!> Mid is the middle on the current part
!> Warning : only works on northern hemisphere !
!> Returns only the position on the previous lat
!-----------------------------------------------------------------------------
function left_north_cell_northern(cur, mid) result (j)
  integer, intent(in) :: cur, mid
  integer :: j
  if (cur == 1) then
    j = cur
  else if (cur < mid + 1) then
    j = cur - 1
  else 
    j = cur - 2
  end if
end function

!-----------------------------------------------------------------------------
!> Fast computation of the rightmost north neighbour of the cell number cur
!> Returns only the position on the previous lat
!> Warning : only works on northern hemisphere !
!> @param mid : the middle on the current part
!> @param last : the number of cells on the current latitude band
!-----------------------------------------------------------------------------
function right_north_cell_northern(cur, mid, last) result (j)
  integer, intent(in) :: cur, mid, last
  integer :: j
  if (cur < mid) then
    j = cur
  else if (cur < last) then
    j = cur - 1
  else 
    j = cur - 2
  end if

end function

end module
