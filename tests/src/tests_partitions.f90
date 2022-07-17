!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Tests_partitions
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Performs all cells tests.
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Test the partitioning
! Must run with 1 proc only
!-------------------------------------------------------------------------------

module Tests_partitions
use Tests_output
use Pangolin_run
use Partitioning_class
use Band_grid_class

contains

!-------------------------------------------------------------------------------
! Start all the tests
!-------------------------------------------------------------------------------
function run_partitions_tests() result(nb_fails)
  integer :: nb_fails, nb_parts

  nb_parts = get_nb_parts_Partitioning(list_parts) 

  nb_fails = sphere_is_partitioned(nb_parts)
  nb_fails = no_partition_isolated(nb_parts) + nb_fails
  nb_fails = reflexivity(nb_parts) + nb_fails

end function

!-------------------------------------------------------------------------------
! Check all partitions have at least one neighbour 
!-------------------------------------------------------------------------------
function no_partition_isolated(nb_parts) result(nb_fails)
  type(Partition), pointer :: cur
  integer, intent(in) :: nb_parts
  integer :: nb, i
  integer :: nb_fails

  nb_fails = 0
  write (*,'(a)', advance="no") "Checking no partition is isolated..."
  if (nb_parts == 1) then
    write (*,*) "(skipped)"
    return
  end if

  do i = 1, nb_parts
    cur => list_parts%parts(i)
    nb = nb_neighbours_Partition(cur)
    if (nb == 0) then
      nb_fails = nb_fails + 1
    end if
  end do
  call check_nb_fails(nb_fails, nb_parts)
end function


!-------------------------------------------------------------------------------
! Check all neighbours of partitions have this partition as their neighbour (and
! in the righ direction !)
!-------------------------------------------------------------------------------
function reflexivity(nb_parts) result(nb_fails)
  type(Partition), pointer :: cur
  integer, intent(in) :: nb_parts
  integer :: nb_fails, i

  nb_fails = 0
  write (*,'(a)', advance="no") "Checking partition relations are reflexive..."
  if (nb_parts == 1) then
    write (*,*) "(skipped)"
    return
  end if

  do i = 1, nb_parts
    cur => list_parts%parts(i)

    call east_west_reflexivity(cur, i, nb_fails)
    call north_reflexivity(cur, i, nb_fails)
    call south_reflexivity(cur, i, nb_fails)

  end do
  call check_nb_fails(nb_fails, nb_parts)
end function

!-------------------------------------------------------------------------------
! Check that for a partition the eart-west neighbour relations is reflexive
!-------------------------------------------------------------------------------
subroutine east_west_reflexivity(cur, i_cur, nb_fails)
  type(Partition), pointer :: cur
  type(Partition), pointer :: prev
  integer :: neighb, test
  integer, intent(in) :: i_cur
  logical :: failed
  integer :: nb_fails

  ! Check east neighbour
  call east_neighbours(cur, neighb)
  if (neighb > -1) then
    ! Beware, partitions numbers start from 0
    prev => list_parts%parts(neighb+1)
    call west_neighbours(prev, test)
    ! Same here
    if (test /= i_cur-1) then
      !write (*,*) 
      !write (*,*) "Failed : east neighbour relation is not reflexive for",cur%id
      !  write (*,*) "at", cur%i_pos, cur%j_pos
      nb_fails = nb_fails + 1
    end if
  else
  end if

  ! Check west neighbour
  call west_neighbours(cur, neighb)
  !if (cur%id == 65) write (*,*) "neighb for 65", neighb
  if (neighb > -1) then
    prev => list_parts%parts(neighb+1)
    call east_neighbours(prev, test)
    if (test /= i_cur-1) then
      !write (*,*) 
      !write (*,*) "Failed : west neighbour relation is not reflexive for", &
      !  cur%id
      !write (*,*) "at", cur%i_pos, cur%j_pos
      failed = .True.
      nb_fails = nb_fails + 1
    end if
  end if
end subroutine

!-------------------------------------------------------------------------------
! Check that for a partition the north neighbour relations is reflexive
!-------------------------------------------------------------------------------
subroutine north_reflexivity(cur, i_cur, nb_fails)
  type(Partition) :: cur
  type(Partition), pointer :: prev
  integer :: neighb, neighb2, test, test2
  integer, intent(in) :: i_cur
  logical :: failed, failed_tmp
  integer :: nb_fails, i, j
  integer :: id_test = -18

  failed_tmp = .True.

  if (cur%id == id_test) then
    write (*,*) 
    write (*,*) "#######################################"
  end if

  ! Check north neighbours : there can be several neighbours
  call north_neighbours(cur, neighb, neighb2)
  if (cur%id == id_test) then
    write (*,*) "id", cur%id, "north neighb", neighb, neighb2
  end if
  ! This allow us to check all the cases : if one of them is -1, there are both
  ! equal  (either positive or -1)
  if (neighb < 0) neighb = neighb2
  if (neighb2 < 0) neighb2 = neighb

  ! If there is at least one north neighbour
  if (neighb > -1) then
    ! Search among all of them
    do i = neighb, neighb2 
      prev => list_parts%parts(i+1)
      call south_neighbours(prev, test, test2)
      if (cur%id == id_test) write (*,*) "trying south neighb for", i,"at", prev%i_pos, "result = ", test, test2
      ! Among all the south neighbours, there must be at least one which is the
      ! current partition

      if (test < 0) test = test2
      if (test2 < 0) test2 = test

      if (neighb > -1) then
        failed_tmp = .True.
        do j = test, test2
          ! i_cur > 0
          if (j == i_cur-1) then 
            failed_tmp = .False.
          end if
        end do
        if (failed_tmp) then
          !write (*,*) 
          !write (*,*) "Failed : north neighbour relation is not reflexive for"&
          !  , cur%id
          !write (*,*) "at", cur%i_pos, cur%j_pos
          !!, cur%i_pos,cur%j_pos
          !write (*,*) "Neighbours are", neighb, neighb2
          nb_fails = nb_fails + 1
        end if
      end if
    end do
  end if
end subroutine

!-------------------------------------------------------------------------------
! Check that for a partition the south neighbour relations is reflexive
!-------------------------------------------------------------------------------
subroutine south_reflexivity(cur, i_cur, nb_fails)
  type(Partition) :: cur
  type(Partition), pointer :: prev
  integer, intent(in) :: i_cur
  integer :: neighb, neighb2, test, test2
  integer :: nb_fails, i, j
  logical :: failed_tmp
  integer :: id_test = -20

  failed_tmp = .True.

  if (cur%id == id_test) then
    write (*,*) 
    write (*,*) "#######################################"
  end if

  ! Check south neighbours : there can be several neighbours
  call south_neighbours(cur, neighb, neighb2)
  if (cur%id == id_test) then
    write (*,*) "cur id", cur%id, "south neighb", neighb, neighb2
  end if
  ! This allow us to check all the cases : if one of them is -1, there are both
  ! equal  (either positive or -1)
  if (neighb < 0) neighb = neighb2
  if (neighb2 < 0) neighb2 = neighb

  ! If there is at least one south neighbour
  if (neighb > -1) then
    ! Search among all of them
    do i = neighb, neighb2
      prev => list_parts%parts(i+1)
      call north_neighbours(prev, test, test2)
      if (cur%id == id_test) then
        write (*,*) "trying north neighb for", i,"at", prev%i_pos, "result = ", test, test2
      end if
      ! Among all the north neighbours, there must be at least one which is the
      ! current partition

      if (test < 0) test = test2
      if (test2 < 0) test2 = test

      if (neighb > -1) then
        failed_tmp = .True.
        do j = test, test2
          ! i_cur > 0
          if (j == i_cur-1) then 
            failed_tmp = .False.
            exit
          end if
        end do
        if (failed_tmp) then
          !write (*,*) 
          !write (*,*) "Failed : south neighbour relation is not reflexive for"&
          !  , cur%id
          !write (*,*) "at", cur%i_pos, cur%j_pos
          nb_fails = nb_fails + 1
        end if
      end if
    end do
  end if
end subroutine


!-------------------------------------------------------------------------------
! Check the grid is a partition of the sphere (for now, a third only).
!-------------------------------------------------------------------------------
function sphere_is_partitioned(nb_parts) result(nb_fails)
  type(Band_grid), pointer :: grid
  type(Partition), pointer :: cur
  integer, intent(in) :: nb_parts
  integer :: i, j
  integer :: i_start, i_end, nb_lon
  double precision :: area, area_tmp, area_ref
  double precision :: lat1, lat2, dlon
  double precision :: coef, prec
  double precision :: PI = 4.D0*DATAN(1.D0)
  integer :: nb_fails
  logical :: test
  coef = pi/180

  area = 0
  write (*, '(a)', advance="no") "Checking all the sphere is partitioned ..."

  do i = 1, nb_parts
    cur => list_parts%parts(i)
    grid => cur%grid
    ! For each latitude band, compute the area

    call interior_lat_indices(i_start, i_end, grid)
    !print *, "i_start, i ned", i_start, i_end, "size", size(grid%dlon)
    lat1 = get_lat_min(grid)
    ! write (*,*) "indice nblat", grid%nb_lat, i_start, i_end

    test = is_on_last_band(cur%i_pos, cur%zone)
    !write (*,*) "i pos", cur%i_pos
    do j = i_start, i_end
      lat2 = lat1 - grid%dlat(1)
      ! We assume lat1 > lat2
      nb_lon = nb_lon_Band_grid(j, grid) 
      dlon = nb_lon*grid%dlon(j)
      area_tmp = (sin(lat1*coef) - sin(lat2*coef))*coef*dlon
      area = area + area_tmp
      lat1 = lat2
    end do
  end do

  area_ref = 4.*pi
  !  write (*,*) "area ", area, area_ref
  prec = get_precision()

  if (abs(area_ref - area) < prec) then
    nb_fails = 0
  else
    nb_fails = 1
  end if
  call check_nb_fails(nb_fails, nb_parts)
end function

!-------------------------------------------------------------------------------
! The precision is a tenth of the minimal dlon
!-------------------------------------------------------------------------------
function get_precision() result(prec)
  integer :: nb_lat2, nb_cells
  double precision :: prec

  nb_lat2 = get_nb_lat2_Global_grid()
  nb_cells = nb_cells_lat_sector(nb_lat2)
  prec = 0.001*(120./nb_cells)
end function

end module
