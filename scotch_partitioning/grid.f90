!*******************************************************************************
! Contains subroutines and function for grid operations
!*******************************************************************************
module grid
  implicit none

contains
  !*******************************************************************************
  ! Compute the neighbors of cell j (j in [1,nb_cells(i)]) at latitude i
  ! Neighbours are numbered from the first cell
  ! 2 neighbours at most at i-1
  ! (no neighbours by default)
  !*******************************************************************************
  subroutine neighbours_prev(j,i,n_zone,nb_lat,nb_lat2,size_max,pos,neighbours)
    integer, intent(in) :: j,i,size_max,nb_lat,nb_lat2
    integer :: neighbours(size_max)
    integer :: pos, n_prev,n_cur,j1,j2, offset,n_zone
    real :: ratio

    if (i > 1) then
      n_prev = nb_cells_lat(i-1,nb_lat,nb_lat2)
      n_cur = nb_cells_lat(i,nb_lat,nb_lat2)

      ! Offset for true numbering (not only at latitude i)
      offset = n_zone*n_prev
      if (i  > 2) then
        offset = offset + 3*nb_cells(i-2,nb_lat,nb_lat2)
      end if

      ! Equator case
      if (i == nb_lat2+1) then
        neighbours(pos) = j + offset
        pos = pos + 3 
        return
      end if

      ratio = real(n_prev)/n_cur
      j1 = floor(ratio*(j-1) + 1)
      ! Approximation error
      if (j == n_cur) then
        j2 = n_prev
      else
        j2 = ceiling(ratio*j)
      end if

      neighbours(pos) = j1 + offset
      ! According to the hemisphere, we may have 2 or 3 previous neighbours
      if (i < nb_lat2+1) then 
        ! No data duplication
        if (j2 /= j1) then
          neighbours(pos + 1) = j2 + offset
        end if
      else
        neighbours(pos + 1) = j1+1 + offset
        if (j2 /= j1+1) then
          neighbours(pos + 2) = j2 + offset
        end if
      end if
    end if
    pos = pos + 3 ! One empty case in the northern hemisphere
  end subroutine

  !*******************************************************************************
  ! Neighbours at current latitude : 3 at most
  ! (no neighbours by default)
  !*******************************************************************************
  subroutine neighbours_cur(j,i,n_zone,nb_lat,nb_lat2,size_max,pos,neighbours)
    integer, intent(in) :: j,i,size_max,nb_lat,nb_lat2
    integer :: neighbours(size_max)
    integer :: pos,n_cur, offset,n_zone,j2
    n_cur = nb_cells_lat(i,nb_lat,nb_lat2)

    ! Offset for true numbering (not only at latitude i)
    offset = 0
    if (i > 1) then
      offset = offset + 3*nb_cells(i-1,nb_lat,nb_lat2) 
    end if

    ! Special case if we are on the first or last cell (0 or 2pi)
    j2 = j + n_zone*n_cur
    if (j2 == 1) then
      neighbours(pos) = 3*n_cur + offset
    else
      neighbours(pos) = j2-1 + offset
    end if
    if (j2 == 3*n_cur) then
      neighbours(pos+1) = 1 + offset
    else
      neighbours(pos+1) = j2+1 + offset
      end if
    pos = pos + 2
  end subroutine

  !*******************************************************************************
  ! Neighbours at next latitude (no neighbours by default)
  !*******************************************************************************
  subroutine neighbours_next(j,i,n_zone,nb_lat,nb_lat2,size_max,pos,neighbours)
    integer, intent(in) :: j,i,size_max,nb_lat,nb_lat2
    integer :: neighbours(size_max)
    integer :: pos,n_next,n_cur,j1,j2, offset,n_zone
    real :: ratio

    if (i < nb_lat) then
      n_cur = nb_cells_lat(i,nb_lat,nb_lat2)
      n_next = nb_cells_lat(i+1,nb_lat,nb_lat2)

      ! Offset for true numbering (not only at latitude i)
      offset = 3*nb_cells(i,nb_lat,nb_lat2) + n_zone*n_next
      ! Equator case
      if (i == nb_lat2) then
        neighbours(pos) = j + offset
        pos = pos + 3 
        return
      end if

      ratio = real(n_next)/n_cur
      j1 = floor(ratio*(j-1) + 1)
      ! Approximation error
      if (j == n_cur) then
        j2 = n_next
      else
        j2 = ceiling(ratio*j)
      end if

      neighbours(pos) = j1 + offset
      ! According to the hemisphere, we may have 2 or 3 next neighbours
      if (i < nb_lat2) then 
        neighbours(pos+1) = j1+1 + offset
        if (j2 /= j1+1) then
          neighbours(pos+2) = j2 + offset
        end if
      else
        if (j2 /= j1) then
          neighbours(pos+1) = j2 + offset
        end if
      end if
    end if
    pos = pos + 3
  end subroutine


  !*******************************************************************************
  ! Returns all the neighbours of a cell j at latitude i on a zone
  ! Result is an array : 
  ! 1..3 : neighbours at i-1
  ! 4..5 : neighbours at i
  ! 5..8 : neighbours at i+1
  !*******************************************************************************
  subroutine find_neighbours(j,i,n_zone,n_cells,n_max,nb_lat,nb_lat2,size_max,neighbours)
    integer, intent(in) :: j,i,size_max,n_cells,nb_lat,nb_lat2
    integer :: neighbours(size_max)
    integer :: pos,n_cur,n_prev,n_next,k,n_max,n_zone
    neighbours(1:size_max) = -1
    pos = 1

    call neighbours_prev(j,i,n_zone,nb_lat,nb_lat2,size_max,pos,neighbours)
    call neighbours_cur(j,i,n_zone,nb_lat,nb_lat2,size_max,pos,neighbours)
    if (i < n_max) then
      call neighbours_next(j,i,n_zone,nb_lat,nb_lat2,size_max,pos,neighbours)
    end if

    ! Check data does not overflow
    do k = 1,size_max
      if (neighbours(k) > 3*n_cells) then
        write (*,*) "A cell is outside the range !",k,i
      end if
    end do
  end subroutine

  !*******************************************************************************
  ! Number of cells at latitude i on a zone
  !*******************************************************************************
  function nb_cells_lat(i,nb_lat,nb_lat2)
    integer, intent(in) :: i,nb_lat,nb_lat2
    integer :: i2,nb_cells_lat
    i2 = i
    if (i > nb_lat2) then
      i2 = nb_lat + 1 - i;
    end if
    nb_cells_lat = 2*i2 - 1;
    return
  end function

  !*******************************************************************************
  ! Total number of cells at latitude i on a zone
  !*******************************************************************************
  recursive function nb_cells(i,nb_lat,nb_lat2) result(n_cells)
    integer i,nb_lat,nb_lat2
    integer i2,n_cells
    if (i > nb_lat2) then
      i2 = nb_lat - i;
      n_cells = 2*nb_cells(nb_lat2,nb_lat,nb_lat2) - nb_cells(i2,nb_lat,nb_lat2);
    else
      n_cells = i*i;
    end if
    return
  end function


  !****************************************************************************
  ! Compute the number of neighbours (positive entries)
  !****************************************************************************
  function nb_neighbours(size_max,neighbours,n_cells)
    integer,allocatable :: neighbours(:)
    integer :: k, nb_neighbours, size_max, n_cells
    nb_neighbours = 0
    do k = 1,size_max
      if (neighbours(k) > -1) then
        nb_neighbours = nb_neighbours + 1
      end if
    end do
    return
  end function

  !****************************************************************************
  ! Compute the number of arcs up to latitude i for all the sphere
  ! We add the number of arcs on each latitude = nb_nodes*2 (previous and next
  ! lat) + 2*nb_cells (current lat)
  ! At the boundary, cells are linked
  !****************************************************************************
  function nb_arcs(i,nb_lat,nb_lat2)
    integer, intent(in) :: i,nb_lat,nb_lat2
    integer :: nb_arcs,k
    nb_arcs = 0
    do k = 1,i-1
      nb_arcs = nb_arcs + 6*nb_nodes(k,nb_lat,nb_lat2) + 6*nb_cells_lat(k,nb_lat,nb_lat2)
    end do
    nb_arcs = nb_arcs + 6*nb_cells_lat(i,nb_lat,nb_lat2)
    return
  end function

  !****************************************************************************
  ! Number of nodes of latitude i (with i in [1,nb_lat]) on a zone
  !****************************************************************************
  function nb_nodes(i,nb_lat,nb_lat2)
    integer,intent(in) :: i,nb_lat,nb_lat2
    integer :: nb_nodes
    if (i < nb_lat2) then
      nb_nodes = 4*i-1
    ! At the Equator, only nb_cells nodes
    else if (i == nb_lat2) then
      nb_nodes = nb_cells_lat(i,nb_lat,nb_lat2)
    else
      nb_nodes = 4*(nb_lat-i)-1
    end if
    ! At the poles
    if (nb_nodes < 0) then
      nb_nodes = 0
    end if
    return
  end function

  !****************************************************************************
  ! Write all the neighbours for the cell i in the array neighbours
  ! Also compute the sum of the neighbours
  !****************************************************************************
  subroutine get_neighbours(id_file,size_max,neighbours,n_cells,sum_neighbours)
    character(len=*),parameter :: format_i='(i7)'
    integer,allocatable :: neighbours(:)
    integer,intent(in) :: size_max,id_file,n_cells
    integer :: sum_neighbours,k,n_neighbours,cur
    n_neighbours = nb_neighbours(size_max,neighbours,n_cells)
    sum_neighbours = sum_neighbours + n_neighbours
    write (id_file,'(i4)',advance='no') n_neighbours
    do k = 1,size_max
      cur = neighbours(k) 
      if (cur > -1) then
        write (id_file,format_i,advance='no') neighbours(k)
      end if
    end do
    write (id_file,*)
  end subroutine

end module
