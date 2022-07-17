!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Global_grid
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> The global grid for the sphere. Features : area-preserving grid.
!> Contains temporarily the data for input/output
!
!-------------------------------------------------------------------------------

module Global_grid_class
use Configuration_class
implicit none

public
type Global_grid
  !> Total number of latitudes from North Pole
  integer :: nb_lat
  !> Number of latitudes from North Pole to the Equator
  integer :: nb_lat2
  !> Number of sectors on a hemisphere
  !> An hemisphere is divided into sectors where the grid is identical
  !> In our case, it will be 3*2 = 6
  integer :: nb_sectors
  !> Height of the cells (set to constant for now)
  double precision, allocatable :: dlat(:)

end type

type (Global_grid), target, private :: g_grid

! MPI type for exchanging data (defined once)
integer, private :: mpi_global_grid

character(*), parameter :: fname_global = "global_grid_class.f90"
contains 

!-----------------------------------------------------------------------------
!> Constructor for global grid (as a global variable)
!-----------------------------------------------------------------------------
subroutine new_Global_grid(nb_lat2)
  integer, intent(in) :: nb_lat2
  integer :: rank, ierr
  integer :: nb_dlat
  nb_dlat = 1

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call check_mpi_error(ierr, "get number processes", "new_Global_grid",&
    fname_global)

  print *, "New global grid with ", nb_lat2, "latitudes on an hemisphere"
  g_grid%nb_lat2 = nb_lat2
  g_grid%nb_lat = 2*g_grid%nb_lat2
  if (has_single_partition()) then
    g_grid%nb_sectors = 1
  else
    g_grid%nb_sectors = 3
  end if

  if (rank == 0) then
    allocate (g_grid%dlat(1))
    g_grid%dlat(1) = 90.d0/g_grid%nb_lat2
  end if
end subroutine

subroutine broadcast_Global_grid()
  integer :: n, rank, ierr


  ! Send dlat size to other process
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call check_mpi_error(ierr, "rank", "broadcast_Global_grid", fname_global)

  if (rank == 0) n = size(g_grid%dlat)
  call mpi_bcast(n, 1 , mpi_integer, 0, mpi_comm_world, ierr)
  if (rank /= 0) allocate(g_grid%dlat(n))

  call new_MPI_Global_grid()
  call mpi_bcast(g_grid, 1, mpi_global_grid, 0, mpi_comm_world, ierr)
  call check_mpi_error(ierr, "broadcast", "broadcast_Global_grid", fname_global)
end subroutine

!-----------------------------------------------------------------------------
!> Derived MPI type for sending/receiving grid
!-----------------------------------------------------------------------------
subroutine new_MPI_Global_grid()
  ! Length, displacement, type
  integer, parameter :: n = 4
  integer :: mlength(n), mtype(n)
  integer(kind=mpi_address_kind) :: mlocation(n)
  integer(kind=mpi_address_kind) :: start
  integer :: i, ierr
  character(*), parameter :: func_name = "new_MPI_Global_grid" 

  ! Get the address and check the error code
  call mpi_get_address(g_grid, start, ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_global)

  call get_address(g_grid%nb_lat, mlocation(1), func_name, fname_global)
  call get_address(g_grid%nb_lat2, mlocation(2), func_name, fname_global)
  call get_address(g_grid%nb_sectors, mlocation(3), func_name, fname_global)
  call get_address(g_grid%dlat, mlocation(4), func_name, fname_global)

  mlength = (/1, 1, 1, size(g_grid%dlat)/)
  mtype = (/mpi_integer, mpi_integer, mpi_integer, mpi_double_precision/)

  do i = 1, size(mlength)
    mlocation(i) = mlocation(i) - start
  end do

  call mpi_type_create_struct(size(mlength), mlength, mlocation, mtype, &
    mpi_global_grid, ierr)
  call check_mpi_error(ierr, "create struct type", func_name, fname_global)
  call mpi_type_commit(mpi_global_grid, ierr)
  call check_mpi_error(ierr, "commit type", func_name, fname_global)

end subroutine


!-----------------------------------------------------------------------------
!> Return the number of cells at latitude 90 - i/2 on a sector with i > 0
!> We assume nb_lat >= nb_lat2
!-----------------------------------------------------------------------------
function nb_cells_lat_sector(i) result(nb)
  integer, intent(in) :: i
  integer :: nb, i2
  i2 = i
  if (i > g_grid%nb_lat2) then
    i2 = g_grid%nb_lat - i + 1
  end if
  nb = 2*i2-1
end function

!-----------------------------------------------------------------------------
!> Return the number of cells at latitude 90 - i/2 on all sectors
!> For faster computing, we suppose there always are 3 sectors
!-----------------------------------------------------------------------------
function nb_cells_lat(i) result(nb)
  integer, intent(in) :: i
  integer :: nb

  nb = 3*nb_cells_lat_sector(i)
end function

!-----------------------------------------------------------------------------
!> Return the total number of cells up to latitude 90 - i/2 on a sector with i > 0
!-----------------------------------------------------------------------------
function sum_nb_cells_sector(i) result(nb)
  integer, intent(in) :: i
  integer :: nb
  if (i < 0) then
    call print_error("Must have a latitude number > 0 in sum_nb_cells_sector", &
      "sum_nb_cells_sector", fname_global);
  else if (i == 0) then
    nb = 0
    return
  end if

  ! On the other hemisphere, we substract from the total
  if (i > g_grid%nb_lat2) then
    nb = 2*(g_grid%nb_lat2)**2 - (g_grid%nb_lat- i)**2
  else
    nb = i**2
  end if
end function

!-----------------------------------------------------------------------------
!> Return the total number of cells up to latitude 90 - i/2 on all sectors
!-----------------------------------------------------------------------------
function sum_nb_cells(i) result(nb)
  integer, intent(in) :: i
  integer :: nb

  nb = 3*sum_nb_cells_sector(i)
end function

function nb_cells_total() result (nb)
  integer :: nb

  nb = 6*g_grid%nb_lat2**2
end function

!-----------------------------------------------------------------------------
!> Return the total number of zonal winds up to latitude 90 - i/2 on all sectors
!-----------------------------------------------------------------------------
function sum_nb_zonal(i) result(nb)
  integer, intent(in) :: i
  integer :: nb

  if (i < 0) then
    call print_error("Must have a latitude number > 0 in sum_nb_zonal", &
      "nb_cells_total_sector", fname_global);
  else if (i == 0) then
    nb = 0
    return
  end if

  ! On the other hemisphere, we substract from the total
  if (i > g_grid%nb_lat2) then
    nb = 6*(g_grid%nb_lat2)**2 + g_grid%nb_lat
    nb = nb - (3*(g_grid%nb_lat- i)**2 + (g_grid%nb_lat- i))
  else
    nb = 3*i**2 + i
  end if
end function

function nb_merid_winds_total() result(nb)
  integer :: nb

  nb = sum_nb_merid(2*g_grid%nb_lat2-1)
end function
 
!-----------------------------------------------------------------------------
!> Return the number of merid winds at latitude 90 - i/2 on all sectors
!-----------------------------------------------------------------------------
function nb_merid_lat(i) result(nb)
  integer, intent(in) :: i
  integer :: i2, nb, nb_lat2

  nb_lat2 = g_grid%nb_lat2
  if (i < nb_lat2 ) then
    nb = 3*(4*i-1)
  ! One-to-one correspondance
  else if (i == nb_lat2) then
    nb = nb_cells_lat(i)
  else
    i2 = 2*nb_lat2-i
    nb = 3*(4*i2-1)
  end if
end function
 
!-----------------------------------------------------------------------------
!> Return the total number of merid winds up to latitude 90 - i/2 on all sectors
!-----------------------------------------------------------------------------
recursive function sum_nb_merid(i) result(nb)
  integer, intent(in) :: i
  integer :: nb, nb_lat2, i2

  if (i < 0) then
    call print_error("Must have a latitude number > 0 in sum_nb_merid", &
      "sum_nb_merid", fname_global);
  else if (i == 0) then
    nb = 0
    return
  end if

  nb_lat2 = g_grid%nb_lat2
  ! On the other hemisphere, we substract from the total
  if (i < nb_lat2) then
    ! Analytical formal = 3*sum(4k-1) for k up to i
    nb = 3*i*(2*i+1)
  ! One to one correspondance at the equator
  else if (i == nb_lat2) then
    i2 = i-1
    nb = 3*i2*(2*i2+1) + nb_cells_lat(i)
    ! No merid winds at south pole
  else if (i == 2*nb_lat2) then
    nb = sum_nb_merid(i-1)
  else
    ! Total number
    i2 = nb_lat2-1
    nb = 2*(3*i2*(2*i2+1)) + nb_cells_lat(nb_lat2)
    ! No merid winds at south pole, so 2*nb_lat2-1
    i2 = 2*nb_lat2-1-i
    nb = nb - 3*i2*(2*i2+1)
  end if
end function


function nb_zonal_winds_total() result (nb)
  integer :: nb

  nb = 6*g_grid%nb_lat2**2 + g_grid%nb_lat
end function

!-----------------------------------------------------------------------------
!> Compute the number of latitude from the number of cells
!-----------------------------------------------------------------------------
function nb_lat2_from_nb_cells(nb) result(nb_lat2)
  integer, intent(in) :: nb
  integer :: nb_lat2

  ! nb_cells = 6*nb_lat2^2 analytically
  nb_lat2 = int(sqrt(nb/6.d0))
  if (6*nb_lat2**2 /= nb) then
    call print_error("Wrong number of cells, should be 6n^2",&
      "nb_lat2_from_nb_cells", fname_global)
  end if
end function

!-----------------------------------------------------------------------------
!> Return the north neighbours for cell (i,j) in global grid
!> @param j_neighb1, j_neighb2 : interval containing the neighbours
!-----------------------------------------------------------------------------
subroutine north_neighbour_cell_Global(j_neighb1, j_neighb2, i, j)
  integer, intent(in) :: i, j
  integer :: j_neighb1, j_neighb2
  integer :: i2
  if (i < g_grid%nb_lat2 + 1) then
    call north_neighbour_cell_northern(j_neighb1, j_neighb2, i, j)
  else if (i == g_grid%nb_lat2 + 1) then
    j_neighb1 = j
    j_neighb2 = j
  else
    i2 = g_grid%nb_lat - i + 1
    call south_neighbour_cell_northern(j_neighb1, j_neighb2, i2, j)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Return the south neighbours for cell (i,j) in global grid
!> @param j_neighb1, j_neighb2 : interval containing the neighbours
!-----------------------------------------------------------------------------
subroutine south_neighbour_cell_Global(j_neighb1, j_neighb2, i, j) 
  integer, intent(in) :: i, j
  integer :: j_neighb1, j_neighb2
  integer :: i2

  if (i < g_grid%nb_lat2) then
    call south_neighbour_cell_northern(j_neighb1, j_neighb2, i, j)
  else if (i == g_grid%nb_lat2) then
    j_neighb1 = j
    j_neighb2 = j
  else
    i2 = g_grid%nb_lat - i + 1
    call north_neighbour_cell_northern(j_neighb1, j_neighb2, i2, j)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Returns cell position in the first sector, along with the initial sector
!> Also return the middle cell
!-----------------------------------------------------------------------------
subroutine switch_first_sector(j2, sector, mid, i, j)
  integer, intent(in) :: i, j
  integer :: j2, sector, mid
  integer :: nb_cells_sector

  nb_cells_sector = nb_cells_lat_sector(i)
  j2 = mod(j, nb_cells_sector)
  if (j2 == 0) j2 = nb_cells_sector
  sector = (j-1)/nb_cells_sector

  mid = (nb_cells_sector + 1)/2
end subroutine

!-----------------------------------------------------------------------------
!> Only works for the northern hemisphere
!-----------------------------------------------------------------------------
subroutine north_neighbour_cell_northern(j_neighb1, j_neighb2, i, j)
  integer, intent(in) :: i, j
  integer :: j_neighb1, j_neighb2
  integer :: j2, sector, mid, offset
  integer :: nb_cells, nb_cells_next


  nb_cells = nb_cells_lat(i)
  nb_cells_next = nb_cells_lat(i-1)
  if (j < 1 .or. j > max(nb_cells, nb_cells_next)) then
    print *, "ij", i, j
    call print_error("cell outside range", "north_neighbour_cell", fname_global)
  end if
  call switch_first_sector(j2, sector, mid, i, j)

  nb_cells = nb_cells_lat_sector(i)
  nb_cells_next = nb_cells_lat_sector(i-1)
  if (j2 == 1) then
    j_neighb1 = 1
    j_neighb2 = j_neighb1
  else if (j2 == nb_cells) then
    j_neighb1 = nb_cells_next
    j_neighb2 = j_neighb1
  else if (j2 == mid) then
    j_neighb1 = j2-1
    j_neighb2 = j_neighb1
  else
    offset = 0
    if (j2 > mid) offset = -1
    j_neighb1 = j2-1 + offset
    j_neighb2 =  j2 + offset
  end if

  ! Going back to complete grid
  j_neighb1 = j_neighb1 + sector*nb_cells_next
  j_neighb2 = j_neighb2 + sector*nb_cells_next

end subroutine

!-------------------------------------------------------------------------------
!> Only works for the northern hemisphere
!-------------------------------------------------------------------------------
subroutine south_neighbour_cell_northern(j_neighb1, j_neighb2, i, j)
  integer, intent(in) :: i, j
  integer :: j_neighb1, j_neighb2
  integer :: nb_cells_sector, j2, mid, sector
  integer :: nb_cells, nb_cells_next
  integer :: offset

  nb_cells = nb_cells_lat(i)
  nb_cells_next = nb_cells_lat(i+1)
  if (j < 1 .or. j > max(nb_cells, nb_cells_next)) then
    print *, "ij", i,j, "n nnext", nb_cells, nb_cells_next
    call print_error("cell outside range", "south neighbour cell", fname_global)
  end if

  call switch_first_sector(j2, sector, mid, i, j)

  if (j2 == nb_cells) then
    j_neighb1 = j2
    j_neighb2 = j_neighb1
  else if (j2 == mid) then
    j_neighb1 = j2
    j_neighb2 = j2 + 2
  else
    offset = 0
    if (j2 > mid) offset = 1
    j_neighb1 = j2 + offset
    j_neighb2 = j2 + 1 + offset
  end if

  ! Going back to complete grid
  nb_cells_sector = nb_cells_lat_sector(i+1)
  j_neighb1 = j_neighb1 + sector*nb_cells_sector
  j_neighb2 = j_neighb2 + sector*nb_cells_sector
end subroutine

!-----------------------------------------------------------------------------
!> Return the width of a cell at latitude 90 - i/2 with i > 0
!-----------------------------------------------------------------------------
function get_dlon_Global_grid(i) result (dlon)
  integer, intent(in) :: i
  double precision dlon

  dlon = 120.d0/nb_cells_lat_sector(i)
end function

!-----------------------------------------------------------------------------
!> Returns the number of latitude on an hemisphere
!-----------------------------------------------------------------------------
function get_nb_lat2_Global_grid() result (nb)
  integer :: nb
  nb = g_grid%nb_lat2
end function

!-----------------------------------------------------------------------------
!> Returns the number of all latitudes
!-----------------------------------------------------------------------------
function get_nb_lat_Global_grid() result (nb)
  integer :: nb
  nb = g_grid%nb_lat
end function

!-----------------------------------------------------------------------------
!> Returns the number of all latitudes
!-----------------------------------------------------------------------------
function get_nb_sectors() result (nb)
  integer :: nb
  nb = g_grid%nb_sectors
end function

!-----------------------------------------------------------------------------
!> Associate a pointer to the global dlat
!> Returns the number of all latitudes
!-----------------------------------------------------------------------------
function  get_dlat_Global_grid()  result(dlat)
  double precision :: dlat(1)

  dlat = g_grid%dlat
end function

!-----------------------------------------------------------------------------
!> Accessor : associates a pointer to the class.
!-----------------------------------------------------------------------------
subroutine get_Global_grid(ptr)
  type (Global_grid), pointer :: ptr
  ptr => g_grid
end subroutine

!-----------------------------------------------------------------------------
!> Get number of latitude
!-----------------------------------------------------------------------------
function nb_lat_Global_grid() result (nb)
  integer :: nb
  nb = g_grid%nb_lat
end function

!-----------------------------------------------------------------------------
!> Get half number of latitude
!-----------------------------------------------------------------------------
function nb_lat2_Global_grid() result (nb)
  integer :: nb
  nb = g_grid%nb_lat2
end function

!-----------------------------------------------------------------------------
!> Free the temporary grid
!-----------------------------------------------------------------------------
subroutine free_Global_grid()

  if (allocated(g_grid%dlat)) then
    deallocate (g_grid%dlat)
  end if
end subroutine

end module
