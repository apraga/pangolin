!------------------------------------------------------------------------------- ! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: List_Tracers
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Contains a list of tracers corresponding to one chemical scheme.
!> It contains only tracer ratio. Contrary to other class, there is only one
!> instance of the list, and it is only used in the module.
!
!-------------------------------------------------------------------------------
module List_Tracers_class
use Data_check
use Message
use Parameters

implicit none

type List_Tracers
  !! Everything is private by default
  !private

  !> Format of the arrays : (cell, tracer) as fortran stores by column
  !> And we browse all of the grid for a single tracer

  !> We use 2D-arrays as each tracer has the same number of cells.
  !> Slope used for advection
  double precision, allocatable :: slope(:,:)
  !> Used for meridional advection
  double precision, allocatable :: gradient(:,:)
  !> Tracer ratio  (with ghost cells)
  double precision, allocatable :: ratio(:,:)
  !> Tracer fluxes used for advection
  double precision, allocatable :: F_merid(:,:)
  double precision, allocatable :: F_zonal(:,:)
  !> A list of MPI types, one for each neighbour. Here we have interior cells for
  !> all tracers. Same type for slope, gradient and ratio
  integer, allocatable :: mpi_interior_tracers(:)
  !> A list of MPI type, one for each neighbour. Here we have ghost cells for
  !> all tracers.
  integer, allocatable :: mpi_ghost_tracers(:)

#ifdef SYMMETRIC_SPLIT
  !> Temporary buffers
  double precision, allocatable :: buffer(:,:,:)
#endif
end type

! For output_
character(*), parameter :: fname_tracers = "list_tracers_class.F90"



contains 

!-------------------------------------------------------------------------------
!> Assume each tracer is defined on the same number of cells
!> Need the size of slope/gradient/ratio, F_zonal, F_merid
!-------------------------------------------------------------------------------
subroutine new_List_Tracers(this, size_tracer, size_zonal, size_merid)
  type(List_Tracers) :: this
  integer, intent(in) :: size_tracer, size_zonal, size_merid
  character(*), parameter :: func_name ="new_List_Tracers" 

  call check_allocated(this%ratio, "ratio", func_name, fname_tracers)
  allocate(this%ratio(size_tracer, nb_tracers))

  call check_allocated(this%slope, "slope", func_name, fname_tracers)
  allocate(this%slope(size_tracer, nb_tracers))

  call check_allocated(this%gradient, "gradient", func_name, fname_tracers)
  allocate(this%gradient(size_tracer, nb_tracers))

  call check_allocated(this%F_zonal, "F_zonal", func_name, fname_tracers)
  allocate(this%F_zonal(size_zonal, nb_tracers))

  call check_allocated(this%F_merid, "F_merid", func_name, fname_tracers)
  allocate(this%F_merid(size_merid, nb_tracers))

#ifdef SYMMETRIC_SPLIT
  allocate(this%buffer(2,size_tracer, nb_tracers))
#endif

end subroutine

subroutine init_tracer_ratio_all(this)
  type(List_Tracers) :: this
  this%ratio = 0.
end subroutine

subroutine init_tracer_data(this)
  type(List_Tracers) :: this

  this%slope = 0.
  this%gradient = 0.
  this%F_zonal = 0.
  this%F_merid = 0.
end subroutine

subroutine unset_gradients_all(this)
  type(List_Tracers) :: this

  this%gradient = UNDEFINED
end subroutine

#ifdef SYMMETRIC_SPLIT
!> Save all ratio in a temporary array
subroutine cp_ratio_to_buffer_all(this, id)
  type(List_Tracers) :: this
  integer, intent(in) :: id

  this%buffer(id,:,:) = this%ratio
end subroutine

!> Save all ratio in a temporary array
!> Assume buffer already allocated
subroutine cp_buffer_to_ratio_all(this, id)
  type(List_Tracers) :: this
  integer, intent(in) :: id

  this%ratio = this%buffer(id,:,:)
end subroutine

!> Average the ratio with the one temporarily stored
subroutine average_with_buffer_all(this, id)
  type(List_Tracers) :: this
  integer :: size_tracer, nb_tracers
  integer, intent(in) :: id
  double precision, pointer :: p(:,:)

  this%ratio = 0.5*(this%ratio + this%buffer(id,:,:))
end subroutine
#endif


!> Create one MPI type identical for all tracers

!-------------------------------------------------------------------------------
!> Creates a list of MPI types for interior cells, one for each partition 
!> neighbours, and for all tracers
!-------------------------------------------------------------------------------
subroutine new_MPI_interior_tracers(this, mpi_interior_cells)
  type(List_Tracers) :: this
  integer, intent(in) :: mpi_interior_cells(:) 
  integer :: n, j, k, ierr
  character(*), parameter :: func_name = "new_MPI_interior_tracers"
  integer :: rank
  integer, allocatable :: mtype(:), mlength(:)
  integer(kind=mpi_address_kind), allocatable :: mlocation(:)

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  n = size(mpi_interior_cells)
  allocate (this%mpi_interior_tracers(n))
  this%mpi_interior_tracers = MPI_DATATYPE_NULL
  allocate (mlocation(nb_tracers))
  allocate (mlength(nb_tracers))
  allocate (mtype(nb_tracers))

  ! Create a list of contiguous type (same for each tracer)
  do k = 1, n

    do j = 1, nb_tracers
      call mpi_get_address(this%ratio(1,j), mlocation(j), ierr)
    end do
    mlocation = mlocation - mlocation(1)

    mlength = 1
    mtype = mpi_interior_cells(k)
    call mpi_type_create_struct(nb_tracers, mlength, mlocation, mtype, &
        this%mpi_interior_tracers(k), ierr)
    call check_mpi_error(ierr, "create vector type", func_name, fname_tracers)
    
    call mpi_type_commit(this%mpi_interior_tracers(k), ierr)
    call check_mpi_error(ierr, "create interior type", func_name, fname_tracers)
  end do
  deallocate(mlocation)
  deallocate(mlength)
  deallocate(mtype)

end subroutine

!-------------------------------------------------------------------------------
!> Creates a list of MPI types for ghost cells, one for each partition neighbour.
!-------------------------------------------------------------------------------
subroutine new_MPI_ghost_tracers(this, mpi_ghost_cells)
  type(List_Tracers) :: this
  integer, intent(in) :: mpi_ghost_cells(:) 
  integer :: n, j, k, ierr
  integer :: offset
  character(*), parameter :: func_name = "new_MPI_ghost_tracers"
  integer, allocatable :: mtype(:), mlength(:)
  integer(kind=mpi_address_kind), allocatable :: mlocation(:)

  n = size(mpi_ghost_cells)
  allocate (this%mpi_ghost_tracers(n))
  allocate (mlocation(nb_tracers))
  allocate (mlength(nb_tracers))
  allocate (mtype(nb_tracers))

  do k = 1, n
    ! Data is stored is a 2d array so we need to replicate the type for ghost 
    ! cells on each column

    do j = 1, nb_tracers
      call mpi_get_address(this%ratio(1,j), mlocation(j), ierr)
    end do
    mlocation = mlocation - mlocation(1)

    mlength = 1
    mtype = mpi_ghost_cells(k)
    call mpi_type_create_struct(nb_tracers, mlength, mlocation, mtype, &
        this%mpi_ghost_tracers(k), ierr)
    call check_mpi_error(ierr, "create vector type", func_name, fname_tracers)
    call mpi_type_commit(this%mpi_ghost_tracers(k), ierr)
    call check_mpi_error(ierr, "create ghost type", func_name, fname_tracers)
  end do
  deallocate (mlocation)
  deallocate(mlength)
  deallocate(mtype)
end subroutine

!-------------------------------------------------------------------------------
!> Start a sending request for all tracers data
!> @param i_pos : position in the MPI type array
!> @param k : partition neighbour
!> @param l : position in the request array
!-------------------------------------------------------------------------------
subroutine start_send_data(this, i_pos, datatype, k, request, l)
  type(List_Tracers) :: this
  integer, intent(in) :: i_pos, datatype
  integer, intent(in) :: k
  integer :: request(:), ierr, l
  integer :: mpitype, tsize
  double precision, pointer :: ptr(:,:)

  mpitype = get_tracer_mpi_interior(this, i_pos)
#ifdef DEBUG
  if (l < 1 .or. l > size(request)) then
    call print_error("outside request array","start_receive_data", fname_tracers)
  end if
#endif

  call associate_pointer(this, ptr, datatype)
  call mpi_isend(ptr, 1, mpitype, k, 0, mpi_comm_world, request(l), ierr)
  l = l + 1
end subroutine

!> Associate the pointer to the data for a given tracer
subroutine associate_tracer_pointer(this, ptr, datatype, tracer)
  type(List_Tracers), target :: this
  integer, intent(in) :: datatype, tracer
  double precision, pointer :: ptr(:)

  if (datatype == IS_RATIO) then
    ptr => this%ratio(:,tracer)
  else if (datatype == IS_SLOPE) then
    ptr => this%slope(:,tracer)
  else if (datatype == IS_GRADIENT) then
    ptr => this%gradient(:,tracer)
  else
    call print_error("Wrong datatype", "associate_pointer", fname_tracers)
  end if

end subroutine

subroutine associate_pointer(this, ptr, datatype)
  type(List_Tracers), target :: this
  integer, intent(in) :: datatype
  double precision, pointer :: ptr(:,:)

  if (datatype == IS_RATIO) then
    ptr => this%ratio
  else if (datatype == IS_SLOPE) then
    ptr => this%slope
  else if (datatype == IS_GRADIENT) then
    ptr => this%gradient
  else
    call print_error("Wrong datatype", "associate_pointer", fname_tracers)
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Start a receiving request for all tracers
!> @param i_pos : position in the MPI type array
!> @param k : partition neighbour
!> @param l : position in the request array
!-------------------------------------------------------------------------------
subroutine start_receive_data(this, i_pos, datatype, k, request, l)
  type(List_Tracers) :: this
  integer, intent(in) :: i_pos, datatype
  integer, intent(in) :: k
  integer :: request(:), ierr, l
  integer :: mpitype, tsize
  double precision, pointer :: ptr(:,:)

  mpitype = get_tracer_mpi_ghost(this, i_pos)
#ifdef DEBUG
  if (l < 1 .or. l > size(request)) then
    call print_error("Outside request array","start_receive_data", fname_tracers)
  end if
#endif

  call associate_pointer(this, ptr, datatype)
  call mpi_irecv(ptr, 1, mpitype, k, 0, mpi_comm_world, request(l), ierr)
  l = l + 1
end subroutine

!-------------------------------------------------------------------------------
!>  Check a tracer mpi type contains only cells inside the array
!-------------------------------------------------------------------------------
subroutine check_mpi_types(this, offsets, blocklens, inside)
  type(List_Tracers) :: this
  integer :: blocklens(:), offsets(:)
  integer, intent(in) :: inside
  integer :: i, cur, ratio_size

  ratio_size = size(this%ratio(:,1))

  do i = 1, size(offsets)
    cur = offsets(i) + blocklens(i)

    if (cur > ratio_size .or. offsets(i) < 0) then
      write (*,*) "offset len",offsets(i),blocklens(i), "i", i
      print *, "total size", ratio_size
      if (inside == 1) then
        call print_error("Error in interior MPI type, cell outside the array", &
          "check_mpi_types", fname_tracers)
      else
        call print_error("Error in ghost MPI type, cell outside the array", &
          "check_mpi_types", fname_tracers)
      end if
    end if
  end do

end subroutine


!-------------------------------------------------------------------------------
!> Getter function. Do not use them directly. Band_grid_class has better getter
!> based on grid indices.
!-------------------------------------------------------------------------------
function get_tracer_ratio(this, pos, tracer) result(val)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, tracer
  double precision :: val

  val = this%ratio(pos, tracer)
end function

function get_tracer_slope(this, pos, tracer) result(val)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, tracer
  double precision :: val

  val = this%slope(pos, tracer)
end function

function get_tracer_gradient(this, pos, tracer) result(val)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, tracer
  double precision :: val

  val = this%gradient(pos, tracer)
end function

function get_ratio_size(this, tracer) result(nb)
  type(List_Tracers) :: this
  integer, intent(in) :: tracer
  integer :: nb
  
  nb = size(this%ratio(:,tracer))
end function

function get_tracer_zonal_flux(this, pos, tracer) result(val)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, tracer
  double precision :: val

  val = this%F_zonal(pos, tracer)
end function

function get_tracer_merid_flux(this, pos, tracer) result(val)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, tracer
  double precision :: val

  val = this%F_merid(pos, tracer)
end function

!-------------------------------------------------------------------------------
!> Same as getters.
!-------------------------------------------------------------------------------
subroutine set_tracer_ratio(this, val, pos, tracer)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, tracer
  double precision :: val

  this%ratio(pos, tracer) = val
end subroutine

subroutine set_tracer_gradient(this, val, pos, tracer)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, tracer
  double precision :: val

  this%gradient(pos, tracer) = val
end subroutine

subroutine set_tracer_slope(this, val, pos, tracer)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, tracer
  double precision :: val

  this%slope(pos, tracer) = val
end subroutine

!-------------------------------------------------------------------------------
!> We can use these function directly (no referencing needed)
!> @param id : tracer id
!-------------------------------------------------------------------------------
subroutine set_tracer_zonal_flux(this, val, pos, id)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, id
  double precision :: val

  this%F_zonal(pos, id) = val
end subroutine

subroutine set_tracer_merid_flux(this, val, pos, id)
  type(List_Tracers) :: this
  integer, intent(in) :: pos, id
  double precision :: val

  this%F_merid(pos, id) = val
end subroutine

subroutine get_tracer_ratio_array(this, ptr)
  type(List_Tracers), target :: this
  double precision, pointer :: ptr(:,:)

  ptr => this%ratio
end subroutine

function get_tracer_mpi_interior(this, k) result(val)
  type(List_Tracers) :: this
  integer :: val, k

#ifdef DEBUG
  if (k < 1 .or. k > size(this%mpi_interior_tracers)) then
    call print_error("indice outside mpi type array", &
      "get_tracer_mpi_interior", fname_tracers)
  end if
#endif
  val = this%mpi_interior_tracers(k)
end function

function get_tracer_mpi_ghost(this, k) result(val)
  type(List_Tracers) :: this
  integer :: val, k

#ifdef DEBUG
  if (k < 1 .or. k > size(this%mpi_ghost_tracers)) then
    call print_error("indice outside mpi type array", &
      "get_tracer_mpi_ghost", fname_tracers)
  end if
#endif
  val = this%mpi_ghost_tracers(k)
end function



subroutine free_List_Tracers(this)
  type(List_Tracers) :: this
  integer :: i, ierr

  call free_array(this%ratio, "ratio", fname_tracers)
  call free_array(this%slope, "slope", fname_tracers)
  call free_array(this%gradient, "gradient", fname_tracers)
  call free_array(this%F_zonal, "F_zonal", fname_tracers)
  call free_array(this%F_merid, "F_merid", fname_tracers)

#ifdef SYMMETRIC_SPLIT
  if (allocated(this%buffer)) deallocate(this%buffer)
#endif

  if (allocated(this%mpi_interior_tracers) .and. &
  allocated(this%mpi_ghost_tracers)) then
    do i = 1, size(this%mpi_interior_tracers)
      call mpi_type_free(this%mpi_interior_tracers(i), ierr)
      call check_error(ierr,"trying to free interior MPI type", &
        "free_List_Tracers", fname_tracers)

      call mpi_type_free(this%mpi_ghost_tracers(i), ierr)
      call check_error(ierr,"trying to free interior MPI type", &
        "free_List_Tracers", fname_tracers)
    end do
    call free_array(this%mpi_interior_tracers, "mpi_interior_tracers",&
      fname_tracers)
    call free_array(this%mpi_ghost_tracers, "mpi_ghost_tracers",&
      fname_tracers)
  end if

end subroutine


end module
