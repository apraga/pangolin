!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Data_check
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Data checks (arrays, parameters etc)
!
!-------------------------------------------------------------------------------
module Data_check

use Message
use Parameters

#ifdef HDF5
! Reduce compilation time, as we only need one variable
use hdf5, only: hsize_t
#endif

implicit none

#ifdef WITH_PGF90
! Non standard library for PGF. Allow us to read command line arguments
include 'lib3f.h'
#endif


!-------------------------------------------------------------------------------
!> Generic function for checking arrays of int or double precision
!-------------------------------------------------------------------------------
interface check_allocated
  module procedure check_allocated_int, check_allocated_double, &
  check_allocated_double2d
end interface

!-------------------------------------------------------------------------------
!> Generic function for freeing int/double precision arrays
!-------------------------------------------------------------------------------
interface free_array

#ifdef HDF5
  module procedure free_array_int, free_array_double, free_array_double2d,&
  free_array_int_hdf5
#else
  module procedure free_array_int, free_array_double, free_array_double2d
#endif
end interface

!-------------------------------------------------------------------------------
!> Generic function for getting the address of int/double precision arrays
!-------------------------------------------------------------------------------
interface get_address
  module procedure get_address_double, get_address_real, get_address_int, &
    get_address_long_int, get_address_double_array, get_address_int_array, &
    get_address_char

end interface


contains 

!-----------------------------------------------------------------------------
!> Check the arguments are valid
!> Arguments : 
!>   --config=filename (anther config file)
!>   --no-run (only writes the partitioning, no simulation)
!-----------------------------------------------------------------------------
subroutine check_arguments(configfile)
  character(*) :: configfile
  integer :: nb_args, i, rank, ierr
  character(LINE_WIDTH) :: arg

  ! Check optional arguments
  nb_args = iargc()
  !write (*,*) "nb args", nb_args

  if (nb_args > 2) then
    call print_error("Too much arguments (1 needed)", "check_arguments", &
      "simulation.f90")
  end if

  ! Default values
  configfile = "config"
  ! Read arguments
  do i = 1, nb_args
    call getarg(i, arg)
    if (trim(arg(1:8)) == "--config") then
      configfile = arg(10:)
      elseif (trim(arg(1:8)) == "--no-run") then
      NO_RUN = .True.
      call mpi_comm_rank(mpi_comm_world, rank, ierr) 
      if (rank == 0) then
        call print_warning("Simulation will not run, checks disabled")
      end if
    else
      call print_error("Argument not recognized : "// trim(arg) // "/" &
        // "Usage : simulation [--config=filename] [--no-run]", &
        "check_arguments", "simulation.f90")
    end if
  end do

end subroutine


!-------------------------------------------------------------------------------
!> Manage MPI errors. We just print a custom error message and abort.
!> We only do that in debug mode 
!> @param ierr : return code
!> @param descr : small description of the work being done
!> @param func_name : function name
!> @param fname : file name
!-------------------------------------------------------------------------------
subroutine check_mpi_error(ierr, descr, func_name, fname)
  integer, intent(in) :: ierr
  character(*) :: descr, func_name, fname
  character(LINE_WIDTH) :: error_string
  integer :: length, res_err

#ifdef DEBUG
  if (ierr /= mpi_success) then
    call mpi_error_string(ierr, error_string, length, res_err)
    call print_error(trim(error_string)//" for operation "//descr, func_name, fname)
    call mpi_abort(mpi_comm_world, res_err)
  end if
#endif

end subroutine


!-------------------------------------------------------------------------------
!> Check if the array of integer is already allocated
!> @param func_name : function name
!> @param fname : file name
!-------------------------------------------------------------------------------
subroutine check_allocated_int(cur, dname, func_name, fname)
  integer, allocatable, intent(in) :: cur(:)
  character(*) :: dname, func_name, fname

  if (allocated(cur)) then
    call print_error(dname//" already allocated", func_name, fname)
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Check if the array of integer is already allocated
!> @param func_name : function name
!> @param fname : file name
!-------------------------------------------------------------------------------
subroutine check_allocated_double(cur, dname, func_name, fname)
  double precision, allocatable, intent(in) :: cur(:)
  character(*) :: dname, func_name, fname

  if (allocated(cur)) then
    call print_error(dname//" already allocated", func_name, fname)
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Check if 2d integer array is already allocated
!-------------------------------------------------------------------------------
subroutine check_allocated_double2d(cur, dname, func_name, fname)
  double precision, allocatable, intent(in) :: cur(:,:)
  character(*) :: dname, func_name, fname

  if (allocated(cur)) then
    call print_error(dname//" already allocated", func_name, fname)
  end if
end subroutine



!-------------------------------------------------------------------------------
!> Check if there is an error
!> @param ierr : error code to be checked
!> @param mesg : message to be printed
!> @param func_name : function name
!> @param fname : file name
!-------------------------------------------------------------------------------
subroutine check_error(ierr, mesg, func_name, fname)
  integer, intent(in) :: ierr
  character(*) :: mesg, func_name, fname

  if (ierr /= 0) then
    call print_error(mesg, func_name, fname)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Get address for MPI of double precision
!-----------------------------------------------------------------------------
subroutine get_address_double(val, loc, dname, fname)
  double precision :: val
  character(*) :: dname, fname
  integer(kind=mpi_address_kind) :: loc
  integer :: ierr

  call mpi_get_address(val, loc, ierr)
  call check_mpi_error(ierr, "get address", dname, fname)

end subroutine


!-----------------------------------------------------------------------------
!> Get address for MPI of double precision array
!-----------------------------------------------------------------------------
subroutine get_address_double_array(val, loc, dname, fname)
  double precision, allocatable :: val(:)
  character(*) :: dname, fname
  integer(kind=mpi_address_kind) :: loc
  integer :: ierr

  call mpi_get_address(val, loc, ierr)
  call check_mpi_error(ierr, "get address", dname, fname)

end subroutine

!-----------------------------------------------------------------------------
!> Get address for MPI of real
!-----------------------------------------------------------------------------
subroutine get_address_real(val, loc, dname, fname)
  real :: val
  character(*) :: dname, fname
  integer(kind=mpi_address_kind) :: loc
  integer :: ierr

  call mpi_get_address(val, loc, ierr)
  call check_mpi_error(ierr, "get address", dname, fname)

end subroutine


!-----------------------------------------------------------------------------
!> Get address for MPI of string
!-----------------------------------------------------------------------------
subroutine get_address_char(val, loc, dname, fname)
  character(*) :: val
  character(*) :: dname, fname
  integer(kind=mpi_address_kind) :: loc
  integer :: ierr

  call mpi_get_address(val, loc, ierr)
  call check_mpi_error(ierr, "get address", dname, fname)

end subroutine

!-----------------------------------------------------------------------------
!> Get address for MPI of int
!-----------------------------------------------------------------------------
subroutine get_address_int(val, loc, dname, fname)
  integer :: val
  character(*) :: dname, fname
  integer(kind=mpi_address_kind) :: loc
  integer :: ierr

  call mpi_get_address(val, loc, ierr)
  call check_mpi_error(ierr, "get address", dname, fname)

end subroutine

subroutine get_address_long_int(val, loc, dname, fname)
  integer(kind=8) :: val
  character(*) :: dname, fname
  integer(kind=mpi_address_kind) :: loc
  integer :: ierr

  call mpi_get_address(val, loc, ierr)
  call check_mpi_error(ierr, "get address", dname, fname)

end subroutine

!> Get address for MPI of int array
subroutine get_address_int_array(val, loc, dname, fname)
  integer :: val(:)
  character(*) :: dname, fname
  integer(kind=mpi_address_kind) :: loc
  integer :: ierr

  call mpi_get_address(val, loc, ierr)
  call check_mpi_error(ierr, "get address", dname, fname)

end subroutine

!> Destructor for array of int
subroutine free_array_int(cur, dname, fname)
  integer, allocatable :: cur(:)
  character(*) :: dname, fname
  integer :: ierr

  if (allocated(cur)) then
    deallocate (cur, stat=ierr)
    call check_error(ierr,"trying to free "//dname,"free_array", fname)
  end if

end subroutine

#ifdef HDF5
!> Destructor for array of int (HDF5 format). Avoid to link HDF5 to much of the
!> code
subroutine free_array_int_hdf5(cur, dname, fname)
  integer(hsize_t), allocatable :: cur(:)
  character(*) :: dname, fname
  integer :: ierr

  if (allocated(cur)) then
    deallocate (cur, stat=ierr)
    call check_error(ierr,"trying to free "//dname,"free_array", fname)
  end if

end subroutine
#endif

!> Destructor for array of double precision
subroutine free_array_double(cur, dname, fname)
  double precision, allocatable :: cur(:)
  character(*) :: dname, fname
  integer :: ierr

  if (allocated(cur)) then
    deallocate (cur, stat=ierr)
    call check_error(ierr,"trying to free "//dname,"free_array", fname)
  end if
end subroutine

subroutine free_array_double2d(cur, dname, fname)
  double precision, allocatable :: cur(:,:)
  character(*) :: dname, fname
  integer :: ierr

  if (allocated(cur)) then
    deallocate (cur, stat=ierr)
    call check_error(ierr,"trying to free "//dname,"free_array", fname)
  end if
end subroutine

!> Destructior for an array of MPI types 
subroutine free_array_mpi_types(cur, dname, fname)
  character(*), intent(in) :: dname, fname
  integer, allocatable :: cur(:)
  integer :: i, ierr

  if (allocated(cur)) then
    do i = 1, size(cur)
      call mpi_type_free(cur(i), ierr)
      call check_error(ierr,"trying to free "//dname,  "free_array_mpi_types", &
        fname)
    end do
    deallocate(cur)
  end if

end subroutine


end module
