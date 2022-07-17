!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Profiling
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Set of functions to log different times used in profiling.
!
!-------------------------------------------------------------------------------

module Profiling
  use Parameters
  use Configuration_class
  ! Contains mpi functions
  use Message

  double precision :: times(5) = 0.
  double precision :: tmp(5)
  integer :: ZONAL_TOTAL = 1
  integer :: ZONAL_ADVEC = 2
  integer :: MERID_TOTAL = 3
  integer :: MERID_ADVEC = 4
  integer :: ADVECTION = 5
!  integer :: WAIT_GHOST = 6
!  integer :: WAIT_SLOPE = 7
contains 

  subroutine start_timer(ttype)
    integer, intent(in) :: ttype
        tmp(ttype) = mpi_wtime()
  end subroutine

  ! Appends time difference to array
  subroutine stop_timer(ttype)
    integer, intent(in) :: ttype

    times(ttype) = mpi_wtime() - tmp(ttype) + times(ttype)

  end subroutine

  function get_timer(ttype) result(t)
    double precision :: t
    integer, intent(in) :: ttype

    t = times(ttype)
  end function

  function output_filename(rank, nb_procs) result (fname)
    integer, intent(in) :: rank, nb_procs
    character(LINE_WIDTH) :: fname
    character(8) :: tmp

    write (tmp, '(i8)') nb_procs
    fname = trim(get_output_dir())//"/profiling_"//adjustl(trim(tmp))
    
  end function

  !-------------------------------------------------------------------------------
  !> Gather all data to master proc which write it to a file
  !> May hang with open MPI on some PC
  !-------------------------------------------------------------------------------
  subroutine write_profiling(rank)
    character(LINE_WIDTH) :: fname
    double precision, allocatable :: alltimes(:)
    integer, intent(in) :: rank
    integer :: id, nb_procs
    integer :: i, j, n, ierr

    n = size(times)

    ! FIXME Every process allocate the array to avoid error from Intel MPI debug
    ! mode
    call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
    allocate(alltimes(nb_procs*n))

    call mpi_gather(times, n, mpi_double_precision, alltimes, n, &
      mpi_double_precision, 0, mpi_comm_world, ierr)

!print *, "rank", rank, "times", times
    if (rank == 0) then
      id = 2
      fname = output_filename(rank, nb_procs)
      open(unit=id, file=trim(fname), action="write")

      write (id, '(a)') "# proc zonal_total zonal_advec merid_total merid_advec advection "
      do i = 1, nb_procs
        write (id, '(i5, 6'//DBLE_FORMAT_S//')') i-1, (alltimes((i-1)*n + j), j=1, n)
      end do

      close(id)
    end if

  end subroutine
end module

