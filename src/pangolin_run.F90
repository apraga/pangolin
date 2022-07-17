!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Pangolin_run
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> A module which starts the simulation
!
! REVISION HISTORY:
!  24 Avril 2012 - Initial Version
!-------------------------------------------------------------------------------
module Pangolin_run
use Configuration_class
use Partitioning_class
use Advection_2d
!use Message
implicit none

type(Partitioning), target, save :: list_parts

contains 

!-------------------------------------------------------------------------------
!> Initialize everything : configuration, initialization, timers
!-------------------------------------------------------------------------------
subroutine initialize_all()
  integer :: rank, nb_procs, ierr, rc
  character(LINE_WIDTH) :: configfile
  integer :: nb_cells, hdferr

  ! Initialize MPI
  call mpi_init(ierr)

  if (ierr .ne. mpi_success) then
    print *,'Error starting MPI program. Terminating.'
    call mpi_abort(mpi_comm_world, rc, ierr)
  end if

  ! MPI errors are managed manually. It is the role of the programmer to manage it.
  call mpi_errhandler_set(mpi_comm_world, mpi_errors_return, ierr);

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)

  if (rank == 0) call open_logfile()

  ! Check command-line arguments
  call check_arguments(configfile)

#ifdef HDF5
  call h5open_f(hdferr)
#endif
  ! Initialization of the partitioning and the grid
  call new_Configuration(rank, configfile)

  call set_parameters()
  if (rank == 0) call print_step("Initializing simulation done")

  call new_Partitioning(list_parts)
  if (rank == 0) call print_step("Partitioning done")

  ! Start run only if we want to
  if (nb_procs == 1) then
    if (get_total_nb_partitions_Configuration() > 1) then
      call print_warning("The sequential version will not run for several"//&
        "partitions")
      NO_RUN = .True.
    end if
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Finalize everything : configuration, initialization, timers
!-------------------------------------------------------------------------------
subroutine finalize_all()
  integer :: rank, ierr, hdferr

  call mpi_comm_rank(mpi_comm_world, rank, ierr) 
  if (rank == 0) call print_step("Cleaning...")

  ! Wait until all process are here (even idle ones)
  call write_profiling(rank)

  call free_Partitioning(list_parts)
  call free_Configuration()
  call free_Global_grid()

#ifdef HDF5
  ! Close hDF5 library
  call h5close_f(hdferr)
#endif

  if (rank == 0) call close_logfile()

  call mpi_finalize(ierr)
end subroutine

!-----------------------------------------------------------------------------
!> Start the simulation (advection, chemistry)
!-----------------------------------------------------------------------------
subroutine start_simulation()
  integer :: rank, i, ierr
  integer :: hdferr

  call mpi_comm_rank(mpi_comm_world, rank, ierr) 


  if (NO_RUN) then
    if (rank == 0) call print_step("Skipping simulation...")
    return
  end if

  if (rank == 0) call print_step("Starting simulation...")
  call start_advection(list_parts, rank)
  if (rank == 0) call print_step("Simulation done")

end subroutine

end module
