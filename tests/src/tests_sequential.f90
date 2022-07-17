!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Tests_sequential
!
!> @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Contains all the tests.
!> Must run with 1 proc only.
!
!-------------------------------------------------------------------------------

program tests_sequential
use Tests_partitions
use Tests_cells
use Tests_borders
use Pangolin_run


character(*), parameter :: configfile = "config"
integer :: nb_procs, rank, ierr, rc
integer :: nb_fails

call initialize_all()

! Partition tests
nb_fails = run_partitions_tests()
! Cells tests
nb_fails = run_cells_tests() + nb_fails
! Borders tests
nb_fails = run_borders_tests() + nb_fails

call finalize_all()


end program
