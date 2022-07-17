!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Tests_parallel
!
!> @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Contains all the tests.
!> Can use several processors
!
!-------------------------------------------------------------------------------

program tests_parallel
use Message
use Tests_advection_2d

integer :: nb_fails

call initialize_all()
! 2d advection tests
nb_fails = 0
call run_advection_2d_tests(nb_fails)

call finalize_all()

!contains

end program
