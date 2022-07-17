!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! PROGRAM: Pangolin
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Main program.
!
!-------------------------------------------------------------------------------

program Pangolin
use Pangolin_run
!use Message
implicit none

call initialize_all()
call start_simulation()
call finalize_all()

contains


end program
