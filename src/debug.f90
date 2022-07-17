!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! PROGRAM: Debug
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Contains debugging subroutines
!
!-------------------------------------------------------------------------------


module Debug

use, intrinsic :: iso_c_binding, only: c_int

implicit none

interface

  !-----------------------------------------------------------------------------
  ! C binding for sleep function
  !-----------------------------------------------------------------------------
  subroutine sleep_f (seconds) bind ( C, name="sleep" )
    use, intrinsic :: iso_c_binding, only: c_int
    integer (c_int), intent (in), VALUE :: seconds
  end subroutine

  !-----------------------------------------------------------------------------
  ! C binding for getpid function
  !-----------------------------------------------------------------------------
  function getpid_f () bind ( C, name="getpid" )
    use, intrinsic :: iso_c_binding, only: c_int
    integer (c_int) :: getpid_f
  end function

end interface

end module 
