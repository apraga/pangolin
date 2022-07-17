!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Parameters
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Defines some constants for all the code. There are the only variables in
!> upper case.
!
!-------------------------------------------------------------------------------

module Parameters

! Desactive most checks and run
logical :: NO_RUN = .False.

!> Fixed value
integer :: NB_TRACERS = 1

!> Allow the use of buffers, which are needed for symmetric operation splitting
logical :: SYMMETRIC_SPLIT = .False.

! Format for double precision
character(*), parameter :: DBLE_FORMAT = 'F23.16'
! In scientific notation. X adds a blank space for readability
character(*), parameter :: DBLE_FORMAT_S = '(X,E23.16)'
double precision :: DBLE_PREC
! Maximum line width for ASCII. Must be constant (so a parameter)
integer, parameter :: LINE_WIDTH = 300
! Non-initialized value
double precision, parameter :: UNDEFINED = -9999.

! Earth radius
double precision, parameter :: RADIUS = 6.3172e6
double precision :: PI = 4.D0*atan(1.D0)

! For larger integer needed to store date
integer, parameter :: k12 = selected_int_kind(12)

enum, bind(c)
! Used for reading/writing data
enumerator :: IS_RATIO, IS_ZWINDS, IS_MWINDS, IS_PARTITIONING
! Used for sending/receiving data
enumerator :: IS_SLOPE, IS_GRADIENT
! Used for neighbour relations
enumerator :: NORTH, EAST, SOUTH, WEST
! Used for selecting the proper wind array
enumerator :: IS_PREV, IS_CUR, IS_NEXT
end enum

contains 

subroutine set_parameters()
  call set_precision()
end subroutine

subroutine set_precision()
  double precision :: test
  double precision :: cur_prec = precision(test)
  ! 100 times the precision
  DBLE_PREC = 100.d0*exp(-log(10.)*cur_prec)

end subroutine

end module

