!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Output
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> A module for writing out data files
!
! REVISION HISTORY:
!  24 Avril 2012 - Initial Version
!-------------------------------------------------------------------------------
module Output

implicit none

contains

!-------------------------------------------------------------------------------
!> Small subroutine which open sequentally a file, write the rank and 2 doubles 
!> and closes it.
!-------------------------------------------------------------------------------
subroutine write_wrapper_double(val1, val2, rank, id_file, fname, fformat, acces)
  character(*), intent(in) :: fname, fformat
  character(*), intent(in) :: acces
  double precision, intent(in) :: val1, val2
  integer, intent(in) :: id_file, rank

  open(unit=id_file, file=fname, access=acces, action="write")
  write (id_file, fformat) rank, val1, val2
  close (id_file)

end subroutine

!-------------------------------------------------------------------------------
!> Small subroutine which open sequentally a file, write the rank, 1 double and 
!> 1 integer and closes it.
!-------------------------------------------------------------------------------
subroutine write_wrapper(val1, val2, rank, id_file, fname, fformat, acces)
  character(*), intent(in) :: fname, fformat
  character(*), intent(in) :: acces
  double precision, intent(in) :: val1
  integer, intent(in) :: val2
  integer, intent(in) :: id_file, rank

  open(unit=id_file, file=fname, access=acces, action="write")
  write (id_file, fformat) rank, val1, val2
  close (id_file)

end subroutine



end module
