!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Tests output
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Formatting for the tests output
!
!-------------------------------------------------------------------------------


module Tests_output

contains

subroutine check_nb_fails(nb_fails, nb_parts)
  integer, intent(in) :: nb_fails, nb_parts
  if (nb_fails > 0) then
    write (*,'(i6,a)', advance="no") nb_fails," fails"
  else
    write (*,'(a)', advance="no") "OK"
  end if
  write (*,'(a,i6,a)') " (", nb_parts, " partitions)"

end subroutine

end module


