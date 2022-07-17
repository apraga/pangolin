!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Message
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> A module for printing error, warnings etc to stdout.
!> Errors and warnings are duplicated in a log file, written only by the root
!> process. The file is opened at the beginning of the simulation and close at the
!> end.
!
!-------------------------------------------------------------------------------
module Message

#ifdef WITH_IFORT
use Ifport
#endif

implicit none
include "mpif.h"

character(80) :: logfile
integer :: logid = 99
contains

!-----------------------------------------------------------------------------
!> Print message when a step is validated
!-----------------------------------------------------------------------------
subroutine print_step(mesg)
  character(*) :: mesg

  call print_mesg("* "//mesg)
     
end subroutine

!> Print any kind of string to stdout and log
recursive subroutine print_mesg(mesg, mesg_int)
  character(*) :: mesg
  character(8) :: str
  integer, optional :: mesg_int
  integer :: rank, ierr

  if (present(mesg_int)) then
    write (str, '(i8)') mesg_int
    call print_mesg(mesg//str)
    return
  endif

  write (*, '(a)') mesg
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  if (rank == 0) then
    call print_to_log(mesg)
  end if

end subroutine

subroutine set_logfile()
  character(20) :: id
  integer :: pid
  
  pid = getpid()
  write (id, *), pid
  logfile = "pangolin.log_"//adjustl(id)
  print *, "logfile", logfile
end subroutine

! Called only by root
subroutine open_logfile()
  call set_logfile()

  open(unit=logid, file=trim(logfile), action="write")
end subroutine

subroutine close_logfile()
  close(unit=logid)
end subroutine


!-----------------------------------------------------------------------------
!> Print error message and exit
!> @param mesg : message to be printed
!> @param func_name : function name
!> @param fname : file name
!-----------------------------------------------------------------------------
subroutine print_error(mesg, func_name, fname)
  character(*) :: mesg, func_name, fname
  integer :: rank, ierr

  write (*,*) ""
  write (*,'(5a)') "Error in ",func_name," (",fname,") :"
  write (*,*) "    ",mesg

  ! Write a log too
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  if (rank == 0) then
    call print_to_log("")
    call print_to_log("Error in "//func_name//" ("//fname//") :")
    call print_to_log(mesg)
  end if

  ! Nasty stop
  call mpi_abort(ierr)

  ! This hangs for maste-only subroutines
  !call mpi_finalize()
  !call exit(1)
end subroutine

!> Assume the file is opened and the process is root
subroutine print_to_log(mesg)
  character(*) :: mesg
  write (logid, '(a)') mesg
end subroutine

!-----------------------------------------------------------------------------
!> Print a warning
!> @param mesg : message to be printed
!-----------------------------------------------------------------------------
subroutine print_warning(mesg)
  character(*) :: mesg
  integer :: rank, ierr

  write (*,'(2a)') "Warning : ",mesg

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  if (rank == 0) then
    call print_to_log("Warning : "//mesg)
  end if
end subroutine

end module
