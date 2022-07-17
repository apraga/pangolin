!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-----------------------------------------------------------------------------
!
! MODULE: Plot
!
!> @author
!> Alexis Praga
!
! DESCRIPTION: 
!> Write data for later plotting
!> Must run with 1 proc only
!
!-------------------------------------------------------------------------------

program plot

use Initialization_class
use Partition_class

character(*), parameter :: configfile = "config"
integer :: nb_procs, ierr

call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
if (nb_procs > 1) then
  call print_error("Only one proc for the tests", "main", "plotting.f90")
end if

call new_Configuration(0, configfile)
call new_Initialization(0)

call last_band_size("../output/area.dat")

call mpi_finalize(ierr)

contains

!-------------------------------------------------------------------------------
!> Print last band min and max area for plotting. Also print square area
!-------------------------------------------------------------------------------
subroutine last_band_size(fname)
  integer :: nb_parts, i
  character(*) :: fname
  type(Partition), pointer :: cur
  integer :: id, vol_min, vol_max, vol, vol_max_last
  double precision :: ratio

  id = 2
  call get_nb_partitions_Initialization(nb_parts) 

  ! It is the caller responsability to have an empty file the first call
  write (*,*) "Check the file is empty the first time !"
  open(unit=id, file=fname, access="append", action="write")
  write (id, '(i5)', advance = "no") nb_parts

  call get_partition_Initialization(1, cur)
  vol_max = partition_volume(cur)
  vol_min = vol_max 
  vol_max_last = 0
  do i = 2, nb_parts
    call get_partition_Initialization(i, cur)
    vol = partition_volume(cur)

    vol_max = max(vol, vol_max)
    vol_min = min(vol, vol_min)
    if (is_on_last_band(cur)) then
      vol_max_last = max(vol, vol_max_last)
    end if

  end do
  write (*,*) "vol max last band", vol_max, vol_max_last
  ratio= double precision(vol_max)/vol_min
  !write (id, '(f15.3)') ratio
  write (id, '(2i10, f15.3)') vol_max, vol_max_last, ratio
  close(id)
end subroutine

end program
