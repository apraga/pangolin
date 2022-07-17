!*******************************************************************************
! Partition the reduced grid :
! - write graph and geometry files
! - partition the graph
! - write output
! This is the main program, subroutines are defined elsewhere
!*******************************************************************************
program partitioning
  use grid
  use output
  character(len=*),parameter :: filename='mesh_init'
  integer nb_lat,nb_lat2
  nb_lat2 = 90
  nb_lat = 2*nb_lat2
  n_max = nb_lat

  ! Writing files
  call write_graph(1, filename, n_max, nb_lat, nb_lat2)
  call write_geom(2, filename, n_max, nb_lat, nb_lat2)

  call partition(filename, n_max, nb_lat, nb_lat2)
end program

!*****************************************************************************
! Initialize partitioning strategy
!*****************************************************************************
subroutine init_strategy(stradat,ierr)
  include "scotchf.h"
  character(1000) :: strategy
  integer :: ierr,id_file
  doubleprecision stradat(scotch_graphdim) 

  strategy = "r{sep=f}"
  call scotchfstratinit(stradat(1), ierr)
  call scotchfstratgraphmap(stradat(1), strategy,ierr) 

  ! Save for checking
  id_file = 3
  open(id_file,file="test.strat")
  call scotchfstratsave(stradat(1), fnum(id_file), ierr)
  close(id_file)
end subroutine

!*****************************************************************************
! Partitioning
!*****************************************************************************
subroutine partition(filename,n_max, nb_lat, nb_lat2)
  use grid
  use output
  include "scotchf.h"

  character(len=*) :: filename
  integer :: n_max,nb_lat,nb_lat2

  ! Declaration for Scotch
  doubleprecision grafdat(scotch_graphdim) ! Graph data
  doubleprecision archdat(scotch_graphdim) ! Architecture
  doubleprecision stradat(scotch_graphdim) ! Strategy
  doubleprecision mapdat(scotch_graphdim) ! Mapping structure
  integer :: nb_parts, ierr,retval,id_file,id_file2
  integer,allocatable :: parttab(:)
  integer :: view_stat = 2

  ! Init graph 
  call scotchfgraphinit(grafdat(1), retval)

  ! Load graph
  id_file = 1
  open(id_file,file=filename//".grf")
  call scotchfgraphload(grafdat(1), fnum(id_file), 1, 0, ierr)
  call scotchfgraphcheck(grafdat(1), ierr)
  close(id_file)

  n_cells = 3*nb_cells(n_max, nb_lat, nb_lat2)
  allocate (parttab(n_cells))

  ! Partition
  nb_parts = 72

  ! Init partitionnig strategy
  call init_strategy(stradat, ierr)

  if (view_stat == 1) then
    ! Simple call
    call scotchfgraphpart(grafdat(1), nb_parts, stradat, parttab, ierr)
  else
    ! Use mapping instead (for stat)
    ! Initialization : architecture, mapping
    call scotchfarchinit(archdat(1), ierr)
    call scotchfarchcmplt(archdat(1), nb_parts)

    call scotchfgraphmapinit(grafdat(1), mapdat(1), archdat(1), parttab, ierr)
    !call scotchfgraphmap(grafdat(1), archdat(1), stradat(1),parttab, ierr)
    call scotchfgraphmapcompute(grafdat(1), mapdat(1), stradat(1),ierr)

    id_file = 2
    open(id_file,file=filename//".stat")
    call scotchfgraphmapview(grafdat(1),mapdat(1),fnum(id_file),ierr)
    close(id_file)

    id_file = 3
    open(id_file,file=filename//".map")
    call scotchfgraphmapsave(grafdat(1),mapdat(1),fnum(id_file),ierr)
    close(id_file)

  end if

  ! Write ps
  call write_ps(filename, parttab, n_max, nb_lat, nb_lat2)

  ! Free memory
  call scotchfgraphexit (grafdat(1))
  call scotchfarchexit(archdat(1));
  deallocate (parttab)
end subroutine
