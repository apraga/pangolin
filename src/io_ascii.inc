!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: IO
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> A set of I/O functions for ASCII format : reading, writing files. Used by
!> the partitioning class (must be included)
!
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> Sequential write for all tracers concentration or winds
!> @param dtype : IS_RATIO, IS_ZWINDS, IS_MWINDS, IS_PARTITIONING
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!> @param tracer : tracer id
!-----------------------------------------------------------------------------
subroutine write_ascii_sequential_wrapper(dtype, iter, this, array)
  type(Partitioning), target :: this
  type(Partition), pointer :: part
  integer, intent(in) :: dtype, iter
  integer, intent(in), optional :: array
  integer :: tracer, k

  do k = 1, get_nb_parts_Partitioning(this) 
    part => this%parts(k)
    if (dtype == IS_RATIO) then
      ! Print all tracer
      do tracer = 1, nb_tracers
        call write_ascii_sequential(dtype, iter, part, array, tracer)
      end do
    else
      call write_ascii_sequential(dtype, iter, part, array)
    end if
  end do

end subroutine

!-----------------------------------------------------------------------------
!> Sequential write for  a tracer concentration or winds. 
!> Naive way : we pass a token allowing to write. 
!> @param dtype : IS_RATIO, IS_ZWINDS, IS_MWINDS, IS_PARTITIONING
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!> @param tracer : tracer id
!-----------------------------------------------------------------------------
subroutine write_ascii_sequential(dtype, iter, part, array, tracer)
  type(Partition) :: part
  integer, intent(in) :: dtype, iter
  integer, intent(in), optional :: array, tracer
  integer :: nb_procs, ierr, rank
  character(line_width) :: fname
  character(line_width) :: fname_tmp

  fname = get_output_filename(dtype, iter, tracer)

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  if (rank == 0) call print_writing_mesg(dtype)

  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  call check_mpi_error(ierr, "get number processes", "write_ascii_sequential", &
    fname_list)

  if (rank == 0) then
    ! Write to a temporary file except in the monoproc version
    fname_tmp = trim(fname)
    if (nb_procs > 1) fname_tmp = trim(fname)//".tmp"
    call read_write_data_file_seq("write", part, "", trim(fname_tmp), dtype, &
      array, tracer)

    if (nb_procs > 1) call send_token(1, "write_ascii_sequential")

    !else if (rank < 3) then
  else
    ! Read from temporary file and write to final
    call receive_token(rank-1, "write_ascii_sequential")
    fname_tmp = trim(fname)//".tmp"

    call read_write_data_file_seq("readwrite", part, trim(fname_tmp), trim(fname),&
      dtype, array, tracer)
    ! Output become temporary file if we are not the last
    if (rank < nb_procs - 1) then
      ! if (rank < 2) then
      ierr = system("cp "//trim(fname)//" "//trim(fname_tmp))
      call send_token(rank+1, "write_ascii_sequential")

    else
      ierr = system("rm -f "//trim(fname_tmp))
    end if
  end if
end subroutine

!!-----------------------------------------------------------------------------
!!> Sequential write for concentration or winds. Naive way : we pass a token allowing
!!> to write. Old version : each process write its own data
!!> @param dtype : IS_RATIO, IS_ZWINDS, IS_MWINDS
!!> @param iter : current iteration
!!-----------------------------------------------------------------------------
!subroutine write_ascii_sequential_old(dtype, iter, distrib, tracer)
!  type(Partitioning) :: distrib
!  integer, intent(in) :: dtype, iter
!  integer, intent(in), optional :: tracer 
!  integer :: nb_procs, ierr
!  integer :: rank, nb
!
!  call mpi_comm_rank(mpi_comm_world, rank, ierr) 
!
!  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
!  call check_mpi_error(ierr, "get number processes", "write_ascii_sequential", &
!    fname_list)
!
!  if (rank == 0) then
!    call write_data(dtype, iter, "sequential", distrib, tracer)
!    if (nb_procs > 1) call send_token(1, "write_ascii_sequential")
!
!  else 
!    call receive_token(rank-1, "write_ascii_sequential")
!
!    call write_data(dtype, iter, "append", distrib, tracer)
!    if (rank < nb_procs - 1) call send_token(rank+1, "write_ascii_sequential")
!  end if
!end subroutine

!-----------------------------------------------------------------------------
!> Each process write all of its information sequentially
!> Write either concentration or winds according to the dtype
!> cell or the neighbours.
!> @param filename : name of the file without the directory
!> @param dtype : IS_RATIO, IS_ZWINDS, IS_MWINDS, IS_PARTITIONING
!-----------------------------------------------------------------------------
subroutine write_data_sequential(dtype, filename, acces, this, tracer)
  type(Partitioning) :: this
  integer, intent(in) :: dtype, tracer
  integer :: output, id, i
  character(*) :: acces
  character(line_width) :: filename
  logical :: write_id
  integer :: nb

  id = 2

  open(unit=id, file=filename, access=acces, action="write")

  nb = get_nb_parts_Partitioning(this)
  do i=1, nb
    if (dtype == IS_RATIO) then
      call write_concentration(this%parts(i), id, tracer)
    else  if (dtype == IS_ZWINDS) then
      call write_zonal_winds(this%parts(i), id)
    else if (dtype == IS_MWINDS) then
      call write_merid_winds(this%parts(i), id)
    else if (dtype == IS_PARTITIONING) then
      call write_partition_id(this%parts(i), id)
    end if
  end do
  close(id)

end subroutine


!-------------------------------------------------------------------------------
!> Send a token for sequential read/write to another process
!-------------------------------------------------------------------------------
subroutine send_token(dest, fname)
  integer, intent(in) :: dest
  character :: dummy
  character(*), intent(in) :: fname
  integer :: ierr, rank

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_ssend(dummy, 1, mpi_character, dest, 0, mpi_comm_world, ierr)
  call check_mpi_error(ierr, "sending token", fname, &
    fname_list)

end subroutine

!-------------------------------------------------------------------------------
!> Receive a token for sequential read/write from another process
!-------------------------------------------------------------------------------
subroutine receive_token(source, fname)
  integer, intent(in) :: source
  character :: dummy
  integer :: status(mpi_status_size)
  character(*), intent(in) :: fname
  integer :: ierr, rank

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_recv(dummy, 1, mpi_character, source, 0, mpi_comm_world, status, ierr)
  call check_mpi_error(ierr, "receiving token", fname, &
    fname_list)

end subroutine

!-----------------------------------------------------------------------------
!> Read a line in the file, and return the value.
!> Assumes we have 3 floats and the first is the value
!-----------------------------------------------------------------------------
function read_value_from_line(id, io) result(val)
  character(line_width) :: string
  integer, intent(in) :: id
  integer :: io, j
  double precision :: tmp(3), val

  read (id, '(a)', iostat=io) string
  val = undefined
  if (io > 0) then
    call print_error("Error in reading file", "read_value_from_line", fname_list)
  else if (io == 0) then
    read (string, '(3'//DBLE_FORMAT//')', iostat=io) (tmp(j),j=1,3)
    val = tmp(1)
  end if
end function

!-----------------------------------------------------------------------------
!> Check file access
!-----------------------------------------------------------------------------
subroutine check_read_write_access(faccess, io, filename)
  character(*), intent(in) :: faccess, filename
  integer, intent(in) :: io
  character(line_width) :: message

  if (faccess /= "read" .and. faccess /= "write" .and. &
    faccess /= "readwrite") then
  call print_error("Wrong argument : should be read or write", &
    "read_data_file", fname_list)
end if

if (io > 0) then
  if (faccess == "read") then
    message = "Cannot open "//filename//" , check the number of file descriptors"
  else
    message = "Write must be sequential"
  end if
  call print_error(trim(message), "read_data_file", fname_list)
end if
end subroutine

!-----------------------------------------------------------------------------
!> Open files for reading or writing data
!-----------------------------------------------------------------------------
subroutine open_data_files_readwrite(faccess, file_in, id_in, file_out, id_out)
  character(*), intent(in) :: file_in, file_out, faccess
  integer :: id_in, id_out 
  integer :: io_in, io_out
  logical :: f_exists

  id_in = int(UNDEFINED)
  id_out = int(UNDEFINED)

  if (len(file_in) > 0) then
    id_in = 2
    inquire(file=file_in, exist=f_exists )
    if (.not. f_exists) then
      call print_error(file_in//" does not exists", "open_data_files_readwrite",&
      fname_list)
    end if

    open(unit=id_in, file=file_in, action="read", iostat=io_in)
    call check_read_write_access(faccess, io_in, file_in)
  end if

  ! If we write, we may need a temporary file for reading already existing data
  if (len(file_out) > 0) then
    id_out = 3
    open(unit=id_out, file=file_out, action="write", iostat=io_out)
    call check_read_write_access(faccess, io_out, file_out)
  end if
end subroutine

! TODO for winds too
subroutine check_nb_lines(file_in, rtype)
  character(*), intent(in) :: file_in
  integer, intent(in) :: rtype
  character(line_width) :: string
  integer :: nb_lines, nb_lat2, io, id

  id = 2
  open(unit=id, file=file_in, action="read")
  nb_lines = 0
  do 
    read (id, '(a)', iostat=io) string
    if (io > 0) then
      call print_error("Error in reading file", "check_nb_lines", fname_list)
    else if (io < 0) then
      exit
    end if
    nb_lines = nb_lines + 1
  end do
  close(id)

  if (rtype == IS_RATIO) then
    nb_lat2 = get_nb_lat2_Configuration()
    if (nb_lines /= 6*nb_lat2*nb_lat2) then
      call print_error("Wrong number of lines in "//file_in, "check_nb_lines",&
        fname_list)
    end if
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Read or write one data type from a file. Strategy : we search the correct
!> position in the global grid and read/write there.
!> It is the responsability of the caller to ensure we can access the file. For
!> writing, it means the file must be locked for 1 process.
!> Assume partitions on a band have the same length !
!> @param rtype : integer defining the type of data to read
!> @param faccess : "read" or "write"
!> @param file_in :  file to be read (if not empty, so in readwrite or read)
!> @param file_out :  file to be write (if not empty so in readwrite or write)
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!> @param tracer : tracer id for reading/writing ratio
!-----------------------------------------------------------------------------
subroutine read_write_data_file_seq(faccess, part, file_in, file_out, rtype, &
    array, tracer)
  type (Partition), target :: part
  type (Band_grid), pointer :: grid
  character(*), intent(in) :: file_in, file_out, faccess
  integer, intent(in) :: rtype
  integer, intent(in), optional :: tracer, array
  integer :: i, j, k, first_lat
  integer :: i_loc, nb_lat_tot
  integer :: id_in, id_out
  integer :: rank, ierr, nb_lat

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  if (rtype == IS_RATIO .and. .not. present(tracer)) then
    call print_error("Need a tracer for writing ratio data", &
    "read_write_data_file_seq", fname_list)
  end if

#ifdef DEBUG
  if (faccess == "read" .and. rank == 0) call check_nb_lines(file_in, rtype)
#endif

  call open_data_files_readwrite(faccess, file_in, id_in, file_out, id_out)

  ! Assume there are 1 partition per processus
  grid => part%grid

  call set_read_write_sections(i_loc, first_lat, nb_lat, faccess, rtype, grid)

  ! Read until we get to the interesting part (in latitude)
  do i = 1, first_lat-1
    call skip_line(id_in, id_out, i, faccess, part, rtype)
  end do

  call read_write_content(id_in, id_out, i_loc, first_lat, nb_lat, &
    part, faccess, rtype, tracer, array)

  ! Skip to the end when writing (as we must create all entries)
  ! TODO optimize it : only the first process need that
  if (faccess /= "read") then
    nb_lat_tot = get_nb_lat_Global_grid()
    if (rtype == IS_MWINDS) nb_lat_tot = nb_lat_tot-1

    do i = first_lat + nb_lat, nb_lat_tot
      call skip_line(id_in, id_out, i, faccess, part, rtype)
    end do
  end if

  if (len(file_in) > 0) close(id_in)
  if (len(file_out) > 0) close(id_out)
end subroutine

!-------------------------------------------------------------------------------
!> Find partition start and end in a file containing global data (winds or ratio)
!> @param i_loc : first latitude indice in the partition
!> @param first_lat : first latitude globally
!> @param nb_lat : nb of latitudes to use
!-------------------------------------------------------------------------------
subroutine set_read_write_sections(i_loc, first_lat, nb_lat, faccess, rtype, grid)
  type (Band_grid) :: grid
  character(*), intent(in) :: faccess
  integer, intent(in) :: rtype
  integer :: i_loc, first_lat, nb_lat
  double precision :: lat_min

  ! We read ghost cells for ratio and meridional winds 
  if (faccess == "read" .and. rtype /= IS_ZWINDS) then
    i_loc = 1
    nb_lat = grid%nb_lat
    lat_min = grid%lat_min

    if (rtype == IS_MWINDS) nb_lat = nb_lat-1
    ! but use only interior cells for reading zonal or writing 
  else
    i_loc = first_interior_lat(grid)
    nb_lat = get_nb_lat_Band_grid(grid)
    lat_min = get_lat_min(grid)

    ! Last interior latitude, except for last latitude on the global grid
    if (rtype == IS_MWINDS) then
      if (.not. has_south_ghost_cells(grid)) nb_lat = nb_lat-1
    end if
  end if

  first_lat = first_global_lat(lat_min, grid)
end subroutine

!--------------------------------------------------------------------------------
!> Skip a complete latitude line in a data file
!--------------------------------------------------------------------------------
subroutine skip_line(id_in, id_out, i, faccess, part, rtype)
  integer, intent(in) :: id_in, id_out, i, rtype
  character(*) :: faccess
  type (Partition) :: part
  integer :: j

  do j = 1, nb_cells_lat(i)
    call skip_cell(id_in, id_out, i, j, faccess, .True., rtype, part) 
  end do
end subroutine

!--------------------------------------------------------------------------------
!> Do the actual readwrite for the partition
!> @param i_loc : position on the grid
!> @param id_in, id_out : file handles for in and out files (can be null)
!> @param first_lat, nb_lat : first lat on the global grid
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!--------------------------------------------------------------------------------
subroutine read_write_content(id_in, id_out, i_loc, first_lat, nb_lat, &
    part, faccess, rtype, tracer, array)
  type (Partition), target :: part
  type (Band_grid), pointer :: grid
  integer, intent(in) :: id_in, id_out, rtype
  character(*), intent(in) :: faccess
  integer, intent(in), optional :: tracer, array
  integer :: i_loc, j_loc
  integer :: i, j, first_lat, nb_lat
  integer :: j_min, j_max, k
  integer :: nb_ghosts_w, nb_ghosts_e
  double precision :: lon_min
  logical :: merid_f, merid_l
  logical :: write_merid

  ! We read east/west ghost cells only for reading meridional winds
  ! The write is done only for interior cells. This suppose there is no hole in
  ! the grid
  merid_f = (rtype == IS_MWINDS .and. .not. is_first_on_band(part) &
    .and.  faccess == "read")
  merid_l = (rtype == IS_MWINDS .and. .not. is_last_on_band(part) &
    .and.  faccess == "read")

  k = 1
  write_merid = (faccess /= "read" .and. rtype == IS_MWINDS) 
  if (write_merid) call skip_winds_north_ghost(k, part)

  grid => part%grid
  ! i is the global latitude, i_loc the latitude in the grid
  do i = first_lat, first_lat + nb_lat -1

    ! Works also for N/S ghost cells
    j_loc = first_interior_lon(i_loc, grid)
    if (merid_f) j_loc = 1

    lon_min = cell_center_lon(i_loc, j_loc, grid)

    j_min = int(lon_min / grid%dlon(i_loc) + 1)
    j_max = j_min + nb_lon_Band_grid(i_loc, grid) - 1

    if (merid_f) j_max = j_max + grid%nb_ghosts_west(i_loc)
    if (merid_l) j_max = j_max + grid%nb_ghosts_east(i_loc)

    ! We skip ghost cells which are linked to interior cells
    if (write_merid) call skip_winds_west_ghost(i_loc, k, part)

    ! Read until we get to the interesting part (in longitude)
    do j = 1, nb_cells_lat(i)
      if (j >= j_min .and. j <= j_max) then
        call read_write_cell(id_in, id_out, i_loc, j_loc, i, j, k, part, &
          faccess, rtype, tracer, array)
        j_loc = j_loc + 1

      else
        call skip_cell(id_in, id_out, i, j, faccess, .False., rtype, part)
      end if
    end do

    ! Skip ghost cells 
    if (write_merid) call skip_winds_east_ghost(i_loc, k, part)

    i_loc = i_loc + 1
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Read the cell data if we are in reading mode, write it if we are in readwrite
!> or write
!> @param i, j : local position on the grid
!> @param i_glob, j_glob : global position on the latitude line
!-------------------------------------------------------------------------------
subroutine read_write_cell(id_in, id_out, i, j, i_glob, j_glob, k, part, &
    faccess, rtype, tracer, array)
  type (Partition) :: part
  character(*), intent(in) :: faccess
  integer, intent(in) :: i, j, i_glob, j_glob, rtype
  integer, intent(in) :: id_in, id_out
  integer, intent(in), optional :: tracer, array
  integer :: k, io
  logical :: skip

  if (faccess == "read") then
    call read_cell(id_in, i, j, i_glob, j_glob, k, part, rtype, tracer, array)
  else 
    call write_cell(id_in, id_out, i, j, j_glob, k, part, tracer, rtype)
  end if

end subroutine

!-----------------------------------------------------------------------------
!> Skip data corresponding to one cell
!> @param i, j : global position on the grid
!> @param line : true if we are skipping a complete line
!-----------------------------------------------------------------------------
subroutine skip_cell(id_in, id_out, i, j, faccess, line, rtype, part) 
  type (Partition) :: part
  character(*), intent(in) :: faccess
  integer, intent(in) :: i, j
  integer, intent(in) :: rtype, id_in, id_out
  logical :: line
  integer :: k, j1, j2

  ! Write dummy data
  if (faccess == "write") then
    call skip_cell_write(id_out, i, j, rtype, part)
    ! Copy data from first file to second
  else if (faccess == "readwrite") then
    call copy_cell(id_in, id_out, i, j, rtype, part)
    ! Just read
  else
    call skip_cell_read(id_in, i, j, line, rtype, part)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Skip data corresponding to one cell (write dummy values)
!> @param i, j : global position on the grid
!-----------------------------------------------------------------------------
subroutine skip_cell_write(id, i, j, rtype, part)
  type (Partition) :: part
  integer, intent(in) :: i, j
  integer, intent(in) :: rtype, id
  integer :: k, j1, j2, io

  if (rtype == IS_RATIO) then
    call write_dummy_value(id, io)

  else if (rtype == IS_ZWINDS) then
    ! Write east and west value on first cell
    if (j == 1) then
      call write_dummy_value(id, io)
    end if

    call write_dummy_value(id, io)

  else
    ! Skip south neighbours
    call south_neighbour_cell_Global(j1, j2, i, j)
    do k = j1, j2
      call write_dummy_value(id, io)
    end do
  end if
end subroutine


subroutine write_dummy_value(id, io)
  integer, intent(in) :: id, io
  character(*), parameter :: fformat = '(3'//DBLE_FORMAT//')'

  write (id, fformat) UNDEFINED, UNDEFINED, UNDEFINED
end subroutine

subroutine copy_line(id_in, id_out)
  integer, intent(in) :: id_in, id_out
  character(line_width) :: line
  integer :: io

  read (id_in, '(a)', iostat=io) line
  write (id_out, '(a)', iostat=io) trim(line)

end subroutine

!-----------------------------------------------------------------------------
!> Copy data from one file to another
!> @param i, j : global position on the grid
!-----------------------------------------------------------------------------
subroutine copy_cell(id_in, id_out, i, j, rtype, part)
  type (Partition) :: part
  integer, intent(in) :: i, j
  integer, intent(in) :: rtype, id_in, id_out
  integer :: k, j1, j2
  double precision :: val

  if (rtype == IS_RATIO) then
    call copy_line(id_in, id_out)

  else if (rtype == IS_ZWINDS) then
    ! Write east and west value on first cell
    if (j == 1) then
      call copy_line(id_in, id_out)
    end if
    call copy_line(id_in, id_out)

  else
    ! Skip south neighbours
    call south_neighbour_cell_Global(j1, j2, i, j)
    do k = j1, j2
      call copy_line(id_in, id_out)
    end do
  end if
end subroutine


!-----------------------------------------------------------------------------
!> Skip data corresponding to one cell (go to the next cell)
!> @param i, j : global position on the grid
!-----------------------------------------------------------------------------
subroutine skip_cell_read(id, i, j, line, rtype, part) 
  type (Partition) :: part
  integer, intent(in) :: i, j
  integer, intent(in) :: rtype, id
  integer :: k, j1, j2, io
  logical :: line
  double precision :: val

  if (rtype == IS_RATIO) then
    val = read_value_from_line(id, io)

  else if (rtype == IS_ZWINDS) then
    ! Only skip west value
    val = read_value_from_line(id, io)
    if (line .and. j == 1)  then
      val = read_value_from_line(id, io)
    end if

  else
    ! Skip south neighbours
    call south_neighbour_cell_Global(j1, j2, i, j)
    do k = j1, j2
      !if (is_interior_cell(i+1, k, part%grid)) then
      val = read_value_from_line(id, io)
      !end if
    end do
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Read cell data and store it
!> @param i, j : local position in the partition
!> @param j_glob : global position in the partition
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!-----------------------------------------------------------------------------
subroutine read_cell(id, i, j, i_glob, j_glob, k_wind, part, rtype, tracer, array)
  type (Partition) :: part
  integer, intent(in) :: i, j, i_glob, j_glob 
  integer, intent(in) :: rtype, tracer, id
  integer, intent(in), optional :: array
  integer :: k_wind
  logical :: skip

  if (rtype == IS_RATIO) then
    call read_cell_ratio(id, i, j, tracer, part%grid)
  else if (rtype == IS_ZWINDS) then
    call read_cell_zonal_winds(id, i, j, k_wind, array, part%grid)
  else
    call read_cell_merid_winds(id, i, j, i_glob, j_glob , k_wind, array, part)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Set cell ratio from file
!-----------------------------------------------------------------------------
subroutine read_cell_ratio(id, i, j, tracer, grid) 
  type(Band_grid) :: grid
  integer, intent(in) :: i, j, tracer, id
  integer :: io
  double precision :: val
  integer :: rank, ierr

  val = read_value_from_line(id, io)
  ! EOF
  if (io < 0) return
  call set_cell_ratio(val, i, j, tracer, grid) 
end subroutine


!-----------------------------------------------------------------------------
!> Set cell zonal winds from file
!> @param i, j : local cell position
!> @param k : current position in the wind array
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!-----------------------------------------------------------------------------
subroutine read_cell_zonal_winds(id, i, j, k, array, grid) 
  type(Band_grid) :: grid
  integer, intent(in) :: i, j, id, array
  integer :: io, k, offset
  double precision :: val

  ! For the first cell read twice (east and west border)
  if (j == first_interior_lon(i, grid)) then
    val = read_value_from_line(id, io)
    if (io < 0) return
    ! Convert from m/s to degree/s
    !val = val*180./(RADIUS*PI)

    call set_zonal_winds(val, k, array, grid) 
    k = k + 1
  end if

  val = read_value_from_line(id, io)
  ! Convert from m/s to degree/s
  !val = val*180./(RADIUS*PI)
  if (io < 0) return
  call set_zonal_winds(val, k, array, grid) 
  k = k + 1
end subroutine

!-----------------------------------------------------------------------------
!> Set cell south meridional winds from file. We need a Partition for neighbour
!> computation.
!> Assume i < nb_lat of the current partition, so the winds can be computed
!> @param i, j : cell position in the partition
!> @param i_glob, j_glob : cell position on the global grid
!> @param k : current position in the wind array
!> @param array : IS_PREV, IS_CUR, IS_NEXT for correct array
!-----------------------------------------------------------------------------
subroutine read_cell_merid_winds(id, i, j, i_glob, j_glob, k, array, part) 
  type(Partition) :: part
  integer, intent(in) :: i, j, id, array
  integer, intent(in) :: i_glob, j_glob
  integer :: io, k, l, p
  integer :: j_start, j_end
  integer :: j_start_loc, j_end_loc
  integer :: first, last
  double precision :: val

  integer :: rank, ierr, nb_procs
  call mpi_comm_rank(mpi_comm_world, rank, ierr)

  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)

  ! Global neighbour as we read all the cell neighbour
  call south_neighbour_cell_Global(j_start, j_end, i_glob, j_glob)
  ! Need a local index for neighbour
  call south_neighbour_cell_Partition(j_start_loc, j_end_loc, i, j, part, i+1)
  p = j_start_loc

  first = first_cell_position_global(i+1, part)
  last = last_cell_position_global(i+1, part)

  do l = j_start, j_end
    val = read_value_from_line(id, io)
    ! Convert from m/s to degree/s
    !val = val*180./(RADIUS*PI)

    if (io < 0) return

    if (l < first .or. l > last) cycle

    ! Only store fluxes to interior cells 
    if (to_or_from_interior_cell(p, i, j, i+1, part%grid)) then
      call set_merid_winds_Partition(val, k, array, part)
      k = k + 1
    end if
    p = p + 1

  end do

end subroutine

!-----------------------------------------------------------------------------
!> Read cell data and store it
!> @param i, j : local position in the partition
!> @param k_wind : current position in the winds array
!-----------------------------------------------------------------------------
subroutine write_cell(id_in, id_out, i, j, j_glob, k, part, tracer, rtype)
  type (Partition) :: part
  integer, intent(in) :: i, j, j_glob, tracer
  integer, intent(in) :: rtype, id_in, id_out
  integer :: j_start, j_end
  integer :: k, l, io
  double precision :: val

  ! For each write in the output file, we must skip a line in the input file
  if (rtype == IS_RATIO) then
    call write_cell_ratio(id_out, i, j, tracer, part%grid)
    if (id_in /= UNDEFINED) val = read_value_from_line(id_in, io)

  else if (rtype == IS_ZWINDS) then
    ! Write only east border, except for the first cell on the latitude line
    ! This means we may have some discontinuites between the partitions, but we
    ! assume we don't
    if (j == first_interior_lon(i, part%grid)) then
      if (j_glob == 1) then
        call write_cell_zonal_wind(id_out, i, j, k, part%grid, is_west=.True.)
        if (id_in /= UNDEFINED) val = read_value_from_line(id_in, io)
      end if
      ! Skip it anyway
      k= k+ 1
    end if

    ! East border
    call write_cell_zonal_wind(id_out, i, j, k, part%grid, is_west=.False.)
    if (id_in /= UNDEFINED) val = read_value_from_line(id_in, io)
    k= k+ 1

  else
    call south_neighbour_cell_Partition(j_start, j_end, i, j, part, i+1)
    do l = j_start, j_end
      ! Only write interior cells, as the other should have been done by other
      ! partitions
      if (is_interior_cell(i, j, part%grid)) then
        call write_cell_merid_wind(id_out, i, j, l, k, part%grid)
        k= k+ 1
        if (id_in /= UNDEFINED) val = read_value_from_line(id_in, io)
      end if
    end do

  end if
end subroutine


!-----------------------------------------------------------------------------
!> Write cell ratio to file
!-----------------------------------------------------------------------------
subroutine write_cell_ratio(id, i, j, tracer, grid)
  type(Band_grid) :: grid
  integer, intent(in) :: i, j, id, tracer
  double precision :: val
  double precision :: lat, lon
  integer :: rank, ierr
  character(*), parameter :: fformat = '(3'//DBLE_FORMAT//')'

  val = get_cell_ratio(i, j, tracer, grid) 
  lat = cell_center_lat(i, grid)
  lon = cell_center_lon(i, j, grid)

  write (id, fformat) val, lat, lon

end subroutine

!-----------------------------------------------------------------------------
!> Write cell zonal wind at the position k
!> @param is_west : true if we print a west border, false for an east border
!-----------------------------------------------------------------------------
subroutine write_cell_zonal_wind(id, i, j, k, grid, is_west)
  type(Band_grid) :: grid
  integer, intent(in) :: i, j, k, id
  logical, intent(in) :: is_west
  double precision :: val
  double precision :: lat, lon
  integer :: rank, ierr
  character(*), parameter :: fformat = '(3'//DBLE_FORMAT//')'

  val = grid%zonal_winds(k)
  ! Convert from degre/s to m/s
  !val = val*(RADIUS*PI)/180.

  lat = cell_center_lat(i, grid)
  if (is_west) then
    lon = cell_west_lon(i, j, grid)
  else
    lon = cell_east_lon(i, j, grid)
  end if

  write (id, fformat) val, lat, lon

end subroutine

!-----------------------------------------------------------------------------
!> Write cell meridional wind between (i,j) and (i+1,j_neighb).
!> @param k : wind position in the array
!> @param j_neighb : cell neighbour at (i+1)
!-----------------------------------------------------------------------------
subroutine write_cell_merid_wind(id, i, j, j_neighb,  k, grid)
  type(Band_grid) :: grid
  integer, intent(in) :: i, j, j_neighb, k, id
  double precision :: val
  double precision :: lat, lon
  integer :: rank, ierr
  character(*), parameter :: fformat = '(3'//DBLE_FORMAT//')'

  val = grid%merid_winds(k)
  ! Convert from degre/s to m/s
  !val = val*(RADIUS*PI)/180.

  lat = cell_south_lat(i, grid)
  lon = cell_interface_middle(i, j, i+1, j_neighb, grid)

  write (id, fformat) val, lat, lon

end subroutine

!-----------------------------------------------------------------------------
!> Write the partitioning: 
!> @param grid : partition file for cells. i-th line is the partition number 
!> of the i-th cell (on the global grid) 
!> @param uwinds : partition file for zonal winds. i-th line is the partition
!> number of the i-th wind (on the global grid) 
!> @param vwinds : partition file for merid winds. Same as above.
!-----------------------------------------------------------------------------
subroutine write_Partitioning()
end subroutine

!-----------------------------------------------------------------------------
!> Write a header with the number of cells in each process
!> @param nb_cells :: array to write
!-----------------------------------------------------------------------------
subroutine write_header(id, nb_cells)
  integer, intent(in) :: id, nb_cells(:)
  integer :: i

  write (id, '(a)', advance="no") "# "
  do i = 1, size(nb_cells)
    write (id, '(i10)', advance="no") nb_cells(i)
  end do
  write (id, *) 

end subroutine

!-----------------------------------------------------------------------------
!> Returns the number of lines in a data file, ignoring blank lines and 
!> comments. Also check the files have 3 floats on each line
!-----------------------------------------------------------------------------
function nb_correct_lines(filename) result(nb)
  integer :: nb, id
  character(*) :: filename
  double precision :: tmp(3)
  character(line_width) :: string
  integer :: io, io2, j
  logical :: ok

  id = 3
  open (unit=id, file=filename, iostat=io)
  if (io > 0) then
    call print_error("Error opening "//filename, "nb_lines", fname_list)
  end if

  nb = 0
  do 
    read (id, '(a)', iostat=io) string
    if (io > 0) then
      call print_error("Error in reading file", "nb_lines", fname_list)
      ! EOF
    else if (io < 0) then
      exit
    end if

    ! Skip comment and blank lines
    ok = scan(string, '#') == 0 .and. string /= ""

    if (ok) then
      read (string, '(3'//DBLE_FORMAT//')', iostat=io2) (tmp(j), j = 1, 3)
      if (io2 > 0) then
        print *, "Line should have 3 floats : ", string
        call print_error("Wrong file format", "nb_lines", fname_list)
      end if
      nb = nb + 1
    end if
  end do
  close(id)

end function

!-----------------------------------------------------------------------------
!> No interpolation. 
!-----------------------------------------------------------------------------
subroutine write_zonal_winds(this, id_file)
  type(Partition), target :: this
  type(Band_grid), pointer :: grid
  integer, intent(in) :: id_file
  integer :: i, j, k
  character(10) :: fformat
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: l, j_n1, j_n2
  double precision :: u, v, lat, lon

  grid => this%grid

  fformat = '(3'//DBLE_FORMAT//')'

  call interior_lat_indices(i_start, i_end, grid)

  k = 1
  ! Write u-component
  do i = i_start, i_end

    call interior_lon_indices(j_start, j_end, grid, i)
    lat = cell_center_lat(i, grid)
    lon = cell_west_lon(i, j_start, grid)
    write (id_file,fformat) grid%zonal_winds(k), lat, lon
    k = k + 1

    do j = j_start, j_end
      lon = cell_east_lon(i, j, grid)
      write (id_file,fformat) grid%zonal_winds(k), lat, lon
      k = k + 1
    end do
  end do

end subroutine

!-----------------------------------------------------------------------------
!> No interpolation. 
!-----------------------------------------------------------------------------
subroutine write_merid_winds(this, id_file)
  type(Partition), target :: this
  type(Band_grid), pointer :: grid
  integer, intent(in) :: id_file
  integer :: i, j, k
  character(10) :: fformat
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: l, j_n1, j_n2
  double precision :: u, v, lat, lon

  grid => this%grid
  fformat = '(3'//DBLE_FORMAT//')'

  call interior_lat_indices(i_start, i_end, grid)

  ! Write v-component
  k = 1
  if (.not. has_south_ghost_cells(grid)) i_end = i_end - 1
  do i = i_start, i_end

    call interior_lon_indices(j_start, j_end, grid, i)
    lat = cell_south_lat(i, grid)
    do j = j_start, j_end
      call south_neighbour_cell_Partition(j_n1, j_n2, i, j, this, i+1)
      if (j_n1 > 0 .and. j_n2 > 0) then
        do l = j_n1, j_n2
          lon = cell_interface_middle(i, j, i+1, l, grid)
          write (id_file,fformat) grid%merid_winds(k), lat, lon
          k = k + 1
        end do
      end if
    end do
  end do
end subroutine

!!-----------------------------------------------------------------------------
!!> Sequential write for times of each process.  Naive way : we pass a 
!!> token allowing to write.
!!> Write in a file total time and advection, in an other computing advection 
!!> time and number of cells
!!> @param t_tot : total time
!!> @param t_simu : simulation time (advection)
!!> @param t_advec : advection time (computing only)
!!> @param nb_cells : total number of cells 
!!-----------------------------------------------------------------------------
!subroutine write_times_sequential(t_tot, t_simu, t_advec, nb_cells, rank)
!  double precision, intent(in) :: t_tot, t_simu, t_advec
!  integer, intent(in) :: rank, nb_cells
!  type(Configuration) :: config
!  character(7) :: nb_procs2
!  character(*), parameter :: fformat = '(i7, 2'//DBLE_FORMAT//')'
!  character(*), parameter :: fformat2 = '(i7, f'//DBLE_FORMAT//', i7)'
!  character(line_width) :: fname, filename
!  character(line_width) :: fname2, filename2
!  integer :: nb_procs, dummy
!  integer :: id_file, id_file2
!  integer :: ierr, status(mpi_status_size)
!
!  ! Write all the times in a single file
!  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
!  call check_mpi_error(ierr, "get number processes", "write_times_sequential",&
!    fname_list)
!
!  write (nb_procs2, "(i7)") nb_procs
!
!  config = get_Configuration()
!  filename = "totaltime_"//trim(adjustl(nb_procs2))
!  fname = trim(config%output_dir)//"/"//trim(filename)
!
!  ! File format
!  id_file = 2
!
!  ! For second writing
!  id_file2 = 3
!  filename2 = "computingtime_"//trim(adjustl(nb_procs2))
!  fname2 = trim(config%output_dir)//"/"//trim(filename2)
!  id_file2 = 3
!
!  if (rank == 0) then
!    call write_wrapper_double(t_tot, t_simu, rank, id_file, fname, fformat,&
!      "sequential")
!    call write_wrapper(t_advec, nb_cells, rank, id_file2, fname2, fformat2, "sequential")
!    if (nb_procs > 1) then
!      call mpi_send(dummy, 1, mpi_integer, 1, 0, mpi_comm_world, ierr)
!      call check_mpi_error(ierr, "send for data writing", "write_times_sequential",&
!        fname_list)
!
!    end if
!  else
!    call mpi_recv(dummy, 1, mpi_integer, rank-1, 0, mpi_comm_world, status, &
!      ierr)
!    call check_mpi_error(ierr, "send for data writing", "write_times_sequential",&
!      fname_list)
!
!    call write_wrapper_double(t_tot, t_simu, rank, id_file, fname, fformat, &
!      "append")
!    call write_wrapper(t_advec, nb_cells, rank, id_file2, fname2, fformat2, "append")
!    if (rank < nb_procs - 1) then
!      call mpi_send(dummy, 1, mpi_integer, rank+1, 0, mpi_comm_world, ierr)
!      call check_mpi_error(ierr, "receive for data writing", "write_times_sequential",&
!        fname_list)
!    end if
!  end if
!end subroutine
!
!!-------------------------------------------------------------------------------
!!> Small subroutine which open sequentally a file, write the rank and 2 doubles 
!!> and closes it.
!!-------------------------------------------------------------------------------
!subroutine write_wrapper_double(val1, val2, rank, id_file, fname, fformat, acces)
!  character(*), intent(in) :: fname, fformat
!  character(*), intent(in) :: acces
!  double precision, intent(in) :: val1, val2
!  integer, intent(in) :: id_file, rank
!
!  open(unit=id_file, file=fname, access=acces, action="write")
!  write (id_file, fformat) rank, val1, val2
!  close (id_file)
!
!end subroutine
!
!!-------------------------------------------------------------------------------
!!> Small subroutine which open sequentally a file, write the rank, 1 double and 
!!> 1 integer and closes it.
!!-------------------------------------------------------------------------------
!subroutine write_wrapper(val1, val2, rank, id_file, fname, fformat, acces)
!  character(*), intent(in) :: fname, fformat
!  character(*), intent(in) :: acces
!  double precision, intent(in) :: val1
!  integer, intent(in) :: val2
!  integer, intent(in) :: id_file, rank
!
!  open(unit=id_file, file=fname, access=acces, action="write")
!  write (id_file, fformat) rank, val1, val2
!  close (id_file)
!
!end subroutine
!
!!-----------------------------------------------------------------------------
!!> Write a value for each partition and the center of the grid
!!-----------------------------------------------------------------------------
!subroutine write_center_data_Partitioning(this, fname, rank)
!  type(Partitioning), target, intent(in):: this
!  type(Partition), pointer :: cur
!  integer, intent(in) :: rank
!  character(*), intent(in) :: fname
!
!  integer :: i, id
!  character(line_width) :: filename, dirtmp
!  id = 3
!
!  dirtmp = get_output_dir()
!  filename = trim(dirtmp)//"/"//trim(fname) 
!  write (*,*) "file",trim(filename)
!
!  ! The master process is the first to write and must not append
!  if (rank == 0) then
!    open(unit=id, file=trim(filename), action="write")
!  else
!    open(unit=id, file=trim(filename), position="append", action="write")
!  end if
!  do i=1, sum(this%nb_parts)
!    cur => this%parts(i)
!    write (id, '(f10.3)', advance="no") cur%grid%ratio(1) 
!    write (id, '(i5)') rank
!  end do
!  close(id)
!end subroutine
!
!
!!-----------------------------------------------------------------------------
!!> Write partitioning for testing : only the neighbours
!!-----------------------------------------------------------------------------
!subroutine write_list_neighbours_id(this, filename)
!  type(Partitioning) :: this
!  type(Partition) :: cur
!  integer :: i, id, k_start, k_end
!  integer :: j, neighb
!  character(4) :: isint = "(i3)"
!  character(*) :: filename
!  character(line_width) :: fname
!  type(Configuration) :: config
!  id = 2
!
!  config = get_Configuration()
!  fname = trim(config%output_dir)//"/"//filename
!
!  write (*,*) filename
!  open(unit=id, file=fname, action="write")
!  do i=1, sum(this%nb_parts)
!    cur = this%parts(i)
!
!    write (id, isint, advance="no") i
!    call west_neighbours(cur, neighb)
!    write (id, isint, advance="no") neighb
!    call east_neighbours(cur, neighb)
!    write (id, isint, advance="no") neighb
!
!    call north_neighbours(cur, k_start, k_end)
!    do j = k_start, k_end
!      write (id, isint, advance="no") neighb
!    end do
!
!    call south_neighbours(cur, k_start, k_end)
!    do j = k_start, k_end
!      write (id, isint, advance="no") neighb
!    end do
!    write (id,*)
!
!  end do
!  close(id)
!end subroutine
!
!
!
!!-----------------------------------------------------------------------------
!!> Write a link between a partition and its neighbours
!!-----------------------------------------------------------------------------
!subroutine write_neighbours_Partitioning(this, id_file)
!  type(Partitioning) :: this
!  type(Partition) :: cur
!  integer :: id
!  integer, intent(in) :: id_file
!  integer :: i, neighb, j
!  integer :: k_start, k_end
!  character(line_width) :: filename, odir
!  id = 3
!
!  odir = get_output_dir()
!  filename = trim(odir)//"/neighbours.dat"
!  open(unit=id,file=filename,action="write")
!  ! Write the links with its neighbours
!
!  do i=1, sum(this%nb_parts)
!    cur = this%parts(i)
!    call write_center_Partitition(cur, id_file)
!
!    call west_neighbours(cur, neighb)
!    call write_neighbour_center(this, neighb,  id_file)
!    call east_neighbours(cur, neighb)
!    call write_neighbour_center(this, neighb, id_file)
!
!    call north_neighbours(cur, k_start, k_end)
!    do j = k_start, k_end
!      call write_neighbour_center(this, j, id_file)
!    end do
!
!    call south_neighbours(cur, k_start, k_end)
!    do j = k_start, k_end
!      call write_neighbour_center(this, j,  id_file)
!    end do
!
!    write (id_file,*)
!  end do
!  close(id)
!end subroutine
!
!!-----------------------------------------------------------------------------
!!> Write the center of a partition's grid from its id
!!-----------------------------------------------------------------------------
!subroutine write_neighbour_center(this, neighb, id_file)
!  type(Partitioning) :: this
!  integer, intent(in) :: id_file, neighb
!
!  if (neighb > sum(this%nb_parts)) then
!    call print_error("Wrong neighbour id : more than number of partitions", &
!      "write_neighbour_center", fname_list)
!  end if
!
!  if (neighb > 0) then
!    write (id_file, '(a)',advance="no") " "
!    call write_center_Partitition(this%parts(neighb), id_file)
!  end if
!end subroutine
!
