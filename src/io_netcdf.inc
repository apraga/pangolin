!-------------------------------------------------------------------------------
!> Read ratio data from the netcdf format.
!> The partition has alread computed the start, count and stride of the data
!> Then the I/O is easy, and possibly parallel
!-------------------------------------------------------------------------------
subroutine read_ratio_netcdf4(grid, fname, tracer) 
  type(Band_Grid) :: grid
  character(*), intent(in):: fname
  character(*), parameter :: func_name = "read_netcdf4"
  character(NF90_MAX_NAME) :: nb_name
  integer, intent(in) :: tracer
  integer :: varid, start
  integer :: nb_cells_id, nb_cells
  integer :: ncid, status
  integer :: i, i_start, i_end, k
  integer :: first(1), nb_lon(1)
  double precision, pointer :: q(:)

  ! Open read only
  status = nf90_open(fname, NF90_NOWRITE, ncid)
  call check_nf90_status(status, func_name)

  ! Get number of cells
  status = nf90_inq_dimid(ncid, "nb_cells", nb_cells_id)
  call check_nf90_status(status, func_name)
  status = nf90_inquire_dimension(ncid, nb_cells_id, nb_name, nb_cells)
  call check_nf90_status(status, func_name)

  if (nb_cells /= nb_cells_total()) then
    call print_error("Wrong number of cells", func_name, fname_list)
  end if
  print *, "nb cells", nb_cells
  allocate (q(nb_cells))

  ! Get the varid of the data variable, based on its name.
  call get_tracer_pointer(grid, q, IS_RATIO, tracer)
  status = nf90_inq_varid(ncid, "q", varid) 
  call check_nf90_status(status, func_name)

  ! Read the data line by line
  k = 1
  start = 1
  call interior_lat_indices(i_start, i_end, grid)
  do i = i_start, i_end
    start = start + get_nb_ghosts_east(i, grid)
    first = (/grid%netcdf_start(k)/)
    nb_lon = (/get_nb_lon(i, grid)/)
    status = nf90_get_var(ncid, varid, q(start:start+nb_lon(1)), first, nb_lon)
    call check_nf90_status(status, func_name)

    ! Go to end of line
    start = start + nb_lon(1) + get_nb_ghosts_east(i, grid)
    k = k+ 1
  end do

  ! Close the file, freeing all resources.
  status = nf90_close(ncid) 
  call check_nf90_status(status, func_name)

end subroutine

!-----------------------------------------------------------------------------
!> Each process write all of its information sequentially
!> Write either concentration or winds according to the dtype
!> cell or the neighbours.
!> @param filename : name of the file without the directory
!> @param dtype : IS_RATIO, IS_ZWINDS, IS_MWINDS, IS_PARTITIONING
!-----------------------------------------------------------------------------
subroutine write_data_netcdf(dtype, iter, this)
  type(Partitioning) :: this
  integer, intent(in) :: dtype, iter
  character(line_width) :: filename
  integer :: output, id, i
  logical :: write_id
  integer :: nb, nb_tracer, tracer

  nb = get_nb_parts_Partitioning(this)
  do i=1, nb
    do id = 1, nb_tracers
      filename = get_output_filename(dtype, iter, id)

      if (dtype == IS_RATIO) then
        call write_ratio_netcdf4(this%parts(i), filename, id)
      else  if (dtype == IS_ZWINDS) then
        print *, "not implemented"
      else if (dtype == IS_MWINDS) then
        print *, "not implemented"
      else if (dtype == IS_PARTITIONING) then
        print *, "not implemented"
      end if
    end do
  end do
  close(id)

end subroutine


!-------------------------------------------------------------------------------
!> Write ratio data from the netcdf format.
!> The partition has alread computed the start, count and stride of the data
!> Then the I/O is easy, and possibly parallel
!-------------------------------------------------------------------------------
subroutine write_ratio_netcdf4(part, fname, tracer) 
  type(Partition), target :: part
  type(Band_Grid), pointer :: grid
  character(*), intent(in):: fname
  character(*), parameter :: func_name = "write_netcdf4"
  character(NF90_MAX_NAME) :: nb_name
  integer, intent(in) :: tracer
  integer :: varid, k, start, start_c
  integer :: i, i_start, i_end
  integer :: nb_cells_id, nb_cells
  integer :: ncid, status
  integer :: ids(2), varids(4)
  double precision, pointer :: q(:), lat(:), lon(:)

  grid => part%grid

  status = nf90_create(trim(fname), IOR(NF90_NETCDF4, NF90_MPIIO), ncid, &
    comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
  call check_nf90_status(status, func_name)

  call write_netcdf_header(ncid, ids, varids)

  ! Data
  ! lon, lat, levels
  status = nf90_def_var(ncid, "q", NF90_DOUBLE, ids(1), varids(3))
  call check_nf90_status(status, func_name)
  !
  status = nf90_enddef(ncid)
  call check_nf90_status(status, func_name)

  ! Get the pointers
  call get_tracer_pointer(part%grid, q, IS_RATIO, tracer)
  lat => part%grid%center_lat
  lon => part%grid%center_lon

  ! Read the data line by line
  k = 1
  start = 1
  start_c = 1
  call interior_lat_indices(i_start, i_end, grid)
  do i = i_start, i_end
    start = start + get_nb_ghosts_east(i, grid)

    call write_netcdf_line(i, k, ncid, varids, lat, lon, q, start, start_c, grid)
  end do

  status = nf90_close(ncid)
  call check_nf90_status(status, func_name)
end subroutine

!-------------------------------------------------------------------------------
!> Write one line of data for cells coordinates and ratio
!> As input, we must have the pointers
!-------------------------------------------------------------------------------
subroutine write_netcdf_line(i, k, ncid, varids, lat, lon, q, start, start_c,&
    grid)
  type (Band_Grid) :: grid
  integer :: first(1), nb_lon(1)
  integer :: i, ncid, status
  integer, intent(in) :: varids(4)
  integer, intent(inout) :: start, start_c, k
  double precision, intent(in) :: q(:), lat(:), lon(:)
  character(*), parameter :: func_name = "write_netcdf_line"

  first = (/grid%netcdf_start(k)/)
  nb_lon = (/get_nb_lon(i, grid)/)

!  if ( i <= 2) then
!    print *, "coords", start_c, start_c + nb_lon(1), i, grid%nb_lon(i)
!  print *, "lon", lon(start_c:start_c + nb_lon(1)-1), "to", first, nb_lon
!end if
    !nb_ghost_cells_lat(i, grid), get_nb_lon(i, grid)
  !print *, "lon", lon(start_c:start_c + nb_lon(1))
  ! Write cells positions
  status = nf90_put_var(ncid, varids(1), lat(start_c:start_c+nb_lon(1)-1), &
    first, nb_lon) 
  call check_nf90_status(status, func_name)
  status = nf90_put_var(ncid, varids(2), lon(start_c:start_c+nb_lon(1)-1), &
    first, nb_lon) 
  call check_nf90_status(status, func_name)

  ! Write ratio
  status = nf90_put_var(ncid, varids(3), q(start:start+nb_lon(1)), first, nb_lon)
  call check_nf90_status(status, func_name)

  ! Go to end of line
  start = start + nb_lon(1) + get_nb_ghosts_east(i, grid)
  start_c = start_c + nb_lon(1)
  k = k+ 1

end subroutine

subroutine write_netcdf_header(ncid, ids, varids)
  integer :: ids(2), varids(4)
  integer :: status, dimids(3)
  integer :: ncid, i, nb_cells
  character(80) :: names(2)
  integer :: nb_lat2, nb_elts(2)
  character(80)  :: dimensions(2)
  character(*), parameter :: func_name = "write_netcdf4_header"

  ! Define metadata
  dimensions = (/ "nb_cells", "lev     "/)
  nb_cells = nb_cells_total()
  nb_elts = (/ nb_cells, 1 /)
  do i = 1, size(dimensions)
    status = nf90_def_dim(ncid, trim(dimensions(i)), nb_elts(i), ids(i))
    call check_nf90_status(status, func_name)
  end do

  ! Coordinates
  names = (/ "center_lat", "center_lon"/)
  do i = 1, size(names)
    status = nf90_def_var(ncid, trim(names(i)), NF90_DOUBLE, ids(1), varids(i))
    call check_nf90_status(status, func_name)
  end do

  ! Units
  status = nf90_put_att(ncid, varids(1), "units", "degrees_north")
  call check_nf90_status(status, func_name)
  status = nf90_put_att(ncid, varids(2), "units", "degrees_east")
  call check_nf90_status(status, func_name)
end subroutine


subroutine check_nf90_status(status, func_name)
  integer, intent ( in) :: status
  character(*), intent(in) :: func_name

  if(status /= nf90_noerr) then 
    call print_error(trim(nf90_strerror(status)), func_name, fname_list)
  end if

end subroutine
