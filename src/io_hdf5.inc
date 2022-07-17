!-------------------------------------------------------------------------------
!> Read file in parallel with HDF5 : we set the file access as parallel and use
!> hyperslab to read from it
!-------------------------------------------------------------------------------
subroutine read_ratio_hdf5(grid, fname, tracer) 
  type(Band_Grid) :: grid
  character(*), intent(in):: fname
  character(*), parameter :: func_name = "read_hdf5"
  integer, intent(in) :: tracer
  integer(hsize_t) :: dims(1)
  integer(hid_t) :: plist_id, memspace
  integer :: error, rank, ierr
  integer :: i, i_start, i_end, k
  integer :: first(1), nb_lon(1)
  double precision, pointer :: q(:)
  integer(hid_t) :: file_id, dataspace, dset_id ! handles
  ! Access type, overridden for parallel I/O

  call set_file_access_read(file_id, fname, plist_id)
  ! Open dataset
  call h5dopen_f(file_id, "ratio", dset_id, error)

  call set_dataset_access(plist_id, error)

  ! Custom grid selection to read from
  ! Get dataset's dataspace identifier.
  call h5dget_space_f(dset_id, dataspace, error)

  call check_file_size(fname, dataspace, IS_RATIO)

  call create_ratio_hyperslabs_file(dataspace, grid)
  ! Only read interior cells
  call get_tracer_pointer(grid, q, IS_RATIO, tracer)
  call create_ratio_hyperslabs_mem(memspace, q, dset_id, grid)

  ! Read data with hyperslabs
  dims = size(q)
  q = UNDEFINED

#ifdef PARALLEL_HDF5
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, q, dims, error, memspace, &
    dataspace, xfer_prp = plist_id )
#else
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, q, dims, error, memspace, &
    dataspace)
#endif

  call check_hdf5_error(error, "Fails to read ratio")

  call mpi_comm_rank(mpi_comm_world, rank, ierr) 

  call close_hdf5_read(dataspace, dset_id, plist_id, file_id)

end subroutine

subroutine read_winds_hdf5(part, fname, rtype, array)
  type (Partition), target :: part
  character(*), intent(in):: fname
  integer, intent(in) :: rtype, array

  if (rtype == IS_ZWINDS) then
    call read_zonal_winds_hdf5(part, fname, array)
  else if (rtype == IS_MWINDS) then
    call read_merid_winds_hdf5(part, fname, array)
  else
    call print_error("Wrong wind type", "read_winds_hdf6", fname_list)
  end if

end subroutine

subroutine read_zonal_winds_hdf5(part, fname, array)
  type (Partition), target :: part
  character(*), intent(in):: fname
  integer, intent(in) :: array
  integer(hid_t) :: plist_id, memspace
  integer :: error
  integer(hid_t) :: file_id, dataspace, dset_id ! handles
  double precision, pointer :: winds(:)
  integer(hsize_t) :: dims(1)

  call set_file_access_read(file_id, fname, plist_id)
  !print *, "reading", fname
  call h5dopen_f(file_id, "u", dset_id, error)

  call set_dataset_access(plist_id, error)

  ! Get dataset's dataspace identifier.
  call h5dget_space_f(dset_id, dataspace, error)

  call check_file_size(fname, dataspace, IS_ZWINDS)

  ! Custom grid selection to read from
  call create_zonal_hyperslabs_file(dataspace, part%grid)

  winds => part%grid%zonal_winds
  dims = size(winds)

  ! We need to create a memory dataspace, even if we take all memory (segfault
  ! otherwise)
  call h5screate_simple_f(1, dims, memspace, error)

#ifdef PARALLEL_HDF5
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, winds, dims, error, memspace, &
    dataspace, xfer_prp = plist_id)
#else
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, winds, dims, error, memspace, &
    dataspace)
#endif
  call check_hdf5_error(error, "Fails to read zonal winds")

  ! Close datasets, property list and file
  call close_hdf5_read(dataspace, dset_id, plist_id, file_id)

end subroutine

subroutine read_merid_winds_hdf5(part, fname, array)
  type (Partition), target :: part
  character(*), intent(in):: fname
  integer, intent(in) :: array
  integer(hid_t) :: plist_id, memspace
  integer :: error
  integer(hid_t) :: file_id, dataspace, dset_id ! handles
  double precision, pointer :: winds(:)
  integer(hsize_t) :: dims(1), n1, n2

  call set_file_access_read(file_id, fname, plist_id)
  !print *, "reading", fname
  call h5dopen_f(file_id, "v", dset_id, error)

  call set_dataset_access(plist_id, error)

  ! Get dataset's dataspace identifier.
  call h5dget_space_f(dset_id, dataspace, error)

  call check_file_size(fname, dataspace, IS_MWINDS)

  ! Custom grid selection to read from
  call create_merid_hyperslabs_file(dataspace, part%grid, .True.)

  winds => part%grid%merid_winds
  dims = size(winds)

  ! We need to create a memory dataspace, even if we take all memory (segfault
  ! otherwise)
  call h5screate_simple_f(1, dims, memspace, error)
  !call h5sget_select_npoints_f(memspace, n1, error) 
  !call h5sget_select_npoints_f(dataspace, n2, error) 
  !print *, "nb points", n1, n2

#ifdef PARALLEL_HDF5
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, winds, dims, error, memspace,&
    dataspace, xfer_prp = plist_id)
#else
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, winds, dims, error, memspace,&
    dataspace)
#endif
  call check_hdf5_error(error, "Fails to read merid winds")

  ! Close datasets, property list and file
  call close_hdf5_read(dataspace, dset_id, plist_id, file_id)

end subroutine

subroutine check_hdf5_error(error, msg)
  integer, intent(in) :: error
  character(*), intent(in) :: msg

  if (error < 0) then
    call print_error(msg, "check_hdf5_error", fname_list)
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Close everything opened by hdf5
!-------------------------------------------------------------------------------
subroutine close_hdf5_read(dataspace, dset_id, plist_id, file_id)
  integer(hid_t) :: dataspace, file_id, plist_id, dset_id
  integer :: error

  call h5sclose_f(dataspace , error)
  call h5dclose_f(dset_id , error)
#ifdef PARALLEL_HDF5
  call h5pclose_f(plist_id, error)
#endif
  call h5fclose_f(file_id, error)
end subroutine

!-------------------------------------------------------------------------------
!> Close everything opened by hdf5
!-------------------------------------------------------------------------------
subroutine close_hdf5_write(dataspace, dset_id, plist_id, file_id)
  integer(hid_t) :: dataspace, dset_id(:), plist_id, file_id
  integer :: error, i

  call h5sclose_f(dataspace, error)
  ! Close datasets
  do i = 1, size(dset_id)
    call h5dclose_f(dset_id(i) , error)
  end do
#ifdef PARALLEL_HDF5
  call h5pclose_f(plist_id, error)
#endif
  call h5fclose_f(file_id, error)

end subroutine


!-------------------------------------------------------------------------------
!> Access the file and return the property list for dataset access
!> @param plist_id : property list for data set write
!-------------------------------------------------------------------------------
subroutine set_file_access_read(file_id, fname, plist_id)
  character(*), intent(in):: fname
  integer(hid_t), intent(inout) :: file_id, plist_id
  integer :: error

#ifdef PARALLEL_HDF5
  ! Setup file access property list with parallel I/O access.
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
#endif

  ! Open file for sequential or parallel access
#ifdef PARALLEL_HDF5
  call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
  call h5pclose_f(plist_id, error)
#else
  call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
#endif

#ifdef PARALLEL_HDF5
  ! Create property list for independent read
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
#endif


end subroutine

!-------------------------------------------------------------------------------
!> Check the number of cells in ratio file
!-------------------------------------------------------------------------------
subroutine check_file_size(fname, dataspace, rtype)
  character(*), intent(in):: fname
  character(5) :: ftype
  integer(hid_t) :: dataspace
  integer, intent(in) :: rtype
  integer :: error, nb, nb_ref

  nb = file_size_from_dataspace(dataspace)
  if (rtype == IS_RATIO) then
    ftype = "ratio"
    nb_ref = nb_cells_total()
  else if (rtype == IS_ZWINDS) then
    ftype = "zonal"
    nb_ref = nb_zonal_winds_total()
  else if (rtype == IS_MWINDS) then
    ftype = "merid"
    nb_ref = nb_merid_winds_total()
  else
    call print_error("Unknown type", "check_file_size", fname_list)
  end if

  if (nb /= nb_ref) then
    call print_error("Wrong "//ftype//" size", "check_file_size", fname_list)
  end if
end subroutine

function file_size_from_dataspace(dataspace) result (nb)
  integer(hsize_t) :: dims_array(1), maxdims(1)
  integer(hid_t) :: dataspace
  integer :: nb, error

  call h5sget_simple_extent_dims_f(dataspace, dims_array, maxdims, error)
  nb = dims_array(1)
end function

function file_size(fname) result (nb)
  character(*), intent(in):: fname
  integer(hsize_t) :: dims_array(1), maxdims(1)
  integer(hid_t) :: dataspace, dset_id, file_id
  integer :: nb, error

  call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
  call h5dopen_f(file_id, "lat", dset_id, error)
  call h5dget_space_f(dset_id, dataspace, error)
  nb = file_size_from_dataspace(dataspace)

  call h5sclose_f(dataspace, error)
  call h5dclose_f(dset_id , error)
  call h5fclose_f(file_id, error)
end function

!-------------------------------------------------------------------------------
!> Define a union of selection to read (or write) ratio from the global grid
!> One latitude line is a hyperslab
!-------------------------------------------------------------------------------
subroutine create_ratio_hyperslabs_file(dataspace, grid)
  type(Band_Grid) :: grid
  integer(hid_t), intent(inout) :: dataspace
  integer(hsize_t) :: nb_lon(1), offset(1)
  integer :: error, i, k
  integer :: i_start, i_end

  call interior_lat_indices(i_start, i_end, grid)

  nb_lon = nb_lon_Band_grid(i_start, grid)
  k = 1
  offset = grid%ratio_offset(k)
  call H5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, nb_lon, error)
  k = k + 1

  do i = i_start+1, i_end
    nb_lon = nb_lon_Band_grid(i, grid)
    offset = grid%ratio_offset(k)

    call H5Sselect_hyperslab_f(dataspace, H5S_SELECT_OR_F, offset, nb_lon, &
      error)
    k = k + 1
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Define a union of selection to read (or write) ratio from the global grid
!> One latitude line is a hyperslab
!-------------------------------------------------------------------------------
subroutine create_zonal_hyperslabs_file(dataspace, grid)
  type(Band_Grid) :: grid
  integer(hid_t), intent(inout) :: dataspace
  integer(hsize_t) :: nb_lon(1), offset(1)
  integer :: error, i, k, nb_cells
  integer :: i_start, i_end

  call interior_lat_indices(i_start, i_end, grid)

  nb_lon = nb_lon_Band_grid(i_start, grid)+1
  k = 1
  offset = grid%zonal_offset(k)
  call H5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, nb_lon, error)
  k = k + 1

  do i = i_start+1, i_end
    nb_lon = nb_lon_Band_grid(i, grid)+1
    offset = grid%zonal_offset(k)

    call H5Sselect_hyperslab_f(dataspace, H5S_SELECT_OR_F, offset, nb_lon, &
      error)
    k = k + 1
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Define a union of selection to read (or write) ratio from the global grid
!> One latitude line is a hyperslab
!> @param is_reading : true if we read from file, false if we write to it
!-------------------------------------------------------------------------------
subroutine create_merid_hyperslabs_file(dataspace, grid, is_reading)
  type(Band_Grid) :: grid
  logical, intent(in) :: is_reading
  integer(hid_t), intent(inout) :: dataspace
  integer(hsize_t) :: length(1), offset(1)
  integer :: error, i, k, nb_cells
  integer :: i_start, i_end, total
  integer :: n

  ! Read north interface but do not write it
  if (is_reading) then
    k = 1
  else
    call interior_lat_indices(i_start, i_end, grid)
    k = i_start
  end if

  length = grid%merid_length(k)
  offset = grid%merid_offset(k)
  call H5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, length, error)

  do i = k+1, size(grid%merid_offset)
    length = grid%merid_length(i)
    offset = grid%merid_offset(i)
    !    if (length(1)+offset(1) > 243) then
    !      print *, "merid length", grid%merid_length(i), i
    !      print *, "merid offset", grid%merid_offset(i), i
    !
    !    end if

    call H5Sselect_hyperslab_f(dataspace, H5S_SELECT_OR_F, offset, length, &
      error)
  end do

end subroutine

!> Define a chunk of memory to write from and set dims_src
subroutine create_merid_hyperslabs_mem(memspace, dims_src, grid)
  type(Band_Grid) :: grid
  integer(hid_t), intent(inout) :: memspace
  integer(hsize_t), intent(inout) :: dims_src(1)
  integer(hsize_t) :: length(1), offset(1)
  integer :: error

  ! Only skip the north interface
  offset = grid%merid_length(1)
  length = sum(grid%merid_length) - offset(1)
  !print *, "offset length mem", offset, length, size(grid%merid_offset)
  call H5Sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset, length, error)
  dims_src = length

end subroutine



!-------------------------------------------------------------------------------
!> Define a union of selection to read (or write) ratio into memory
!> One latitude line is a hyperslab
!-------------------------------------------------------------------------------
subroutine create_ratio_hyperslabs_mem(memspace, q, dset_id, grid)
  type(Band_Grid) :: grid
  integer(hid_t), intent(in) :: dset_id
  integer(hid_t), intent(inout) :: memspace
  integer(hsize_t) :: nb_lon(1), dims_mem(1)
  integer :: error, i, sum_prev
  integer :: i_start, i_end, rank, ierr
  integer(hsize_t) :: offset(1)
  double precision, pointer, intent(in) :: q(:)

  ! Create memory dataspace.
  dims_mem = size(q)
  call h5screate_simple_f(1, dims_mem, memspace, error)

  call mpi_comm_rank(mpi_comm_world, rank, ierr) 

  call interior_lat_indices(i_start, i_end, grid)
  nb_lon = nb_lon_Band_grid(i_start, grid)
  sum_prev = get_nb_ghosts_north(grid)
  offset = sum_prev + get_nb_ghosts_west(i_start, grid)

  ! Select hyperslab in memory
  call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset, nb_lon, &
  error)
!  if (rank == 5) then
!    print *, "nb ghost", sum_prev
!    print *, "i=", i_start, "offset, nb lon", offset, nb_lon, "total", grid%nb_lon(i_start), &
!    "east", get_nb_ghosts_west(i_start, grid)
!    end if
  sum_prev = sum_prev + grid%nb_lon(i_start)

  do i = i_start+1, i_end
    nb_lon = nb_lon_Band_grid(i, grid)
    offset = sum_prev + get_nb_ghosts_west(i, grid)
    call H5Sselect_hyperslab_f(memspace, H5S_SELECT_OR_F, offset, nb_lon, &
      error)
!  if (rank == 5) then
!    print *, "i=", i, "offset, nb lon", offset, nb_lon, "total", grid%nb_lon(i), &
!    "east", get_nb_ghosts_west(i, grid)
!  end if

    sum_prev = sum_prev + grid%nb_lon(i)
  end do
end subroutine

subroutine write_data_hdf5(dtype, iter, this)
  type(Partitioning) :: this
  integer, intent(in) :: dtype, iter
  character(line_width) :: filename
  integer :: output, id, i
  logical :: write_id
  integer :: nb, nb_tracer, tracer

  nb = get_nb_parts_Partitioning(this)
  do i=1, nb

    ! Write all tracers
    if (dtype == IS_RATIO) then
      do id = 1, nb_tracers
        filename = get_output_filename(dtype, iter, id)
        call write_ratio_hdf5(this%parts(i), filename, id)
      end do

    else  if (dtype == IS_ZWINDS) then
      filename = get_output_filename(dtype, iter)
      call write_zonal_winds_hdf5(this%parts(i), filename)

    else if (dtype == IS_MWINDS) then
      filename = get_output_filename(dtype, iter)
      call write_merid_winds_hdf5(this%parts(i), filename)
    else if (dtype == IS_PARTITIONING) then
      print *, "not implemented"
    end if
  end do

end subroutine

!-------------------------------------------------------------------------------
!> Create the file 
!-------------------------------------------------------------------------------
subroutine set_file_access_write(fname, file_id, plist_id)
  character(*), intent(in):: fname
  integer(hid_t), intent(inout) :: file_id, plist_id
  integer :: error

#ifdef PARALLEL_HDF5
  ! Setup file access property list with parallel I/O access.
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! Create the file for both case (parallel and sequntial)
  call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, error,&
    access_prp = plist_id)
  call h5pclose_f(plist_id, error)

#else
  print *, "trying to create", trim(fname)
  call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, error)
#endif

end subroutine

subroutine set_dataset_access(plist_id, error)
  integer(hid_t), intent(inout) :: plist_id
  integer, intent(inout) :: error

#ifdef PARALLEL_HDF5
  ! Create property list for independent write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
#endif
end subroutine

!-------------------------------------------------------------------------------
!> Write file in parallel with HDF5 : we set the file access as parallel and use
!> hyperslab to write in it
!-------------------------------------------------------------------------------
subroutine write_ratio_hdf5(part, fname, tracer) 
  type(Partition), target :: part
  type(Band_Grid), pointer :: grid
  character(*), intent(in) :: fname
  integer, intent(in) :: tracer
  integer :: error, i
  integer(hsize_t) :: dims_dest(1), dims_src(1)
  integer(hid_t) :: plist_id 
  integer(hid_t) :: dset_id(3) ! handles
  integer(hid_t) :: file_id, dataspace
  integer(hid_t) :: memspace, memspace_coord
  double precision, pointer :: q(:), lat(:), lon(:)
  integer :: rank, ierr


  call set_file_access_write(fname, file_id, plist_id)

  ! Get the pointers
  grid => part%grid
  call get_tracer_pointer(grid, q, IS_RATIO, tracer)
  lat => grid%center_lat
  lon => grid%center_lon
  dims_dest = nb_cells_total()
  dims_src = size(lat)

  ! One dataspace for the file, one for the coordinates in memory
  call h5screate_simple_f(1, dims_dest, dataspace, error)
  call h5screate_simple_f(1, dims_src, memspace_coord, error)

  ! Create dataset for ratio and cells coordinates (array of double)
  call h5dcreate_f(file_id, "ratio", H5T_NATIVE_DOUBLE, dataspace, dset_id(1), &
    error)
  call h5dcreate_f(file_id, "lat", H5T_NATIVE_DOUBLE, dataspace, dset_id(2), &
    error)
  call h5dcreate_f(file_id, "lon", H5T_NATIVE_DOUBLE, dataspace, dset_id(3), &
    error)

  call set_dataset_access(plist_id, error)

  ! Custom grid selection to write to
  call create_ratio_hyperslabs_file(dataspace, grid)
  ! Custom memory selection for ratio to read from
  call create_ratio_hyperslabs_mem(memspace, q, dset_id(1), grid)

  ! Write the data to the dataset by hyperslabs
  ! Works for both parallel and sequential I/O
  call hdf5_dataset_write(dset_id(1), q, dims_src, error, memspace, &
    dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write q in ratio file")
  call hdf5_dataset_write(dset_id(2), lat, dims_src , error, &
    memspace_coord, dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write lat in ratio file")
  call hdf5_dataset_write(dset_id(3), lon, dims_src, error, &
    memspace_coord, dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write lon in ratio file")

  ! Close everything
  call close_hdf5_write(dataspace, dset_id, plist_id, file_id)

end subroutine

!-------------------------------------------------------------------------------
!> Write file in parallel with HDF5 : we set the file access as parallel and use
!> hyperslab to write in it
!-------------------------------------------------------------------------------
subroutine write_zonal_winds_hdf5(part, fname) 
  type(Partition), target :: part
  type(Band_Grid), pointer :: grid
  character(*), intent(in) :: fname
  integer :: error, i
  integer(hsize_t) :: dims_src(1), dims_dest(1)
  integer(hid_t) :: plist_id 
  integer(hid_t) :: dset_id(3) ! handles
  integer(hid_t) :: file_id, dataspace
  integer(hid_t) :: memspace, memspace_coord
  double precision, pointer :: u(:), lat(:), lon(:)
  integer :: rank, ierr

  call set_file_access_write(fname, file_id, plist_id)

  ! Get the pointers
  grid => part%grid
  call get_winds_pointer(grid, u, IS_ZWINDS, IS_CUR)
  lat => grid%u_lat
  lon => grid%u_lon
  dims_src = size(u)
  dims_dest = nb_zonal_winds_total()

  ! One dataspace for the file, one for the coordinates in memory
  call h5screate_simple_f(1, dims_dest, dataspace, error)
  !call h5screate_simple_f(1, dims_src, memspace_coord, error)

  ! Create dataset for ratio and cells coordinates (array of double)
  call h5dcreate_f(file_id, "u", H5T_NATIVE_DOUBLE, dataspace, dset_id(1), &
    error)
  call h5dcreate_f(file_id, "lat", H5T_NATIVE_DOUBLE, dataspace, dset_id(2), &
    error)
  call h5dcreate_f(file_id, "lon", H5T_NATIVE_DOUBLE, dataspace, dset_id(3), &
    error)

  call set_dataset_access(plist_id, error)

  ! Custom grid selection to write to
  call create_zonal_hyperslabs_file(dataspace, grid)
  ! We need to create a memory dataspace, even if we take all memory (segfault
  ! otherwise)
  call h5screate_simple_f(1, dims_src, memspace, error)

  ! Write the data to the dataset by hyperslabs
  call hdf5_dataset_write(dset_id(1), u, dims_src, error, memspace, &
    dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write u in zonal winds file")
  call hdf5_dataset_write(dset_id(2), lat, dims_src , error,  &
    memspace, dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write lat in zonal winds file")
  call hdf5_dataset_write(dset_id(3), lon, dims_src, error, &
    memspace, dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write lon in zonal winds file")

  ! Close everything
  call close_hdf5_write(dataspace, dset_id, plist_id, file_id)

end subroutine

!-------------------------------------------------------------------------------
!> Write file in parallel with HDF5 : we set the file access as parallel and use
!> hyperslab to write in it
!-------------------------------------------------------------------------------
subroutine write_merid_winds_hdf5(part, fname) 
  type(Partition), target :: part
  type(Band_Grid), pointer :: grid
  character(*), intent(in) :: fname
  integer :: error, i
  integer(hsize_t) :: dims_mem(1), dims_dest(1), dims_src(1)
  integer(hid_t) :: plist_id 
  integer(hid_t) :: dset_id(3) ! handles
  integer(hid_t) :: file_id, dataspace
  integer(hid_t) :: memspace, memspace_coord
  double precision, pointer :: v(:), lat(:), lon(:)
  integer :: rank, ierr
  integer(hsize_t) :: n1, n2

  call set_file_access_write(fname, file_id, plist_id)

  ! Get the pointers
  grid => part%grid
  call get_winds_pointer(grid, v, IS_MWINDS, IS_CUR)
  lat => grid%v_lat
  lon => grid%v_lon
  dims_mem = size(v)
  dims_dest = nb_merid_winds_total()

  ! One dataspace for the file, one for the coordinates in memory
  call h5screate_simple_f(1, dims_dest, dataspace, error)

  ! Create dataset for ratio and cells coordinates (array of double)
  call h5dcreate_f(file_id, "v", H5T_NATIVE_DOUBLE, dataspace, dset_id(1), &
    error)
  call h5dcreate_f(file_id, "lat", H5T_NATIVE_DOUBLE, dataspace, dset_id(2), &
    error)
  call h5dcreate_f(file_id, "lon", H5T_NATIVE_DOUBLE, dataspace, dset_id(3), &
    error)

  call set_dataset_access(plist_id, error)

  ! Custom grid selection to write to
  call create_merid_hyperslabs_file(dataspace, grid, .False.)

  call h5screate_simple_f(1, dims_mem, memspace, error)
  ! We need to create a memory dataspace to avoid writing winds at the north
  ! interface
  if (has_north_ghost_cells(grid)) then
    call create_merid_hyperslabs_mem(memspace, dims_src, grid)
  end if

  !call h5sget_select_npoints_f(memspace, n1, error) 
  !call h5sget_select_npoints_f(dataspace, n2, error) 
 
  ! Write the data to the dataset by hyperslabs
  call hdf5_dataset_write(dset_id(1), v, dims_src, error, memspace, &
    dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write v in merid winds file")

  call hdf5_dataset_write(dset_id(2), lat, dims_src , error,  &
    memspace, dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write lat in merid winds file")

  call hdf5_dataset_write(dset_id(3), lon, dims_src, error, &
    memspace, dataspace, plist_id)
  call check_hdf5_error(error, "Fails to write lon in merid winds file")

  ! Close everything
  call close_hdf5_write(dataspace, dset_id, plist_id, file_id)

end subroutine

! Wrapper for parallel or sequential dataset write
! Write only double
! Return the error code as argument
subroutine hdf5_dataset_write(dset_id, src, dims_src, error, memspace, &
    dataspace, plist_id)
  double precision, intent(in) :: src(:)
  integer, intent(inout) :: error
  integer(hsize_t) :: dims_src(1)
  integer(hid_t), intent(in) :: plist_id, dset_id
  integer(hid_t), intent(in) :: memspace, dataspace
 
#ifdef PARALLEL_HDF5
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, src, dims_src , error,  &
    memspace, dataspace, xfer_prp = plist_id)
#else
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, src, dims_src , error,  &
    memspace, dataspace)
#endif

end subroutine
