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
!> A set of I/O functions : reading, writing files. Used by the partitioning 
!> class (must be included)
!> Mostly an interface to ASCII or netcdf writing
!
!-------------------------------------------------------------------------------

#ifdef HDF5
#include "io_hdf5.inc"
#else
#include "io_ascii.inc"
#endif


subroutine print_writing_mesg(rtype)
  integer, intent(in) :: rtype
  if (rtype == IS_RATIO) then
    call print_mesg("Writing ratio (sequential)")
  else if (rtype == IS_ZWINDS) then
    call print_mesg("Writing zonal winds (sequential)")
  else
    call print_mesg("Writing meridial winds (sequential)")
  end if
end subroutine


subroutine free_Partitioning(this)
  type(Partitioning) :: this
  integer :: k, ierr
  integer :: rank

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call check_mpi_error(ierr, "get number processes", "free_Partitioning", &
    fname_list)

  do k=1, sum(this%nb_parts)
    call free_Partition(this%parts(k))
  end do

  if (allocated(this%parts)) then
    deallocate(this%parts, stat=ierr)
    call check_error(ierr, "Trying to free partitions list",&
      "free_Partitioning", fname_list)
  end if

  if (allocated(this%nb_parts)) then
    deallocate(this%nb_parts, stat=ierr)
    call check_error(ierr, "Trying to free number of partitions",&
      "free_Partitioning", fname_list)
  end if
end subroutine

!-----------------------------------------------------------------------------
!> Each process write all of its information sequentially
!> Write either concentration or winds according to the dtype
!> cell or the neighbours.
!> @param filename : name of the file without the directory
!> @param dtype : IS_RATIO, IS_ZWINDS, IS_MWINDS, IS_PARTITIONING
!-----------------------------------------------------------------------------
subroutine write_data(dtype, iter, this)
  type(Partitioning) :: this
  integer, intent(in) :: dtype, iter
  integer :: output, id, i
  logical :: write_id


#ifdef HDF5
  call write_data_hdf5(dtype, iter, this)
#else
  call write_ascii_sequential_wrapper(dtype, iter, this)
#endif

end subroutine


!-----------------------------------------------------------------------------
!> Read the number of cells in concentration files and deduces the number of
!> latitudes
!-----------------------------------------------------------------------------
function nb_lat2_from_file(iter) result(nb)
  integer :: nb_l, nb, iter
  character(line_width) :: filename

  filename = trim(get_input_filename(IS_RATIO, iter, 1))

#ifdef HDF5
  nb_l = file_size(trim(filename))
#else
  nb_l = nb_correct_lines(trim(filename))
  print *, "nb correct lines" , nb_l
#endif

  nb = nb_lat2_from_nb_cells(nb_l)
end function

!-----------------------------------------------------------------------------
!> Read ratio
!> @param iter :: iteration (start from 0)
!-----------------------------------------------------------------------------
subroutine read_ratio(iter, this)
  type (Partitioning) :: this
  character(line_width) :: fname
  integer, intent(in) :: iter
  integer :: rank, ierr, id
  integer :: i, K

  call mpi_comm_rank(mpi_comm_world, rank, ierr)

  do k = 1, get_nb_parts_Partitioning(this)
    do id = 1, nb_tracers
      fname = trim(get_input_filename(IS_RATIO, iter, id))
      if (rank == 0) then
        call print_mesg("Reading ratio (sequential)"//trim(fname))
      end if

#ifdef HDF5
      call read_ratio_hdf5(this%parts(k)%grid, trim(fname), id) 
#else
      call read_write_data_file_seq("read", this%parts(k), trim(fname),&
       "", IS_RATIO, 0, id)
#endif
    end do
  end do

end subroutine

!> Read winds at two different times for future interpolation
!> Disabled if iter_s = iter_e
subroutine read_winds(iter_s, iter_e, part)
  type (Partition):: part
  character(line_width) :: fname
  integer :: rank, ierr, i
  integer, intent(in) :: iter_s, iter_e

  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  if (rank == 0) then
    if (iter_s == iter_e) then
      call print_mesg("Reading winds, no interpolation")
    else
      call print_mesg("Reading winds, with interpolation")
    end if
  end if

  ! No interpolation
  if (iter_s == iter_e) then
    fname = trim(get_input_filename(IS_ZWINDS, iter_s, 1))
    call read_winds_wrapper(part, trim(fname), IS_ZWINDS, IS_CUR)

    fname = trim(get_input_filename(IS_MWINDS, iter_s, 1))
    call read_winds_wrapper(part, trim(fname), IS_MWINDS, IS_CUR)

  else
    fname = trim(get_input_filename(IS_ZWINDS, iter_s, 1))
    call read_winds_wrapper(part, trim(fname), IS_ZWINDS, IS_PREV)
    fname = trim(get_input_filename(IS_ZWINDS, iter_e, 1))
    call read_winds_wrapper(part, trim(fname), IS_ZWINDS, IS_NEXT)

    fname = trim(get_input_filename(IS_MWINDS, iter_s, 1))
    call read_winds_wrapper(part, trim(fname), IS_MWINDS, IS_PREV)
    fname = trim(get_input_filename(IS_MWINDS, iter_e, 1))
    call read_winds_wrapper(part, trim(fname), IS_MWINDS, IS_NEXT)

    ! Don't forget to set current winds to prev
    call set_winds_to_prev(part)
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Wrapper to read winds either for ascii or HDF5
!> @param rtype : IS_ZWINDS, IS_MWINDS
!> @param array : IS_CUR, IS_PREV, IS_NEXT
!-------------------------------------------------------------------------------
subroutine read_winds_wrapper(part, fname, rtype, array)
  type (Partition) :: part
  character(*), intent(in) :: fname
  integer, intent(in) :: rtype, array
 
#ifdef HDF5
  call read_winds_hdf5(part, trim(fname), rtype, array)
#else
  call read_write_data_file_seq("read", part, trim(fname), "", rtype, array)
#endif

end subroutine

!-------------------------------------------------------------------------------
!> Update winds between t=(iter-1)*dt and t=iter=dt
!> Either analytically, or by reading data. If T_winds is different from
!> T_output, we interpolate linearly
!> A null period means winds are constant
!-------------------------------------------------------------------------------
subroutine update_winds(iter, dt, this)
  type (Partitioning), target :: this
  type (Partition), pointer :: part
  integer :: ierr, iter
  integer :: T, k
  double precision :: dt

  do k = 1, get_nb_parts_Partitioning(this)
    part => this%parts(k)

    if (is_data_analytical()) then
      call cv_winds((iter-1)*dt, iter*dt, part)
    else
      ! Constant winds, only read once
      if (get_T_winds() == 0) then
        if (iter == 1) call read_winds(iter, iter, part)
        return
      end if

      if (must_readwinds(iter)) then
        ! T winds must be a multiple of dt
        T = int(get_T_winds()*60/dt)
        call read_winds(iter, iter+T, part)
      else 
        call interpolate_winds(iter, part)
      end if
    end if
  end do 

  ! Sync processes
  call mpi_barrier(mpi_comm_world, ierr)
end subroutine

subroutine set_winds_to_prev(part)
  type (Partition) :: part

  part%grid%zonal_winds = part%grid%zonal_winds_prev 
  part%grid%merid_winds = part%grid%merid_winds_prev 
end subroutine

!-------------------------------------------------------------------------------
!> Winds are generated at time k*dt and we want them at time (iter-1/2)*dt
!-------------------------------------------------------------------------------
subroutine interpolate_winds(iter, part)
  type (Partition),target :: part
  type (Band_grid), pointer :: grid
  integer, intent(in) :: iter
  double precision :: alpha
  integer :: T, t_prev
  double precision :: dt
  character(10) :: str_t1, str_t2, str_a

  ! Everything is in minutes
  dt = get_dt()
  T = get_T_winds()
  ! Time is actually 1 iteration behind
  t_prev = int((iter-1)*dt/T)*T
  alpha = ((iter-0.5)*dt - t_prev)/T
  grid =>part%grid 

  write (str_t1, '(i10)') t_prev
  write (str_t2, '(i10)') t_prev+T
  write (str_a, '(F10.5)') alpha
  call print_mesg("interpolating between "//str_t1//" and "//str_t2//", alpha"//str_a)

  grid%zonal_winds = (1.-alpha)*grid%zonal_winds_prev +  &
    alpha*grid%zonal_winds_next 
  grid%merid_winds = (1.-alpha)*grid%merid_winds_prev +  &
    alpha*grid%merid_winds_next 

end subroutine

!-------------------------------------------------------------------------------
!> Skip winds on the first ( north ghost) latitude line. Update indice.
!-------------------------------------------------------------------------------
subroutine skip_winds_north_ghost(k, part)
  type (Partition) :: part
  integer :: k, j, l
  integer :: j_start, j_end

  if (.not. has_north_ghost_cells(part%grid)) return

  do j = 1, get_nb_lon(1, part%grid)
    call south_neighbour_cell_Partition(j_start, j_end, 1, j, part, 2)
    do l = j_start, j_end
      if (to_or_from_interior_cell(l, 1, j, 2, part%grid))  then
        k = k + 1
      end if
    end do
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Skip winds on west ghost cells. Update indice.
!-------------------------------------------------------------------------------
subroutine skip_winds_west_ghost(i, k, part)
  type (Partition) :: part
  integer :: k, i, j, l
  integer :: j_start, j_end

  if (is_first_on_band(part)) return
  if (is_north_south_ghost(i, part%grid)) return

  do j = 1, get_nb_ghosts_west(i, part%grid)
    call south_neighbour_cell_Partition(j_start, j_end, i, j, part, i+1)
    do l = j_start, j_end
      if (to_or_from_interior_cell(l, i, j, i+1, part%grid))  then
        k = k + 1
      end if
    end do
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Skip winds on east ghost cells. Update indice.
!-------------------------------------------------------------------------------
subroutine skip_winds_east_ghost(i, k, part)
  type (Partition) :: part
  integer :: k, j, i, nb, nb_lon
  integer :: j_start, j_end, l

  if (is_last_on_band(part)) return
  if (is_north_south_ghost(i, part%grid)) return

  nb = get_nb_ghosts_east(i, part%grid)
  nb_lon = part%grid%nb_lon(i)
  do j = nb_lon - nb+1, nb_lon
    call south_neighbour_cell_Partition(j_start, j_end, i, j, part, i+1)
    do l = j_start, j_end
      if (to_or_from_interior_cell(l, i, j, i+1, part%grid))  then
        k = k + 1
      end if
    end do
  end do
end subroutine

subroutine skip_merid_winds(i, j, k, part)
  type (Partition) :: part
  integer :: k, j, i
  integer :: j_start, j_end

  call south_neighbour_cell_Partition(j_start, j_end, i, j, part, i+1)
  k = k + j_end - j_start + 1

end subroutine

!-----------------------------------------------------------------------------
!> Each process write all of its information sequentially
!-----------------------------------------------------------------------------
subroutine write_partitioning_sequential(dtype, filename, acces, this, tracer)
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
  print *, "Not done"
!  do i=1, nb
!    if (dtype == IS_RATIO) then
!      call write_concentration(this%parts(i), id, tracer)
!    else  if (dtype == IS_ZWINDS) then
!      call write_zonal_winds(this%parts(i), id)
!    else if (dtype == IS_MWINDS) then
!      call write_merid_winds(this%parts(i), id)
!    else if (dtype == IS_PARTITIONING) then
!      call write_partition_id(this%parts(i), id)
!    end if
!  end do
  close(id)

end subroutine



