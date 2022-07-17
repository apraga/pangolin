! Small test to write data according to pangolin format

program write_pango
use netcdf
implicit none
include 'mpif.h'
 
character(*), parameter :: fname = "test_pango.nc"
integer :: status, varids(4)
integer :: ncid, ids(3), i
integer :: dimids(3)
character(80) :: names(2)
integer :: nb_lat2, nb_elts(2)
character(80)  :: dimensions(2)
integer :: rank, nb_procs, ierr

call MPI_Init(ierr)
call MPI_Comm_rank(mpi_comm_world, rank, ierr)
call MPI_Comm_size(mpi_comm_world, nb_procs, ierr)

print *, "Writing ", trim(fname)
status = nf90_create(trim(fname), IOR(NF90_NETCDF4, NF90_MPIIO), ncid, &
!status = nf90_create(trim(fname), IOR(NF90_CLASSIC_MODEL, NF90_MPIIO), ncid, &
  comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)


call check_netcdf_err(status)

! Define metadata
nb_lat2 = 5
dimensions = (/ "nb_cells", "lev     "/)
nb_elts = (/ 6*nb_lat2**2, 1 /)
do i = 1, size(dimensions)
  status = nf90_def_dim(ncid, trim(dimensions(i)), nb_elts(i), ids(i))
  call check_netcdf_err(status)
end do

! Coordinates
names = (/ "center_lat", "center_lon"/)
do i = 1, size(names)
  status = nf90_def_var(ncid, trim(names(i)), NF90_DOUBLE, ids(1), varids(i))
  call check_netcdf_err(status)
end do

! Units
status = nf90_put_att(ncid, varids(1), "units", "degrees_north")
call check_netcdf_err(status)
status = nf90_put_att(ncid, varids(2), "units", "degrees_east")
call check_netcdf_err(status)

! Data
! lon, lat, levels
status = nf90_def_var(ncid, "q", NF90_DOUBLE, ids(1), varids(3))
call check_netcdf_err(status)
!
status = nf90_enddef(ncid)
call check_netcdf_err(status)

call set_data(ncid, varids, nb_lat2, nb_elts(1))

status = nf90_close(ncid)
call check_netcdf_err(status)

call MPI_Finalize(ierr)

contains

subroutine set_data(ncid, varids, nb_lat2, nb_cells)
  integer, intent(in) :: varids(:), ncid
  integer, intent(in) :: nb_lat2, nb_cells
  double precision, allocatable :: lat(:), lon(:), q(:)
  double precision :: dlat, dlon
  integer :: i, j, k, nb_lon
  double precision, parameter :: PI= 3.14159265359
  double precision :: coef = PI/180.
  integer :: rank, ierr
  integer :: start(1), stride(1)
  integer :: k1, k2

  print *, "nb_cells", nb_cells
  allocate(lat(nb_cells))
  allocate(lon(nb_cells))
  allocate(q(nb_cells))

  k = 1
  call MPI_Comm_rank(mpi_comm_world, rank, ierr)
  dlat = 90./nb_lat2
  do i = 1, 2*nb_lat2
    nb_lon = 2*i-1
    if (i > nb_lat2) then
      nb_lon = 2*(2*nb_lat2 - i+1) - 1
    end if
    dlon = 120.d0/nb_lon
    do j = 1, 3*(nb_lon)
      lat(k) = 90. -(i-0.5)*dlat

      lon(k) = (j-0.5)*dlon
      if (lon(k) < 9 .and. lon(k) > 8) then
        print *, "dlon", dlon, "result", lon(k)
      end if
      q(k) = cos(2*lat(k)*coef)*cos(4*lon(k)*coef)
      k = k + 1
    end do
  end do

  status = nf90_put_var(ncid, varids(1), lat, start=(/1/), count=(/nb_cells/))
  call check_netcdf_err(status)
  status = nf90_put_var(ncid, varids(2), lon)
  call check_netcdf_err(status)
  !! dummy levels
  status = nf90_put_var(ncid, varids(3), q)
  call check_netcdf_err(status)

  deallocate(lat)
  deallocate(lon)
  deallocate(q)

end subroutine


subroutine check_netcdf_err(status)
  integer :: status
  if (status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    call exit(1)
  end if
end subroutine


end program
