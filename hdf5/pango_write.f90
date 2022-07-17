! ************************************************************
!
!  This example shows how to read and write data to a
!  dataset.  The program first writes integers to a dataset
!  with dataspace dimensions of DIM0xDIM1, then closes the
!  file.  Next, it reopens the file, reads back the data, and
!  outputs it to the screen.
!
!  This file is intended for use with HDF5 Library verion 1.8
!
! ************************************************************

program main

use hdf5

implicit none

character(23), parameter :: filename = "ratio_1_201301010000.h5"
integer :: hdferr
integer :: nb_lat2 = 80
integer(HID_T) :: file, space, dset(3) ! handles
integer(HSIZE_T) :: dims(1)
double precision , allocatable :: lat(:), lon(:), q(:)
integer :: i, j, k, nb_lon
double precision :: coef = 4*atan(1.d0)/180.
double precision :: dlat, dlon
!
! Initialize FORTRAN interface.
!
call h5open_f(hdferr)

dims = 6*nb_lat2**2

allocate (q(dims(1)))
allocate (lat(dims(1)))
allocate (lon(dims(1)))
k = 1
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
    q(k) = cos(2*lat(k)*coef)*cos(4*lon(k)*coef)
    k = k + 1
  end do
end do

! Create a new file (ecrase older version)
call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)
! Create dataspace "space"
call h5screate_simple_f(1, dims, space, hdferr)
! Create dataset for ratio and cells coordinates (array of double)
call h5dcreate_f(file, "ratio", H5T_NATIVE_DOUBLE, space, dset(1), hdferr)
call h5dcreate_f(file, "center_lat", H5T_NATIVE_DOUBLE, space, dset(2), hdferr)
call h5dcreate_f(file, "center_lon", H5T_NATIVE_DOUBLE, space, dset(3), hdferr)

! Write the data to the dataset.
call h5dwrite_f(dset(1), H5T_NATIVE_DOUBLE, q, dims, hdferr)
call h5dwrite_f(dset(2), H5T_NATIVE_DOUBLE, lat, dims, hdferr)
call h5dwrite_f(dset(3), H5T_NATIVE_DOUBLE, lon, dims, hdferr)
! Close and release resources.
do i = 1, size(dset)
  call h5dclose_f(dset(i) , hdferr)
end do
call h5fclose_f(file , hdferr)

!Try to read the file
q = -9999999.
dims = -1
print *, "file", filename
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
call h5dopen_f (file, "ratio", dset(1), hdferr)
call h5dread_f(dset(1), H5T_NATIVE_DOUBLE, q, dims, hdferr)
print *, "dims", dims
!print *, q

call h5dclose_f(dset(1) , hdferr)
call h5fclose_f(file , hdferr)
call h5close_f(hdferr)
end program
