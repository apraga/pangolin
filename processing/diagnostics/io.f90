! Read data for pangolin on ASCII format

module IO
use hdf5

implicit none

character(*), parameter :: fformat = "f23.16"
! Needs a folder, an initial data file and a final one
character(80) :: folder, init_file, final_file, out_file ! 
character(200) :: init_file_c, final_file_c ! Final file names
logical :: is_hdf5
integer :: nb_lat2

enum, bind(C) 
  enumerator :: QSTART=1, QEND=2
end enum

contains

! Read filenames from namelist and set them
subroutine get_filenames(fname)
  character(*), intent(in) :: fname
  integer :: n1, n2
  logical :: is_hdf5_1, is_hdf5_2

  namelist /files/ folder, init_file, final_file, out_file, nb_lat2
  open(8,file=fname, status='OLD', recl=80, delim='APOSTROPHE')
  read(8, nml=files)
  close(8)

  init_file_c = adjustr(trim(folder))//adjustr(trim(init_file))
  final_file_c = adjustr(trim(folder))//adjustr(trim(final_file))
  n1 = len(trim(init_file_c))
  n2 = len(trim(final_file_c))

  is_hdf5_1 =  (init_file_c(n1-2:n1) .eq. ".h5")
  is_hdf5_2 =  (final_file_c(n2-2:n2) .eq. ".h5")
  if ((is_hdf5_1 .and. .not. is_hdf5_2) .or. &
    (.not. is_hdf5_2 .and. is_hdf5_2)) then
    print *, "Both files must be on the same format (HDF5 or ascii)"
    call exit(1)
  end if
  is_hdf5 = (is_hdf5_1 .and. is_hdf5_2)
  if (is_hdf5) then
    print *, "HDF5 format"
  else 
    print *, "Ascii format"
  endif

end subroutine

function get_nb_cells() result(nb_cells)
  integer :: nb_cells
  nb_cells = 6*nb_lat2**2
end function

function get_resolution() result(res)
  double precision :: res
  res = 120./(2*nb_lat2-1)
end function

subroutine get_ratio_wrapper(q, qtype)
  double precision, allocatable :: q(:)
  integer :: qtype

  if (qtype == QSTART) then
    if (is_hdf5) then
      call get_ratio_hdf5(q, trim(init_file_c))
    else
      call get_ratio_ascii(q, trim(init_file_c))
    end if
    print *, "init file", trim(init_file_c)

  else if (qtype == QEND) then
    if (is_hdf5) then
    call get_ratio_hdf5(q, trim(final_file_c))
  else
    call get_ratio_ascii(q, trim(final_file_c))
  end if
    print *, "final file", trim(final_file_c)
  else
    print *, "Wrong type"
    call exit(1)
  end if
end subroutine

! Get ratio from file
subroutine get_ratio_ascii(q, fname)
  double precision, allocatable :: q(:)
  character(*), intent(in) :: fname
  integer :: k, io
  integer :: fid=2
  double precision :: a, b, c

  k = 1
  open(unit=fid, file=fname)
  do 
    read(fid, '(3'//fformat//')', iostat=io) a, b, c
    if (io < 0) then
      exit
    else if (io > 0) then
      print *, "Error reading file"
      call exit(1)
    else
      q(k) = a
      k = k + 1
    end if
  end do
  close(fid)
end subroutine

subroutine get_ratio_hdf5(q, fname)
  double precision, allocatable :: q(:)
  character(*), intent(in) :: fname
  integer(hid_t) :: dataspace, file_id, dset_id
  integer :: error
  integer(hsize_t) :: dims(1)

  call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
  ! Open dataset
  call h5dopen_f(file_id, "ratio", dset_id, error)

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, q, dims, error)
  call h5sclose_f(dataspace , error)
  call h5dclose_f(dset_id , error)
  call h5fclose_f(file_id, error)

end subroutine

! Number of lines from initial file
function get_nb_lines() result(nb_lines)
  character(200) :: fname
  integer :: nb_lines
 
  if (is_hdf5) then
    nb_lines = get_nb_lines_hdf5()
  else
    nb_lines = get_nb_lines_ascii()
  end if

end function

function get_nb_lines_hdf5() result(nb_lines)
  character(200) :: fname
  integer :: nb_lines , io
  integer :: fid=2
  logical :: file_e
  double precision :: a, b, c

  call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
  call h5dopen_f(file_id, "lat", dset_id, error)
  call h5dget_space_f(dset_id, dataspace, error)

  call h5sget_simple_extent_dims_f(dataspace, dims_array, maxdims, error)
  nb_lines = dims_array(1)

  call h5sclose_f(dataspace, error)
  call h5dclose_f(dset_id , error)
  call h5fclose_f(file_id, error)

end function


function get_nb_lines_ascii() result(nb_lines)
  character(200) :: fname
  integer :: nb_lines , io
  integer :: fid=2
  logical :: file_e
  double precision :: a, b, c

  nb_lines = 0
  inquire( file=trim(init_file_c), exist=file_e )
  if (.not. file_e) then
    print * , "File "//trim(init_file_c)// " not found."
    call exit(1)
  end if

  open(unit=fid, file=trim(init_file_c), action="read")
  do 
    read(fid, '(3'//fformat//')', iostat=io) a, b, c
    if (io < 0) then
      exit
    else if (io > 0) then
      print *, "Error reading file"
      call exit(1)
    else
      nb_lines = nb_lines + 1
    end if
  end do
  close(fid)
end function


end module

