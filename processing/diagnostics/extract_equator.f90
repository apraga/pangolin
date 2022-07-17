! Extract data at the Equator
program extract_equator

implicit none
integer :: nb_lines, nb_lat
integer :: nb_cells
real :: dlat
character(200) :: input, output
character(*), parameter :: fformat = "f23.16"

namelist /nfile/ input, output

! Read namelist
open(8,file='config_extract', status='OLD', recl=80, delim='APOSTROPHE')
read(8, nml=nfile)
close(8)

! Do stuff
nb_lines = get_nb_lines(input)
nb_lat  =sqrt(real(nb_lines)/6)
dlat = 90./nb_lat
call get_band(adjustr(trim(input)), adjustr(trim(output)), dlat)


contains

! Extract the cells at the latitude closest to the Equator (northern 
! hemisphere) to output file
subroutine get_band(input, output, dlat)
  character(*), intent(in) :: input, output
  integer :: k, io
  integer :: fid=2
  integer :: fid2=3
  real :: a, b, c, dlat

  print *, "file", input
  print *, "out", output
  k = 1
  open(unit=fid, file=input)
  open(unit=fid2, file=output, action="write")
  do 
    read(fid, '(3'//fformat//')', iostat=io) a, b, c
    if (io < 0) then
      exit
    else if (io > 0) then
      print *, "Error reading file"
      call exit(1)
    else
      ! Read or count data
      if (b < dlat .and. b > 0) then
        write (fid2, *) a, b, c
      end if
    end if
  end do
  close(fid)
  close(fid2)
end subroutine

function get_nb_lines(fname) result(nb_lines)
  character(*), intent(in) :: fname
  integer :: nb_lines, io
  integer :: fid=2
  logical :: file_e
  real :: a, b, c

  nb_lines = 0
  inquire( file=fname, exist=file_e )
  if (.not. file_e) then
    print * , "File "//fname// " not found."
    call exit(1)
  end if

  open(unit=fid, file=fname, action="read")
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


end program
