! Compute filament diagnostic for pangolin, using Lauritzen module
! Files are read from a namelist

program pango_diagnostics
use diag
use io
implicit none

integer :: nb_args
character(80):: arg
! Initial and final concentration, wit cells areas
double precision, allocatable :: q0(:), q(:), areas(:)

! Check optional arguments
nb_args = iargc()
!write (*,*) "nb args", nb_args
if (nb_args > 2) then
  print *, "Needs only one argument"
  call exit(1)
end if

call set_arrays()

call getarg(1, arg)
if (trim(arg(1:11)) == "--filaments") then
  call filaments()
else if (trim(arg(1:5)) == "--err") then
  call error()
else
  print *, "unrecognized options",  trim(arg(1:11))
  print *, "Usage : ./pango_diagnostics [--filaments] [--err]"
  call exit(1)
end if

contains 

subroutine error()
  double precision :: l2, loo, l_quad

  if (size(q) /= size(areas)) then
    print *, "Error size differs between areas and q"
    call exit(1)
  end if

  l2 = error_2()
  loo = error_oo()
  l_quad = error_quad()

  print *, "l2, loo, l_quad", l2, loo, l_quad
  print *, "max diff", maxval(abs(q-q0)), "max", maxval(abs(q0))

  open(unit = 2, file=trim(out_file), status='replace')
  write (2, *) "#resolution nb_cells l2 loo lquad"
  write (2, *) get_resolution(), get_nb_cells(), l2, loo, l_quad
  close(2)
end subroutine

function error_oo() result(loo)
  double precision :: loo
  integer :: i

  loo = maxval(abs(q-q0))/maxval(abs(q0))
end function

function error_2() result(l2)
  double precision :: err_up, err_down, l2
  integer :: i

  err_up = 0
  err_down = 0
  do i = 1, size(q)
    err_up = err_up + (q(i) - q0(i))**2*areas(i)
    err_down = err_down +(q0(i))**2*areas(i)
  end do
  l2 = sqrt(err_up /err_down)
end function

! Weighted quadratic error
function error_quad() result(l_quad)
  double precision :: l_quad
  integer :: i

  l_quad = 0
  do i = 1, size(q)
    l_quad = l_quad + (q(i) - q0(i))**2*areas(i)
  end do
end function


! Read concentrations and set cell areas
subroutine set_arrays()
  integer :: nb_lines

  call get_filenames("config")
  nb_lines = get_nb_lines()

  allocate(q0(nb_lines))
  allocate(q(nb_lines))
  allocate(areas(nb_lines))

  call get_ratio_wrapper(q0, 1)
  call generate_areas(areas, nb_lines)
  call get_ratio_wrapper(q, 2)

end subroutine

subroutine free_arrays()
  double precision, allocatable :: q(:), areas(:), q0(:)

  deallocate(q)
  deallocate(q0)
  deallocate(areas)

end subroutine


! Filament diagnostics
subroutine filaments()
  !call filament_diag(K,f1,dA,fila_t0,linit)
  double precision, allocatable :: filts_init(:)
  integer :: nb_lines

  ! Read config from files
  print *, "Filaments diagnostic"

  allocate(filts_init(size(q)))

  call filament_diag(size(q0), q0, areas, filts_init, .True., out_file)
  call filament_diag(size(q), q, areas, filts_init, .False., out_file)

  deallocate(filts_init)
end subroutine

! In degrees
subroutine generate_areas(areas, nb_lines)
  integer :: i, nb_lat2, j, nb_lon, k
  integer, intent(in) :: nb_lines
  double precision, allocatable :: areas(:)
  double precision :: lat, dlat
  double precision :: dlon, cos_lat
  double precision :: PI = 4.D0*atan(1.D0) 

  nb_lat2 = sqrt(nb_lines/6.)
  dlat = 0.5*PI/nb_lat2
  print *, "nblat2", nb_lat2
  k = 1
  do i = 1, 2*nb_lat2
    nb_lon = nb_cells(i, nb_lat2)
    dlon = 2*PI/nb_lon
    do j = 1, nb_lon
      ! Special case at the poles : triangular area
      if (i == 1 .or. i == 2*nb_lat2) then
        areas(k) = 2*PI*(1-cos(dlat))/3.
      ! Exact formula
      else
        areas(k) = dlon*(-cos(i*dlat) + cos((i-1)*dlat))
      end if
      k = k + 1
    end do
  end do
end subroutine

function nb_cells(i, nb_lat2) result(n)
  integer, intent(in) :: i, nb_lat2
  integer :: n, i2

  i2 = i
  if (i > nb_lat2) i2 = 2*nb_lat2 -i + 1
  n = 3*(2*i2-1)
end function


end program
