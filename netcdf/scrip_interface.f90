!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This routine is the driver for computing the addresses and weights 
!     for interpolating between two grids on a sphere.
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: scrip.f,v 1.6 2001/08/21 21:06:44 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!***********************************************************************

!-------------------------------------------------------------------------------
! This interface is then SCRIP main program rewritten with some added
! subroutines
!-------------------------------------------------------------------------------

module ScripInterface
use kinds_mod                  ! module defining data types
use constants                  ! module for common constants
use iounits                    ! I/O unit manager
use timers                     ! CPU timers
use grids                      ! module with grid information
use remap_vars                 ! common remapping variables
use remap_conservative         ! routines for conservative remap
use remap_distance_weight      ! routines for dist-weight remap
use remap_bilinear             ! routines for bilinear interp
use remap_bicubic              ! routines for bicubic  interp
use remap_write                ! routines for remap output
use netcdf

implicit none

contains 

!> Performs a scrip interpolation between 2 grid files.
! Mappings are written in two files.
subroutine scrip_interpol(grid1_file, grid2_file, interp_file1, interp_file2, &
    map_method, normalize_opt)!, output_opt)

  character (char_len) :: &
    grid1_file,  & ! filename of grid file containing grid1
    grid2_file,   & ! filename of grid file containing grid2
    interp_file1, & ! filename for output remap data (map1)
    interp_file2, & ! filename for output remap data (map2)
    map1_name,    & ! name for mapping from grid1 to grid2
    map2_name,    & ! name for mapping from grid2 to grid1
    map_method,   & ! choice for mapping method
    normalize_opt,& ! option for normalizing weights
    output_opt    ! option for output conventions

  integer (kind=int_kind) ::  nmap          ! number of mappings to compute (1 or 2)
  integer :: n

  luse_grid1_area = .false.
  luse_grid2_area = .false.
  ! 1 for grid1 -> grid2, 2 for grid1 <-> grid2
  num_maps      = 2
  map_type      = 1
  map1_name     = 'unknown'
  map2_name     = 'unknown'
  output_opt    = 'scrip'
  restrict_type = 'latitude'
  num_srch_bins = 90

  !-----------------------------------------------------------------------
  !
  !     initialize timers
  !
  !-----------------------------------------------------------------------

  call timers_init
  do n=1,max_timers
    call timer_clear(n)
  end do

  select case(map_method)
  case ('conservative')
    map_type = map_type_conserv
    luse_grid_centers = .false.
  case ('bilinear')
    map_type = map_type_bilinear
    luse_grid_centers = .true.
  case ('bicubic')
    map_type = map_type_bicubic
    luse_grid_centers = .true.
  case ('distwgt')
    map_type = map_type_distwgt
    luse_grid_centers = .true.
  case default
    stop 'unknown mapping method'
  end select

  select case(normalize_opt(1:4))
  case ('none')
    norm_opt = norm_opt_none
  case ('frac')
    norm_opt = norm_opt_frcarea
  case ('dest')
    norm_opt = norm_opt_dstarea
  case default
    stop 'unknown normalization option'
  end select

  !-----------------------------------------------------------------------
  !
  !     initialize grid information for both grids
  !
  !-----------------------------------------------------------------------

  call grid_init(grid1_file, grid2_file)

  write(stdout, *) ' Computing remappings between: ',grid1_name
  write(stdout, *) '                          and  ',grid2_name

  !-----------------------------------------------------------------------
  !
  !     initialize some remapping variables.
  !
  !-----------------------------------------------------------------------

  call init_remap_vars

  !-----------------------------------------------------------------------
  !
  !     call appropriate interpolation setup routine based on type of
  !     remapping requested.
  !
  !-----------------------------------------------------------------------

  select case(map_type)
  case(map_type_conserv)
    call remap_conserv
  case(map_type_bilinear)
    call remap_bilin
  case(map_type_distwgt)
    call remap_distwgt
  case(map_type_bicubic)
    call remap_bicub
  case default
    stop 'Invalid Map Type'
  end select

  !-----------------------------------------------------------------------
  !
  !     reduce size of remapping arrays and then write remapping info
  !     to a file.
  !
  !-----------------------------------------------------------------------

  if (num_links_map1 /= max_links_map1) then
    call resize_remap_vars(1, num_links_map1-max_links_map1)
  endif
  if ((num_maps > 1) .and. (num_links_map2 /= max_links_map2)) then
    call resize_remap_vars(2, num_links_map2-max_links_map2)
  endif

  call write_remap(map1_name, map2_name, interp_file1, interp_file2, output_opt)

  !-----------------------------------------------------------------------

end subroutine

!-----------------------------------------------------------------------
! Interpolate between two uni dimensional arrays using the mapping stored in the netcdf file
! interp_f1 (created by create_mapping)
!-----------------------------------------------------------------------
subroutine interpolate_array(temp2, temp1, interp_f1, normalize_opt)
  character(char_len), intent(in) :: interp_f1
  character(char_len) :: normalize_opt
  double precision :: temp2(:), temp1(:)
  integer, allocatable :: src_address(:), dst_address(:)
  real, allocatable :: dst_area(:), dst_frac(:)
  double precision, allocatable :: remap_matrix(:,:)
  integer :: num_links, num_wgts
  integer :: n, addr_dest, addr_src

  print *, "Extracting mapping"
  call extract_mapping(src_address, dst_address, dst_area, &
    dst_frac, remap_matrix, num_links,  num_wgts, interp_f1)

  print *, "Interpolating data"
  select case(normalize_opt)
  case ('frac')
    do n = 1, num_links
      addr_dest = dst_address(n)
      addr_src = src_address(n)
      temp2(addr_dest) = temp2(addr_dest) + remap_matrix(1,n)*temp1(addr_src)
    end do
  case ('destarea')
    do n = 1, num_links
      addr_dest = dst_address(n)
      addr_src = src_address(n)
      temp2(addr_dest) = temp2(addr_dest) + &
        remap_matrix(1,n)*temp1(addr_src)/dst_frac(addr_dest)
    end do
  case ('none')
    do n = 1, num_links
      addr_dest = dst_address(n)
      addr_src = src_address(n)
      temp2(addr_dest) = temp2(addr_dest) + &
        remap_matrix(1,n)*temp1(addr_src)/(dst_area(addr_dest)*dst_frac(addr_dest))
      !print *, dst_area(addr_dest)
    end do
  case default
    print *, "unknown remapping",  normalize_opt
  end select
end subroutine


!-----------------------------------------------------------------------
! Extract mappings from netcdf file
! Arrays are allocated here, do not forget to free them afterward
!-----------------------------------------------------------------------
subroutine extract_mapping(src_address, dst_address, dst_area,&
    dst_frac, remap_matrix, num_links, num_wgts, fname)
  character(char_len) :: fname, names(5)
  integer :: num_links, num_wgts
  integer, allocatable :: src_address(:), dst_address(:)
  real, allocatable :: dst_area(:), dst_frac(:)
  integer :: id, ncid, status
  integer :: id_links, id_wgts, i
  double precision, allocatable :: remap_matrix(:,:)

  status = nf90_open(trim(fname), NF90_NOWRITE, ncid)
  call check_netcdf_err(status)

  ! Get dimensions
  status = nf90_inq_dimid(ncid, "num_links", id_links)
  call check_netcdf_err(status)
  status = nf90_inquire_dimension(ncid, id_links, len=num_links)
  call check_netcdf_err(status)

  status = nf90_inq_dimid(ncid, "num_wgts", id_wgts)
  call check_netcdf_err(status)
  status = nf90_inquire_dimension(ncid, id_wgts, len=num_wgts)
  call check_netcdf_err(status)

  ! Get variables
  allocate (src_address(num_links))
  allocate (dst_address(num_links))
  allocate (dst_area(num_links))
  allocate (dst_frac(num_links))
  allocate (remap_matrix(num_wgts, num_links))

  names = (/ &
    "src_address  ", &
    "dst_address  ", &
    "dst_grid_area", &
    "dst_grid_frac", &
    "remap_matrix "/)
  do i = 1, size(names)
    status = nf90_inq_varid(ncid, trim(names(i)), id)
    call check_netcdf_err(status)
    if (i == 1) then
      status = nf90_get_var(ncid, id, src_address)
    else if (i == 2) then
      status = nf90_get_var(ncid, id, dst_address)
    else 
      status = nf90_get_var(ncid, id, remap_matrix)
    end if
    call check_netcdf_err(status)
  end do

end subroutine

!-------------------------------------------------------------------------------
! Small wrapper to check netcdf status
!-------------------------------------------------------------------------------
subroutine check_netcdf_err(status)
  integer :: status
  if (status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    call exit(1)
  end if
end subroutine

end module 
