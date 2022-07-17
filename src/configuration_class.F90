!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Simulation config
!
!> @author
!>  Alexis Praga
!
! DESCRIPTION: 
!>  Contains the global simulation parameters, given by the user in configuration 
!>  file. This module is mostly a parser with some checks
!
!-------------------------------------------------------------------------------

module Configuration_class
use Data_check
use Parameters

#ifdef WITH_IFORT
use Ifport
#endif

implicit none


public
type Configuration
  ! Default values are set here
  ! Chemical reations : "yes" for activating (not implemented)
  character(LINE_WIDTH) :: add_chemic = "no"
  ! Advection type : zonal, meridional, 2d
  character(LINE_WIDTH) :: advection = "2d"
  ! CFL condition
  double precision :: CFL = 0.
  ! Timestep in minutes
  double precision :: dt = 0.
  ! Input file for numerical ratio (relative to root directory)
  character(LINE_WIDTH) :: input_dir = ""
  ! Input directory containing all needed winds (relative to root directory)
  character(LINE_WIDTH) :: input_winds_dir = ""
  ! Total number of bands of partitions on a zone. Set later (during partitioning)
  integer :: nb_bands(6) = 0
  ! Number of bands containing resized partitions (starting from the equator)
  ! Needed for finding neighbours at the equator
  integer :: nb_bands_resized(6) = 0
  ! Number of bands containing square (or rectangular) partitions on a zone. Set later 
  ! (during partitioning)
  integer :: nb_bands_square(6) = 0
  ! Number of latitudes on an hemisphere
  integer :: nb_lat2 = 90
  ! Number of partitions on each sector. Default : the first elements
  ! contains all the number of partitions, until set by the partitioning.
  integer :: nb_partitions(6) = 0
  ! Output data directory and files
  character(LINE_WIDTH) :: output_dir = "output"
  character(LINE_WIDTH) :: output_winds_dir = "output"
  ! Analytical test case : meridional, zonal, solid_rot
  character(LINE_WIDTH) :: test_case = ""
  ! Starting and ending time (for finding winds file)
  integer(kind=k12) :: t_start = UNDEFINED
  integer(kind=k12):: t_end = UNDEFINED
  ! Period for reading winds, format [DD[HH]]MM
  integer(kind=4) :: T_winds = UNDEFINED
  ! Period for output, format HHMM
  integer(kind=4) :: T_output = UNDEFINED

end type

! Here we use a type defined only once in all the program
type(Configuration), private, save :: config

! MPI type for exchanging data (defined once)
integer, private :: mpi_configuration
! MPI type for 8 bytes integer
integer, private :: mpi_longint

! For output
character(*), parameter :: fname_config = "configuration_class.F90"

contains

!-----------------------------------------------------------------------------
!> Constructor for a type defined once (global configuration)
!-----------------------------------------------------------------------------
subroutine new_Configuration(rank, configfile)
  integer :: rank
  character(*) :: configfile
  logical :: file_e

  ! Master process does the init
  if (rank == 0) then

    ! Check the file exists
    inquire( file=trim(configfile), exist=file_e )
    if (.not. file_e) then
      call print_error("Configuration file "//trim(configfile)// &
        " not found.", "new_Configuration", fname_config)
    end if

    open(2, file=trim(configfile))
    call parse_config(2)
    close(2)

    if (.not. NO_RUN) call check_config()
#ifdef VERBOSE
    call print_Configuration()
#endif
    call print_step("Reading configuration file "//trim(configfile)//" done")

  end if

  ! Initialize mpi types
  mpi_configuration = mpi_datatype_null
  call new_mpi_longint()

  ! Send it to other processes
  call broadcast_Configuration()

  call set_times()
end subroutine

subroutine set_times()
  integer :: d1, d2, t1, t2
  d1 = config%t_start/10000
  t1 = mod(config%t_start, 10000)

  d2 = config%t_end/10000
  t2 = mod(config%t_end, 10000)
  ! Split times into a date and time to be more portable
  call set_startend_times(d1, t1, d2, t2)
end subroutine

subroutine new_mpi_longint()
  integer :: ierr

  call mpi_type_contiguous(2, mpi_integer, mpi_longint, ierr)
  call check_mpi_error(ierr, "create type", "new_mpi_longint", fname_config)
  call mpi_type_commit(mpi_longint, ierr)
  call check_mpi_error(ierr, "commit type", "new_mpi_longint", fname_config)


end subroutine

!-------------------------------------------------------------------------------
!> Broadcast the config from master process. Must be done after the partitioning
!> as it can affect the number of bands and partitions.
!-------------------------------------------------------------------------------
subroutine broadcast_Configuration()
  integer :: nb_procs, ierr

  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  call check_mpi_error(ierr, "get number processes", "broadcast_Configuration",&
    fname_config)

  ! Then broadcast it
  if (nb_procs > 1) then
    call new_MPI_Configuration()

    call mpi_bcast(config, 1, mpi_configuration, 0, mpi_comm_world, ierr)
    call check_mpi_error(ierr, "broadcast", "broadcast_Configuration", &
      fname_config)
  end if

  !if (rank == 0) call print_step("* Creating configuration done")
  !call print_step("* Creating configuration done")
end subroutine

!-------------------------------------------------------------------------------
!> Update the number of partitions and bands after the partitioning from the 
!> master process (only for parallel configuration)
!> The caller must check nb_procs > 1
!-------------------------------------------------------------------------------
subroutine update_Configuration()
  integer :: mpi_nb_partitions
  integer :: ierr

  ! Then broadcast it
  mpi_nb_partitions = new_MPI_nb_partitions()
  call mpi_bcast(config, 1, mpi_nb_partitions, 0, mpi_comm_world, ierr)
  call check_mpi_error(ierr, "broadcast config", "update_Configuration", &
    fname_config)

end subroutine

!-----------------------------------------------------------------------------
!> Create a MPI user-defined type for all the configuration
!-----------------------------------------------------------------------------
subroutine new_MPI_Configuration()
  ! Length, displacement, type
  integer, parameter :: n = 18
  integer :: mlength(n), mtype(n)
  integer(kind=mpi_address_kind) :: mlocation(n)
  integer(kind=mpi_address_kind) :: start
  integer :: i, ierr
  character(*), parameter :: func_name = "new_MPI_Configuration" 

  ! Get the address and check the error code
  call mpi_get_address(config, start, ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_config)

  call get_address(config%add_chemic, mlocation(1), func_name, fname_config)
  call get_address(config%advection, mlocation(2), func_name, fname_config)
  call get_address(config%CFL, mlocation(3), func_name, fname_config)
  call get_address(config%dt, mlocation(4), func_name, fname_config)
  call get_address(config%input_dir, mlocation(5), func_name, fname_config)
  call get_address(config%input_winds_dir, mlocation(6), func_name, fname_config)
  call get_address(config%nb_bands, mlocation(7), func_name, fname_config)
  call get_address(config%nb_bands_square, mlocation(8), func_name, fname_config)
  call get_address(config%nb_bands_resized, mlocation(9), func_name, fname_config)
  call get_address(config%nb_partitions, mlocation(10), func_name, fname_config)
  call get_address(config%nb_lat2, mlocation(11), func_name, fname_config)
  call get_address(config%output_dir, mlocation(12), func_name, fname_config)
  call get_address(config%output_winds_dir, mlocation(13), func_name, fname_config)
  call get_address(config%test_case, mlocation(14), func_name, fname_config)
  call get_address(config%t_start, mlocation(15), func_name, fname_config)
  call get_address(config%t_end, mlocation(16), func_name, fname_config)
  call get_address(config%T_winds, mlocation(17), func_name, fname_config)
  call get_address(config%T_output, mlocation(18), func_name, fname_config)

  mlength = (/LINE_WIDTH, LINE_WIDTH, &
    1, 1, &
    LINE_WIDTH, LINE_WIDTH, &
    6, 6, 6, 6, 1, &
    LINE_WIDTH,  LINE_WIDTH, LINE_WIDTH, &
    1, 1, 1, 1/)
  mtype = (/mpi_character, mpi_character, &
    mpi_double_precision, mpi_double_precision, &
    mpi_character, mpi_character, &
    mpi_integer,  mpi_integer, mpi_integer, mpi_integer,  mpi_integer, &
    mpi_character, mpi_character, mpi_character, &
    mpi_longint, mpi_longint, mpi_integer, mpi_integer/)

  do i = 1, size(mlength)
    mlocation(i) = mlocation(i) - start!mlocation(1)
  end do
  !mlocation(1) = 0

  call mpi_type_create_struct(size(mlength), mlength, mlocation, mtype, &
    mpi_configuration, ierr)
  call check_mpi_error(ierr, "create struct type", func_name, fname_config)
  call mpi_type_commit(mpi_configuration, ierr)
  call check_mpi_error(ierr, "commit type", func_name, fname_config)

end subroutine

!-----------------------------------------------------------------------------
!> Create a MPI user-defined type for number of partitions and bands.
!-----------------------------------------------------------------------------
function new_MPI_nb_partitions() result(mpi_nb_partitions)
  ! Length, displacement, type
  integer :: mlength(4), mtype(4)
  integer(kind=mpi_address_kind) :: mlocation(4)
  integer(kind=mpi_address_kind) :: start
  integer :: i, ierr
  integer :: mpi_nb_partitions
  character(*), parameter :: func_name = "new_MPI_nb_partitions" 

  call mpi_get_address(config, start, ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_config)

  call mpi_get_address(config%nb_bands, mlocation(1), ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_config)

  call mpi_get_address(config%nb_bands_square, mlocation(2), ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_config)

  call mpi_get_address(config%nb_bands_resized, mlocation(3), ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_config)

  call mpi_get_address(config%nb_partitions, mlocation(4), ierr)
  call check_mpi_error(ierr, "get address", func_name, fname_config)


  mlength = size(config%nb_bands)
  mtype = mpi_integer

  ! Type is used on an instance of configuration
  !do i = 2, size(mlength)
  do i = 1, size(mlength)
    mlocation(i) = mlocation(i) - start !mlocation(1)
  end do
  !mlocation(1) = 0

  call mpi_type_create_struct(size(mlength), mlength, mlocation, mtype, &
    mpi_nb_partitions, ierr)
  call check_mpi_error(ierr, "create struct type", func_name, fname_config)
  call mpi_type_commit(mpi_nb_partitions, ierr)
  call check_mpi_error(ierr, "commit type", func_name, fname_config)
end function


!-----------------------------------------------------------------------------
!> Parser for configuration file (allows custom error)
!-----------------------------------------------------------------------------
subroutine parse_config(id)
  integer,intent(in) :: id
  integer :: stat
  character(LINE_WIDTH) :: dummy
  character(LINE_WIDTH) :: list_args(16)

  call init_list_args(list_args)

  do while(.true.)
    read(id,'(a)',iostat=stat) dummy
    if (stat < 0) then 
      exit
    end if
    call extract_and_set_val(dummy, list_args)
  end do
end subroutine

!-----------------------------------------------------------------------------
!> Set the value given in the string dummy on the form arg = val
!-----------------------------------------------------------------------------
subroutine extract_and_set_val(dummy, list_args)
  character(*) :: dummy
  character(256) :: arg, val
  character(*) :: list_args(:)
  integer :: pos_eq, pos_comm

  pos_eq = index(dummy, '=')
  if (pos_eq < 1) return
  ! Remove comments
  pos_comm = index(dummy, '!')
  if (pos_comm > 0) return

  val = adjustl(adjustr(dummy(pos_eq+1:)))
  arg = adjustl(adjustr(dummy(1:pos_eq-1)))
  call set_arg(trim(arg), trim(val), list_args)
end subroutine

!-----------------------------------------------------------------------------
! Search argument in the list
!-----------------------------------------------------------------------------
function search_arg(arg,list_args) result (pos)
  character(*), intent(in) :: arg, list_args(:)
  integer :: pos
  integer :: i
  pos = -1

  ! Search in the list of arguments
  !print *, "arg", trim(arg)
  do i = 1,size(list_args)
    ! Don't forget to trim
    if (trim(arg) == trim(list_args(i))) then
      pos = i
      exit
    end if
  end do

  if (pos < 0) then 
    call print_error("parameter "//trim(arg)//" does not exist.", &
      "search_arg", fname_config)
  end if
end function


!-----------------------------------------------------------------------------
! Init the list of arguments
!-----------------------------------------------------------------------------
subroutine init_list_args(list_args)
  character(LINE_WIDTH), intent(inout) :: list_args(:)
  list_args(1) = "add_chemic"
  list_args(2) = "advection"
  list_args(3) = "CFL"
  list_args(4) = "dt"
  list_args(5) = "input_dir"
  list_args(6) = "input_winds_dir"
  list_args(7) = "nb_lat2"
  list_args(8) = "nb_partitions"
  list_args(9) = "output_dir"
  list_args(10) = "output_winds_dir"
  list_args(11) = "test_case"
  list_args(12) = "t_start"
  list_args(13) = "t_end"
  list_args(14) = "T_winds"
  list_args(15) = "T_output"

end subroutine

!-----------------------------------------------------------------------------
!> Check user configuration
!-----------------------------------------------------------------------------
subroutine check_config()

  call check_advection()
  call check_chemistry()

  ! Numerical or analytical data
  if (len(trim(config%test_case)) > 0) then
    call print_warning("Analytical case (for debug)")
  end if

  !  cond = len(trim(config%input_ratio)) > 0 .and. len(trim(config%test_case)) > 0
  !  if (cond ) then
  !    call print_error("Data is numerical, cannot use analytical test case."&
  !      //" Change input_ratio or test_case.", "check_config", fname_config)
  !  end if

  ! Analytical test case
  !call check_analytical()
  call check_timestep()
  call check_nb_procs()
  call check_directories()
  call check_times()
  call check_T_output()

end subroutine

subroutine check_analytical()
  logical :: cond

  if (len(trim(config%test_case)) > 0) then
    cond = trim(config%test_case) == "zonal" .or. &
      trim(config%test_case) == "meridional" .or. &
      trim(config%test_case) == "solid_rot"
    if (.not. cond) then
      call print_error("Analytical test case not implemented", "check_config",&
        fname_config)
    end if
  end if
end subroutine

subroutine check_chemistry()
  if (config%add_chemic /= "no") then
    call print_error("No chemical reactions implemented yet.", "check_config",&
      fname_config)
  end if
end subroutine

subroutine check_T_output()
  if (config%T_output < 1 .and. config%T_output /= UNDEFINED) then
    call print_error("T_output must be > 0", "check_T_output", fname_config)
  end if
  if (config%T_output == UNDEFINED) then
    call print_warning("T_output not initialized (or very large)")
  end if
end subroutine


subroutine check_times()
  if (config%t_start == UNDEFINED) then
    call print_error("Undefined t_start", "check_times", fname_config)
  end if

  if (config%t_end == UNDEFINED) then
    call print_error("Undefined t_end", "check_times", fname_config)
  end if

  if (config%T_winds == UNDEFINED) then
    call print_error("Undefined T_winds", "check_times", fname_config)
  end if

  if (config%t_end < config%t_start) then
    call print_error("t_end must be >= t_start", "check_times", fname_config)
  end if
end subroutine


subroutine check_timestep()
  if (config%dt < 0.) then
    call print_error("Negative timestep", "check_config", fname_config)
  end if

  if (config%dt == 0. .and. config%CFL == 0.)  then
    call print_error("You must set either dt or the CFL", "check_timestep",&
      fname_config)
  end if

  if (config%dt == 0.)  then
    call print_warning("Timestep will be set automatically by CFL")
  else if (config%CFL == 0.)  then
    call print_warning("CFL set from timestep")
  end if
end subroutine

subroutine check_advection()
  logical :: cond

  cond = trim(config%advection) == "zonal" .or. &
    trim(config%advection) == "meridional" .or. &
    trim(config%advection) == "2d"

  if (trim(config%advection) == "zonal") then
    call print_warning("Zonal advection only")
  else if (trim(config%advection) == "meridional") then
    call print_warning("Meridional advection only")
  end if

  if (.not. cond) then
    call print_error("Wrong advection type : should be zonal, meridional, 2d.",&
      "check_config", fname_config)
  end if
end subroutine

subroutine check_nb_procs()
  integer :: nb_parts
  character(7) :: nb
  integer :: nb_procs, ierr

  ! Check we have one processor for each partitions
  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  call check_mpi_error(ierr, "get number processes", "check_config", fname_config)
  nb_parts = sum(config%nb_partitions)

  call print_mesg("Nb processes =", nb_parts)
  call print_mesg("Nb partitions =", nb_procs)

  write (nb, '(i7)') nb_parts
  if (nb_procs > 1 .and. nb_parts /= nb_procs) then
    call print_error("The parallel version must have 1 process per partition."&
      //" Use "//nb//" processes", "check_config", fname_config)
  end if

  ! Check the number of processes is a multiple of 3 (or 1)
  if (nb_parts /= 1 .and. modulo(nb_parts, 3) /= 0) then
    call print_error("Number of partitions must be 1 or a multiple of 3.", &
      "check_config", fname_config)
  end if
end subroutine

subroutine check_directories()
  call check_directory(trim(config%input_dir), .True.)
  call check_directory(trim(config%input_winds_dir), .True.)
  call check_directory(trim(config%output_dir), .False.)
  call check_directory(trim(config%output_winds_dir), .False.)
end subroutine

!-----------------------------------------------------------------------------
! Check if the file exist when the data is numerical. Strings are already 
! trimmed
!-----------------------------------------------------------------------------
subroutine check_num_data(cur_file)
  character(*), intent(in) :: cur_file
  logical :: dir_e

  if (len(cur_file) > 0) then
    ! Check input files if they exists
    inquire( file=cur_file, exist=dir_e )
    if ( dir_e ) then
#ifdef VERBOSE
      write(*,'(3a)'), "File ",cur_file," exists : OK"
#endif
    else
      call print_error("File "//cur_file//" does not exist","check_num_data", &
        fname_config)
    end if
  end if
end subroutine

!-----------------------------------------------------------------------------
! Check the directory exists, otherwise create it if create = True
! String is already trimmed
!-----------------------------------------------------------------------------
subroutine check_directory(str, is_input)
  character(*), intent(in) :: str
  logical :: dir_e, is_input
  integer :: ierr

#ifdef WITH_IFORT
  inquire(directory=str, exist=dir_e )
#else
  inquire(file=str, exist=dir_e )
#endif

  if (is_input) then
    if (.not. dir_e) then
      call print_error("Directory "//str//" does not exist", "check_directory",&
        fname_config) 
    else
      return
    end if
  end if

  if ( dir_e ) then
#ifdef VERBOSE
    call print_warning("Directory "//str//" exists")
#endif
  else
    call print_warning("Directory "//str//" does not exist, creating it")
    ierr = system('mkdir '//str)
    if (ierr /= 0) then
      call print_error("Failed to create directory", "check_directory",&
        fname_config)
    end if
  end if
end subroutine

!-----------------------------------------------------------------------------
! Search and set arg
!-----------------------------------------------------------------------------
subroutine set_arg(arg, val, list_args)
  character(*), intent(in) :: arg, list_args(:), val
  integer :: pos = -1

  pos = search_arg(arg,list_args)

  ! Select data type for casting
  select case(list_args(pos))
  case ("add_chemic")
    config%add_chemic = trim(val)
  case ("advection")
    config%advection = trim(val)
  case ("CFL")
    read (val, *) config%CFL
  case ("dt")
    read (val, *) config%dt
  case ("input_dir")
    config%input_dir = trim(val)
  case ("input_winds_dir")
    config%input_winds_dir = trim(val)
  case ("nb_lat2")
    read (val, *) config%nb_lat2
  case ("nb_partitions")
    ! Total number of partitions is stored in the first element, and split later
    read (val, *) config%nb_partitions(1) 
  case ("output_dir")
    config%output_dir = trim(val)
  case ("output_winds_dir")
    config%output_winds_dir = trim(val)
  case ("test_case")
    config%test_case = trim(val)
  case ("t_start")
    read (val, *) config%t_start
  case ("t_end")
    read (val, *) config%t_end
  case ("T_winds")
    call store_period(val, config%T_winds)
  case ("T_output")
    call store_period(val, config%T_output)
  case default 
    call print_error("Parameter "//trim(val)//"not found in configuration", &
      "set_arg", fname_config)
  end select
end subroutine

!-------------------------------------------------------------------------------
!> Split period from HHMM format to minutes only. Minutes can be a float
!-------------------------------------------------------------------------------
subroutine store_period(val, T)
  character(*) :: val
  integer(kind=4) :: T, n
  integer :: mm, hh, dd

  n = len(val)
  if (n < 2 .and. T > 0) then
    call print_error("Period format is [DD[HH]]MM", "store_period",&
      fname_config)
  end if

  read (val(n-1:n), *), mm
  hh = 0
  if (n > 3) read (val(n-3:n-2), *), hh
  dd = 0
  if (n > 5) read (val(n-5:n-4), *), dd
  read (val(n-1:n), *), mm
  T = dd*24*60 + hh*60 + mm
end subroutine

!-----------------------------------------------------------------------------
!> Compute the number of bands maximizing the number of square partition_
!-----------------------------------------------------------------------------
subroutine compute_nb_bands_square_Configuration(zone, nb_parts)
  integer, intent(in) :: zone, nb_parts

  config%nb_bands_square(zone) = int(sqrt(dble(nb_parts)))
end subroutine

!-----------------------------------------------------------------------------
!> Compute the number of all bands (may be changed later)
!> @param zone : zone indice
!> @param nb_parts : number of partitions on this zone
!-----------------------------------------------------------------------------
subroutine compute_nb_bands_Configuration(zone, nb_parts) 
  integer, intent(in) :: zone, nb_parts
  integer :: nb_bands

  call compute_nb_bands_square_Configuration(zone, nb_parts)

  nb_bands = config%nb_bands_square(zone)
  ! If the number of partitions is not a square
  if (nb_bands ** 2 < nb_parts) nb_bands = nb_bands + 1
  config%nb_bands(zone) = nb_bands
  !write (*,*) "nb bands", config%nb_bands
end subroutine

!-----------------------------------------------------------------------------
!> Set the total number of bands
!-----------------------------------------------------------------------------
subroutine set_nb_bands_Configuration(nb)
  integer, intent(in) :: nb(:)

  config%nb_bands = nb
end subroutine

!-----------------------------------------------------------------------------
!> Set the number of bands with square partitions
!-----------------------------------------------------------------------------
subroutine set_nb_bands_square_Configuration(nb)
  integer, intent(in) :: nb(:)

  config%nb_bands_square = nb
end subroutine


!-----------------------------------------------------------------------------
!> Set the number of bands on a zone
!> @param zone : zone indice
!-----------------------------------------------------------------------------
subroutine set_nb_bands_Configuration_zone(nb, zone)
  integer, intent(in) :: nb, zone

  config%nb_bands(zone) = nb
end subroutine

!-----------------------------------------------------------------------------
!> Set the number of bands with square partitions on a zone
!> @param zone : zone indice
!-----------------------------------------------------------------------------
subroutine set_nb_bands_square_Configuration_zone(nb, zone)
  integer, intent(in) :: nb, zone

  config%nb_bands_square(zone) = nb
end subroutine

!-----------------------------------------------------------------------------
!> Set the number of resized bands
!-----------------------------------------------------------------------------
subroutine set_nb_bands_resized_Configuration(nb)
  integer :: nb(:)

  config%nb_bands_resized = nb
end subroutine

!-----------------------------------------------------------------------------
!> Set the number of resized bands
!> @param zone : zone indice
!-----------------------------------------------------------------------------
subroutine set_nb_bands_resized_Configuration_zone(nb, zone)
  integer :: nb
  integer, intent(in) :: zone

  config%nb_bands_resized(zone) = nb
end subroutine

subroutine set_dt(dt)
  double precision :: dt
  config%dt = dt
end subroutine

subroutine set_cfl(cfl)
  double precision :: cfl
  config%CFL = cfl
end subroutine
!-----------------------------------------------------------------------------
!> Accessor
!-----------------------------------------------------------------------------
function get_Configuration()
  type(Configuration) :: get_Configuration
  get_Configuration = config
end function

function get_CFL() result(cfl)
  double precision :: cfl
  cfl = config%cfl
end function

!> Result in minutes
function get_dt() result(dt)
  double precision :: dt
  dt = config%dt
end function

!> Write at the end of the timestep, so when iter*dt is a multiple of T
function must_write(iter) result (res)
  integer, intent(in) :: iter
  logical :: res
  integer :: r

  call remainder_end_t(r, iter, config%T_output)

  res = (config%T_output > 0 .and. r == 0)
end function

function must_write_final(iter) result (res)
  integer, intent(in) :: iter
  logical :: res
  integer :: r

  call remainder_end_t(r, iter+1, config%T_output)

  res = (config%T_output > 0  .and. r /= 0)
end function


function must_readwinds(iter) result (res)
  integer, intent(in) :: iter
  logical :: res
  integer :: r

  call remainder_end_t(r, iter-1, config%T_winds)
  res = (config%T_winds > 0 .and. r == 0)
end function

subroutine print_iterations(nb_iter)
  double precision :: dt, cfl
  integer, intent(in) :: nb_iter
  integer :: year1, year2
  integer :: month1, month2
  integer :: day1, day2
  integer :: hour1, hour2
  integer :: minute1, minute2
  character(50) :: str1, str2
  real :: nb_days


  dt = config%dt*60
  cfl = config%cfl
  ! C functions
  call extract_date_wrapper(config%t_start, year1, month1, day1, hour1, minute1)
  call extract_date_wrapper(config%t_end, year2, month2, day2, hour2, minute2)

  ! Must cast everything to string, ugly

  write (str1, '(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') &
      year1, "/",month1, "/",day1, " ", hour1, ":",minute1
  write (str2, '(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') &
      year2, "/",month2, "/",day2, " ", hour2, ":",minute2

  call print_mesg("Run from "//trim(str1)//" to "//trim(str2))
  nb_days = dt*nb_iter/(24*3600)

  write (str1, '(a,f10.3,a)') "Duration", nb_days , " days "
  write (str2, '(a,i7,a,f7.2,a,f10.5,a)') "(nb iter=", nb_iter, ", dt=", &
    dt/60, "min, CFL=",cfl,")"
  call print_mesg(trim(str1)//trim(str2))
end subroutine

subroutine extract_date_wrapper(t, y, m, d, h, minute)
  integer(kind=k12), intent(in) :: t
  integer :: y, m, d, h, minute
  integer :: d1, t1

  d1 = t / 10000
  t1 = mod(t, 10000)
  call extract_date(d1, t1, y, m, d, h, minute)
end subroutine

!-----------------------------------------------------------------------------
!> Get the number of resized bands
!> @param zone : zone indice
!-----------------------------------------------------------------------------
function get_nb_bands_resized_Configuration(zone) result (nb)
  integer :: nb
  integer, intent(in) :: zone

  nb = config%nb_bands_resized(zone)
end function

!-----------------------------------------------------------------------------
!> Get the first and last iterations
!> @param zone : zone indice
!-----------------------------------------------------------------------------
subroutine get_iterations(tstart, tend)
  integer(kind=k12) :: tstart, tend

  tstart = config%t_start
  tend = config%t_end
end subroutine

function get_output_dir() result(fname)
  character(LINE_WIDTH) :: fname
  fname = config%output_dir
end function

function get_t_start() result(t)
  integer(kind=k12) :: t
  t = config%t_start
end function


!-----------------------------------------------------------------------------
!> Returns either the concentration or wind data filename 
!> Data is written at the end of the timestep, so iter*dt
!> @param dtype : IS_RATIO, IS_ZWINDS, IS_MWINDS
!-----------------------------------------------------------------------------
function get_output_filename(dtype, iter, tracer) result(fname)
  integer, intent(in) :: dtype, iter
  integer, intent(in), optional :: tracer
  integer :: date, time
  character(LINE_WIDTH) :: fname
  character(12) :: t_str

  call get_time(date, time, iter+1)
  write (t_str, '(i8, i4.4)') date, time

  fname = config%output_winds_dir
  if (dtype == IS_RATIO) fname = config%output_dir

  if (dtype == IS_RATIO) then
    call get_ratio_filename(fname, tracer)
  else if (dtype == IS_ZWINDS) then
    fname = trim(fname) // "/u"
  else if (dtype == IS_MWINDS) then
    fname = trim(fname) // "/v"
  else if (dtype == IS_PARTITIONING) then
    call print_warning("Default partitioning output file")
    fname = trim(fname) // "/partitioning"
  else
    call print_error("Wrong data name", "get_output_filename",&
      fname_config)
    return
  end if
  ! All files contains the iteration 

#ifdef HDF5
  fname = trim(fname) // "_"//t_str//".h5"
#else
  fname = trim(fname) // "_"//t_str//".dat"
#endif
end function

!-------------------------------------------------------------------------------
!> Appends tracer id to filename for I/O
!-------------------------------------------------------------------------------
subroutine get_ratio_filename(fname, tracer) 
  integer, intent(in) :: tracer
  character(LINE_WIDTH) :: fname
  character(4) :: tracer_str

  fname = trim(fname) // "/ratio"

  write (tracer_str, '(i4)') tracer
  fname = trim(fname) // "_" // adjustl(tracer_str)
end subroutine


!-----------------------------------------------------------------------------
!> Returns either the concentration or wind data filename 
!> @param dtype : IS_RATIO, IS_ZWINDS, IS_MWINDS
!-----------------------------------------------------------------------------
function get_input_filename(dtype, iter, tracer) result(fname)
  integer, intent(in) :: dtype, iter, tracer
  integer :: date, time
  logical :: file_e
  character(LINE_WIDTH) :: fname
  character(12) :: t_str

#ifdef DEBUG
  if (iter < 1) then
    call print_error("iteration must be > 0", "get_input_filename", &
      fname_config)
  end if
#endif

  call get_time(date, time, iter)
  write (t_str, '(i8, i4.4)') date, time

  fname = config%input_winds_dir

  if (dtype == IS_RATIO) then
    fname = config%input_dir
    call get_ratio_filename(fname, tracer)
  else if (dtype == IS_ZWINDS) then
    fname = trim(fname) // "/u"
  else if (dtype == IS_MWINDS) then
    fname = trim(fname) // "/v"
  else if (dtype == IS_PARTITIONING) then
    call print_warning("Default partitioning output file")
    fname = trim(fname) // "/partitioning"
  else
    call print_error("Wrong data name", "get_input_filename",&
      fname_config)
    return
  end if

  ! All files contains the iteration 
#ifdef HDF5
  fname = trim(fname) // "_"//t_str//".h5"
#else
  fname = trim(fname) // "_"//t_str//".dat"
#endif

  ! Check if file exists
  inquire( file=trim(fname), exist=file_e )
  if (.not. file_e) then
    call print_error("File "//trim(fname)// " not found.", &
      "get_input_filename", fname_config)
  end if

end function



!-----------------------------------------------------------------------------
!> Get total number partitions
!-----------------------------------------------------------------------------
function get_total_nb_partitions_Configuration() result(nb)
  integer :: nb
  nb = sum(config%nb_partitions)
end function

!-----------------------------------------------------------------------------
!> Get total number partitions on a sector
!-----------------------------------------------------------------------------
function get_total_nb_partitions_sector_Configuration() result(nb)
  integer :: nb
  nb = config%nb_partitions(1) + config%nb_partitions(4)
end function


!-----------------------------------------------------------------------------
!> Get total number partitions on a zone
!-----------------------------------------------------------------------------
function get_nb_partitions_Configuration(zone) result(nb)
  integer, intent(in) :: zone
  integer :: nb
  nb = config%nb_partitions(zone)
end function

!-----------------------------------------------------------------------------
!> Get total number partitions
!-----------------------------------------------------------------------------
function get_total_nb_partitions() result(nb)
  integer :: nb
  nb = sum(config%nb_partitions)
end function


!-----------------------------------------------------------------------------
!> Get number of latitudes on a hemisphere
!-----------------------------------------------------------------------------
function get_nb_lat2_Configuration() result(nb)
  integer :: nb
  nb = config%nb_lat2
end function


!-----------------------------------------------------------------------------
!> Get the number of bands with square partitions
!> @param zone : zone indice
!-----------------------------------------------------------------------------
function get_nb_bands_square_Configuration(zone) result (nb)
  integer, intent(in) :: zone
  integer :: nb
  nb = config%nb_bands_square(zone)
end function

!-----------------------------------------------------------------------------
!> Get the total number of bands on a zone.
!> @param zone : zone indice
!-----------------------------------------------------------------------------
function get_nb_bands_Configuration(zone) result (nb)
  integer, intent(in) :: zone
  integer :: nb
  nb = config%nb_bands(zone)
end function

!-----------------------------------------------------------------------------
!> Get the total number of bands on all sectors
!-----------------------------------------------------------------------------
function get_total_nb_bands_Configuration() result (nb)
  integer :: nb
  nb = sum(config%nb_bands)
end function


!-----------------------------------------------------------------------------
!> Get the total number of bands on 2 sectors (north + south)
!-----------------------------------------------------------------------------
function get_total_nb_bands_Configuration_sector() result (nb)
  integer :: nb
  nb = config%nb_bands(1) + config%nb_bands(4)
end function

!-----------------------------------------------------------------------------
!> Returns true if there is only one partition
!-----------------------------------------------------------------------------
function has_single_partition() result (res)
  logical :: res
  res = (sum(config%nb_partitions) == 1)
end function

!-----------------------------------------------------------------------------
!> Returns true if there is only 3 partitions (one for each sector)
!-----------------------------------------------------------------------------
function has_three_partitions() result (res)
  logical :: res
  res = (sum(config%nb_partitions) == 3)
end function

!-----------------------------------------------------------------------------
!> Check if the concentration data is set analytically
!-----------------------------------------------------------------------------
function is_data_analytical() result(res)
  logical :: res
  !res = .not. is_data_numerical()
  res = (len(trim(config%test_case)) > 0)
end function

!-----------------------------------------------------------------------------
!> Check if the concentration data is set analytically
!-----------------------------------------------------------------------------
function is_data_numerical() result(res)
  logical :: res
  !res = (len(trim(config%input_dir)) > 0)
  res = .not. is_data_analytical()
end function

!-----------------------------------------------------------------------------
!> Meridional advection if numerical case or analytical case (except for zonal
!> test case)
!-----------------------------------------------------------------------------
function is_advection_meridional_only() result(res)
  logical :: res
  res = trim(config%advection) == "meridional"
end function

!-----------------------------------------------------------------------------
!> Zonal advection if numerical case or analytical case (except for meridional
!> test case)
!-----------------------------------------------------------------------------
function is_advection_zonal_only() result(res)
  logical :: res
  res = trim(config%advection) == "zonal"
end function

!-----------------------------------------------------------------------------
!> Set total number partitions in a zone (with idle process)
!-----------------------------------------------------------------------------
subroutine set_nb_partitions_Configuration_zone(nb, zone)
  integer, intent(in) :: nb, zone

  config%nb_partitions(zone) = nb
end subroutine


!-----------------------------------------------------------------------------
!> Set total number partitions (with idle process)
!-----------------------------------------------------------------------------
subroutine set_nb_partitions_Configuration(nb)
  integer, intent(in) :: nb(:)
  integer :: nb_size

  ! Special case : 3 partitions
  nb_size = size(nb)
  config%nb_partitions(1:nb_size) = nb
end subroutine

!-----------------------------------------------------------------------------
!> Get analytical test case
!-----------------------------------------------------------------------------
function get_test_case() result (str)
  character(LINE_WIDTH) :: str

  str = config%test_case
end function

!> Result in minutes
function get_T_winds() result (t)
  integer :: t
  t = config%T_winds
end function


!-----------------------------------------------------------------------------
!> Print
!-----------------------------------------------------------------------------
subroutine print_Configuration()
  write (*,'(a)') "Parameters = "
  write (*,'(a,i10)') "   Nb partitions = ",config%nb_partitions
  write (*,'(a,i10,a,i10)') "   Times = ",config%t_start,"..",config%t_end
  write (*,'(2a)') "   Add chemical reactions = ", trim(config%add_chemic)
  write (*,'(2a)') "   Advection type = ", trim(config%advection)
  write (*,'(2a)') "   Input ratio directory = ",trim(config%input_dir)
  write (*,'(2a)') "   Input winds directory = ", trim(config%input_winds_dir)
  write (*,'(2a)') "   Output winds directory = ", trim(config%output_winds_dir)
  write (*,'(2a)') "   Output directory = ",trim(config%output_dir)
end subroutine

!-----------------------------------------------------------------------------
!> Destructor
!-----------------------------------------------------------------------------
subroutine free_Configuration()
  integer :: ierr

  if (mpi_configuration /= mpi_datatype_null) then
    call mpi_type_free(mpi_configuration, ierr)
    call check_error(ierr,"trying to free config",  "free_Configuration", &
      fname_config)
  end if

  if (mpi_longint /= mpi_datatype_null) then
    call mpi_type_free(mpi_longint, ierr)
    call check_mpi_error(ierr, "trying to free mpi long int type", &
      "free_Configuration", fname_config)
  end if

end subroutine

end module
