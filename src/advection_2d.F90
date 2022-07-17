!-------------------------------------------------------------------------------
! CERFACS, Aviation and Environment
!-------------------------------------------------------------------------------
!
! MODULE: Advection_2d
!
!>  @author
!> Alexis Praga
!
! DESCRIPTION: 
!> A module for 2D-advection.
!
! REVISION HISTORY:
!  24 Avril 2012 - Initial Version
!-------------------------------------------------------------------------------
module Advection_2d
use Partitioning_class
use Profiling

!#define NOLIMITER
implicit none
double precision :: dt = 0

logical :: can_do_zonal, can_do_merid
double precision, allocatable :: mass_init(:), mass_final(:) 
character(*), parameter :: fname_advection = "advection_2d.F90"

contains

!-----------------------------------------------------------------------------
!> Allocate and init data
!-----------------------------------------------------------------------------
subroutine init_advection(distrib, nb)
  type (Partitioning), target :: distrib
  type (Partition), pointer :: part
  double precision :: dt, cfl
  integer, intent(in) :: nb
  integer :: i

  if (nb > 1) then
    call print_error("Should have 1 process per partition", "init_advection", &
      "advection_2d.f90")
  end if

  do i = 1, nb
    part => distrib%parts(i)
    call init_advection_data(part%grid)
    call init_requests(part)

    call count_halo_interior(part)
  end do

  dt = get_dt()
  cfl = get_cfl()
  if (dt ==  0) then
    call set_estimation(cfl, distrib, .True.)
  else
    call set_estimation(dt*60, distrib, .False.)
    ! Must be done for all process
    call set_timestep_date(dt*60)
  end if


  allocate(mass_init(NB_TRACERS))
  allocate(mass_final(NB_TRACERS))
  can_do_zonal = .not. is_advection_meridional_only()
  can_do_merid = .not. is_advection_zonal_only()

end subroutine


!-----------------------------------------------------------------------------
!> Start an 2D-advection
!-----------------------------------------------------------------------------
subroutine start_advection(distrib, rank)
  type(Partitioning), target :: distrib
  type(Partition), pointer :: part
  integer, intent(in) :: rank
  integer :: nb_iter, i, k, ierr, nb_parts
  double precision :: t2, t1

  nb_parts = get_nb_parts_Partitioning(distrib)
  call init_advection(distrib, nb_parts)

  dt = get_dt()*60
  call get_nb_iterations(nb_iter, dt)
  if (rank == 0) call print_iterations(nb_iter)

  ! Init mass
  call mass_tracer_Partitioning(mass_init, distrib)
 
  call write_data(IS_RATIO, 0, distrib)
  call write_data(IS_ZWINDS, 0, distrib)
  call write_data(IS_MWINDS, 0, distrib)

  ! Synchronize all processes
  call mpi_barrier(mpi_comm_world, ierr)

  do k = 1, nb_iter
    call print_iteration(k, rank)
    ! Winds are truly read here
    ! Processes are sync'd after that !
    call update_winds(k, dt, distrib)

    do i = 1, nb_parts
      part => distrib%parts(i)

      call start_timer(ADVECTION)
      call advection_wrapper(part, k)

      call stop_timer(ADVECTION)

#ifdef DEBUG
      call check_data(part)
#endif
    end do

    if (must_write(k) ) then
      call write_data(IS_RATIO, k, distrib)
      call print_iteration(k, rank, .True.)
      call new_mass(rank, distrib)
    end if

  end do

  ! Write concentration if we are not on a multiple of the period
  if (.not. must_write(nb_iter)) then
    call new_mass(rank, distrib)
    call write_data(IS_RATIO, nb_iter, distrib)

  end if

end subroutine

subroutine advection_wrapper(part, k)
  type (Partition) :: part
  integer, intent(in) :: k

  ! select the splittinng strategy
  call basic_splitting(part)

  !call alternate_directions(part, k)
  !call symmetric_splitting(part)
end subroutine

!> Basic operation splitting : zonal then advection
subroutine basic_splitting(part)
  type (Partition) :: part

  call zonal_advection(part, .True.)
  call meridional_advection(part, .False.)
end subroutine

!> Alternate operation splitting : on each timestep, try a direction different 
!> of the previous
subroutine alternate_directions(part, k)
  type (Partition) :: part
  integer, intent(in) :: k

  if (mod(k, 2) == 1) then
    call zonal_advection(part, .True.)
    call meridional_advection(part, .False.)
  else
    call meridional_advection(part, .True.)
    call zonal_advection(part, .False.)
  end if
end subroutine

#ifdef SYMMETRIC_SPLIT
!> Take the mean value of boh directions in splitting
subroutine symmetric_splitting(part)
  type (Partition) :: part

  call cp_ratio_to_buffer(part%grid, 1)
  call zonal_advection(part, .True.)
  call meridional_advection(part, .False.)

  call cp_ratio_to_buffer(part%grid, 2)

  ! restore initial condition
  call cp_buffer_to_ratio(part%grid, 1)
  call meridional_advection(part, .True.)
  call zonal_advection(part, .False.)
  call average_with_buffer(part%grid, 2)
end subroutine
#endif

subroutine new_mass(rank, distrib)
  type (Partitioning) :: distrib
  integer, intent(in) :: rank
  integer :: tracer
  double precision :: m_f

  call mass_tracer_Partitioning(mass_final, distrib)
  if (rank /= 0) return
  do tracer = 1, NB_TRACERS
    !write (*,*) "Tracer", tracer
    call print_mass(mass_init(tracer), mass_final(tracer), rank)
  end do

end subroutine

subroutine print_mass(m_i, m_f, rank)
  integer, intent(in) :: rank
  double precision, intent(in) :: m_i, m_f
  character(80) :: str

  write (str, '(a, 2'//DBLE_FORMAT//', '//DBLE_FORMAT_S//')')  "Mass : ", &
    m_i, m_f, m_f - m_i
  call print_mesg(str)
end subroutine


!> Write to log file for master rank and eventually to stdout
subroutine print_iteration(k, rank, to_screen)
  integer, intent(in) :: k, rank
  logical, optional :: to_screen
  character(6) :: i

  if (rank == 0) then
    write(i, '(i6)'), k

    if (present(to_screen) .and. to_screen) then
      print *, "-----------------------------------------------"
      print *, "Iter"//i
    else
      call print_to_log("-----------------------------------------------")
      call print_to_log("Iter"//i)

    end if
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Set simulation time step with 2D CFL condition, or CFL from timestep.
!> @param input : either timestep or CFL
!-------------------------------------------------------------------------------
subroutine set_estimation(input, distrib, is_timestep)
  type(Partitioning), target :: distrib
  type(Partition), pointer :: part
  double precision, intent(in) :: input
  logical, intent(in) :: is_timestep
  integer :: k, nb, ierr
  double precision :: estim_loc, estim, dummy
  integer :: rank

  if (is_timestep) then
    dummy = 9999999.
  else
    dummy = -9999999.
  end if
  dt = dummy

  nb = get_nb_parts_Partitioning(distrib)
  do k = 1, nb
    part => distrib%parts(k)

    estim_loc = find_local_estim(input, part, dummy, is_timestep)
  end do

  ! Then find the min in every process
  if (is_timestep) then
    call mpi_allreduce(estim_loc, estim, 1, mpi_double_precision, mpi_min, mpi_comm_world, ierr)
  else
    call mpi_allreduce(estim_loc, estim, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierr)
  end if

  if (estim == dummy) then
    if (is_timestep) then
      call print_error("Failed to compute timestep", "set_estimation", &
        fname_advection)
    else
      call print_error("Failed to compute CFL", "set_estimation", &
        fname_advection)
    end if
  else
    if (is_timestep) then
      call set_dt(estim)
      ! Must be done for all process, as C functions do not have yet acess to the times 
      call set_timestep_date(estim*60)
    else
      call set_cfl(estim)
    end if
  end if
end subroutine

!> Winds are in degree/s and we return a timestep in minutes
function find_local_estim(input, part, dummy, is_timestep) result(estim_loc)
  type(Partitioning), target :: distrib
  type(Partition), pointer :: part
  type(Band_Grid), pointer :: grid
  double precision, intent(in) :: input, dummy
  logical, intent(in) :: is_timestep

  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: i, j, k, nb
  integer :: ierr, k_z, k_m
  double precision :: estim_loc, dlon, dlat
  double precision :: cur
  double precision :: lat, coslat

  estim_loc = dummy
  grid => part%grid

  ! Find local dt
  call interior_lat_indices(i_start, i_end, grid)
  dlat = get_dlat(grid)
  k_z = 1
  k_m = 1

  do i = i_start, i_end-1
    call interior_lon_indices(j_start, j_end, grid, i)
    lat = cell_south_lat(i, grid)
    coslat = cos(lat*pi/180.d0)
    dlon = get_dlon(i, grid)*coslat

    do j = j_start, j_end
      call update_local_estim(estim_loc, input, i, j, k_z, k_m, dlat, dlon, &
        part, is_timestep) 
    end do
    ! Last zonal wind on band
    k_z = k_z + 1
  end do

  ! Convert to minutes
  if (is_timestep) estim_loc = estim_loc/60

end function

!-------------------------------------------------------------------------------
!> Compute the timestep or CFL inside the cell (i,j) according to the type of
!> advection
!> @param input : either given cfl or timestep
!> dlon : true distance on the sphere
!-------------------------------------------------------------------------------
subroutine update_local_estim(estim, input, i, j, k_z, k_m, dlat, dlon, part, &
    is_timestep) 
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical :: is_timestep
  integer :: k_z, k_m
  integer :: j_n1, j_n2, l
  integer, intent(in) :: i, j
  double precision, intent(in) :: dlat, input, dlon
  double precision :: estim, dt_loc, dlon_tmp
  double precision :: val, cur_m, cur_z

  grid => part%grid
  ! u*dt /dlon < 1
  if (.not. is_advection_meridional_only()) then
    val = abs(grid%zonal_winds(k_z))

    if (is_timestep) then
      cur_z = input*dlon/val
    else
      cur_z = val*input/dlon
    end if
  end if

  ! (v1 dlon1 + v2 dlon2)*dt / (dlat*dlon) < 1
  if (.not. is_advection_meridional_only()) then
    call south_neighbour_cell_Partition(j_n1, j_n2, i, j, part, i+1)
    cur_m = 0
    do l = j_n1, j_n2
      val = abs(grid%merid_winds(k_m))
      dlon_tmp = cell_interface_length(i, j, i+1, l, part%grid)

      cur_m = cur_m + val*dlon_tmp
      k_m = k_m + 1
    end do
  end if

  if (is_timestep) then
    cur_m = input*dlat*dlon/cur_m
  else
    cur_m = cur_m*input/(dlat*dlon)
    !cur_m = cur_m*input
  end if

  k_z = k_z + 1
  if (is_timestep) then
    estim = min(estim, min(cur_m, cur_z))
  else
    estim = max(estim, max(cur_m, cur_z))
  end if
  !  end if

end subroutine


!-------------------------------------------------------------------------------
!> Zonal advection on a partition, west to east.
!-------------------------------------------------------------------------------
subroutine zonal_advection(part, first_split)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: first_split
  double precision :: t1, t2
  integer :: l, nb_neighb

  if (.not. can_do_zonal) return

  grid => part%grid

  call start_timer(ZONAL_TOTAL)

  ! Special case : a single partition on a band does not receive or send (as it
  ! is alone)
  if (get_total_nb_partitions() == 1) then
    call start_timer(ZONAL_ADVEC)
    call zonal_advection_sequential(part, first_split)
    call stop_timer(ZONAL_ADVEC)
  else

    ! First send/receive ghost cells
    call start_mpi_exchange(part, IS_RATIO, zonal=.True., merid=.False.)

    ! Compute slope inside
    call start_timer(ZONAL_ADVEC)
    call set_zonal_slope(grid, interior=.True.)
    call stop_timer(ZONAL_ADVEC)

    call wait_for_all(part%requests)

    ! Compute slope on the border
    call start_timer(ZONAL_ADVEC)
    call set_zonal_slope(grid, interior=.False.)
    call stop_timer(ZONAL_ADVEC)

    ! Then we can send the slope on the border
    call start_mpi_exchange(part, IS_SLOPE, zonal=.True.,merid=.False.)

    ! Cover communication
    call start_timer(ZONAL_ADVEC)
    call zonal_tracer_fluxes(part, interior=.True.)
    call stop_timer(ZONAL_ADVEC)

    ! Wait for all requests
    call wait_for_all(part%requests)

    ! Finish fluxes computation
    call start_timer(ZONAL_ADVEC)
    call zonal_tracer_fluxes(part, interior=.False.)
    call stop_timer(ZONAL_ADVEC)

    ! Update tracer ratio
    call start_timer(ZONAL_ADVEC)
    call zonal_advection_nocomm(grid, first_split)
    call stop_timer(ZONAL_ADVEC)

  end if
  call stop_timer(ZONAL_TOTAL)

end subroutine



!-------------------------------------------------------------------------------
!> Meridional advection on a partition, north to south.
!-------------------------------------------------------------------------------
subroutine meridional_advection(part, first_split)
  type (Partition) :: part
  integer :: nb_neighb, l
  logical, intent(in) :: first_split
  double precision :: t1, t2
  integer :: i_start, i_end, i
  integer :: j_start, j_end, j, k_f

  if (.not. can_do_merid) return

  ! Special case : a single partition on a band does not receive or send (as it
  ! is alone)
  call start_timer(MERID_TOTAL)

  if (get_total_nb_partitions() == 1) then
    call start_timer(MERID_ADVEC)
    call merid_advection_sequential(part, first_split)
    call stop_timer(MERID_ADVEC)
  else

    call unset_gradients(part%grid)

    ! First send/receive ghost cells
    call start_mpi_exchange(part, IS_RATIO, zonal=.True., merid=.True.)
    ! Compute gradient inside (needed for neighbours interpolation)
    call start_timer(MERID_ADVEC)
    call set_zonal_gradient(part, interior=.True.)
    call stop_timer(MERID_ADVEC)

    call wait_for_all(part%requests)

    call start_timer(MERID_ADVEC)
    call set_zonal_gradient(part, interior=.False.)
    call stop_timer(MERID_ADVEC)

    ! Send it
    call start_mpi_exchange(part, IS_GRADIENT, zonal=.True., merid=.True.)
    ! Cover communications by computing slope on the border
    call start_timer(MERID_ADVEC)
    call set_merid_slope(part, interior=.True.)
    call stop_timer(MERID_ADVEC)

    call wait_for_all(part%requests)

    call start_timer(MERID_ADVEC)
    call set_merid_slope(part, interior=.False.)
    call stop_timer(MERID_ADVEC)

    ! Then we can send the slope on the border
    call start_mpi_exchange(part, IS_SLOPE, zonal=.True., merid=.True.)
    ! Cover communication
    call start_timer(MERID_ADVEC)
    call merid_tracer_fluxes(part, interior=.True.)
    call stop_timer(MERID_ADVEC)

    call wait_for_all(part%requests)

    ! Finish fluxes computation
    call start_timer(MERID_ADVEC)
    call merid_tracer_fluxes(part, interior=.False.)
    call stop_timer(MERID_ADVEC)

    ! Update tracer ratio
    call start_timer(MERID_ADVEC)
    call merid_advection_nocomm(part, first_split)
    call stop_timer(MERID_ADVEC)

  end if
  call stop_timer(MERID_TOTAL)

end subroutine

!-------------------------------------------------------------------------------
!> Sequential zonaladvection
!> Not optimized as we do not need to store the fluxes
!-------------------------------------------------------------------------------
subroutine zonal_advection_sequential(part, first_split)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: first_split
  integer :: i_start, i_end, i
  integer :: j_start, j_end, j
  integer :: k_f, tracer
  double precision :: cur, val

  grid => part%grid
  do tracer = 1, NB_TRACERS
    ! Set all slopes
    call interior_lat_indices(i_start, i_end, grid)
    do i = i_start, i_end
      call interior_lon_indices(j_start, j_end, grid, i)

      do j = j_start, j_end
        cur = zonal_slope(i, j, tracer, grid)
        call set_cell_slope(cur, i, j, tracer, grid)
      end do
    end do

    ! Set all fluxes
    k_f = 1

    do i = i_start, i_end
      call interior_lon_indices(j_start, j_end, grid, i)

      ! Periodic boundaries conditions
      val = zonal_tracer_flux(i, j_end, j_start, k_f, tracer, part)
      call set_zonal_flux(grid, val, k_f, tracer)
      k_f = k_f + 1

      do j = j_start, j_end-1
        val = zonal_tracer_flux(i, j, j+1, k_f, tracer, part)
        call set_zonal_flux(grid, val, k_f, tracer)
        k_f = k_f + 1
      end do

      val = zonal_tracer_flux(i, j_end, j_start, k_f, tracer, part)
      call set_zonal_flux(grid, val, k_f, tracer)
      k_f = k_f + 1
    end do
  end do

  call zonal_advection_nocomm(grid, first_split)
  !call print_border_seq(grid)
  !call check_zonal(part)
end subroutine

!-------------------------------------------------------------------------------
!> Sequential meridional advection
!> Not optimized as we do not need to store the fluxes
!-------------------------------------------------------------------------------
subroutine merid_advection_sequential(part, first_split)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: first_split
  integer :: i_start, i_end, i
  integer :: j_start, j_end, j
  integer :: k_f, tracer
  double precision :: cur

  grid => part%grid

  ! Set all slopes
  do tracer = 1, NB_TRACERS
    call interior_lat_indices(i_start, i_end, grid)
    do i = i_start, i_end
      call interior_lon_indices(j_start, j_end, grid, i)

      ! No slope at the poles (done inside merid slope)
      do j = j_start, j_end
        cur = merid_slope(i, j, tracer, part)
        call set_cell_slope(cur, i, j, tracer, grid)
      end do
    end do

    ! Set all fluxes
    k_f = 1

    do i = i_start, i_end-1
      call interior_lon_indices(j_start, j_end, grid, i)

      do j = j_start, j_end
        call merid_tracer_cell_fluxes(i, j, k_f, tracer, part)
      end do
    end do
  end do

  call merid_advection_nocomm(part, first_split)
  !call check_merid(part)
end subroutine


!-------------------------------------------------------------------------------
!> Zonal advection for all cells (except ghost cells) 
!> Fluxes must be computed
!> Assumes i_start <= i_end and j_start <= j_end 
!-------------------------------------------------------------------------------
subroutine zonal_advection_nocomm(grid, first_split)
  type (Band_Grid) :: grid
  logical, intent(in) :: first_split
  integer :: i, j, tracer
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: k_f, ierr, rank
  double precision :: m_i
  logical :: cond

  call interior_lat_indices(i_start, i_end, grid)
  do tracer = 1, NB_TRACERS
    ! Indices for the fluxes
    k_f = 1
    do i = i_start, i_end
      call interior_lon_indices(j_start, j_end, grid, i)

      ! Cell area : we use a mass density of 1 for the air, so the air mass is 
      ! actually the cell area (which is constant for all j)
      m_i = cell_area(i, j_start, grid)

      do j = j_start, j_end
        !cond = abs(cell_center_lat(i, grid) - 1.68) < 0.01 .and. &
        !abs(cell_center_lon(i, j, grid) - 0.38) < 0.01
        cond = abs(cell_center_lat(i, grid) - 1.68) < 0.01 .and. &
        abs(cell_center_lon(i, j, grid)) > 360-0.8

        call mpi_comm_rank(mpi_comm_world, rank, ierr) 
        !cond = (rank == 4 .and. i == 2 .and. j == 2) .or. &
        !(rank == 0 .and. i == 79 .and. j == 1)
        !if (cond) then
        !  print *, "coucou", i, j
        !end if
!        if (cond) then
!        print *, "cell", cell_center_lat(i, grid), cell_center_lon(i, j, grid)
!        print *, "data", get_cell_ratio(i, j, 1, grid)
     !   print *, "ij", i, j, rank
      !end if
        ! Update using winds at k_f and k_f+1
        call update_tracer_zonal(i, j, k_f, grid, m_i, tracer, first_split)
     !   if (cond) then
     !   print *, "data after", get_cell_ratio(i, j, 1, grid)
     ! end if

        k_f = k_f + 1
      end do
      k_f = k_f + 1
    end do
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Meridional advection (except for ghost cells) 
!> Need fluxes to be computed
!-------------------------------------------------------------------------------
subroutine merid_advection_nocomm(part, first_split)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: first_split
  integer :: i, j, tracer
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: j_p1, j_p2, j_n1, j_n2
  integer :: k_p, k_n, nb_lon
  double precision :: m_i

  grid => part%grid
  call interior_lat_indices(i_start, i_end, grid)

  do tracer = 1, NB_TRACERS
    k_p = 1
    k_n = 1

    ! Indices for the prev and next fluxes
    call starting_winds_indices(k_p, k_n, i_start, part)

    do i = i_start, i_end
      call interior_lon_indices(j_start, j_end, grid, i)
      ! Browse all cells, except ghost on the boundaries
      !call lon_indices(j_start, j_end, part, i)

      m_i = cell_area(i, j_start, grid)
      call skip_first_merid_winds(k_n, k_p, i, part)
      do j = j_start, j_end
        call update_tracer_merid(i, j, k_p, k_n, tracer, m_i, part, first_split)
      end do
      call skip_last_merid_winds(k_n, k_p, i, part)

    end do
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Update the tracer ratio and air mass for zonal advection
!> according to a Van Leer finite-volume scheme
!> @param m_i : the cell area is equal to the mass
!> @param k_f : position in the flux array of the west flux
!-------------------------------------------------------------------------------
subroutine update_tracer_zonal(i, j, k_f, grid, m_i, tracer, first_split)
  type (Band_Grid) :: grid
  logical, intent(in) :: first_split
  integer, intent(in) :: i, j, k_f, tracer
  integer :: rank
  double precision :: U_prev, U_next
  double precision :: F_prev, F_next
  double precision, intent(in) :: m_i
  double precision :: q_new, q_cur
  double precision :: m_air, m_tracer
  double precision :: m_air_prev, ierr
  logical :: cond

  U_prev = zonal_air_flux(k_f, grid)
  U_next =  zonal_air_flux(k_f+1, grid)
  ! Indices are not the same as for the cells centers
  F_prev = get_zonal_flux(grid, k_f, tracer)
  F_next = get_zonal_flux(grid, k_f+1, tracer)

  if (first_split) then
    m_air_prev = m_i
  else
    m_air_prev = get_cell_air_mass(i, j, grid)
  end if

  m_air = m_air_prev + U_prev - U_next
  q_cur = get_cell_ratio(i, j, tracer, grid)
  m_tracer = m_air_prev*q_cur + F_prev - F_next

  q_new = m_tracer/m_air

!  call mpi_comm_rank(mpi_comm_world, rank, ierr) 
!  cond = (rank == 4 .and. i == 2 .and. j == 2) .or. &
!  (rank == 0 .and. i == 79 .and. j == 1)
!  if (cond) then
!    print *, "ij", i, j, "old q", q_cur, "new q", q_new
!    print *, "air flux", U_prev, U_next
!    print *, "trace flux", F_prev, F_next
!  end if 

#ifdef DEBUG
  call check_positivity(m_air, m_tracer)
  call check_cell_ratio("Zonal", q_new, i, j)
#endif

  if (abs(m_air) < DBLE_PREC) q_new = 0

  call set_cell_ratio(q_new, i, j, tracer, grid)
  call set_cell_air_mass(m_air, i, j, grid)
end subroutine

! FIXME avoid spamming when errors
subroutine check_positivity(m_air, m_tracer)
  double precision :: m_air, m_tracer
  character(10) :: str

  if (m_air < DBLE_PREC) then
    write (str, '(F10.5)') m_air
    call print_warning("Air mass negative "//str)
  end if
  if (m_tracer < DBLE_PREC) then
    write (str, '(F10.5)') m_tracer
    call print_warning("Tracer mass negative"//str)
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Must be positive and inferior to 1
!-------------------------------------------------------------------------------
subroutine check_cell_ratio(mesg, q_new, i, j)
  integer, intent(in) :: i, j
  double precision, intent(in) :: q_new
  character(*) :: mesg

  if (q_new > 1.d0 + 1e-13) then
    print *, mesg//" q > 1.", q_new, "at", i, j
  end if

  if (q_new < -1e-13) then
    print *, mesg//" q < 0.", q_new, "at", i, j
  end if

end subroutine

!-------------------------------------------------------------------------------
!> Update the tracer ratio and air mass for meridional advection
!> according to a Van Leer finite-volume scheme. 
!> We use the mass set in the zonal advection.
!> @param k_p, k_n : position in the flux array of the prev (next) flux
!> @param m_i : air mass at the beginning of timestep
!-------------------------------------------------------------------------------
subroutine update_tracer_merid(i, j, k_p, k_n, tracer, m_i, part, first_split)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: i, j, tracer
  logical, intent(in) :: first_split

  integer :: j_n1, j_n2, j_p1, j_p2
  integer :: k_n, k_p, l
  double precision, intent(in) :: m_i
  double precision :: U_prev, U_next
  double precision :: F_prev, F_next
  double precision :: q_new, q_cur
  double precision :: m_air, m_tracer
  double precision :: m_air_prev, tmp
  integer :: j_first, j_last

  grid => part%grid

  ! Can have several neighbours so we sum
  U_prev = 0
  F_prev = 0
  call north_neighbour_cell_Partition(j_p1, j_p2, i, j, part, i-1)
  if (j_p1 > 0 .and. j_p2 > 0) then
    do l = j_p1, j_p2

      if (to_or_from_interior_cell(l, i, j, i-1, part%grid)) then
        U_prev = U_prev + merid_air_flux(i, j, i-1, l, k_p, grid)
        F_prev = F_prev + get_merid_flux(grid, k_p, tracer)
        k_p = k_p + 1
      end if
    end do

  end if

  U_next = 0
  F_next = 0
  call south_neighbour_cell_Partition(j_n1, j_n2, i, j, part, i+1)
  if (j_n1 > 0 .and. j_n2 > 0) then
    do l = j_n1, j_n2
      if (to_or_from_interior_cell(l, i, j, i+1, part%grid)) then
        U_next = U_next + merid_air_flux(i, j, i+1, l, k_n, grid)
        F_next = F_next + get_merid_flux(grid, k_n, tracer)
        k_n = k_n + 1
      end if

    end do
  end if

  if (first_split) then
    m_air_prev = m_i
  else
    m_air_prev = get_cell_air_mass(i, j, grid)
  end if

  ! Input : positive winds means toward the south
  m_air = m_air_prev + U_prev - U_next
  q_cur = get_cell_ratio(i, j, tracer, grid)

  m_tracer = m_air_prev*q_cur + F_prev - F_next
  q_new = m_tracer/m_air

#ifdef DEBUG
  call check_positivity(m_air, m_tracer)
  call check_cell_ratio("Meridional", q_new, i, j)
#endif

  ! If null air mass
  if (abs(m_air) < DBLE_PREC) q_new = 0

  call set_cell_ratio(q_new, i, j, tracer, grid)

  ! FIXME
  ! All tracers needs initial mass, so set it at the then
  ! FIXME 
  ! No need to update if meridional advec is the last step of dimensional
  ! splitting
  if (tracer == NB_TRACERS) then
    call set_cell_air_mass(m_air, i, j, grid)
  end if
end subroutine


!-------------------------------------------------------------------------------
!> Compute the tracer transfers and store them
!> @param k : the position in the slope array
!-------------------------------------------------------------------------------
subroutine zonal_tracer_fluxes(part, interior)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer :: i, j, k_f, tracer
  logical, intent(in) :: interior
  integer :: i_start, i_end
  integer :: j_start, j_end
  double precision :: val

  grid => part%grid

  do tracer = 1, NB_TRACERS
    call interior_lat_indices(i_start, i_end, grid)
    ! Indice for position in fluxes array 
    k_f = 1
    do i = i_start, i_end
      call interior_lon_indices(j_start, j_end, grid, i)

      ! Beware, we are on cells boundaries now
      ! Store flux between j and j+1 into position j (k_f actually)
      if (interior) then

        k_f = k_f + 1
        do j = j_start, j_end-1
          val = zonal_tracer_flux(i, j, j+1, k_f, tracer, part)
          call set_zonal_flux(grid, val, k_f, tracer)
          k_f = k_f + 1
        end do

      else

        val = zonal_tracer_flux(i, j_start-1, j_start, k_f, tracer, part)
        call set_zonal_flux(grid, val, k_f, tracer)
        k_f = k_f + (j_end - j_start + 1)
        val = zonal_tracer_flux(i, j_end, j_end+1, k_f, tracer, part)
        call set_zonal_flux(grid, val, k_f, tracer)

      end if
      k_f = k_f + 1

    end do
  end do
end subroutine

! Count cells
function count_zonal_tracer_fluxes(part, interior) result (nb)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: interior
  integer :: i_start, i_end, i
  integer :: j_start, j_end, j
  integer :: nb
  nb = 0

  grid => part%grid
  call interior_lat_indices(i_start, i_end, grid)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, grid, i)
    if (interior) then
      nb = nb + j_end-1-j_start+1
    else
      nb = nb + 3
    end if
  end do
end function


!-------------------------------------------------------------------------------
!> Update array indice so we skip the cell south fluxes. Assume there are 
!> neighbours.
!-------------------------------------------------------------------------------
subroutine skip_merid_flux(k_f, i, j, part)
  type (Partition) :: part
  integer, intent(in) :: i, j
  integer :: j1, j2, k_f, l

  call south_neighbour_cell_Partition(j1, j2, i, j, part, i+1)
  do l = j1, j2
    if (to_or_from_interior_cell(l, i, j, i+1, part%grid)) then
      k_f = k_f + 1
    end if
  end do

end subroutine

!-------------------------------------------------------------------------------
!> Compute the zonal gradients for meridional advection
!-------------------------------------------------------------------------------
subroutine set_zonal_gradient(part, interior)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: interior
  integer :: i_start, i_end

  grid => part%grid

  ! Does not compute gradient on north/south ghost cells, they are received
  ! later
  call interior_lat_indices(i_start, i_end, grid)
  !i_start = 1
  !i_end = grid%nb_lat

  if (interior) then
    call zonal_gradient_interior(i_start, i_end, grid)
  else
    call zonal_gradient_border(i_start, i_end, grid)
  end if
end subroutine


subroutine zonal_gradient_interior(i_start, i_end, grid)
  type (Band_Grid):: grid
  integer, intent(in) :: i_start, i_end
  integer :: j_start, j_end
  integer :: i, j, tracer
  double precision :: cur

  do tracer = 1, NB_TRACERS
    do i = i_start+1, i_end-1

      call interior_lon_indices(j_start, j_end, grid, i)

      do j =j_start, j_end
        cur = zonal_gradient(i, j, tracer, grid)
        call set_cell_gradient(cur, i, j, tracer, grid)
      end do
    end do
  end do
end subroutine

! Count interior cells
function count_zonal_gradient_interior(grid) result(nb)
  type (Band_Grid) :: grid
  integer :: i_start, i_end, i
  integer :: j_start, j_end, j
  integer :: nb
  nb = 0

  call interior_lat_indices(i_start, i_end, grid)
  do i = i_start+1, i_end-1

    call interior_lon_indices(j_start, j_end, grid, i)

    nb = nb + j_end - j_start + 1
  end do

end function

subroutine zonal_gradient_border(i_start, i_end, grid)
  type (Band_Grid) :: grid
  integer, intent(in) :: i_start, i_end
  integer :: j_start, j_end
  integer :: i, j, tracer
  double precision :: cur

  do tracer = 1, NB_TRACERS
    do i = i_start, i_end
      call interior_lon_indices(j_start, j_end, grid, i)

      if (i > i_start .and. i < i_end) then
        cur = zonal_gradient(i, j_start, tracer, grid)
        call set_cell_gradient(cur, i, j_start, tracer, grid)
        cur = zonal_gradient(i, j_end, tracer, grid)
        call set_cell_gradient(cur, i, j_end, tracer, grid)

      else

        do j = j_start, j_end
          cur = zonal_gradient(i, j, tracer, grid)
          call set_cell_gradient(cur, i, j, tracer, grid)
        end do
      end if
    end do
  end do
end subroutine

! Count cells
function count_zonal_gradient_border(grid) result(nb)
  type (Band_Grid) :: grid
  integer :: i_start, i_end, i
  integer :: j_start, j_end, j
  integer :: nb
  nb = 0

  call interior_lat_indices(i_start, i_end, grid)
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, grid, i)

    if (i > i_start .and. i < i_end) then
      nb = nb + 2
    else
      nb = nb + j_end-j_start+1
    end if
  end do

end function

!-------------------------------------------------------------------------------
!> Compute the tracer transfers and store them
!-------------------------------------------------------------------------------
subroutine merid_tracer_fluxes(part, interior)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: interior
  integer :: i_start, i_end

  grid => part%grid

  call interior_lat_indices(i_start, i_end, grid)
  if (has_north_ghost_cells(grid)) i_start = i_start - 1
  if (.not. has_south_ghost_cells(grid)) i_end = i_end - 1

  if (interior) then
    call merid_fluxes_interior(i_start, i_end, part)
  else
    call merid_fluxes_border(i_start, i_end, part)
  end if
end subroutine

!-------------------------------------------------------------------------------
!> Compute the fluxes at the interior. The lat indices must be precomputed.
!-------------------------------------------------------------------------------
subroutine merid_fluxes_border(i_start, i_end, part) 
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: i_start, i_end
  integer ::  g_start, g_end
  integer :: int_start, int_end
  integer :: i, j, k_f, tracer

  ! Indice for fluxes array 
  do tracer = 1, NB_TRACERS
    k_f = 1

    do i = i_start, i_end
      call lon_indices(g_start, g_end, part, i)

      call  merid_interior(int_start, int_end, i, i_end, g_start, g_end, part%grid)

      do j = g_start, g_end
        ! Do north and south border
        if (i == i_start .or. i == i_end) then
          call merid_tracer_cell_fluxes(i, j, k_f, tracer, part)
        else
          if (j >= int_start .and. j <= int_end) then
            call skip_merid_flux(k_f, i, j, part)
          else
            call merid_tracer_cell_fluxes(i, j, k_f, tracer, part)
          end if
        end if
      end do
    end do
  end do

end subroutine

function count_merid_tracer_fluxes(part, interior) result(nb)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: interior
  integer :: i_start, i_end, nb

  grid => part%grid

  call interior_lat_indices(i_start, i_end, grid)
  if (has_north_ghost_cells(grid)) i_start = i_start - 1
  if (.not. has_south_ghost_cells(grid)) i_end = i_end - 1

  if (interior) then
    nb = count_merid_fluxes_interior(i_start, i_end, part)
  else
    nb = count_merid_fluxes_border(i_start, i_end, part)
  end if
end function

! Count cells
function count_merid_fluxes_border(i_start, i_end, part) result(nb)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: i_start, i_end
  integer ::  g_start, g_end, j, i
  integer :: int_start, int_end
  integer :: nb
  nb = 0

  do i = i_start, i_end
    call lon_indices(g_start, g_end, part, i)
    call  merid_interior(int_start, int_end, i, i_end, g_start, g_end, part%grid)

    do j = g_start, g_end
      if (i == i_start .or. i == i_end) then
        nb = nb + 1
      else
        if (.not. (j >= int_start .and. j <= int_end)) then
          nb = nb + 1
        end if
      end if
    end do
  end do

end function


!-------------------------------------------------------------------------------
!> Compute the indices for interior cells in respect to meridional fluxes 
!> (toward the south only)
!-------------------------------------------------------------------------------
subroutine merid_interior(int_start, int_end, i, i_end, g_start, g_end, grid)
  type (Band_Grid) :: grid
  integer, intent(in) :: i, i_end
  integer, intent(in) :: g_start, g_end
  integer :: int_start, int_end, j_start, j_end
  integer :: nb_ghostw, nb_ghoste

  call interior_lon_indices(j_start, j_end, grid, i)
  ! We browse ghost cells for fluxes toward interior cells
  ! and the border is of size nb_ghosts
  ! This also deals with the special case of sector ghost cells

  nb_ghoste = get_nb_ghosts_east(i, grid)
  nb_ghostw = get_nb_ghosts_west(i, grid)

  ! Special case slope may be > 1 
  if (i < i_end) then
    nb_ghoste = max(nb_ghoste, get_nb_ghosts_east(i+1, grid))
    nb_ghostw = max(nb_ghostw, get_nb_ghosts_west(i+1, grid))
  end if
  int_start = j_start+nb_ghostw
  int_end = j_end - nb_ghoste

end subroutine

!-------------------------------------------------------------------------------
!> Compute the fluxes at the interior. The lat indices must be precomputed.
!-------------------------------------------------------------------------------
subroutine merid_fluxes_interior(i_start, i_end, part) 
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: i_start, i_end
  integer :: int_start, int_end, g_start, g_end
  integer :: i, j, k_f, tracer

  ! Indice for fluxes array 
  grid => part%grid

  do tracer = 1, NB_TRACERS
    k_f = 1
    do i = i_start, i_end
      call lon_indices(g_start, g_end, part, i)
      call  merid_interior(int_start, int_end, i, i_end, g_start, g_end, part%grid)

      do j = g_start, g_end
        ! Skip north and south
        if (i == i_start .or. i == i_end) then
          call skip_merid_flux(k_f, i, j, part)

        else
          if (j >= int_start .and. j <= int_end) then
            call merid_tracer_cell_fluxes(i, j, k_f, tracer, part)
          else
            call skip_merid_flux(k_f, i, j, part)
          end if
        end if
      end do
    end do
  end do

end subroutine

! Count cells
function count_merid_fluxes_interior(i_start, i_end, part) result(nb)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: i_start, i_end
  integer :: int_start, int_end, g_start, g_end
  integer :: i, j, k_f
  integer :: nb
  nb = 0

  ! Indice for fluxes array 
  grid => part%grid

  do i = i_start, i_end
    call lon_indices(g_start, g_end, part, i)
    call  merid_interior(int_start, int_end, i, i_end, g_start, g_end, part%grid)

    do j = g_start, g_end
      ! Skip north and south
      if (i /= i_start .and.i /= i_end) then
        if (j >= int_start .and. j <= int_end) then
          nb = nb + 1
        end if
      end if
    end do
  end do

end function

!-------------------------------------------------------------------------------
!> Set the meridional south fluxes of cells at (i,j)
!> @param k_f : position in the flux array
!> @param tracer : tracer id
!-------------------------------------------------------------------------------
subroutine merid_tracer_cell_fluxes(i, j, k_f, tracer, part)
  type (Partition) :: part
  integer, intent(in) :: i, j, tracer
  integer :: j_neigh1, j_neigh2, l
  integer :: k_f
  double precision :: val

  call south_neighbour_cell_Partition(j_neigh1, j_neigh2, i, j, part, i+1)
  ! Only compute fluxes to interior cells
  do l = j_neigh1, j_neigh2
    if (to_or_from_interior_cell(l, i, j, i+1, part%grid)) then
      val = merid_tracer_flux(i, j, l, k_f, tracer, part)
      call set_merid_flux(part%grid, val, k_f, tracer)
      k_f = k_f + 1
    end if
  end do

end subroutine

!-------------------------------------------------------------------------------
!> Compute tracer flux between cell (i,j) and cell (i, j_neighb)
!> @param k : position in the winds array
!-------------------------------------------------------------------------------
function zonal_tracer_flux(i, j, j_neighb, k, tracer, part) result(flux)
  type (Partition) :: part
  integer, intent(in) :: i, j, j_neighb, tracer
  integer :: k
  double precision :: dlon, lat, q_var
  double precision :: U, q_hat, flux
!  integer :: rank, ierr
!  logical :: cond

  U = zonal_air_flux(k, part%grid)

  lat = cell_center_lat(i, part%grid) 
  dlon = get_dlon(i, part%grid)*cos(lat*pi/180d0)
  q_var = part%grid%zonal_winds(k)*dt/dlon
  q_hat = compute_q_hat(q_var, i, j, i, j_neighb, tracer, part)

  flux = U*q_hat

!  call mpi_comm_rank(mpi_comm_world, rank, ierr) 
!  cond = (rank == 4 .and. i == 2 .and. j_neighb == 2) .or. &
!  (rank == 0 .and. i == 79 .and. j_neighb == 1)
!  if (cond) then
!    print *, "q var", q_var, "q hat", q_hat
!    print *, "cur ratio", get_cell_ratio(i, j, tracer, part%grid), "at", i, j
!    print *, "cur slope", get_cell_slope(i, j, tracer, part%grid)
! 
!  print *, "zonal flux", flux
!end if
 
end function

!-------------------------------------------------------------------------------
!> Compute tracer flux between cell (i,j) and (i+1, j_neighb)
!> @param k : position in the winds array
!> @param tracer : tracer id
!-------------------------------------------------------------------------------
function merid_tracer_flux(i, j, j_neighb, k, tracer, part) result(flux)
  type (Partition) :: part
  integer, intent(in) :: i, j, j_neighb,tracer
  integer :: k
  double precision :: U, q_hat, flux
  double precision :: q_var, dlat
  integer :: rank, ierr, nb_procs

  U = merid_air_flux(i, j, i+1, j_neighb, k, part%grid)

  ! TODO optimize
  dlat = get_dlat(part%grid)

  q_var = part%grid%merid_winds(k)*dt/dlat
  q_hat = compute_q_hat(q_var, i, j, i+1, j_neighb, tracer, part)

  flux = U*q_hat
end function

!-------------------------------------------------------------------------------
!> Compute air transfer at interface k between cell (i,j) and (i_neighb, j_neighb)
!-------------------------------------------------------------------------------

function merid_air_flux(i, j, i_neighb, j_neighb, k, grid) result (U)
  type (Band_Grid) :: grid
  integer, intent(in) :: i, j, i_neighb, j_neighb, k
  double precision :: U, dlon

  ! Interface length
  dlon = cell_interface_length(i, j, i_neighb, j_neighb, grid)
  U = grid%merid_winds(k)*dt*dlon

end function


!-------------------------------------------------------------------------------
!> Compute air transfer at interface k
!-------------------------------------------------------------------------------
function zonal_air_flux(k, grid) result (U)
  type (Band_Grid) :: grid
  integer, intent(in) :: k
  double precision :: U, dlat

  dlat = get_dlat(grid)
  U = grid%zonal_winds(k)*dt*dlat
  !  if (grid%zonal_winds(k) == UNDEFINED) then 
  !    print *, "undefined wind"
  !  end if


end function

!-------------------------------------------------------------------------------
!> Compute q hat with a second order scheme between (i,j) and (i_neighb,j_neighb)
!> Works for zonal and meridional.
!> @param q_var : proportion of tracer going out/in. It is better than to
!> compute U/m_i as we can simply the fraction
!> @param tracer : tracer id
!-------------------------------------------------------------------------------
function compute_q_hat(q_var, i, j, i_neighb, j_neighb, tracer, part) &
    result (q_hat)
  type (Partition) :: part
  integer, intent(in) :: i, j, tracer
  integer, intent(in) :: i_neighb, j_neighb
  double precision, intent(in) :: q_var
  double precision :: q_hat, dq, q

  ! q_var = U/mi so the sign of U is the sign of q_var
  if (q_var > 0) then
    q = get_cell_ratio(i, j, tracer, part%grid)

    dq = get_cell_slope(i, j, tracer, part%grid)
    q_hat = q + 0.5d0*(1d0 - q_var)*dq

  else
    q = get_cell_ratio(i_neighb, j_neighb, tracer, part%grid)
    dq = get_cell_slope(i_neighb, j_neighb, tracer, part%grid)

    q_hat = q - 0.5d0*(1d0 + q_var)*dq
  end if
end function

!-------------------------------------------------------------------------------
!> Set zonal slope for zonal border cells or the interior (the rest)
!-------------------------------------------------------------------------------
subroutine set_zonal_slope(grid, interior)
  type (Band_Grid) :: grid
  logical, intent(in) :: interior
  integer :: i, j, tracer
  integer :: i_start, i_end
  integer :: j_start, j_end
  double precision :: cur

  do tracer = 1, NB_TRACERS
    call interior_lat_indices(i_start, i_end, grid)

#ifdef HAS_OPENMP
    !$OMP PARALLEL DO PRIVATE(j_start, j_end, j)
#endif
    do i = i_start, i_end
      call interior_lon_indices(j_start, j_end, grid, i)
      !  id = OMP_GET_THREAD_NUM()

      if (interior) then
        ! Skip first and last cell
        do j = j_start+1, j_end-1
          cur = zonal_slope(i, j, tracer, grid)
                    call set_cell_slope(cur, i, j, tracer, grid)
        end do
      else
        cur = zonal_slope(i, j_start, tracer, grid)
        call set_cell_slope(cur, i, j_start, tracer, grid)
        cur = zonal_slope(i, j_end, tracer, grid)
        call set_cell_slope(cur, i, j_end, tracer, grid)

      end if
    end do
#ifdef HAS_OPENMP
    !$OMP END PARALLEL DO
#endif
  end do
end subroutine

!-------------------------------------------------------------------------------
!> Set meridional slope for meridional border cells or the interior (the rest)
!> The border has 2 layers in this case due to neighbour interpolation
!-------------------------------------------------------------------------------
subroutine set_merid_slope(part, interior)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: interior
  integer :: i, j, k, tracer
  integer :: i_start, i_end, int_start, int_end
  integer :: g_start, g_end
  integer :: j_start, j_end
  integer :: step, n_west, n_east
  double precision :: cur

  grid => part%grid
  call interior_lat_indices(i_start, i_end, grid)

  do tracer = 1, NB_TRACERS
    ! First and last latitude lines
    if (.not. interior) then
      do k = i_start, i_end
        if (k > i_start .and. k < i_end) cycle

        call interior_lon_indices(j_start, j_end, grid, k)
        do j = j_start, j_end
          cur = merid_slope(k, j, tracer, part)
          call set_cell_slope(cur, k, j, tracer, grid)
        end do
      end do
    end if

    do i = i_start+1, i_end-1
      call interior_lon_indices(j_start, j_end, grid, i)

      call  merid_interior(int_start, int_end, i, i_end, g_start, g_end, grid)

      if (interior) then
        ! Skip first and last cells as we may need to compute the slope for
        ! neighbours. This depends on the number of ghost cells on previous and
        ! next line
        do j = int_start, int_end
          cur = merid_slope(i, j, tracer, part)
          call set_cell_slope(cur, i, j, tracer, grid)
        end do
      else
        do j = j_start, int_start-1
          cur = merid_slope(i, j, tracer, part)
          call set_cell_slope(cur, i, j, tracer, grid)
        end do
        do j = int_end+1, j_end
          cur = merid_slope(i, j, tracer, part)
          call set_cell_slope(cur, i, j, tracer, grid)
        end do
      end if
    end do
  end do
end subroutine

! Count cells
function count_merid_slope(part, interior) result(nb)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  logical, intent(in) :: interior
  integer :: i_start, i_end, i, k
  integer :: j_start, j_end, j
  integer :: nb
  nb = 0

  grid => part%grid
  call interior_lat_indices(i_start, i_end, grid)

  if (.not. interior) then
    do k = i_start, i_end
      if (k > i_start+2 .and. k < i_end-2) cycle
      call interior_lon_indices(j_start, j_end, grid, k)
      nb = nb + j_end - j_start + 1
    end do
  end if

  do i = i_start+2, i_end-2
    call interior_lon_indices(j_start, j_end, grid, i)

    if (interior) then
      nb = nb + j_end-2 - (j_start+2) + 1
    else
      nb = nb +4
    end if
  end do
end function


!-------------------------------------------------------------------------------
!> Compute zonal slope at i,j with limitation
!-------------------------------------------------------------------------------
function zonal_slope(i, j, tracer, grid) result (res)
  type (Band_Grid) :: grid
  integer, intent(in) :: i, j, tracer
  integer :: j_prev, j_next
  double precision :: q_prev, q_cur, q_next
  double precision :: res

  j_prev = j-1
  j_next = j+1
  if (has_single_partition()) then
    if (j_prev < 1) j_prev = grid%nb_lon(i)
    if (j_next > grid%nb_lon(i)) j_next = 1
  end if

  !call start_timer(ZONAL_ADVEC)
  q_prev = get_cell_ratio(i, j_prev, tracer, grid)
  !call stop_timer(ZONAL_ADVEC)
  q_cur = get_cell_ratio(i, j, tracer, grid)
  q_next = get_cell_ratio(i, j_next, tracer, grid)

  res = slope_limitation(q_prev, q_cur, q_next)

end function

!-------------------------------------------------------------------------------
!> Compute meridional slope at i,j with limitation
!-------------------------------------------------------------------------------
function merid_slope(i, j, tracer, part) result (res)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: tracer
  integer :: i, j
  double precision :: q_prev, q_cur, q_next
  double precision :: res, q_interp

  grid => part%grid

  ! No slope at the poles
  res = 0.
  if (is_north_pole(i, part) .or. is_south_pole(i, part)) return

  ! Interpolates prev and next value
  if (i == 1) then
    q_prev = 0
  else
    q_prev = ratio_merid_interpolation(i, j, i-1, tracer, part)
  end if

  if (i == grid%nb_lat) then
    q_next = 0
  else
    q_next = ratio_merid_interpolation(i, j, i+1, tracer, part)
  end if

  q_cur = get_cell_ratio(i, j, tracer, grid)

  res = slope_limitation(q_prev, q_cur, q_next)
end function

! For testing
!function interpolation_linear(i_neighb, j1, j2, i, j, grid) result (q_interp)
!  type (Band_Grid) :: grid
!  integer, intent(in) :: j1, j2, i_neighb, i,j
!  double precision :: q_next, lon1, lon2
!  double precision :: coslat, q_prev, slope, mid_cur
!  double precision :: q_interp, lat
!
!  lat = cell_center_lat(i_neighb, grid)
!  coslat = cos(lat*pi/180.d0)
!  q_next = get_cell_ratio(i_neighb, j2, grid)
!  q_prev = get_cell_ratio(i_neighb, j1, grid) 
!  slope = q_next - q_prev
!  mid_cur = cell_center_lon(i, j, grid)
!  lon1 = cell_center_lon(i_neighb, j1, grid)
!  lon2 = cell_center_lon(i_neighb, j2, grid)
!  if (lon1 == lon2) return
!  slope = slope / ((lon2 - lon1)*coslat)
!  q_interp = q_prev  + (mid_cur - lon1)*coslat*slope
!end function

!-------------------------------------------------------------------------------
!> Piecewise interpolation of neighbouring cells concentrations
!> Better than a linear interpolation from the cell concentrations directly.
!-------------------------------------------------------------------------------
function ratio_merid_interpolation(i, j, i_neighb, tracer, part) result(q_interp)
  type (Partition), target :: part
  type (Band_Grid), pointer :: grid
  integer, intent(in) :: i, j, i_neighb, tracer
  integer :: j1, j2, j_neighb
  double precision :: slope, dx, q_interp
  double precision :: mid_cur, mid_neighb
  double precision :: q_prev, coslat, lat
  integer :: nb_procs, ierr, j_prev, rank

  grid => part%grid
  if (i_neighb < i) then
    call north_neighbour_cell_Partition(j1, j2, i, j, part, i_neighb)
  else
    call south_neighbour_cell_Partition(j1, j2, i, j, part, i_neighb)
  end if

  if (j1 < 0 .or. j2 < 0) then
    call print_error("Negative neighbours", "ratio_merid_interpolation", &
      fname_advection)
  end if

  ! Select the closest cell for interpolation
  ! Even for 3 neighbours, no need to use more than the first two neighbours
  mid_cur = cell_center_lon(i, j, grid)
  j_neighb = j1
  if (mid_cur > cell_east_lon(i_neighb, j_neighb, grid) + DBLE_PREC) then
    j_neighb = j_neighb+1
  end if

  if (has_single_partition()) then
    slope = zonal_gradient(i_neighb, j_neighb, tracer, grid)

  else
    slope = get_cell_gradient(i_neighb, j_neighb, tracer, grid)
  end if

  ! Interpolation is done from the cell center
  mid_neighb = cell_center_lon(i_neighb, j_neighb, grid)

  q_prev = get_cell_ratio(i_neighb, j_neighb, tracer, grid)
  lat = cell_center_lat(i_neighb, grid)
  coslat = cos(lat*pi/180.d0)

  ! Negative for the right neigbour
  dx = (mid_cur - mid_neighb)*coslat
  ! Piecewise interpolation from the gradient
  q_interp = q_prev + dx*slope


end function

!-------------------------------------------------------------------------------
!> Compute zonal gradient, used for meridional advection (in neighbour
!> interpolation)
!> interpolation
!> @param i, j : current cell
!-------------------------------------------------------------------------------
function zonal_gradient(i, j, tracer, grid) result (slope)
  type (Band_Grid) :: grid
  integer, intent(in) :: i, j, tracer
  integer :: j_prev, j_next, nb_lon
  double precision :: q_prev, q_next, lat, dlon
  double precision :: slope, coslat
  integer :: nb_procs, ierr, jtmp

  ! Use ghost cells if more than 1 partition
  j_prev = j - 1
  j_next = j + 1

  ! Boundaries for single partition
  if (has_single_partition()) then
    nb_lon = grid%nb_lon(i)
    if (j_prev < 1) j_prev = nb_lon
    if (j_next > nb_lon) j_next = 1
  end if

  q_prev = get_cell_ratio(i, j_prev, tracer, grid)
  q_next = get_cell_ratio(i, j_next, tracer, grid)

  lat = cell_center_lat(i, grid)
  coslat = cos(lat*pi/180.d0)
  dlon = grid%dlon(i)*coslat
  slope = 0.5d0*(q_next - q_prev)/dlon

end function

!-------------------------------------------------------------------------------
!> Compute slope with limitation (van Leer)
!-------------------------------------------------------------------------------
function slope_limitation(q_prev, q_cur, q_next) result (res)
  double precision :: q_prev, q_cur, q_next
  double precision :: res
  integer :: sign1, sign2, sign3

#ifdef NOLIMITER
  ! No limiter for testing
  res = 0.5*(q_next - q_prev)
#else

  res = 0.
  sign1 = int(sign(1.d0, q_next - q_cur))
  sign2 = int(sign(1.d0, q_next - q_prev))
  sign3 = int(sign(1.d0, q_cur - q_prev))
  ! q_cur must be between q_prev and q_next
  if (sign1 == sign2 .and. sign2 == sign3) then
    res = min(0.5d0*abs(q_next-q_prev), 2d0*abs(q_next-q_cur))
    res = min(res, 2d0*abs(q_cur - q_prev))
    res = sign(1.d0, q_next - q_cur)*res
  end if
#endif
end function

!-------------------------------------------------------------------------------
!> Meridional advection for border cells (which do need ghost cells)
!-------------------------------------------------------------------------------
subroutine meridional_advection_interior(part)
  type (Partition) :: part
  integer :: i
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: l

  call interior_lat_indices(i_start, i_end, part%grid)

  if (i_start > i_end) return
  l = 0
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, part%grid, i)
    ! Interior only
    j_start = j_start + 1
    j_end = j_end - 1

    l = l + j_end - j_start + 1
    if (j_start > j_end) return
  end do

end subroutine

!-------------------------------------------------------------------------------
!> Meridional advection for border cells (which do need ghost cells)
!-------------------------------------------------------------------------------
subroutine meridional_advection_border(part)
  type (Partition) :: part
  integer :: i
  integer :: i_start, i_end
  integer :: j_start, j_end


  call interior_lat_indices(i_start, i_end, part%grid)

  if (i_start > i_end) return
  do i = i_start, i_end
    call interior_lon_indices(j_start, j_end, part%grid, i)
    ! Interior only
    j_start = j_start + 1
    j_end = j_end - 1

    if (j_start > j_end) return
  end do

end subroutine


!-----------------------------------------------------------------------------
!> Sleep instead of chemistry
!-----------------------------------------------------------------------------
subroutine chemistry(distrib, nb)
  integer, intent(in) :: nb
  type(Partition), pointer :: part
  type(Partitioning), pointer :: distrib
  double precision :: t_advec
  integer :: i, j, k
  integer :: i_start, i_end
  integer :: j_start, j_end
  integer :: ierr

  do k = 1, nb
    part => distrib%parts(k)
    i_start = get_first_lat_Partition(part)
    i_end = get_last_lat_Partition(part)
    do i = i_start, i_end
      j_start = get_first_lon_Partition(i, part)
      j_end = get_last_lon_Partition(i, part)
      do j = j_start, j_end
#ifdef WITH_PGF90
        call system("usleep 12500")
#else
        ierr = system("usleep 12500")
#endif
      end do
    end do
  end do
end subroutine


!-----------------------------------------------------------------------------
!> Wait for a request to finish (blocking)
!-----------------------------------------------------------------------------
subroutine wait_request(request)
  integer :: request
  integer :: ierr, status(mpi_status_size)

  if (request /= mpi_request_null) then
    call mpi_wait(request, status, ierr)
    !write (*,*) "request done"
  end if
end subroutine

!> Count interior and border cells
!> Gather all data to master proc which write it to a file
subroutine count_halo_interior(part)
  type (Partition) :: part
  integer :: nb(9) = 0.
  integer, allocatable :: allnb(:)
  integer :: rank, ierr
  character(LINE_WIDTH) :: fname
  integer :: id, nb_procs
  integer :: i, j, n


  nb(1) = count_zonal_tracer_fluxes(part, interior=.True.)
  nb(2) = count_zonal_tracer_fluxes(part, interior=.False.)
  nb(3) = count_zonal_gradient_interior(part%grid)
  nb(4) = count_zonal_gradient_border(part%grid)
  nb(5) = count_merid_tracer_fluxes(part, interior=.True.)
  nb(6) = count_merid_tracer_fluxes(part, interior=.False.) 
  nb(7) = count_merid_slope(part, interior=.True.)
  nb(8) = count_merid_slope(part, interior=.False.)
  nb(9) = nb_cells_Band_grid(part%grid)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)

  n = size(nb)

  ! Intel does not want the array to be allocated only for master proc
  !if (rank == 0) then
    call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
    allocate(allnb(nb_procs*n))
  !end if

  call mpi_gather(nb, n, mpi_integer, allnb, n, mpi_integer,&
    0, mpi_comm_world, ierr)

  if (rank == 0) then
    id = 2
    fname = nbcells_filename(rank, nb_procs)
    open(unit=id, file=trim(fname), action="write")

    write (id, '(a)') "#proc int_flux_z ghost_flux_z int_grad_m ghost_grad_m"//&
      "int_flux_m ghost_flux_m int_slope_m ghost_slope_m"
    do i = 1, nb_procs
      write (id, '(i5, 9i8)') i-1, (allnb((i-1)*n + j), j=1, n)
    end do

    close(id)
    deallocate(allnb)
  end if
end subroutine

function nbcells_filename(rank, nb_procs) result (fname)
  integer, intent(in) :: rank, nb_procs
  character(LINE_WIDTH) :: fname
  character(8) :: tmp

  write (tmp, '(i8)') nb_procs
  fname = trim(get_output_dir())//"/nb_cells"//adjustl(trim(tmp))

end function


end module
