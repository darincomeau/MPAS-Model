!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_bergs
!
!> \brief   Iceberg model
!> \author  Darin Comeau, LANL
!> \date    29 Sept 2017
!> \details This module is intended to interface the iceberg model with 
!> MPAS-seaice.
!
!-----------------------------------------------------------------------

module seaice_bergs

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_log, only: mpas_log_write

   use mpas_stream_manager
   use mpas_io_units

   implicit none

   private
   save

   public :: &
      seaice_bergs_init, &
      seaice_run_berg_predynamics, &
      seaice_run_berg_dynamics, &
      seaice_run_berg_postdynamics

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_bergs_init
!
!> \brief   Initialize iceberg model         
!> \author  Darin Comeau, LANL
!> \date    19 May 2017
!> \details This routine is intended to initialize the iceberg model if
!> icebergs are turned on based on the namelist option config_use_bergs.
!
!-----------------------------------------------------------------------

  subroutine seaice_bergs_init(domain)

    use seaice_berg_calving, only: &
         seaice_init_berg_calving_fluxes

    use seaice_berg_velocity_solver, only: &
         seaice_init_berg_velocity_solver

    use seaice_berg_advection, only: &
         seaice_init_berg_tracers

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    logical, pointer :: &
         config_use_bergs, & ! flag to turn on / off iceberg model
         config_do_restart   ! flag if restart

    call MPAS_pool_get_config(domain % configs, "config_use_bergs", config_use_bergs)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)


    if (config_use_bergs) then

       call mpas_timer_start("Berg initialization")

       ! initial run
       if (.not. config_do_restart) then

          ! initialize calving fluxes
          call seaice_init_berg_calving_fluxes(domain)

       endif

       ! initialize velocity solver
       call seaice_init_berg_velocity_solver(domain)

       ! initialize tracers
       call seaice_init_berg_tracers(domain)

       call mpas_timer_stop("Berg initialization")

    endif

  end subroutine seaice_bergs_init

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_run_berg_predynamics
!
!> \brief   Wrapper for pre-dynamics routines
!> \author  Darin Comeau, LANL
!> \date    07 Feb 2018
!> \details This subroutine is a wrapper for the iceberg timestep before
!> dynamics.
!
!-----------------------------------------------------------------------

  subroutine seaice_run_berg_predynamics(domain)

    use seaice_berg_state, only: &
         seaice_berg_init_displaced_area

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    logical, pointer :: &
         config_berg_update_open_water ! flag to turn on / off bergs affecting open water

    call MPAS_pool_get_config(domain % configs, "config_berg_update_open_water", config_berg_update_open_water)

    ! reset berg displaced area for timestep
    if (config_berg_update_open_water) call seaice_berg_init_displaced_area(domain)

  end subroutine seaice_run_berg_predynamics

! !-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_run_berg_dynamics
!
!> \brief   Wrapper for iceberg dynamics
!> \author  Darin Comeau, LANL
!> \date    29 Sept 2017
!> \details This subroutine is a wrapper for the iceberg velocity solver
!> over a sea ice time step.
!
!-----------------------------------------------------------------------

  subroutine seaice_run_berg_dynamics(domain, clock, ierr)

    use seaice_berg_velocity_solver, only: &
         seaice_run_berg_velocity_solver

    use seaice_berg_advection, only: &
         seaice_run_berg_advection

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(inout) :: &
         ierr  !< Input/Output:

    logical, pointer :: &
         config_use_berg_velocity_solver, & ! indicates if berg velocity solver is used
         config_use_berg_advection, &
         config_berg_advection_outside_velocity_subcycling

    integer, pointer :: &
         config_berg_velocity_subcycle_number ! berg dynamics subcycling number

    integer :: &
         iBergSubcycle

    ! determine if velocity solver switched on
    call MPAS_pool_get_config(domain % configs, "config_use_berg_velocity_solver", config_use_berg_velocity_solver)
    call MPAS_pool_get_config(domain % configs, "config_use_berg_advection", config_use_berg_advection)
    call MPAS_pool_get_config(domain % configs, "config_berg_velocity_subcycle_number", config_berg_velocity_subcycle_number)
    call MPAS_pool_get_config(domain % configs, "config_berg_advection_outside_velocity_subcycling", config_berg_advection_outside_velocity_subcycling)

    if (config_berg_advection_outside_velocity_subcycling) then ! advection outside subcycling

       do iBergSubcycle = 1, config_berg_velocity_subcycle_number

          ! berg velocity solver
          call mpas_timer_start("Berg velocity solver")
          if (config_use_berg_velocity_solver) call seaice_run_berg_velocity_solver(domain, clock)
          call mpas_timer_stop("Berg velocity solver")

       enddo ! iBergSubcycle

       ! advection by incremental remapping
       call mpas_timer_start("Berg advection")
       if (config_use_berg_advection) call seaice_run_berg_advection(domain, clock, ierr)
       call mpas_timer_stop("Berg advection")

    else ! advection inside subcycling

       do iBergSubcycle = 1, config_berg_velocity_subcycle_number

          ! berg velocity solver
          call mpas_timer_start("Berg velocity solver")
          if (config_use_berg_velocity_solver) call seaice_run_berg_velocity_solver(domain, clock)
          call mpas_timer_stop("Berg velocity solver")

          ! advection by incremental remapping
          call mpas_timer_start("Berg advection")
          if (config_use_berg_advection) call seaice_run_berg_advection(domain, clock, ierr)
          call mpas_timer_stop("Berg advection")

       enddo ! iBergSubcycle

    endif ! config_berg_advection_outside_velocity_subcycling

  end subroutine seaice_run_berg_dynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_run_berg_postdynamics
!
!> \brief   Wrapper for post-dynamics routines
!> \author  Darin Comeau, LANL
!> \date    29 Sept 2017
!> \details This subroutine is a wrapper for the iceberg velocity solver
!> over a sea ice time step.
!
!-----------------------------------------------------------------------

  subroutine seaice_run_berg_postdynamics(domain, clock, itimestep)

    use seaice_berg_decay, only: &
         seaice_run_berg_decay

    use seaice_berg_state, only: &
         seaice_berg_update_mass_and_area, &
         seaice_berg_accumulate_calving_mass, &
         seaice_berg_release_calving_mass, &
         seaice_berg_update_open_water_area

    use seaice_berg_velocity_solver, only: &
         seaice_berg_forcing_for_ice

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(in) :: &
         itimestep !< Input:

    logical, pointer :: &
         config_use_berg_decay, &         ! flag to turn on / off berg decay parameterization
         config_berg_update_open_water, & ! flag to turn on / off bergs affecting open water
         config_berg_seaice_interaction   ! flag to turn on / off bergs affecting open water

    call MPAS_pool_get_config(domain % configs, "config_use_berg_decay", config_use_berg_decay)
    call MPAS_pool_get_config(domain % configs, "config_berg_update_open_water", config_berg_update_open_water)
    call MPAS_pool_get_config(domain % configs, "config_berg_seaice_interaction", config_berg_seaice_interaction)

    ! berg mass and area state variables update
    call mpas_timer_start("Berg update mass")
    call seaice_berg_update_mass_and_area(domain)
    call mpas_timer_stop("Berg update mass")

    ! berg decay solver
    call mpas_timer_start("Berg decay")
    if (config_use_berg_decay) call seaice_run_berg_decay(domain, clock)
    call mpas_timer_stop("Berg decay")

    ! berg calving
    call mpas_timer_start("Berg calving")
    call seaice_berg_accumulate_calving_mass(domain)
    call seaice_berg_release_calving_mass(domain, itimestep)
    call mpas_timer_stop("Berg calving")

    ! berg sea ice interaction
    call mpas_timer_start("Berg sea ice interaction")
    if (config_berg_seaice_interaction) call seaice_berg_forcing_for_ice(domain)
    call mpas_timer_stop("Berg sea ice interaction")

    ! update open water
    if (config_berg_update_open_water) call seaice_berg_update_open_water_area(domain)

  end subroutine seaice_run_berg_postdynamics

!-----------------------------------------------------------------------

end module seaice_bergs
