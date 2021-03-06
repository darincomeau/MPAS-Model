!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_berg_advection
!
!> \brief   Iceberg model advection
!> \author  Darin Comeau, LANL
!> \date    29 Nov 2017
!> \details This module is intended to advect the iceberg model.
!
!-----------------------------------------------------------------------

module seaice_berg_advection

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
      seaice_init_berg_tracers, &
      seaice_run_berg_advection

contains

!-----------------------------------------------------------------------
! Set up iceberg tracers
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_berg_tracers
!
!> \brief   Initialize iceberg tracers
!> \author  Darin Comeau, LANL
!> \date    29 Nov 2017
!> \details This routine is intended to construct iceberg tracers.
!
!-----------------------------------------------------------------------

  subroutine seaice_init_berg_tracers(domain)

    use seaice_advection_incremental_remap_tracers, only: &
         seaice_add_berg_tracers_to_linked_list, &
         bergTracersHeadArray

    type(domain_type), intent(inout) :: &
         domain

    integer, pointer :: &
         nBergCategories

    integer :: iCategory

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBergCategories", nBergCategories)

    allocate(bergTracersHeadArray(nBergCategories))

    do iCategory = 1, nBergCategories

       call seaice_add_berg_tracers_to_linked_list(domain, &
                                                   bergTracersHeadArray(iCategory) % tracerTypePtr, &
                                                   iCategory)

    enddo

  end subroutine seaice_init_berg_tracers

!-----------------------------------------------------------------------
! Iceberg advection
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_run_berg_advection
!
!> \brief   Advect icebergs
!> \author  Darin Comeau, LANL
!> \date    29 Nov 2017
!> \details This routine is intended to advect icebergs by incremental
!> remapping (currently only option).
!
!-----------------------------------------------------------------------

  subroutine seaice_run_berg_advection(domain, clock, ierr)

    use seaice_advection_incremental_remap, only: &
         seaice_run_advection_incremental_remap

    use seaice_advection_incremental_remap_tracers, only: &
         bergTracersHeadArray

    type(domain_type), intent(inout) :: &
         domain

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(inout) :: &
         ierr  !< Input/Output:

    integer, pointer :: &
         nBergCategories

    integer :: iCategory

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nBergCategories", nBergCategories)

    do iCategory = 1, nBergCategories

       call seaice_run_advection_incremental_remap(domain, &
                                                   bergTracersHeadArray(iCategory) % tracerTypePtr, &
                                                   clock, &
                                                   ierr, &
                                                   updateHaloInitIn = .true., &
                                                   updateHaloFinalIn = .true.)

    enddo

  end subroutine seaice_run_berg_advection

!-----------------------------------------------------------------------

end module seaice_berg_advection
