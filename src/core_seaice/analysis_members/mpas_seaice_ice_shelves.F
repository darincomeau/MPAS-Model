! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_ice_shelves
!
!> \brief MPAS sea ice analysis mode member: ice_shelves
!> \author Adrian K. Turner
!> \date   9th April 2016
!> \details
!>  MPAS sea ice analysis mode member: ice_shelves
!>
!-----------------------------------------------------------------------

module seaice_ice_shelves

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: seaice_bootstrap_ice_shelves, &
             seaice_init_ice_shelves, &
             seaice_precompute_ice_shelves, &
             seaice_compute_ice_shelves, &
             seaice_restart_ice_shelves, &
             seaice_finalize_ice_shelves

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_ice_shelves
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th April 2016
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_ice_shelves(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_bootstrap_ice_shelves!}}}

!***********************************************************************
!
!  routine seaice_init_ice_shelves
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th April 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_ice_shelves(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_init_ice_shelves!}}}

!***********************************************************************
!
!  routine seaice_precompute_ice_shelves
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th April 2016
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_ice_shelves(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_precompute_ice_shelves!}}}

!***********************************************************************
!
!  routine seaice_compute_ice_shelves
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th April 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_ice_shelves(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           iceShelvesAMPool, &
           oceanCouplingPool, &
           tracersAggregatePool, &
           meshPool

      real(kind=RKIND), pointer :: &
           iceAreaOverIceShelves

      real(kind=RKIND), dimension(:), pointer :: &
           iceAreaCell, &
           areaCell

      integer, dimension(:), pointer :: &
           landIceMask

      real(kind=RKIND) :: &
           iceAreaOverIceShelvesThisProc

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell

      iceAreaOverIceShelvesThisProc = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
         call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

         call MPAS_pool_get_array(oceanCouplingPool, "landIceMask", landIceMask)
         call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         do iCell = 1, nCellsSolve

            if (landIceMask(iCell) == 1) then

               iceAreaOverIceShelvesThisProc = iceAreaOverIceShelvesThisProc + iceAreaCell(iCell) * areaCell(iCell)

            endif ! landIceMask

         enddo ! iCell

         block => block % next
      enddo

      ! sum over all processors
      call MPAS_pool_get_subpool(domain % blocklist % structs, "iceShelvesAM", iceShelvesAMPool)

      call MPAS_pool_get_array(iceShelvesAMPool, "iceAreaOverIceShelves", iceAreaOverIceShelves)

      call MPAS_dmpar_sum_real(domain % dminfo, iceAreaOverIceShelvesThisProc, iceAreaOverIceShelves)

   end subroutine seaice_compute_ice_shelves!}}}

!***********************************************************************
!
!  routine seaice_restart_ice_shelves
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th April 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_ice_shelves(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_restart_ice_shelves!}}}

!***********************************************************************
!
!  routine seaice_finalize_ice_shelves
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th April 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_ice_shelves(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_finalize_ice_shelves!}}}

!-----------------------------------------------------------------------

end module seaice_ice_shelves

! vim: foldmethod=marker
