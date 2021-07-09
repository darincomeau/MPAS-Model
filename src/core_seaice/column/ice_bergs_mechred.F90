!  SVN:$Id: $
!=======================================================================

! Ridging of sea ice by advancing ice shelves or icebergs
! based on ice_mechred.F90
!
! authors: Elizabeth C. Hunke, LANL

      module ice_bergs_mechred

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, c10, p5, &
          puny, Lfresh, rhoi, rhos
      use ice_itd, only: column_sum, &
                         column_conservation_check
      use ice_mechred, only: asum_ridging, ridge_itd, ridge_shift, ridge_prep
      use ice_colpkg_tracers, only: nt_qice, nt_qsno, nt_fbri, nt_sice
      use ice_warnings, only: add_warning

      implicit none
      save

      private
      public :: ridge_ice_by_bergs

      logical (kind=log_kind), parameter :: &
         l_conservation_check = .false.  ! if true, check conservation
!          l_conservation_check = .true.   ! if true, check conservation

!=======================================================================

      contains

!=======================================================================

! Compute changes in the sea ice thickness distribution due to icebergs.
! Based on ice_mechred.F90
!
! NOTE: This subroutine operates over a single block.
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!   2008: Elizabeth Hunke modified for icebergs.

      subroutine ridge_ice_by_bergs (dt,           dice,        &
                                     ncat,         n_aero,      &
                                     nilyr,        nslyr,       &
                                     ntrcr,        hin_max,     &
                                     aicen,        btrcrn,      &
                                     vicen,        vsnon,       &
                                     aice0,        bice,        &
                                     trcr_depend,  trcr_base ,  &
                                     n_trcr_strata,             &
                                     nt_strata,    l_stop,      &
                                     stop_label,                &
                                     krdg_partic,  krdg_redist, &
                                     mu_rdg,       tr_brine,    &
                                     dardg1dt,     dardg2dt,    &
                                     dvirdgdt,     opening,     &
                                     fpond,                     &
                                     fresh,        fhocn,       &
                                     faero_ocn,                 &
                                     aparticn,     krdgn,       &
                                     aredistn,     vredistn,    &
                                     dardg1ndt,    dardg2ndt,   &
                                     dvirdgndt,                 &
                                     araftn,       vraftn)

!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt     , & ! time step
         dice   , & ! fractional area to be ridged (displaced)
         bice   , & ! fractional berg (and shelf) area
         mu_rdg     ! gives e-folding scale of ridged ice (m^.5)

      integer (kind=int_kind), intent(in) :: &
         ncat   , & ! number of thickness categories
         n_aero , & ! number of aerosol tracers
         nilyr  , & ! number of ice layers
         nslyr  , & ! number of snow layers
         ! nblyr  , & ! number of bio layers
         ntrcr      ! number of tracers in use

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max   ! category limits (m)

      real (kind=dbl_kind), dimension (ncat), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         btrcrn             ! temporary array for tracers
 
      real (kind=dbl_kind), intent(inout) :: & 
         aice0     ! concentration of open water

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      logical (kind=log_kind), intent(out) :: &
         l_stop   ! if true, abort on return

      character (char_len), intent(out) :: &
         stop_label     ! diagnostic information for abort

      integer (kind=int_kind), intent(in) :: &
         krdg_partic  , & ! selects participation function
         krdg_redist      ! selects redistribution function 

      logical (kind=log_kind), intent(in) :: &
         tr_brine       ! if .true., brine height differs from ice thickness

      ! optional history fields
      real (kind=dbl_kind), intent(inout), optional :: &
         dardg1dt  , & ! rate of fractional area loss by ridging ice (1/s)
         dardg2dt  , & ! rate of fractional area gain by new ridges (1/s)
         dvirdgdt  , & ! rate of ice volume ridged (m/s)
         opening   , & ! rate of opening due to divergence/shear (1/s)
         fpond     , & ! fresh water flux to ponds (kg/m^2/s)
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn         ! net heat flux to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         dardg1ndt  , & ! rate of fractional area loss by ridging ice (1/s)
         dardg2ndt  , & ! rate of fractional area gain by new ridges (1/s)
         dvirdgndt  , & ! rate of ice volume ridged (m/s)
         aparticn   , & ! participation function
         krdgn      , & ! mean ridge thickness/thickness of ridging ice
         araftn     , & ! rafting ice area
         vraftn     , & ! rafting ice volume 
         aredistn   , & ! redistribution function: fraction of new ridge area
         vredistn       ! redistribution function: fraction of new ridge volume

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         faero_ocn      ! aerosol flux to ocean (kg/m^2/s)

      ! local variables

      real (kind=dbl_kind), dimension (ncat) :: &
         eicen      , & ! energy of melting for each ice layer (J/m^2)
         esnon      , & ! energy of melting for each snow layer (J/m^2)
         vbrin      , & ! ice volume with defined by brine height (m)
         sicen          ! Bulk salt in h ice (ppt*m)

      ! variables for ridging routines
      real (kind=dbl_kind) :: &
         asum       , & ! sum of ice and open water area
         aksum      , & ! ratio of area removed to area ridged
         msnow_mlt  , & ! mass of snow added to ocean (kg m-2)
         esnow_mlt  , & ! energy needed to melt snow in ocean (J m-2)
         mpond      , & ! mass of pond added to ocean (kg m-2)
         closing_net, & ! net rate at which area is removed    (1/s)
                        ! (ridging ice area - area of new ridges) / dt
         divu_adv   , & ! divu as implied by transport scheme  (1/s)
         opning     , & ! rate of opening due to divergence/shear
                        ! opning is a local variable;
                        ! opening is the history diagnostic variable
         ardg1      , & ! fractional area loss by ridging ice
         ardg2      , & ! fractional area gain by new ridges
         virdg      , & ! ice volume ridged
         aopen          ! area opening due to divergence/shear

      real (kind=dbl_kind), dimension (n_aero) :: &
         maero          ! aerosol mass added to ocean (kg m-2)

      real (kind=dbl_kind), dimension (0:ncat) :: &
         apartic          ! participation function; fraction of ridging
                          ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (ncat) :: &
         hrmin        , & ! minimum ridge thickness
         hrmax        , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp        , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg         , & ! mean ridge thickness/thickness of ridging ice
         ardg1n       , & ! area of ice ridged
         ardg2n       , & ! area of new ridges
         virdgn       , & ! ridging ice volume
         mraftn           ! rafting ice mask

      real (kind=dbl_kind) :: &
         vice_init, vice_final, & ! ice volume summed over categories
         vsno_init, vsno_final, & ! snow volume summed over categories
         eice_init, eice_final, & ! ice energy summed over layers
         vbri_init, vbri_final, & ! ice volume in fbri*vicen summed over categories
         sice_init ,sice_final, & ! ice bulk salinity summed over categories
         esno_init, esno_final    ! snow energy summed over layers

      integer (kind=int_kind), parameter :: &
         nitermax = 20    ! max number of ridging iterations

      integer (kind=int_kind) :: &
         k            , & ! layer index
         n            , & ! thickness category index
         niter        , & ! iteration counter
         it               ! tracer index

      real (kind=dbl_kind) :: &
         dti              ! 1 / dt

      logical (kind=log_kind) :: &
         iterate_ridging, & ! if true, repeat the ridging
         asum_error         ! flag for asum .ne. 1

      character (len=char_len) :: &
         fieldid            ! field identifier

      character(len=char_len_long) :: &
         warning ! warning message

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      msnow_mlt = c0
      esnow_mlt = c0
      maero(:)  = c0
      ardg1     = c0
      ardg2     = c0
      virdg     = c0
      ardg1n(:) = c0
      ardg2n(:) = c0
      virdgn(:) = c0
      mpond     = c0
      aopen     = c0

      divu_adv    = c0
      closing_net = c0
      opning      = c0
      asum = c0

      !-----------------------------------------------------------------
      ! Compute initial values of conserved quantities. 
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         do n = 1, ncat
         eicen(n) = c0
         esnon(n) = c0
         sicen(n) = c0
         vbrin(n) = c0

         do k = 1, nilyr
            eicen(n) = eicen(n) + btrcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo

         do k = 1, nslyr
            esnon(n) = esnon(n) + btrcrn(nt_qsno+k-1,n) &
                     * vsnon(n)/real(nslyr,kind=dbl_kind)
         enddo ! k

         vbrin(n) = vicen(n)
         if (tr_brine) vbrin(n) =  btrcrn(nt_fbri,n) * vicen(n)

         enddo ! n

         call column_sum (ncat, vicen(:), vice_init)
         call column_sum (ncat, vsnon(:), vsno_init)
         call column_sum (ncat, eicen(:), eice_init)
         call column_sum (ncat, esnon(:), esno_init)
         call column_sum (ncat, sicen(:), sice_init)
         call column_sum (ncat, vbrin(:), vbri_init)

      endif ! conservation check

      do niter = 1, nitermax

      !-----------------------------------------------------------------
      ! Compute the thickness distribution of ridging ice
      ! and various quantities associated with the new ridged ice.
      !-----------------------------------------------------------------

         call ridge_itd (ncat,        aice0,       &
                         aicen(:),    vicen(:),    &
                         krdg_partic, krdg_redist, &
                         mu_rdg,                   &
                         aksum,       apartic(:),  &
                         hrmin(:),    hrmax(:),    &
                         hrexp(:),    krdg(:),     &
                         aparticn(:), krdgn(:),    &
                         mraftn(:))

      !-----------------------------------------------------------------
      ! Compute the area opening and closing.
      !-----------------------------------------------------------------

         closing_net = dice/dt  ! net closing rate
         divu_adv = -closing_net
         opning = c0

      !-----------------------------------------------------------------
      ! Redistribute area, volume, and energy.
      !-----------------------------------------------------------------

         call ridge_shift (ntrcr,       dt,            &
                           ncat,        hin_max(:),    &
                           aicen(:),    btrcrn(:,:),   &
                           vicen(:),    vsnon(:),      &
                           aice0,       trcr_depend,   &
                           trcr_base,   n_trcr_strata, &
                           nt_strata,   krdg_redist,   &
                           aksum,       apartic(:),    &
                           hrmin(:),    hrmax(:),      &
                           hrexp(:),    krdg(:),       &
                           closing_net, opning,        &
                           ardg1,       ardg2,         &
                           virdg,       aopen,         &
                           ardg1n(:),   ardg2n(:),     &
                           virdgn(:),                  &
                           nslyr,       n_aero,        &
                           msnow_mlt,   esnow_mlt,     &
                           maero(:),    mpond,         &
                           l_stop,      stop_label,    &
                           aredistn(:), vredistn(:))
         if (l_stop) return

      !-----------------------------------------------------------------
      ! Make sure the new area = 1.  If not (because the closing
      ! and opening rates were reduced above), prepare to ridge again
      ! with new rates.
      !-----------------------------------------------------------------

         ! dcmod - adjust with berg area
         call asum_ridging (ncat, aicen(:), aice0, asum)
         asum = asum + bice

         ! dcmod - if less than 1 due to forced ridging, add missing area back to open water
         if (asum < c1) aice0 = aice0 + (c1 - asum)

         ! check sum again
         call asum_ridging (ncat, aicen(:), aice0, asum)
         asum = asum + bice

         if (abs(asum - c1) < puny) then
            iterate_ridging = .false.
            closing_net = c0
            opning      = c0
         else
            iterate_ridging = .true.
            divu_adv = (c1 - asum) / dt
            closing_net = max(c0, -divu_adv)
            opning = max(c0, divu_adv)
         endif

      !-----------------------------------------------------------------
      ! If done, exit.  If not, prepare to ridge again.
      !-----------------------------------------------------------------

         if (iterate_ridging) then
            write(warning,*) 'Repeat ridging (bergs), niter =', niter
            call add_warning(warning)
            ! ! dcmod debugging messages
            ! write(warning,*) 'area (inc. bergs):', asum
            ! call add_warning(warning)
            ! write(warning,*) 'ice:', asum - bice - aice0
            ! call add_warning(warning)
            ! write(warning,*) 'open water:', aice0
            ! call add_warning(warning)
         else
            ! exit rdg_iteration
            exit
         endif

         if (niter == nitermax) then
            write(warning,*) ' '
            call add_warning(warning)
            write(warning,*) 'Exceeded max number of ridging iterations'
            call add_warning(warning)
            write(warning,*) 'max =',nitermax
            call add_warning(warning)
            l_stop = .true.
            stop_label = "ridge_ice: Exceeded max number of ridging iterations"
            return
         endif

      enddo    ! niter

      !-----------------------------------------------------------------
      ! Compute final values of conserved quantities. 
      ! Check for conservation (allowing for snow thrown into ocean).
      !-----------------------------------------------------------------

      if (l_conservation_check) then

         do n = 1, ncat
         eicen(n) = c0
         esnon(n) = c0
         sicen(n) = c0
         vbrin(n) = c0

         do k = 1, nilyr
            eicen(n) = eicen(n) + btrcrn(nt_qice+k-1,n) &
                     * vicen(n)/real(nilyr,kind=dbl_kind)
         enddo
         do k = 1, nslyr
            esnon(n) = esnon(n) + btrcrn(nt_qsno+k-1,n) &
                     * vsnon(n)/real(nslyr,kind=dbl_kind)
         enddo

         vbrin(n) =  vicen(n)
         if (tr_brine)  vbrin(n) =  btrcrn(nt_fbri,n) * vbrin(n)

         enddo ! ncat

         call column_sum (ncat, vicen(:), vice_final)
         call column_sum (ncat, vsnon(:), vsno_final)
         call column_sum (ncat, eicen(:), eice_final)
         call column_sum (ncat, esnon(:), esno_final)
         call column_sum (ncat, sicen(:), sice_final)
         call column_sum (ncat, vbrin(:), vbri_final)

         vsno_final = vsno_final + msnow_mlt/rhos
         esno_final = esno_final + esnow_mlt

         fieldid = 'vice, ridging'
         call column_conservation_check (fieldid,               &
                                         vice_init, vice_final, &
                                         puny,                  &
                                         l_stop)
         if (l_stop) return

         fieldid = 'vsno, ridging'
         call column_conservation_check (fieldid,               &
                                         vsno_init, vsno_final, &
                                         puny,                  &
                                         l_stop)
         if (l_stop) return

         fieldid = 'eice, ridging'
         call column_conservation_check (fieldid,               &
                                         eice_init, eice_final, &
                                         puny*Lfresh*rhoi,      &
                                         l_stop)
         if (l_stop) return

         fieldid = 'esno, ridging'
         call column_conservation_check (fieldid,               &
                                         esno_init, esno_final, &
                                         puny*Lfresh*rhos,      &
                                         l_stop)
         if (l_stop) return

      endif ! l_conservation_check

      !-----------------------------------------------------------------
      ! Compute ridging diagnostics.
      !-----------------------------------------------------------------

      dti = c1/dt

      if (present(dardg1dt)) then
         dardg1dt = ardg1*dti
      endif
      if (present(dardg2dt)) then
         dardg2dt = ardg2*dti
      endif
      if (present(dvirdgdt)) then
         dvirdgdt = virdg*dti
      endif
      if (present(opening)) then
         opening = aopen*dti
      endif

      if (present(dardg1ndt)) then
         do n = 1, ncat
            dardg1ndt(n) = ardg1n(n)*dti
         enddo
      endif
      if (present(dardg2ndt)) then
         do n = 1, ncat
            dardg2ndt(n) = ardg2n(n)*dti
         enddo
      endif
      if (present(dvirdgndt)) then
         do n = 1, ncat
            dvirdgndt(n) = virdgn(n)*dti
         enddo
      endif
      if (present(araftn)) then
         do n = 1, ncat
            araftn(n) = mraftn(n)*ardg2n(n)
!            araftn(n) = mraftn(n)*ardg1n(n)*p5
         enddo
      endif
      if (present(vraftn)) then
         do n = 1, ncat
            vraftn(n) = mraftn(n)*virdgn(n)
         enddo
      endif

      !-----------------------------------------------------------------
      ! Update fresh water and heat fluxes due to snow melt.
      !-----------------------------------------------------------------

      if (present(fresh)) then
         fresh = fresh + msnow_mlt*dti
      endif
      if (present(fhocn)) then
         fhocn = fhocn + esnow_mlt*dti
      endif
      if (present(faero_ocn)) then
         do it = 1, n_aero
            faero_ocn(it) = faero_ocn(it) + maero(it)*dti
         enddo
      endif
      if (present(fpond)) then
         fpond = fpond - mpond ! units change later
      endif

      !-----------------------------------------------------------------
      ! Check for fractional ice area > 1.
      !-----------------------------------------------------------------

      call asum_ridging (ncat, aicen(:), aice0, asum)
      ! dcmod - adjust with berg area
      asum = asum + bice
      if (abs(asum - c1) > puny) then
         l_stop = .true.

         write(warning,*) ' '
         call add_warning(warning)
         write(warning,*) 'Ridging error (bergs): total area > 1'
         call add_warning(warning)
         write(warning,*) 'area:', asum
         call add_warning(warning)
         write(warning,*) 'n, aicen:'
         call add_warning(warning)
         write(warning,*)  0, aice0
         call add_warning(warning)
         do n = 1, ncat
            write(warning,*) n, aicen(n)
            call add_warning(warning)
         enddo
         return
      endif

      end subroutine ridge_ice_by_bergs

!=======================================================================

      end module ice_bergs_mechred

!=======================================================================
