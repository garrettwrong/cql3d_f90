module achiefn_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use cfpcoefc_mod, only : cfpcoefc
  use cfpgamma_mod, only : cfpgamma
  use diag_mod, only : diag
  use diaggnde_mod, only : diaggnde
  use dskout_mod, only : dskout
  use efield_mod, only : efield
  use exsweep_mod, only : exsweep
  use finit_mod, only : finit
  use impavnc0_mod, only : impavnc0
  use ntdstore_mod, only : ntdstore
  use ntloop_mod, only : ntloop
  use pltmain_mod, only : pltmain
  use pltrun_mod, only : pltrun
  use restvty_mod, only : restvty
  use sourcee_mod, only : sourcee
  use tdoutput_mod, only : tdoutput

  !---END USE

!
!

contains

  subroutine achiefn(kopt)
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      use impavnc0_mod, only : impavnc0
      use pltmain_mod, only : pltmain

      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     kopt=0: This routine advances the  equations in momemtum space
!             over time-step n (or if setup0%cqlpmod.eq.'disabled',
!             transp='enabled',soln_method.eq.'it3drv', then only
!             calc coeffs until at last flux surface.
!     kopt=1: Compute plasma density and energy transfer
!             based on new f (after solution is obtained);
!             when kopt=1, impavnc0 is not called.
!     kopt=2: Calculated h and g functions and associated fluxes,
!             to obtain electrostatic electric field for
!             setup0%cqlpmod="enabled",sbdry="periodic",esfld="enabled".
!             (BH:Evidently, have not fully implemented this yet.
!             No kopt reference is in impanvc0.  It was
!             successfully implemented in the STELLA code (a derivative
!             of Kupfer's FPET code).  See achiefn in stella.f.
!             Also check file Harvey_rf99.pdf, or related paper, from
!             Radiofrequency Power in Plasmas meeting, Anapolis, 1999.)
!     kopt=3: Calculate h and g functions, for time-advancing of electric
!             field according to Ampere-Faraday equations.
!             Iteration on Kupfer's g and h functions is used in
!             obtaining the implicit toroidal electric field.
!             BH131104, with help from YuP.
!
!..................................................................

#ifdef __MPI
!MPI >>>
      include 'mpilib.h'
!MPI <<<
#endif

      character cptline*80
      integer getpid
      real :: cputime


      if (kopt.eq.0  .or.  kopt.eq.2) then

!..................................................
!     bring background profiles up to time step n
!..................................................
!cc      if(nefiter.eq.1) call profiles !-YuP: moved to tdchief

!.....................................................................
!     compute electric field either according to direct specification,
!     or to maintain a target toroidal current.
!.....................................................................
      if (eseswtch.ne."enabled") then
         call efield
      endif

        ! Start time advancement:
      if(nefiter.eq.1) then
      if (kopt.eq.0) then
           n=n+1
           n_(l_)=n ! new time-step for this flux surface
           ! for 2-d (v_par,v_perp) calculation ntloop controls
           ! end of run or restart.
           ! Also updates time.
           call ntloop
      endif
      endif

!..................................................
!     determine coulomb log...
!..................................................
      if(nefiter.eq.1) call cfpgamma

!....................................................
!     Apply additional preloading to f, if n=nfpld
!     (If only initial loading is used, then nfpld=0)
!....................................................
      if (n.eq.nfpld .and. nefiter.eq.1) then
          call finit ! Note: at n=0 finit is called from ainitial.
      endif          ! First call to achiefn is with n=0, then advanced.

!..................................................
!     determine particle sources...
!..................................................
      if (n.ne.1 .and. nefiter.eq.1) call sourcee
      ! Note: at n=0 sourcee is called from ainitial
      ! First call to achiefn is with n=0, then advanced in tdchief.


#ifdef __MPI
!MPI >>>
      if(mpirank.eq.mpiworker) then

!MPI <<<
#endif
!..................................................
!     generate collisional F.P. coefficients...
!..................................................
       call cfpcoefc ! time-consuming calculations
!..................................................................
!     call implicit time advancement routine if implct.eq."enabled"
!..................................................................
      if (implct .eq. "enabled") call impavnc0(kopt)
#ifdef __MPI
!MPI >>>
      endif  ! for if(mpirank.eq.***)
!MPI <<<
#endif

!...........................................................
!     call splitting time advancement if implct = "disabled"
!...........................................................
      if (implct .ne. "enabled") call exsweep

      go to 999 !-> return/end

      endif ! kopt=0 or 2

!..................................................................
!     Simultaneous FP and Ampere-Faraday Eqns solution, kopt=3
!     Following above coding for kopt=0.
!..................................................................
      if (kopt.eq.3) then
         n=n+1  ! n is used by logics in cfpcoefn, impavnc0, etc.
         zdttot=one
         call cfpgamma
         if (adimeth.eq."enabled" .and. transp.eq."enabled" &
              .and. n.ge.nonadi) zdttot=2.0
         timet=timet+zdttot*dtreff
        ! n and time are not really updated yet -
        ! n and time are updated after n_ and time_ are updated.
        call cfpgamma  !determine coulomb log...   But called 6 lines above?
!BH131107:  NEED to check this n=0 issue.
        ! Apply additional preloading to f, if n=nfpld
        ! (If only initial loading is used, then nfpld=0)
        if (n.eq.nfpld) then
               call finit
        endif
        ! Note: at n=0 finit is called from ainitial.
        ! First call to achiefn is with n=0, then advanced.
        if (n.ne.1) call sourcee !determine particle sources...
        ! Note: at n=0 sourcee is called from ainitial
        ! First call to achiefn is with n=0, then advanced in tdchief.

! For ZOW, parallelization (in lr) is done here.
#ifdef __MPI
!MPI >>>
      if(mpirank.eq.mpiworker) then

!MPI <<<
#endif
          !collisional F.P. coefficients:
          call cfpcoefc ! time-consuming calculations
          !call implicit time advancement:
          if (implct .eq. "enabled") call impavnc0(kopt)
#ifdef __MPI
!MPI >>>
      endif  ! for if(mpirank.eq.***)
!MPI <<<
#endif
!BH110309:  Restore n
!BH110309:        n=n-1
      endif  !On kopt.eq.3








      if (kopt.eq.1) then
!..................................................
!     Compute plasma energy and density and energy transfer
!     rate terms...
!     Diagnostics ONLY COMPUTED for cases where transport
!     is not considered (transp="disabled")
!..................................................

      if (transp.eq."disabled" .or. &
        (setup0%cqlpmod.eq."enabled".and.n.ge.nontran)) then
        call diaggnde
!..................................................
!     Compute plasma resistivity for electron runs..
!..................................................
        call restvty
!..................................................
!     Compute various diagnostics...
!..................................................
        call diag
!..................................................
!     Obtain data for time dependent plots
!..................................................
        call ntdstore
      endif  !On transp.eq."disabled" ....

!..................................................
!     plotting logic for 2-D (v-theta) code...
!..................................................
      iplot=0  !Index for nonzero values of increasing nplot(1:nplota)
      nplott=0 !Total nplot(i).ge.0 and .le.nstop
      do i=1,nplota
         if (nplot(i).ge.0 .and. nplot(i).le.nstop) &
             nplott=nplott+1 !YuP:could do this counting in ainsetva?
      enddo
      if (setup0%noplots.ne."enabled1") then
         do i=1,nplott
            if(n.eq.nplot(i)) then
               iplot=i
               tplot(i)=timet ! to be printed in plots (sub.pltrun)
            endif
         enddo
      endif

      if (setup0%lrzmax.eq.1) then
        if (n.eq.nstop .or. iplot.ne.0) call tdoutput(2)
        if (n.eq.nstop .or. iplot.ne.0) call pltmain
        if (n.eq.nstop.and.pltra.eq."enabled") then
           call pltrun
        endif
!cc        if (netcdfnm.ne."disabled" .and. nstop.ne.0) then
!cc           call netcdfrw2(1) !YuP: Why here? Also called from tdchief
!cc        endif
        if(n.eq.nstop) then
           call dskout(l_)
           call pgend
!BH080106           cputime=etime(tarray)
           call cpu_time(cputime)
           write (*,*) 'achiefn/setup0%lrzmax=1: CPU time (seconds)', cputime
           if (ichkpnt.ne."disabled") then
             write(cptline,100) i,ichkpnt
 100         format("chkpnt -p ",i5," -f ",a8)
             print *,cptline
           endif
        endif ! n.eq.nstop
      endif ! setup0%lrzmax=1

!..................................................
!     plotting (and netcdf store) logic for 3-D code ..
!     Used if transp="disabled".
!..................................................
      if (transp .eq. "disabled" .and. setup0%noplots.ne."enabled1") then
!$$$      if (transp .eq. "disabled") then
         if (setup0%lrzmax.gt.1) then
            if (n.ge.nstop) then
               call pltmain
               call netcdfmain
               go to 999 !-> return/end
            endif
            if (iplot.ne.0) call pltmain
         endif
      endif

      go to 999 !-> return/end

      endif ! if kopt=1




 999  return

      end subroutine achiefn


end module achiefn_mod
