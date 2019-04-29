module tdchief_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use achief1_mod, only : achief1
  use achiefn_mod, only : achiefn
  use aclear_mod, only : aclear
  use aindfpa_mod, only : ainadjnl
  use aindfpa_mod, only : ainadjnl_fsetup_setup0
  use aindfpa_mod, only : aindfpa
  use ainplt_mod, only : ainplt
  use ainpltpa_mod, only : ainpltpa
  use ainsetpa_mod, only : ainsetpa
  use ampfar_mod, only : ampfdiff
  use ampfar_mod, only : ampfefldb
  use cfpgamma_mod, only : cfpgamma
  use coefmidt_mod, only : coefmidt
  use coefmidv_mod, only : coefmidv
  use coefstup_mod, only : coefstup
  use coefwti_mod, only : coefwti
  use coefwtj_mod, only : coefwtj
  use diag_mod, only : diag
  use diagentr_mod, only : diagentr_vol
  use diaggnde_mod, only : diaggnde
  use diagimpd_mod, only : diagimpd
  use diagscal_mod, only : diagscal
  use dsk_gr_mod, only : dsk_gr
  use dskout_mod, only : dskout
  use eflditer_mod, only : eflditer
  use esefld_mod, only : efld_cd
  use ntdstore_mod, only : ntdstore
  use pltendn_mod, only : plt_fow_cons
  use pltinit_mod, only : pltinit
  use pltmain_mod, only : pltmain
  use pltrun_mod, only : pltrun
  use profiles_mod, only : profiles
  use r8subs_mod, only : dcopy
  use restvty_mod, only : restvty
  use sigv_mod, only : sigv
  use tddiag_mod, only : tddiag
  use tdinitl_mod, only : tdinitl
  use tdnflxs_mod, only : tdnflxs
  use tdnpadiag_mod, only : tdnpadiag
  use tdoutput_mod, only : tdoutput
  use tdplteq_mod, only : tdplteq
  use tdpltmne_mod, only : tdpltmne
  use tdsxray_mod, only : tdsxray
  use tdtloop_mod, only : tdtloop
  use tdtoaray_mod, only : tdtoaray
  use tdtransp_mod, only : tdtransp
  use tdtrcon_mod, only : tdtrcon
  use tdtrdfus_mod, only : diffus_io
  use tdtrdfus_mod, only : difus_io_scale
  use tdtrfcop_mod, only : tdtrfcop
  use tdtrrsou_mod, only : tdtrrsou
  use tdtrsavf_mod, only : tdtrsavf
  use urfchief_mod, only : urfchief
  use urfwrite_mod, only : urfwrite
  use wparsou_mod, only : wparsou
  use wpavg_mod, only : wpavg
  use wpelecf_mod, only : wpelecf
  use wpsavf_mod, only : wpsavf
  use wptrafx_mod, only : wptrafx
  use wptramu_mod, only : wptramu

  !---END USE

!
!

contains

      subroutine tdchief
      use param_mod
      use comm_mod
      use netcdfrf_mod, only : netcdfrf
      use pltmain_mod, only : pltmain
      use r8subs_mod, only : dcopy
      use ampfar_mod, only : ampfarl, ampfdiff, ampfefldb
      use ampfar_mod, only : ampfinit, ampfsoln
      use ainsetpa_mod, only : ainsetpa
      use ainpltpa_mod, only : ainpltpa
      use ainplt_mod, only : ainplt
      use aindfpa_mod, only : aindfpa, ainadjnl_fsetup_setup0, ainadjnl
      use aclear_mod, only : aclear
      use achiefn_mod, only : achiefn
      use achief1_mod, only : achief1

      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine directs the calculation of CQL3d; it controls
!     input, output, calls to CQL, calls to WD*, and holds the
!     main loops over radius of the toroidal device.
!..................................................................

      include 'name.h'
!MPIINSERT_INCLUDE

      character*8 icall,iplotsxr
      data iflag1/0/

!.......................................................................
!     Open cqlinput NL file and adjust to new setup0/setup structure,
!     if old two setup namelist sections are present [maintaining
!     backwards compatibility].  BH070414.
!.......................................................................
      call ainadjnl(0) ! at mpirank=0 only (but it can change cqlinput)
!MPIINSERT_BARRIER

      !Look for the presence of &FSETUP namelist,
      !adjust cqlinput: rename &FSETUP to &SETUP0 (later, restore cqlinput)
      call ainadjnl_fsetup_setup0(0) !at mpirank=0 (but it can change cqlinput)
      !pause
!MPIINSERT_BARRIER
!.......................................................................
!     Set default values for setup0 namelist
!.......................................................................
      call aindfpa

!..................................................................
!     Read in first namelist, setup0: determines type of run, etc.
!..................................................................
      open(unit=2,file='cqlinput',status='old')
      read(2,setup0)
!MPIINSERT_BARRIER

!..................................................................
!     Set some major input parameters according to
!     first namelist setup
!..................................................................
      call ainsetpa

!.......................................................................
!     Zero/set some arrays
!.......................................................................
      call aclear

      read(2,setup)  ! Gets pltinput variable, for ainplt routine.
      close(2)

      sumdtr=zero

!.................................................................
!     If noplots.ne."enabled1", initialize PGPLOT and plot
!     namelist input and parameters which are set in code before
!     compilation.
!.................................................................
      if (noplots.ne."enabled1") then
        call pltinit   ! Initiates PGPLOT
        call ainplt
        open(unit=2,file='cqlinput',status='old')
        read(2,setup0)         !  re-read 1st (i.e., setup0) namelist
        close(2)
        call ainsetpa          !  re-set according to setup0 nml
        call ainpltpa ! plots out the parameters
      endif


!..................................................................
!     If namelist variable lrzmax.le.1, then run the 2D code and
!     dispense with the "td" module. (Runs thru achief1 and achiefn
!     to a normal exit in ntloop).
!..................................................................
!-YuP: Not used anymore. Use general call to achiefn below, in ll-loop
!-BH:  Except that some lrzmax=1 test cases are still used, and
!-BH:  errors can arise when only only one value of a mesh
!-BH:  is called  for.  For example, tdrmshst fails.
!-BH:  However, additional functionality, such as sigmamod plotting
!-BH:  depdends on lrzmax.ge.4 (BH120308: check this).
      if (lrzmax.le.1) then
         nefiter=1              ! counts iterations; elecfld iterations for
                                ! one flux surface are not functional;
                            ! set nefiter to 1 for logic control in impavnc0
         call achief1           ! YuP: only called during n=0; Why needed?
                                ! BH:  In the past, at least, this call with
                                !      lrzmax=1, time-stepped the soln
                                !      to n=nstop in achiefn.
      endif

!..................................................................
!     call routine which controls array initialization
!     (and sets n=0, n_(1:lrorsa)=0 though call aindflt1).
!..................................................................

      call tdinitl !-> call ainitial
         ! tdinitl-> eqcoord-> eqfndpsi-> eqorbit-> trace flux surf.
         ! solr(l,lr_), solz(l,lr_) are R,Z coords. of flux surface
      dtr0=dtr

!..................................................................
!     Initialize main netCDF write, if netcdfnm.ne."disabled".
!     If netcdfshort.eq.'lngshrtf' determine number of distn saves.
!..................................................................

!MPIINSERT_IF_RANK_EQ_0
      if (netcdfnm.ne."disabled") then
         call netcdfrw2(0)
      endif
      if (netcdfshort.eq.'lngshrtf') then
         isave=0  !Index for nonzero values of increasing nsave(1:nsavea)
         nsavet=0 !Total nsave(i).ge.0 and .le.nstop
         do i=1,nsavea
            if (nsave(i).ge.0 .and. nsave(i).le.nstop) &
                 nsavet=nsavet+1
         enddo
      endif

!MPIINSERT_ENDIF_RANK

!.......................................................................
!     LOOP OVER TIME-STEP n
!.......................................................................
 10   continue  !Loop back point, from end of subroutine

!MPIINSERT_BARRIER
!MPIINSERT_STARTTIMESTEP
      call cpu_time(t_n_start) !-.-.-.-.-.-.-.-.-.-.-.-.-.

!.......................................................................
!     Reset dtr time step if (n+1).eq.nondtr1(i)
!     (achiefn, below, advances n to n+1 near its beginning).
!.......................................................................

      do 15 i=1,ndtr1a
         if ((n+1).eq.nondtr1(i)) then
            dtr=dtr1(i)
            dtreff=dtr
            dttr=dtr*nrstrt
         endif
 15   continue

!.......................................................................
!     Recalculate neutral beam source, if (n+1).eq.nonvphi or noffvphi,
!     or,
!     if at beginning of a beam pulse, and have time-dep background.
!     Set ibeampon= beam pulse on-indicator (used in coefstup)
!         ibeamponp= indicates freya calc carried out for given pulse.
!.......................................................................

      if ((n+1).eq.nonvphi .or. (n+1).eq.noffvphi) then
        write(*,*)'tdchief/nonvphi: call frnfreya, n=',n,'time=',timet
        call frnfreya(frmodp,fr_gyrop,beamplsep,beamponp,beampoffp, &
                      hibrzp,mfm1p,noplots)
      endif

      k=kfrsou  !Only set up for one modulated beam species.
      if(k.gt.0)then
      do is=1,nso
        if ((n+1).ge.nonso(k,is) .and. (n+1).lt.noffso(k,is)) then
!          Check and apply criteria for square wave pulsed beam:
           if (beamplsep.ne.'disabled' .and. k.eq.kfrsou) then
!             Set flag and time indicating time beam first starts
              if (iflag1.eq.0) then
                 iflag1=1
                 timestrt=timet+dtreff !Need to check at time of next
                                       !velocity step.
                 timespn=beamponp+beampoffp
                 timeon=beamponp
              endif
!           Determine if in an on-period of the beam pulse cycle
              iperiod=(timet+dtreff-timestrt)/timespn + em12 !allowing for
                                                    !roundoff near zero.
              time1=iperiod*timespn !Time at beginning of a pulse
              timedf=timet+dtreff-time1
              if (timedf.le.timeon) then
                 ibeampon=1
                 if (ibeamponp.eq.0.and.nbctime.ne.0) then
                 write(*,*)'tdchief/nbctime:frnfreya,n=',n,'time=',timet
                   call frnfreya(frmodp,fr_gyrop,beamplsep,beamponp, &
                         beampoffp,hibrzp,mfm1p,noplots)
                   ibeamponp=1
                 endif
              else
                 ibeampon=0
                 ibeamponp=0
              endif
           endif                !On beamplsep and kfrsou
        endif                   !On (n+1)
      enddo                     !On is
      endif ! k>0

!$$$c    Alternative to above do loop, test coding:
!$$$      ibeampon=1
!$$$      ibeamponp=1


!.......................................................................
!     Obtain transport time-dependent scale factor, for required cases
!.......................................................................
          if (difus_io(1).eq."drrin") then
             if (n_d_rr.ne.0) then
                !Getting time-dep scale factors
                do k=1,n_d_rr   !n_d_rr obtained in diffus_io
                   drrt(k)=difus_io_scale(k,1)
                enddo
             endif
          elseif(difus_io(1).eq."drrdrin") then
             if (n_d_rr.ne.0) then
                !Getting time-dep scale factors
                do k=1,n_d_rr   !n_d_rr obtained in diffus_io
                   drrt(k)=difus_io_scale(k,1)
                   drt(k)=difus_io_scale(k,2)
                enddo
             endif
          endif

!.......................................................................
!     compute source term due to spatial transport operator (ADI method)
!.......................................................................

      if (transp.eq."enabled" .and. adimeth.eq."enabled" .and. &
        (n+1).ge.nonadi) then
        dtreff=dtr / 2.0
        dttr=dtreff
!     radial transport
        if (cqlpmod .ne. "enabled") call tdtrrsou
!     parallel transport
        if (cqlpmod.eq."enabled" .and. n.ge.nontran) call wparsou
      endif

!..................................................................
!     Call urf module..
!..................................................................
      t_urf1=0.
      t_urf2=0.
      t_urf3=0.
      if (urfmod.ne."disabled") then
        call cpu_time(t_urf1) !-.-.-.-.-.-.-.-.-.-.-.-.-.
        call urfchief
        call cpu_time(t_urf2) !-.-.-.-.-.-.-.-.-.-.-.-.-.
!MPIINSERT_IF_RANK_EQ_0
        WRITE(*,*) 'tdchief after urfchief'
!MPIINSERT_ENDIF_RANK
!       Set up netcdf store of rf data.
!MPIINSERT_IF_RANK_EQ_0
        if (netcdfnm.ne."disabled" .and. n.eq.0) then
           do krf=1,mrf !YuP:04-2010: Separate data file for each wave type krf
              call netcdfrf(0,krf) !kopt=0: initialize and write grids,...
              WRITE(*,*) 'after netcdfrf(0,krf)  krf=', krf
           enddo
        endif
!MPIINSERT_ENDIF_RANK
        call cpu_time(t_urf3) !-.-.-.-.-.-.-.-.-.-.-.-.-.
        !YuP[03-2016] Repeat plotting surfaces, but now - with rays
        if (noplots.ne."enabled1" .and. eqmod.eq."enabled" &
            .and. n.eq.0) then
           if(mrfn.gt.0)then ! just in case
              !(mrfn is supposed to be >0 when urfmod.ne."disabled")
              do krf=1,mrfn
               call tdplteq(krf) ! separate page with rays for each krf
              enddo
           endif
        endif
      endif ! urfmod.ne."disabled"

!.......................................................................
!     cqlpmod.eq.enabled:
!     Compute parallel electric field from Poisson equation
!.......................................................................

      if (transp.eq."enabled" .and. cqlpmod.eq."enabled") &
        call wpelecf(11)  ! eleven, not ll

!.......................................................................
!     cqlpmod.eq.enabled:
!     Compute electrostatic parallel electric field from constant flux
!     condition (for CQLP + current drive localized along field line).
!     Uses Kupfer (PoP, 1995) f and g function decompostion in achiefn.
!     This is a nonlocal Efield calc, and therefore all space
!     points must be cycled through.
!     BH131103 Comment added:  This option is not operational.
!              Compare with the stella code, where it is implemented.
!              Also compare with ampfmod coding.
!.......................................................................

      if (cqlpmod.eq."enabled" .and. eseswtch.eq."enabled" .and. &
          sbdry.eq."enabled") then
         kopt=2
         ilend=lrors
         do ll=1,ilend
            call tdnflxs(ll)
            call achiefn(kopt) ! achiefn(2) here (cqlpmod="enabled")
         enddo
         ! XXX this looks like an indexing bug to me
         ! xxx but i used what is defined by  esefld::efld_cd code
         print *, 'XXXCHECK'
         call efld_cd(dz(1:lrors+2,lr_),lrors,vnorm,flux1,flux2,elparnw,flux0)
      endif

!.......................................................................
!     MAIN VEL TIME STEP: loop over each radial position sequentially
!.......................................................................

      ilend=lrors
      if (transp.eq."enabled" .and. cqlpmod.eq."enabled" .and. &
        mod(nummods,10).ge.5 .and. sbdry.ne."periodic" &
        .and. lmidvel.ne.0) ilend=lrors-1

      if (nstop.eq.0) go to 2

!BH131230:  Why won't ampfmod work before MAIN VEL TIME STEP?? nstop=0 issue?
!.......................................................................
!     Computer the time-advanced toroidal electric field from Ampere-
!     Faraday equations using Kupfer h and g functions.
!
!     cqlpmod.ne.enabled has been checked.
!     Initialize electric fields at turnon of ampfmod.
!.......................................................................


      if (ampfmod.eq.'enabled' .and. n.eq.0) then
         !YuP[21-08-2017] added: no time adv. yet, no ampf iterations yet,
         !simply fill-in and save the values at n=0 for plotting.
         !!! Besides, if nonampf=1, we need the zero iteration,
         !!! which is taken from the previous time step,
         !!! and so we need values of elecfldn at n=0 in such case.
         !!! See a note just after call ampfefldb few lines below.
         do niter=0,nampfmax ! or up to nefitera
           do ll=1,lrz
           elecfldn(ll,n,niter)=elecfld(ll)/300.d0 !here: n=0
           enddo
           ! Also need the end (bndry) point:
           elecfldn(lrz+1,n,niter)=elecfldb/300.d0 !here: n=0
           ! Also save ll=0 point (magn.axis point):
           elecfldn(0,n,niter)=elecfld(0)/300.d0   !here: n=0
         enddo
      endif

      if (ampfmod.eq."enabled" .and. n+1.ge.nonampf) then
         ! n is not updated yet. Here n=0,1,2,...,nstop-1.
         ! For example, if nstop=5 and nonampf=5
         ! (meaning: apply ampf calc. at the last step),
         ! we need to start using this part when n+1=nstop=nonampf

         if (n+1.eq.nonampf) then
           call ampfinit ! n is not updated yet, Here n=0,1,...,nstop-1
         endif

         if (kelecg.eq.0) then
            write(*,*)
            WRITE(*,*)'Amp-Faraday eqns require electron gen species'
            write(*,*)
            stop
         endif
!        Find dtr for this time-step
         do i=1,ndtr1a
            if ((n+1).eq.nondtr1(i)) then
               dtr=dtr1(i)
               dtreff=dtr
               dttr=dtr*nrstrt
            endif
         enddo
!        Copy current distribution f into f_
         call dcopy(iyjx2*ngen*lrors,f(0:iyjx2*ngen*lrors-1,0,kelec,1), &
              1,f_(0:iyjx2*ngen*lrors-1,0,kelec,1),1)
!        Bring background profiles up to time step n
         ! No effect if bctime=0 (time-indep. profiles).
         call profiles ! if(ampfmod.eq.'enabled' .and. n+1.ge.nonampf)
                       ! skip elecfld background profiles

!        Get given boundary condition (i.e., edge) elec field, at time
!        advanced position (i.e., elecfldn(ll,n+1,0) at ll=lrz+1 point).
!        If time-dependent (nbctime>0),
!        then it is given by time advanced elecb() or elecin_t(njene,).
         call ampfefldb(n+1,time+dtr)
!        Besides, subr.ampfefldb Sets zero iteration of elecfldn
!        equal to previous time step radial profile, like this:
         ! elecfldn(ll,nn,0)=elecfldn(ll,nn-1,it_prev) ! where nn==n+1 (and ll=0:lrz)
         ! So, if nonampf=1, then we have here n=0, or nn=1,
         ! which means that we need the values of elecfldn at nn-1=0.


!        elecfld(1:lrz) has been set.  If n=1=nonampf, then is given;
!                       If n.ge.nonampf, then set near end of ampfsoln.

!        Determine the time-advanced tor e-field
         kopt=3  !Amp-Faraday option
         ilend=lrz

!        Initialize elecfldn(ll,n+1,it=0) for upcoming ampfsoln
         do ll=1,ilend
            elecfldn(ll,n+1,0)=elecfld(ll)/300.d0
         enddo
         !YuP[21-08-2017] added but then commented out:
         !elecfldn(lrz+1,n+1,0)=elecfldb/300.d0 !Not needed, done in ampfefldb

!        Obtain the h and g functions over all radii and iteratively
!        solve Ampere-Faraday eqns for the tor electric field.
         do it=1,nampfmax
            it_ampf=it

!BH131109: Maybe don't need from bndry to center for fh/fg?   do ll=ilend,1,-1
!BH131109: Backwards is giving bounds problem with abd in impavnc0.
!BH131109: Present system is determining all fh/fg first as function
!BH131109: of radius, then solve AF eqn.  Future solves may require
!BH131109: fix the counter resetting for ll=ilend,1,-1.

            do ll=1,ilend
               call tdnflxs(ll)  !get l_,lr_,..
               call achiefn(kopt) !kopt=3 here !Increments n by 1, for each ll
                                   !Gives fh,fg(,,1,ll) for each ll.
                                   !it_ampf is passed in common block
               ! YuP test/printout
               amp_f_=ampfarl(f_(0:ll,0,kelec,ll),ll)*dtr !ampfarl has 1/dtr factor
               amp_f=ampfarl(f(0:ll,0,kelec,ll),ll)*dtr !ampfarl has 1/dtr factor
               amp_h=ampfarl(fh(0:ll,0,kelec,ll),ll)*dtr
               amp_g=ampfarl(fg(0:ll,0,kelec,ll),ll)*dtr
               write(*,'(a,3i4,3e12.4)') &
                 'after achiefn(3): n,it,ll,integrals f, f-h, g:', &
                  n,it,ll, amp_f, amp_f-amp_h, amp_g
               ! YuP test/printout:
               !Confirmed that f() remains unchanged during
               !iterations in this it loop, and it remains equal to f_()
            enddo
!           Solve for iteratively updated tor electric field
!           n is time advanced value
            call ampfsoln(it,n) !here n=nonampf,...,nstop
!           Presently, ampfdiff does nothing, i.e. gives iflag=1
!           In future, refine this, making it=it(1:lrz).
            itt=it
            call ampfdiff(iflag)  !iflag=0, only if error criteria met.
                                  !Presently a dummy call giving iflag=0
            if (iflag.eq.0) go to 16
         enddo !On it (Note: it is incremented by 1 at the exit of loop)
 16      continue

         !write(*,'(a,3i4)')'tdchief-373: n,itt,it=', n,itt,it

!        Now, setup with time-advanced electric field in elecfld(),
!        then at end of ampfmod interations, call achiefn(0) below.
      endif  !On ampfmod


!...........................................................
!     Loop back point for electric field iteration
!       (if efiter.eq."enabled")
!...........................................................
      nefiter=1 ! counts iterations
      do ll=1,ilend
        call tdnflxs(ll) ! determine l_,lr_, etc.
        nefiter_(l_)=1 ! counts iterations for each flux surface
      enddo

 20   continue  !Loop back point for electric field iteration

      call cpu_time(t_before_soln) !-.-.-.-.-.-.-.-.-.-.-.-.-.
!MPIINSERT_BARRIER

!     Copy current distribution f into f_
      call dcopy(iyjx2*ngen*lrors,f(0:iyjx2*ngen*lrors-1,0,1,1),1, &
           f_(0:iyjx2*ngen*lrors-1,0,1,1),1)

      if (transp.eq."enabled" .and. n.ne.0 .and. adimeth.ne."enabled" &
            .and. soln_method.ne."it3drv" .and. nefiter.eq.1)  then
!.......................................................................
!     With adimeth, one assumes f (in tdtrrsou or wparsou)
!     as the starting point of the new iteration.
!     Thus should not need to redefine here.
!     CQLP case: f defined in wpsavf is assumed to be the good one.
!     soln_method=it3drv:  soln is in f(,,,) on last call on the set of
!     flux surfaces (soln at start of step is in f_(,,,)).
!.......................................................................
        if(cqlpmod.ne."enabled" .and. n.ge.nontran .and. n.lt.nofftran) &
            call dcopy(iyjx2*ngen*lrors,frn_2(0:iyjx2*ngen*lrors-1, &
              0,1,1),1,f(0:iyjx2*ngen*lrors-1,0,1,1),1)
      endif

!..................................................
!     Bring background profiles up to time step n.
!     Exclude updating of electric field elecfld(1:lrz) for
!     ampfmod.eq.enabled .and. n.ge.nonampf [Check this for
!     the eseswtch case also.]
!..................................................

      if(nefiter.eq.1) call profiles

!.......................................................................
!BH131103      do 1 ll=1,ilend   !ilend=lrz for cqlpmod.ne.'enabled'
!
!     Can reverse order of do 1 ll=1,ilend, if needed for purposes
!     of solving Ampere-Faraday eqns starting at given edge voltage.
!     Results are unchanged for simple efield case, BUT additional
!     debugging required for TCV ECCD/radial diffusion test case (CD
!     profiles are similar, but total driven current changes from
!     108 kA to 138 kA.)
!BH131107 Moved the ll=1,ilend reversal to above if(ampfmod.eq.enabled..
!.......................................................................
      do 1 ll=1,ilend   !ilend=lrz for cqlpmod.ne.'enabled'
        !determine local variables depending on flux surface (l_,iy,..)
        call tdnflxs(ll) !-> get l_,lr_,...
        ! Reset time step if (n+1).eq.nondtr1(i). .AND. LRZMAX=1
        do i=1,ndtr1a
           if ((n+1).eq.nondtr1(i)) then
              dtr=dtr1(i)
              dtreff=dtr
              dttr=dtr*nrstrt
           endif
        enddo

!......................................................................
!     If ampfmod.eq.enabled, but n.lt.nonampf, then insert values of
!     elecfldn(1,,) from elecfld [for plotting with .nc].
!......................................................................

      if (ampfmod.eq.'enabled' .and. n+1.lt.nonampf) then
         ! Here, n is NOT advanced yet for the next step.
         ! Here n=0,1,2,...,nstop-1.
         !    n+1=1,2,3,...,nstop
         ! For example, if nstop=5 and nonampf=5
         ! (meaning: apply ampf calc. at the last step),
         ! we need to skip this part
         ! when n+1=nonampf,
         ! but use it when n+1<nonampf
         ! Save into the future-updated n+1= 1,2,3,nonampf-1
!BH170329
         !YuP[21-08-2017] Corrected, for the fields storage BEFORE the ampf is applied:
         do niter=0,nampfmax ! or up to nefitera
           elecfldn(ll,n+1,niter)=elecfld(ll)/300.d0
           ! here: ll=1:lrz. But we also need the end (bndry) point:
           elecfldn(lrz+1,n+1,niter)=elecfldb/300.d0 !here: n+1.lt.nonampf
           ! Can also save ll=0 point (magn.axis point):
           elecfldn(0,n+1,niter)=elecfld(0)/300.d0
         enddo
      endif


!......................................................................
!     TIME ADVANCE ON FLUX SURFACE rz(ll).
!     If cqlpmod.ne.'enabled', transp='enabled', soln_method='it3drv',
!       then following call only adds ll flux surface contributions to
!       the coefficient matrix, and then solves for time advanced distn
!       at last flux surface, ll=lrz.
!......................................................................

!MPIINSERT_MPIWORKER
!    It will insert :
!      if(soln_method.eq.'direct' .and. lrzmax.gt.1) then
!         ! Parallelization for the impavnc0 solver is limited
!         ! to soln_method='direct' (for now)
!         mpiworker= MOD(ll-1,mpisize-1)+1  !1...(mpisize-1)
!      else
!         ! In all other cases, perform calculations
!         ! for all flux surfaces on mpirank=0, then broadcast results
!         mpiworker=0
!      endif

        call achiefn(0)  !--> if implct='enabled', calls impavnc0


!MPIINSERT_SEND_RECV
!     It will send or recv data, but only in case of
!     soln_method='direct' (for now)

 1    continue ! End loop over radius:  New f is obtained for each ll

!MPIINSERT_BARRIER
!MPIINSERT_BCAST_DISTRIBUTION
!MPIINSERT_BCAST_COLL_COEFFS
!MPIINSERT_BCAST_SCAL
!MPIINSERT_BCAST_VELSOU
      call cpu_time(t_after_soln) !-.-.-.-.-.-.-.-.-.-.-.-.-.

      do ll=1,ilend  ! Re-scale f if needed;
        call tdnflxs(ll) ! determine l_,lr_, etc.
        do k=1,ngen  ! Compute density gains and losses, and powers.
           ! For lbdry0='disabled',  Redefine f at v=0 so it is unique:
           ! (For lbdry0='enabled', coeff matrix is set up
           !   to automatically maintain unicity.)
           if (lbdry0.ne."enabled") then !-YuP: moved here from impavnc0
             call dcopy(iyjx2,f(0:iyjx2-1,0,k,l_),1, &
                   fxsp(0:iyjx2-1,0,k,l_),1)
             s=0.
             t=0.
             do 2100 i=1,iy
               s=s+vptb(i,lr_)*cynt2(i,l_)
               t=t+vptb(i,lr_)*cynt2(i,l_)*f(i,1,k,l_)
 2100        continue
             do 2200 i=1,iy
               f(i,1,k,l_)=t/s
 2200        continue
           endif
           !YuP-101227: Skip diagnostics here; It is done after 21_continue.
           !(to skip, set lbdry(1)='fixed'):
           call diagscal(k) !YuP! renorm f if lbdry(k)=scale/consscal
        enddo ! k
      enddo ! ll

      call cpu_time(t_after_diag1) !-.-.-.-.-.-.-.-.-.-.-.-.-.

!..................................................................
!     Check/iterate on electric field to obtain a target current
!..................................................................
      if (efiter.eq."enabled") then
         ! Stop iterations, if beyond iteration limit:
         if (nefiter .gt. nefitera) go to 21
         ! Stop iterations if n<(control turn on step):
         if (n .le. noncntrl) go to 21
         write(*,*)'TDCHIEF/EFLDITER: Time step=',n, &
                   '    Starting iteration #',nefiter
         nefiter_all=0 ! Summ-up nefiter_out (output of eflditer)
         do ll=1,ilend ! Check current for each flux surface
            ! determine local variables depending on flux surface (l_, iy,..)
            call tdnflxs(ll)
            if(nefiter_(l_).ge.1) then
              nefiter_in=nefiter_(l_) ! remember the input value of nefiter.
              ! re-adjust elecfld(l_) if target current is not achieved:
              call eflditer ! nefiter_(l_) is both in and out; now in comm.h !
              ! If target current is achieved, nefiter_(l_)->0 for this flux surface.
              nefiter_out=nefiter_(l_)
              nefiter_all=nefiter_all+nefiter_out
              ! nefiter_all will remain 0 if nefiter_out=0 for each flux surface.
              !!!nefiter_(l_)=nefiter_(l_)+1 ! counts iterations for each flux surf.
              ! Comment the line above to skip iterations for surfaces for which
              ! target current is achieved (otherwise - check them during next iteration)
            endif
         enddo ! ll
         if (nefiter_all.eq.0) go to 21 ! Finished with iterations
         ! Otherwise: Restore distribution function back to value at beginning of
         ! of the time step  (saved near beginning of impanvc/impavnc0).
         do ll=1,ilend
            ! determine local variables depending on flux surface (l_, iy,..)
            call tdnflxs(ll)
            call dcopy(iyjx2,f_(0:iyjx2-1,0,1,l_),1, &
                 f(0:iyjx2-1,0,1,l_),1)
         enddo ! ll
         nefiter=nefiter+1 ! counts iterations
         go to 20 ! Another iteration, using old f() but new elecfld()
      endif

 21   continue

!MPIINSERT_BARRIER

!     Accumulate favg time average, if tavg.ne."disabled"
      if (tavg.ne."disabled") then
         do ii=1,ntavga
         tavg12=tavg2(ii)-tavg1(ii)
         if (tavg12 .gt. 0.d0) then
         if (timet.gt.tavg1(ii) .and. timet.le.tavg2(ii)) then
           !Note: if at the end of run
           ! timet is still less than the very 1st time point
           ! for averaging t=tavg1(1)
           ! (may happen if nstop and/or dtr are too small)
           ! then the averaging process never happens,
           ! and so sumdtr remains zero.
           sumdtr=sumdtr+dtr   !Accumulate dtr, for denom in time-avg
           do ll=1,ilend
           do k=1,ngen
               do i=0,iy
                  do j=0,jx
                    favg(i,j,k,ll)=favg(i,j,k,ll)+dtr*f(i,j,k,ll)
                  enddo
               enddo
           enddo  ! On k
           enddo  ! On ll
         endif  ! On timet
         endif  ! On tavg12>0
         enddo  ! On ii
         ! At the end of run form the distribution averaged over all
         ! "ii-blips" [tavg1(ii);tavg2(ii)]
         if (n .eq. nstop) then
           if(sumdtr.lt.em90)then
             ! YuP[03-01-2016] fix added for sumdtr=0 case
!MPIINSERT_IF_RANK_EQ_0
           WRITE(*,*)'tdchief WARNING: No time averaging of f was done'
           WRITE(*,*)'tdchief WARNING: Probably t=nstop*dtr is less'
           WRITE(*,*)'tdchief WARNING: than tavg1(1).'
           WRITE(*,*)'tdchief WARNING: Instead of favg(), '
           WRITE(*,*)'tdchief WARNING: the last available f() will be'
           WRITE(*,*)'tdchief WARNING: saved into mnemonic.nc file. '
!MPIINSERT_ENDIF_RANK
             do ll=1,ilend
             do k=1,ngen
               do i=0,iy
               do j=0,jx
                  favg(i,j,k,ll)=f(i,j,k,ll)
               enddo
               enddo
             enddo  ! On k
             enddo  ! On ll
           else  ! sumdtr>0
             do ll=1,ilend
             do k=1,ngen
               do i=0,iy
               do j=0,jx
                  favg(i,j,k,ll)=favg(i,j,k,ll)/sumdtr
               enddo
               enddo
             enddo  ! On k
             enddo  ! On ll
           endif ! sumdtr =0 or >0
         endif  ! On n.eq.nstop
      endif  ! On tavg


      do ll=1,ilend  ! perform diagnostics.
        call tdnflxs(ll) ! determine l_,lr_, etc.
        call cfpgamma ! Re-calc. Coul.Log for the new distr.func.
        do k=1,ngen  ! Compute density gains and losses, and powers.
           call coefstup(k) ! To define da...df coeffs, gon(i,j), etc
           call coefmidv(da,1)
           call coefmidv(db,2)
           call coefmidv(dc,3)
           call coefmidt(dd,1)
           call coefmidt(de,2)
           call coefmidt(df,3)
           call coefwtj(k)
           call coefwti(k)
           call diagimpd(k) !->diagentr , get sgain(1:8,k)
        enddo ! k
        call achiefn(1) !-> call diaggnde: <energy>, density, energy transfer
        !  Write out distributions and meshes (controlled by namelist):
        if (n .ge. nstop) call dskout(ll)
        if (n .ge. nstop) call dsk_gr
      enddo ! ll

!MPIINSERT_BARRIER

!......................................................................
!     Integrate power densities over plasma volume
!......................................................................

        call diagentr_vol

      call cpu_time(t_after_diag2) !-.-.-.-.-.-.-.-.-.-.-.-.-.

 2    continue ! if nstop.eq.0

      if (ilend .eq. lrors-1) n_(lrors)=n_(lrors-1)

!.......................................................................
!     Time advance spatial transport equation [except if transport
!     with soln_method=it3drv, not applicable since done in impavnc0].
!.......................................................................

      if (soln_method.ne.'it3drv') then
      if (transp .eq. "enabled" .and. n.ge.nontran .and. n.le.nofftran) &
         then

        if (cqlpmod .ne. "enabled") then
!......................................................................
!     The radial transport (tdtr...) routines follow. Store the results
!     of the velocity space split in fvn. Then call the transport
!     subroutine, tdtransp.
!......................................................................
          call tdtransp
          call tdtrsavf

        else  ! i.e., cqlpmod.eq.enabled

!.......................................................................
!     Parallel transport.
!     First compute E_n+1/2 from Poisson and modify velsou accordingly
!.......................................................................
          call wpelecf(2)
!%OS  call wptrans
!%OS  call wpsavf
!%OS
          dtreff=dtreff/noffso(1,1)
          do 200 in=1,noffso(1,1)
            if (meshy .eq. "fixed_mu") call wptramu
            if (meshy .ne. "fixed_mu") call wptrafx
            call wpsavf
 200      continue
          dtreff=noffso(1,1)*dtreff
!%OS
        endif ! on cqlpmod

      endif  ! on transp, etc.
      endif  ! on soln_method.ne.it3drv

!.......................................................................
!     Plot flag and plot time set:
!.......................................................................

      iplt3d=0
      if (noplots.ne."enabled1") then
         do i=1,nplota
            if(n.eq.nplt3d(i)) then
               iplt3d=i
               tplt3d(i)=timet
            endif
         enddo
      endif

!.......................................................................
!     Diagnostics (if transp=enabled)
!.......................................................................

      do 30 ll=1,lrors
!......................................................................
!     determine local variables depending on flux surface (l_, iy,..)
!......................................................................
        call tdnflxs(ll)

        if (transp.eq."enabled") then
          call diaggnde

!..................................................
!     Compute plasma resistivity for electron runs.
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

!..................................................
!     plotting logic for 3-D code...
!       (iplot set in subroutine achiefn)
!..................................................
          if ((n.ge.nstop .or. iplot.ne.0 .or. iplt3d.ne.0).and. &
              n.ge.nontran .and. n.le.nofftran) then
            call tdtrfcop(1)
            call diaggnde
            if (l_ .eq. lmdpln_) then
              write(*,*)'tdchief before xlndnr: n,l_',n,l_
              do 25 k=1,ngen
                xlndnr(k,lr_)=xlndn(k,lr_)
                energyr(k,lr_)=energy(k,lr_)
                currr(k,lr_)=curr(k,lr_)/3.e9
 25           continue
            endif
            if ((adimeth.eq."enabled" .or. cqlpmod.ne."enabled") .and. &
              n.ge.nontran .and. n.le.nofftran) then
              call tdtrfcop(2)
              call diaggnde
              if (l_ .eq. lmdpln_) then
              write(*,*)'tdchief before xlndnv: n,l_',n,l_
                do 26 k=1,ngen
                  xlndnv(k,lr_)=xlndn(k,lr_)
                  energyv(k,lr_)=energy(k,lr_)
                  currv_(k,lr_)=curr(k,lr_)/3.e9
 26             continue
              endif
            endif
            call tdtrfcop(3)  !Copies all ngen species
            call diaggnde     !Loops over all species
!$$$  Commenting following, in favor of dimensioning pwrrfs with l_,
!$$$    to save results calculated in impavnc0 ==> diagimpd ==> diagentr.
!$$$    BH090806
!$$$            do k=1,ngen
!$$$               call diagentr(3,k) !To get pwrrfs at current l_ surface
!$$$            enddo
            if (noplots.ne."enabled1" .and. &
                (iplot.ne.0 .or. n.ge.nstop) ) then
                !YuP[2017-11-27] Added n.ge.nstop, similar to achiefn
!$$$            if ( iplot.ne.0) then
              call pltmain
            endif
          endif
!.......................................................................
!     write some radially dependent netcdf output, similar to pltmain.
!.......................................................................

          call netcdfmain

        endif !End of transp="enabled" if-block

!......................................................................
!     Store updated information.
!......................................................................
        call tdtoaray
 30   continue !  ll=1,lrors

!.......................................................................
!l    Compute flux surface average current and resistivity for CQLP case
!.......................................................................

      if (transp.eq."enabled" .and. cqlpmod.eq."enabled") call wpavg

!......................................................................
!     Compute conservation constant for transp="enabled"
!......................................................................

      if (transp.eq."enabled" .and. cqlpmod.ne."enabled") then
        call tdtrcon
      endif

!.......................................................................
!     compute some radial diagnostics.
!.......................................................................

!%OS  still to verify and change if cqlpmod=enabled
      call tddiag

!.......................................................................
!     Calc fusion rates, if sigmamod=enabled
!.......................................................................

        if (sigmamod.eq."enabled" .and. msig.gt.0) then
          icall="notfrst"
          call sigv(icall)
        endif ! sigmamod

!.......................................................................
!     Calc xray spectra, if softxry.ne."disabled".and.kelecg.gt.0
!     Spectra is stored in the netCDF output file at each time step.
!     Also, set plot flag  [iplt3d is set above, if plotting this step].
!.......................................................................

      if (lrzmax.lt.3 .and. softxry.ne."disabled") then
         write(*,*)'*******************************************'
         write(*,*)'tdchief:  SXR not computed for lrzmax.lt.3'
         write(*,*)'*******************************************'
      else
         if (softxry.ne."disabled".and.kelecg.gt.0) then
            icall="notfrst"
            ilold=l_
            iplotsxr='no'
            if (iplt3d.ne.0 .or. n.eq.nstop .and. noplots.ne.'enabled1') &
                 iplotsxr='yes'
            call tdsxray(icall,iplotsxr) !Soft Xray diagnostic
            call tdnflxs(ilold)
         endif
      endif


!.......................................................................
!     NPA synthetic diagnostic
!.......................................................................

        if (npa_diag.ne."disabled".and.niong.ne.0) then
           if (n.eq.nstop .or. npa_diag.eq.'ncdf_all') then
              call tdnpadiag(icall)
           endif
        endif

!BH090513:  BOB, check why no call tdnflxs after tdnpadiag??

!.......................................................................
!     plot radial information
!.......................................................................

      if ((iplt3d.ne.0 .or. n.eq.nstop).and.lrzmax .ge. 3) then
        icall="notfrst"
        if (noplots.ne."enabled1") then
          call tdpltmne
!         if equilibrium changed during run:
!%OS      if (eqmod.eq."enabled".and.n.ge.nstop) call tdplteq
            ! Plot the change in Total N of ptcles as a func of rho, for different time.
            call plt_fow_cons ! Can be used for ZOW
        endif
      endif

      if (n.eq.nstop.and.pltra.eq."enabled") then
        call pltrun
      endif

!.......................................................................
!     write radial information
!.......................................................................

      if (iplt3d.ne.0 .or. n.eq.nstop) then
        call tdoutput(2)
      endif

!.......................................................................
!     write netCDF output file, if netcdfnm.ne."disabled"
!.......................................................................
!MPIINSERT_IF_RANK_EQ_0
      if (netcdfnm.ne."disabled" .and. nstop.ne.0) then
         if (netcdfshort.eq.'lngshrtf') then
            isave=0
            do i=1,nsavet
               if(n.eq.nsave(i)) then
                  isave=i
                  tsave(i)=timet
               endif
            enddo
         endif
         call netcdfrw2(1)
         if (urfmod.eq."enabled") then
            do krf=1,mrf !YuP:04-2010: Separate data file for each wave type krf
               call netcdfrf(1,krf)
            enddo
         endif
      endif
!MPIINSERT_ENDIF_RANK
!.......................................................................
!     re-write ray data at last step, if urfmod & urfwrray="enabled"
!.......................................................................

      if (urfmod.ne."disabled" .and. urfwrray.eq."enabled" &
        .and. n.eq.nstop .and. nstop.ne.0) call urfwrite

!.......................................................................
!     Write tavg12, sumdtr, for tavg case
!.......................................................................

      if (tavg.ne."disabled") then
         tavg12=0.d0
         do ii=1,ntavga
            tavg12=tavg12+(tavg2(ii)-tavg1(ii))
         enddo
         write(*,*)
         write(*,*)'tavg.ne.disabled: tavg12 and sumdtr should be close'
         write(*,*)'at end of run: tavg12,sumdtr = ',tavg12,sumdtr
      endif

!.......................................................................
!     check for call exit
!.......................................................................
      call cpu_time(t_end) !-.-.-.-.-.-.-.-.-.-.-.-.-.
!MPIINSERT_BARRIER
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,'(a,f10.3,a,f10.3,a,f10.3)') &
                      'TDCHIEF: tcpu_urfchief=',t_urf2-t_urf1, &
                      '    tcpu_netcdf=',   t_urf3-t_urf2, &
                      '    tcpu_impavnc=', t_after_soln-t_before_soln

         WRITE(*,'(a,2f10.3)') ' tcpu_diagnostics 1 and 2:', &
                                t_after_diag1-t_after_soln, &
                                t_after_diag2-t_after_diag1

         WRITE(*,'(a,i5,a,f10.3)') ' Finished Time Step ======>',  n, &
                                '     tcpu(sec.)==',  t_end-t_n_start
         WRITE(*,'(/)')
!MPIINSERT_ENDIF_RANK

!MPIINSERT_BARRIER
      call tdtloop

      if (n .lt. nstop) go to 10 !-> next time step

!.......................................................................
!     write out radial diffusion coeffs at end of program, if enabled.
!.......................................................................
      if (difus_io(1).eq."drrout") then
         call diffus_io(0)
         call diffus_io(1)
      elseif (difus_io(1).eq."drrdrout") then
         call diffus_io(0)
         call diffus_io(2)
      endif

      return
      end
end module tdchief_mod
