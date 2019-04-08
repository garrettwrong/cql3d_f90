c
c
      subroutine tdchief
      use netcdfrf_mod, only : netcdfrf
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     This routine directs the calculation of CQL3d; it controls
c     input, output, calls to CQL, calls to WD*, and holds the
c     main loops over radius of the toroidal device.
c..................................................................

      include 'param.h'
      include 'comm.h'
      include 'name.h'
CMPIINSERT_INCLUDE

      character*8 icall,iplotsxr
      data iflag1/0/

c.......................................................................
c     Open cqlinput NL file and adjust to new setup0/setup structure,
c     if old two setup namelist sections are present [maintaining
c     backwards compatibility].  BH070414.
c.......................................................................
      call ainadjnl(0) ! at mpirank=0 only (but it can change cqlinput)
CMPIINSERT_BARRIER

      !Look for the presence of &FSETUP namelist,
      !adjust cqlinput: rename &FSETUP to &SETUP0 (later, restore cqlinput)
      call ainadjnl_fsetup_setup0(0) !at mpirank=0 (but it can change cqlinput)
      !pause
CMPIINSERT_BARRIER
c.......................................................................
c     Set default values for setup0 namelist 
c.......................................................................
      call aindfpa
      
c..................................................................
c     Read in first namelist, setup0: determines type of run, etc.
c..................................................................
      open(unit=2,file='cqlinput',status='old')
      read(2,setup0)
CMPIINSERT_BARRIER

c..................................................................
c     Set some major input parameters according to 
c     first namelist setup
c..................................................................
      call ainsetpa

c.......................................................................
c     Zero/set some arrays
c.......................................................................
      call aclear

      read(2,setup)  ! Gets pltinput variable, for ainplt routine.
      close(2)

      sumdtr=zero

c.................................................................
c     If noplots.ne."enabled1", initialize PGPLOT and plot 
c     namelist input and parameters which are set in code before
c     compilation.
c.................................................................
      if (noplots.ne."enabled1") then
        call pltinit   ! Initiates PGPLOT
        call ainplt
        open(unit=2,file='cqlinput',status='old')
        read(2,setup0)         !  re-read 1st (i.e., setup0) namelist
        close(2)
        call ainsetpa          !  re-set according to setup0 nml
        call ainpltpa ! plots out the parameters
      endif


c..................................................................
c     If namelist variable lrzmax.le.1, then run the 2D code and
c     dispense with the "td" module. (Runs thru achief1 and achiefn
c     to a normal exit in ntloop).
c..................................................................
c-YuP: Not used anymore. Use general call to achiefn below, in ll-loop
c-BH:  Except that some lrzmax=1 test cases are still used, and 
c-BH:  errors can arise when only only one value of a mesh
c-BH:  is called  for.  For example, tdrmshst fails.
c-BH:  However, additional functionality, such as sigmamod plotting
c-BH:  depdends on lrzmax.ge.4 (BH120308: check this).
      if (lrzmax.le.1) then
         nefiter=1              ! counts iterations; elecfld iterations for
                                ! one flux surface are not functional; 
                            ! set nefiter to 1 for logic control in impavnc0
         call achief1           ! YuP: only called during n=0; Why needed?
                                ! BH:  In the past, at least, this call with
                                !      lrzmax=1, time-stepped the soln
                                !      to n=nstop in achiefn.
      endif

c..................................................................
c     call routine which controls array initialization
c     (and sets n=0, n_(1:lrorsa)=0 though call aindflt1).
c..................................................................

      call tdinitl !-> call ainitial             
         ! tdinitl-> eqcoord-> eqfndpsi-> eqorbit-> trace flux surf.
         ! solr(l,lr_), solz(l,lr_) are R,Z coords. of flux surface      
      dtr0=dtr

c..................................................................
c     Initialize main netCDF write, if netcdfnm.ne."disabled".
c     If netcdfshort.eq.'lngshrtf' determine number of distn saves.
c..................................................................

CMPIINSERT_IF_RANK_EQ_0
      if (netcdfnm.ne."disabled") then
         call netcdfrw2(0)
      endif
      if (netcdfshort.eq.'lngshrtf') then
         isave=0  !Index for nonzero values of increasing nsave(1:nsavea)
         nsavet=0 !Total nsave(i).ge.0 and .le.nstop
         do i=1,nsavea
            if (nsave(i).ge.0 .and. nsave(i).le.nstop)
     1           nsavet=nsavet+1
         enddo
      endif

CMPIINSERT_ENDIF_RANK

c.......................................................................
c     LOOP OVER TIME-STEP n
c.......................................................................
 10   continue  !Loop back point, from end of subroutine

CMPIINSERT_BARRIER
CMPIINSERT_STARTTIMESTEP
      call cpu_time(t_n_start) !-.-.-.-.-.-.-.-.-.-.-.-.-.

c.......................................................................
c     Reset dtr time step if (n+1).eq.nondtr1(i)
c     (achiefn, below, advances n to n+1 near its beginning). 
c.......................................................................

      do 15 i=1,ndtr1a
         if ((n+1).eq.nondtr1(i)) then
            dtr=dtr1(i)
            dtreff=dtr
            dttr=dtr*nrstrt
         endif
 15   continue

c.......................................................................
c     Recalculate neutral beam source, if (n+1).eq.nonvphi or noffvphi,
c     or, 
c     if at beginning of a beam pulse, and have time-dep background.
c     Set ibeampon= beam pulse on-indicator (used in coefstup)
c         ibeamponp= indicates freya calc carried out for given pulse.
c.......................................................................

      if ((n+1).eq.nonvphi .or. (n+1).eq.noffvphi) then
        write(*,*)'tdchief/nonvphi: call frnfreya, n=',n,'time=',timet
        call frnfreya(frmodp,fr_gyrop,beamplsep,beamponp,beampoffp,
     .                hibrzp,mfm1p,noplots)
      endif

      k=kfrsou  !Only set up for one modulated beam species.
      if(k.gt.0)then
      do is=1,nso
        if ((n+1).ge.nonso(k,is) .and. (n+1).lt.noffso(k,is)) then
c          Check and apply criteria for square wave pulsed beam:
           if (beamplsep.ne.'disabled' .and. k.eq.kfrsou) then
c             Set flag and time indicating time beam first starts
              if (iflag1.eq.0) then
                 iflag1=1
                 timestrt=timet+dtreff !Need to check at time of next
                                       !velocity step.
                 timespn=beamponp+beampoffp
                 timeon=beamponp
              endif
c           Determine if in an on-period of the beam pulse cycle
              iperiod=(timet+dtreff-timestrt)/timespn + em12 !allowing for
                                                    !roundoff near zero.
              time1=iperiod*timespn !Time at beginning of a pulse
              timedf=timet+dtreff-time1
              if (timedf.le.timeon) then
                 ibeampon=1
                 if (ibeamponp.eq.0.and.nbctime.ne.0) then
                 write(*,*)'tdchief/nbctime:frnfreya,n=',n,'time=',timet
                   call frnfreya(frmodp,fr_gyrop,beamplsep,beamponp,
     +                   beampoffp,hibrzp,mfm1p,noplots)
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

c$$$c    Alternative to above do loop, test coding:
c$$$      ibeampon=1
c$$$      ibeamponp=1


c.......................................................................
c     Obtain transport time-dependent scale factor, for required cases
c.......................................................................
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
      
c.......................................................................
c     compute source term due to spatial transport operator (ADI method)
c.......................................................................

      if (transp.eq."enabled" .and. adimeth.eq."enabled" .and. 
     +  (n+1).ge.nonadi) then
        dtreff=dtr / 2.0
        dttr=dtreff
c     radial transport
        if (cqlpmod .ne. "enabled") call tdtrrsou
c     parallel transport
        if (cqlpmod.eq."enabled" .and. n.ge.nontran) call wparsou
      endif

c..................................................................
c     Call urf module..
c..................................................................
      t_urf1=0.
      t_urf2=0.
      t_urf3=0.      
      if (urfmod.ne."disabled") then
        call cpu_time(t_urf1) !-.-.-.-.-.-.-.-.-.-.-.-.-.
        call urfchief
        call cpu_time(t_urf2) !-.-.-.-.-.-.-.-.-.-.-.-.-.
CMPIINSERT_IF_RANK_EQ_0
        WRITE(*,*) 'tdchief after urfchief'
CMPIINSERT_ENDIF_RANK
c       Set up netcdf store of rf data.
CMPIINSERT_IF_RANK_EQ_0
        if (netcdfnm.ne."disabled" .and. n.eq.0) then
           do krf=1,mrf !YuP:04-2010: Separate data file for each wave type krf
              call netcdfrf(0,krf) !kopt=0: initialize and write grids,...
              WRITE(*,*) 'after netcdfrf(0,krf)  krf=', krf
           enddo
        endif
CMPIINSERT_ENDIF_RANK
        call cpu_time(t_urf3) !-.-.-.-.-.-.-.-.-.-.-.-.-.
        !YuP[03-2016] Repeat plotting surfaces, but now - with rays
        if (noplots.ne."enabled1" .and. eqmod.eq."enabled"
     +      .and. n.eq.0) then
           if(mrfn.gt.0)then ! just in case 
              !(mrfn is supposed to be >0 when urfmod.ne."disabled")
              do krf=1,mrfn
               call tdplteq(krf) ! separate page with rays for each krf
              enddo
           endif
        endif
      endif ! urfmod.ne."disabled"

c.......................................................................
c     cqlpmod.eq.enabled:
c     Compute parallel electric field from Poisson equation
c.......................................................................

      if (transp.eq."enabled" .and. cqlpmod.eq."enabled")
     +  call wpelecf(11)  ! eleven, not ll

c.......................................................................
c     cqlpmod.eq.enabled:
c     Compute electrostatic parallel electric field from constant flux
c     condition (for CQLP + current drive localized along field line).
c     Uses Kupfer (PoP, 1995) f and g function decompostion in achiefn.
c     This is a nonlocal Efield calc, and therefore all space 
c     points must be cycled through.
c     BH131103 Comment added:  This option is not operational.
c              Compare with the stella code, where it is implemented.
c              Also compare with ampfmod coding.
c.......................................................................

      if (cqlpmod.eq."enabled" .and. eseswtch.eq."enabled" .and.
     +    sbdry.eq."enabled") then
         kopt=2
         ilend=lrors
         do ll=1,ilend
            call tdnflxs(ll) 
            call achiefn(kopt) ! achiefn(2) here (cqlpmod="enabled")
         enddo
         call efld_cd(dz(1,lr_),lrors,vnorm,flux1,flux2,elparnw,flux0)
      endif

c.......................................................................
c     MAIN VEL TIME STEP: loop over each radial position sequentially
c.......................................................................

      ilend=lrors
      if (transp.eq."enabled" .and. cqlpmod.eq."enabled" .and.
     +  mod(nummods,10).ge.5 .and. sbdry.ne."periodic"
     +  .and. lmidvel.ne.0) ilend=lrors-1

      if (nstop.eq.0) go to 2

cBH131230:  Why won't ampfmod work before MAIN VEL TIME STEP?? nstop=0 issue?
c.......................................................................
c     Computer the time-advanced toroidal electric field from Ampere-
c     Faraday equations using Kupfer h and g functions.
c     
c     cqlpmod.ne.enabled has been checked.
c     Initialize electric fields at turnon of ampfmod.
c.......................................................................


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
c        Find dtr for this time-step
         do i=1,ndtr1a
            if ((n+1).eq.nondtr1(i)) then
               dtr=dtr1(i)
               dtreff=dtr
               dttr=dtr*nrstrt
            endif
         enddo
c        Copy current distribution f into f_
         call dcopy(iyjx2*ngen*lrors,f(0,0,kelec,1),1,f_(0,0,kelec,1),1)
c        Bring background profiles up to time step n
         ! No effect if bctime=0 (time-indep. profiles).
         call profiles ! if(ampfmod.eq.'enabled' .and. n+1.ge.nonampf)
                       ! skip elecfld background profiles
         
c        Get given boundary condition (i.e., edge) elec field, at time
c        advanced position (i.e., elecfldn(ll,n+1,0) at ll=lrz+1 point).
c        If time-dependent (nbctime>0),
c        then it is given by time advanced elecb() or elecin_t(njene,).
         call ampfefldb(n+1,time+dtr)
c        Besides, subr.ampfefldb Sets zero iteration of elecfldn 
c        equal to previous time step radial profile, like this:
         ! elecfldn(ll,nn,0)=elecfldn(ll,nn-1,it_prev) ! where nn==n+1 (and ll=0:lrz)
         ! So, if nonampf=1, then we have here n=0, or nn=1,
         ! which means that we need the values of elecfldn at nn-1=0.
         

c        elecfld(1:lrz) has been set.  If n=1=nonampf, then is given;
c                       If n.ge.nonampf, then set near end of ampfsoln.

c        Determine the time-advanced tor e-field
         kopt=3  !Amp-Faraday option
         ilend=lrz

c        Initialize elecfldn(ll,n+1,it=0) for upcoming ampfsoln
         do ll=1,ilend
            elecfldn(ll,n+1,0)=elecfld(ll)/300.d0
         enddo
         !YuP[21-08-2017] added but then commented out: 
         !elecfldn(lrz+1,n+1,0)=elecfldb/300.d0 !Not needed, done in ampfefldb

c        Obtain the h and g functions over all radii and iteratively
c        solve Ampere-Faraday eqns for the tor electric field.
         do it=1,nampfmax
            it_ampf=it

cBH131109: Maybe don't need from bndry to center for fh/fg?   do ll=ilend,1,-1
cBH131109: Backwards is giving bounds problem with abd in impavnc0.
cBH131109: Present system is determining all fh/fg first as function
cBH131109: of radius, then solve AF eqn.  Future solves may require
cBH131109: fix the counter resetting for ll=ilend,1,-1.

            do ll=1,ilend
               call tdnflxs(ll)  !get l_,lr_,..
               call achiefn(kopt) !kopt=3 here !Increments n by 1, for each ll
                                   !Gives fh,fg(,,1,ll) for each ll.
                                   !it_ampf is passed in common block
               ! YuP test/printout
               amp_f_=ampfarl(f_(0,0,kelec,ll),ll)*dtr !ampfarl has 1/dtr factor
               amp_f=ampfarl(f(0,0,kelec,ll),ll)*dtr !ampfarl has 1/dtr factor
               amp_h=ampfarl(fh(0,0,kelec,ll),ll)*dtr
               amp_g=ampfarl(fg(0,0,kelec,ll),ll)*dtr
               write(*,'(a,3i4,3e12.4)')
     +           'after achiefn(3): n,it,ll,integrals f, f-h, g:',
     +            n,it,ll, amp_f, amp_f-amp_h, amp_g  
               ! YuP test/printout: 
               !Confirmed that f() remains unchanged during
               !iterations in this it loop, and it remains equal to f_()
            enddo
c           Solve for iteratively updated tor electric field
c           n is time advanced value
            call ampfsoln(it,n) !here n=nonampf,...,nstop
c           Presently, ampfdiff does nothing, i.e. gives iflag=1
c           In future, refine this, making it=it(1:lrz).
            itt=it
            call ampfdiff(iflag)  !iflag=0, only if error criteria met.
                                  !Presently a dummy call giving iflag=0
            if (iflag.eq.0) go to 16
         enddo !On it (Note: it is incremented by 1 at the exit of loop)
 16      continue
 
         !write(*,'(a,3i4)')'tdchief-373: n,itt,it=', n,itt,it

c        Now, setup with time-advanced electric field in elecfld(),
c        then at end of ampfmod interations, call achiefn(0) below.
      endif  !On ampfmod


c...........................................................
c     Loop back point for electric field iteration
c       (if efiter.eq."enabled")
c...........................................................
      nefiter=1 ! counts iterations
      do ll=1,ilend
        call tdnflxs(ll) ! determine l_,lr_, etc.
        nefiter_(l_)=1 ! counts iterations for each flux surface
      enddo
      
 20   continue  !Loop back point for electric field iteration

      call cpu_time(t_before_soln) !-.-.-.-.-.-.-.-.-.-.-.-.-.
CMPIINSERT_BARRIER

c     Copy current distribution f into f_
      call dcopy(iyjx2*ngen*lrors,f(0,0,1,1),1,f_(0,0,1,1),1)
      
      if (transp.eq."enabled" .and. n.ne.0 .and. adimeth.ne."enabled"
     +      .and. soln_method.ne."it3drv" .and. nefiter.eq.1)  then
c.......................................................................
c     With adimeth, one assumes f (in tdtrrsou or wparsou) 
c     as the starting point of the new iteration.
c     Thus should not need to redefine here.
c     CQLP case: f defined in wpsavf is assumed to be the good one.
c     soln_method=it3drv:  soln is in f(,,,) on last call on the set of
c     flux surfaces (soln at start of step is in f_(,,,)).
c.......................................................................
        if(cqlpmod.ne."enabled" .and. n.ge.nontran .and. n.lt.nofftran)
     +      call dcopy(iyjx2*ngen*lrors,frn_2(0,0,1,1),1,f(0,0,1,1),1)
      endif
      
c..................................................
c     Bring background profiles up to time step n.
c     Exclude updating of electric field elecfld(1:lrz) for 
c     ampfmod.eq.enabled .and. n.ge.nonampf [Check this for
c     the eseswtch case also.]
c..................................................

      if(nefiter.eq.1) call profiles
                 
c.......................................................................
cBH131103      do 1 ll=1,ilend   !ilend=lrz for cqlpmod.ne.'enabled'
c
c     Can reverse order of do 1 ll=1,ilend, if needed for purposes
c     of solving Ampere-Faraday eqns starting at given edge voltage.
c     Results are unchanged for simple efield case, BUT additional
c     debugging required for TCV ECCD/radial diffusion test case (CD
c     profiles are similar, but total driven current changes from 
c     108 kA to 138 kA.)
cBH131107 Moved the ll=1,ilend reversal to above if(ampfmod.eq.enabled..
c.......................................................................
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

c......................................................................
c     If ampfmod.eq.enabled, but n.lt.nonampf, then insert values of
c     elecfldn(1,,) from elecfld [for plotting with .nc].
c......................................................................

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
cBH170329
         !YuP[21-08-2017] Corrected, for the fields storage BEFORE the ampf is applied:
         do niter=0,nampfmax ! or up to nefitera
           elecfldn(ll,n+1,niter)=elecfld(ll)/300.d0
           ! here: ll=1:lrz. But we also need the end (bndry) point:
           elecfldn(lrz+1,n+1,niter)=elecfldb/300.d0 !here: n+1.lt.nonampf
           ! Can also save ll=0 point (magn.axis point):
           elecfldn(0,n+1,niter)=elecfld(0)/300.d0
         enddo
      endif


c......................................................................
c     TIME ADVANCE ON FLUX SURFACE rz(ll).
c     If cqlpmod.ne.'enabled', transp='enabled', soln_method='it3drv',
c       then following call only adds ll flux surface contributions to
c       the coefficient matrix, and then solves for time advanced distn
c       at last flux surface, ll=lrz.
c......................................................................

CMPIINSERT_MPIWORKER
c    It will insert :
c      if(soln_method.eq.'direct' .and. lrzmax.gt.1) then
c         ! Parallelization for the impavnc0 solver is limited 
c         ! to soln_method='direct' (for now)
c         mpiworker= MOD(ll-1,mpisize-1)+1  !1...(mpisize-1)
c      else
c         ! In all other cases, perform calculations 
c         ! for all flux surfaces on mpirank=0, then broadcast results
c         mpiworker=0
c      endif

        call achiefn(0)  !--> if implct='enabled', calls impavnc0


CMPIINSERT_SEND_RECV
c     It will send or recv data, but only in case of
c     soln_method='direct' (for now)

 1    continue ! End loop over radius:  New f is obtained for each ll
 
CMPIINSERT_BARRIER
CMPIINSERT_BCAST_DISTRIBUTION
CMPIINSERT_BCAST_COLL_COEFFS
CMPIINSERT_BCAST_SCAL
CMPIINSERT_BCAST_VELSOU
      call cpu_time(t_after_soln) !-.-.-.-.-.-.-.-.-.-.-.-.-.

      do ll=1,ilend  ! Re-scale f if needed; 
        call tdnflxs(ll) ! determine l_,lr_, etc.
        do k=1,ngen  ! Compute density gains and losses, and powers.
           ! For lbdry0='disabled',  Redefine f at v=0 so it is unique:
           ! (For lbdry0='enabled', coeff matrix is set up 
           !   to automatically maintain unicity.)
           if (lbdry0.ne."enabled") then !-YuP: moved here from impavnc0
             call dcopy(iyjx2,f(0,0,k,l_),1,fxsp(0,0,k,l_),1)
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

c..................................................................
c     Check/iterate on electric field to obtain a target current
c..................................................................
      if (efiter.eq."enabled") then
         ! Stop iterations, if beyond iteration limit:
         if (nefiter .gt. nefitera) go to 21
         ! Stop iterations if n<(control turn on step):
         if (n. le. noncntrl) go to 21
         write(*,*)'TDCHIEF/EFLDITER: Time step=',n,
     ~             '    Starting iteration #',nefiter
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
            call dcopy(iyjx2,f_(0,0,1,l_),1,f(0,0,1,l_),1)           
         enddo ! ll
         nefiter=nefiter+1 ! counts iterations
         go to 20 ! Another iteration, using old f() but new elecfld()
      endif
      
 21   continue

CMPIINSERT_BARRIER

c     Accumulate favg time average, if tavg.ne."disabled"
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
CMPIINSERT_IF_RANK_EQ_0
           WRITE(*,*)'tdchief WARNING: No time averaging of f was done'
           WRITE(*,*)'tdchief WARNING: Probably t=nstop*dtr is less'
           WRITE(*,*)'tdchief WARNING: than tavg1(1).'
           WRITE(*,*)'tdchief WARNING: Instead of favg(), '
           WRITE(*,*)'tdchief WARNING: the last available f() will be'
           WRITE(*,*)'tdchief WARNING: saved into mnemonic.nc file. '
CMPIINSERT_ENDIF_RANK
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
      
CMPIINSERT_BARRIER
      
c......................................................................
c     Integrate power densities over plasma volume
c......................................................................

        call diagentr_vol
        
      call cpu_time(t_after_diag2) !-.-.-.-.-.-.-.-.-.-.-.-.-.

 2    continue ! if nstop.eq.0

      if (ilend .eq. lrors-1) n_(lrors)=n_(lrors-1)

c.......................................................................
c     Time advance spatial transport equation [except if transport
c     with soln_method=it3drv, not applicable since done in impavnc0].
c.......................................................................

      if (soln_method.ne.'it3drv') then
      if (transp .eq. "enabled" .and. n.ge.nontran .and. n.le.nofftran)
     +   then

        if (cqlpmod .ne. "enabled") then
c......................................................................
c     The radial transport (tdtr...) routines follow. Store the results
c     of the velocity space split in fvn. Then call the transport
c     subroutine, tdtransp.
c......................................................................
          call tdtransp
          call tdtrsavf

        else  ! i.e., cqlpmod.eq.enabled

c.......................................................................
c     Parallel transport.
c     First compute E_n+1/2 from Poisson and modify velsou accordingly
c.......................................................................
          call wpelecf(2)
c%OS  call wptrans
c%OS  call wpsavf
c%OS  
          dtreff=dtreff/noffso(1,1)
          do 200 in=1,noffso(1,1)
            if (meshy .eq. "fixed_mu") call wptramu
            if (meshy .ne. "fixed_mu") call wptrafx
            call wpsavf
 200      continue
          dtreff=noffso(1,1)*dtreff
c%OS  
        endif ! on cqlpmod

      endif  ! on transp, etc.
      endif  ! on soln_method.ne.it3drv

c.......................................................................
c     Plot flag and plot time set:
c.......................................................................

      iplt3d=0
      if (noplots.ne."enabled1") then
         do i=1,nplota
            if(n.eq.nplt3d(i)) then
               iplt3d=i
               tplt3d(i)=timet
            endif
         enddo
      endif

c.......................................................................
c     Diagnostics (if transp=enabled)
c.......................................................................

      do 30 ll=1,lrors
c......................................................................
c     determine local variables depending on flux surface (l_, iy,..)
c......................................................................
        call tdnflxs(ll)
                
        if (transp.eq."enabled") then
          call diaggnde

c..................................................
c     Compute plasma resistivity for electron runs.
c..................................................
          call restvty

c..................................................
c     Compute various diagnostics...
c..................................................
          call diag

c..................................................
c     Obtain data for time dependent plots
c..................................................
          call ntdstore

c..................................................
c     plotting logic for 3-D code...
c       (iplot set in subroutine achiefn)
c..................................................
          if ((n.ge.nstop .or. iplot.ne.0 .or. iplt3d.ne.0).and.
     +        n.ge.nontran .and. n.le.nofftran) then
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
            if ((adimeth.eq."enabled" .or. cqlpmod.ne."enabled") .and.
     +        n.ge.nontran .and. n.le.nofftran) then
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
c$$$  Commenting following, in favor of dimensioning pwrrfs with l_,
c$$$    to save results calculated in impavnc0 ==> diagimpd ==> diagentr.
c$$$    BH090806
c$$$            do k=1,ngen
c$$$               call diagentr(3,k) !To get pwrrfs at current l_ surface
c$$$            enddo
            if (noplots.ne."enabled1" .and. 
     +          (iplot.ne.0 .or. n.ge.nstop) ) then
                !YuP[2017-11-27] Added n.ge.nstop, similar to achiefn
c$$$            if ( iplot.ne.0) then
              call pltmain
            endif
          endif
c.......................................................................
c     write some radially dependent netcdf output, similar to pltmain.
c.......................................................................

          call netcdfmain

        endif !End of transp="enabled" if-block

c......................................................................
c     Store updated information.
c......................................................................
        call tdtoaray
 30   continue !  ll=1,lrors

c.......................................................................
cl    Compute flux surface average current and resistivity for CQLP case
c.......................................................................

      if (transp.eq."enabled" .and. cqlpmod.eq."enabled") call wpavg

c......................................................................
c     Compute conservation constant for transp="enabled"
c......................................................................

      if (transp.eq."enabled" .and. cqlpmod.ne."enabled") then
        call tdtrcon
      endif

c.......................................................................
c     compute some radial diagnostics.
c.......................................................................

c%OS  still to verify and change if cqlpmod=enabled
      call tddiag

c.......................................................................
c     Calc fusion rates, if sigmamod=enabled
c.......................................................................

        if (sigmamod.eq."enabled" .and. msig.gt.0) then
          icall="notfrst"
          call sigv(icall)
        endif ! sigmamod

c.......................................................................
c     Calc xray spectra, if softxry.ne."disabled".and.kelecg.gt.0
c     Spectra is stored in the netCDF output file at each time step.
c     Also, set plot flag  [iplt3d is set above, if plotting this step].
c.......................................................................

      if (lrzmax.lt.3 .and. softxry.ne."disabled") then
         write(*,*)'*******************************************'
         write(*,*)'tdchief:  SXR not computed for lrzmax.lt.3'
         write(*,*)'*******************************************'
      else
         if (softxry.ne."disabled".and.kelecg.gt.0) then
            icall="notfrst"
            ilold=l_
            iplotsxr='no'
            if (iplt3d.ne.0 .or. n.eq.nstop .and. noplots.ne.'enabled1')
     +           iplotsxr='yes'
            call tdsxray(icall,iplotsxr) !Soft Xray diagnostic
            call tdnflxs(ilold)
         endif
      endif


c.......................................................................
c     NPA synthetic diagnostic
c.......................................................................

        if (npa_diag.ne."disabled".and.niong.ne.0) then
           if (n.eq.nstop .or. npa_diag.eq.'ncdf_all') then 
              call tdnpadiag(icall)
           endif   
        endif

cBH090513:  BOB, check why no call tdnflxs after tdnpadiag??

c.......................................................................
c     plot radial information
c.......................................................................

      if ((iplt3d.ne.0 .or. n.eq.nstop).and.lrzmax .ge. 3) then
        icall="notfrst"
        if (noplots.ne."enabled1") then
          call tdpltmne
c         if equilibrium changed during run:
C%OS      if (eqmod.eq."enabled".and.n.ge.nstop) call tdplteq
            ! Plot the change in Total N of ptcles as a func of rho, for different time.
            call plt_fow_cons ! Can be used for ZOW
        endif
      endif

      if (n.eq.nstop.and.pltra.eq."enabled") then
        call pltrun
      endif

c.......................................................................
c     write radial information
c.......................................................................

      if (iplt3d.ne.0 .or. n.eq.nstop) then
        call tdoutput(2)
      endif

c.......................................................................
c     write netCDF output file, if netcdfnm.ne."disabled"
c.......................................................................
CMPIINSERT_IF_RANK_EQ_0
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
CMPIINSERT_ENDIF_RANK
c.......................................................................
c     re-write ray data at last step, if urfmod & urfwrray="enabled"
c.......................................................................

      if (urfmod.ne."disabled" .and. urfwrray.eq."enabled" 
     +  .and. n.eq.nstop .and. nstop.ne.0) call urfwrite

c.......................................................................
c     Write tavg12, sumdtr, for tavg case
c.......................................................................

      if (tavg.ne."disabled") then
         tavg12=0.d0
         do ii=1,ntavga
            tavg12=tavg12+(tavg2(ii)-tavg1(ii))
         enddo
         write(*,*)
         write(*,*)'tavg.ne.disabled: tavg12 and sumdtr should be close'
         write(*,*)'at end of run: tavg12,sumdtr = ',tavg12,sumdtr
      endif

c.......................................................................
c     check for call exit
c.......................................................................
      call cpu_time(t_end) !-.-.-.-.-.-.-.-.-.-.-.-.-.
CMPIINSERT_BARRIER
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,'(a,f10.3,a,f10.3,a,f10.3)')
     +                'TDCHIEF: tcpu_urfchief=',t_urf2-t_urf1,
     +                '    tcpu_netcdf=',   t_urf3-t_urf2,
     +                '    tcpu_impavnc=', t_after_soln-t_before_soln
     
         WRITE(*,'(a,2f10.3)') ' tcpu_diagnostics 1 and 2:',  
     +                          t_after_diag1-t_after_soln,
     +                          t_after_diag2-t_after_diag1
     
         WRITE(*,'(a,i5,a,f10.3)') ' Finished Time Step ======>',  n,
     +                          '     tcpu(sec.)==',  t_end-t_n_start
         WRITE(*,'(/)')
CMPIINSERT_ENDIF_RANK

CMPIINSERT_BARRIER
      call tdtloop

      if (n .lt. nstop) go to 10 !-> next time step

c.......................................................................
c     write out radial diffusion coeffs at end of program, if enabled.
c.......................................................................
      if (difus_io(1).eq."drrout") then
         call diffus_io(0)
         call diffus_io(1)
      elseif (difus_io(1).eq."drrdrout") then
         call diffus_io(0)
         call diffus_io(2)
      endif
      
      return
      end
