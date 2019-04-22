module tdinitl_mod

!
!

contains

      subroutine tdinitl ! called only at n=0
      use param_mod
      use comm_mod
      use tdeqdsk_mod, only : tdeqdsk
      use ampfar_mod, only : ampfalloc
      use ainvnorm_mod, only : ainvnorm
      use ainspec_mod, only : ainspec
      use ainsetva_mod, only : ainsetva
      use ainpla_mod, only : ainpla
      use aingeom_mod, only : aingeom
      use ainitial_mod, only : ainitial
      use aindflt_mod, only : aindflt
      use aindflt1_mod, only : aindflt1
      use ainalloc_mod, only : ainalloc
      use aindfpa_mod , only : ain_transcribe
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!..................................................................
!     This routine initializes the arrays required by CQL3D.
!     It also controls the initialization of distribution
!     functions and quasilinear coefficients.
!..................................................................


      include 'name.h'
!MPIINSERT_INCLUDE

      character*8 icall,iplotsxr
      character*128 filenm ! template for file name with data
      real*8 a_new(njene)  ! working array

      if (n.gt.0) return

!..................................................................
!     set defaults for namelisted quantities and other variables
!..................................................................

      call aindflt
      call eqindflt
      call urfindfl
      call aindflt1

!.......................................................................
!     Read namelist from previous run if rerun case. Thus namelist from
!     this run should only reset some diagnostic or numerical variables.
!     This call will position pointer in file distrfunc to next read f.
!.......................................................................

!BH110201      if (nlrestrt.ne."disabled") call tdreadf(1)
      if (nlrestrt.eq."enabled") call tdreadf(1)

!..................................................................
!     read in namelist input for CQL3D
!..................................................................
      open(unit=2,file='cqlinput',status='old')
        read(2,setup)
        read(2,trsetup)
        read(2,sousetup)
        read(2,eqsetup)
        read(2,rfsetup)
      close(2)

      if (partner.eq."selene") then
        open(unit=18,file='kcqlsel',status='old')
        read(18,2) ncount,noplots
        eqdskalt="enabled"
        close(unit=18)
      endif
 2    format(i5,a8)

!..................................................................
!     Call routine which finds electron and ion species indices.
!..................................................................

      call ainspec

!.......................................................................
!     set variables dependent on input variables
!     (May add an impurity ion, increasing species indices, if
!      iprozeff.ne.'disabled').
!.......................................................................

      call ainsetva

!.......................................................................
!     Allocate arrays and radial transport arrays, if required
!.......................................................................

      call ainalloc
      if (ampfmod.eq."enabled" .and.cqlpmod.ne."enabled") call ampfalloc

      if (transp.eq."enabled" .and. cqlpmod.ne."enabled") call tdtraloc
!%OS  if (transp.eq."enabled") call tdtraloc
      if (transp.eq."enabled" .and. cqlpmod.eq."enabled") call wpalloc

!.......................................................................
!     print namelists
!.......................................................................

      if (nmlstout.eq."enabled") then
         open(unit=2,file='cqlinput',delim='apostrophe',status='old')
         write(6,*)'  In tdinitl: '
         write(6,setup0)
         write(6,setup)
         write(6,trsetup)
         write(6,sousetup)
         write(6,eqsetup)
         write(6,rfsetup)
         close(2)
      elseif (nmlstout.eq."trnscrib") then
         write(6,*)'  In tdinitl: '
         call ain_transcribe("cqlinput")
      else
         write(6,*)
         write(6,*) 'mnemonic = ',mnemonic
         write(6,*)
      endif

!..................................................................
!     Call the initialization routines for the appended modules..
!..................................................................
      call eqinitl
      call urfinitl ! mrfn= is set here also
      call frinitl
      open(unit=2,file='cqlinput',delim='apostrophe',status='old')
      call frset(lrz,noplots,nmlstout)   ! Uses unit 2
      close(2)

      if (machine .ne. "toroidal") call tdwrng(1)
      if (lrzdiff.eq."enabled" .and. frmodp.eq."enabled") &
        call diagwrng(18)

!.....................................................................
!     This routine initializes the normalized theta mesh that is used
!     if the meshes on the various flux surfaces are linked.
!.....................................................................

      call tdtrmuy

!.....................................................................
!     Call routine which initializes certain arrays using input
!     data. Also initialize the radial (rz) mesh and some plasma
!     parameters arrays.
!.....................................................................

      if(lrzmax.gt.1) call tdxinitl

!.....................................................................
!     Call routines to initialize any time-dependent (parabolic)
!     profiles, and to set them up in place of profiles
!     otherwise set up in tdxinitl.
!.....................................................................

      if (nbctime.gt.0) then
         call profiles
      endif

!.....................................................................
!     Determine mesh normalization constant vnorm.
!.....................................................................

      call ainvnorm

!.......................................................................
!l    1.1 Loop over all radial mesh points lr=1,..,lrzmax to determine
!     the plasma and equilibrium parameters on the full radial mesh
!     Also determine the mesh parallel to the magnetic field line
!.......................................................................

      if (cqlpmod.eq."enabled" .and. numclas.eq.1 .and. ls.eq.lsmax)then
        lz=lz/2+1
        lsmax=lsmax/2+1
        ls=ls/2+1
      endif

      do 110 ll=lrzmax,1,-1

!.......................................................................
!     Sets up the parameters on each lrzmax flux surface. Thus does not
!     define lr_=lrindx(l_), but lr_=l_, l_=1,lrzmax
!     (Not like in tdnflxs)
!......................................................................

        l_=ll
        lr_=ll
!     these two indices should not be used in loop 110 => set to -1, to detect
!     out of bound errors
        lmdpln_=-1
        indxlr_=-1

!..................................................................
!     Call an initialization routine which determines flux surface
!     geometry and magnetic field structure.
!..................................................................

        call aingeom !-> eqcoord-> eqfndpsi-> eqorbit-> trace flux surf.
             ! solr(l,lr_), solz(l,lr_) are R,Z coords. of flux surface
!.......................................................................
!     Initialize mesh along the magnetic field line, as well as
!     density and temperature profiles if cqlpmod=enabled
!.......................................................................

        call micxiniz

!..................................................................
!     Copy some radial and non-time dependent diagnostic quantities
!..................................................................

        call tdtoarad

 110  continue ! ll=lrzmax,1,-1

!MPIINSERT_IF_RANK_EQ_0
      do ir=1,lrz
         WRITE(*,'(a,i3,4e13.5)')'ir,rya,rpcon,rmcon,equilpsi=', &
                      ir,rya(ir),rpcon(ir),rmcon(ir),equilpsi(ir)
      enddo
       !pause
!MPIINSERT_ENDIF_RANK



!.......................................................................
!     Determine equilibrium parameters at lower half of cross-section
!     on lrindx(1) flux surface
!.......................................................................

      if (cqlpmod.eq."enabled" .and. numclas.eq.1 .and. ls.eq.lsmax)then
        lz=2*(lz-1)
        lsmax=2*(lsmax-1)
        ls=2*(ls-1)
        call wploweq
      endif

!.....................................................................
!     Redefines mu-mesh at midplane if needed
!.....................................................................

      if (cqlpmod.eq."enabled" .and. meshy.eq."fixed_mu" .and. &
        tfac.lt.0.0) call wptrmuy

!.......................................................................
!l    1.1.2 Initialize some plasma parameters on entire lrzmax mesh
!.......................................................................

      call ainpla

!.....................................................................
!l    1.1.3 First Loop over spatial variable index for which the equations
!     will be solved to compute y-mesh on all l_ indices
!     (this takes micxinit out of ainitial subroutine)
!.....................................................................

      do 113 ll=lrors,1,-1

!......................................................................
!     determine local variables depending on variable index l_
!......................................................................

        call tdnflxs(ll)

!.......................................................................
!     call a routine to determine meshes y, x and related quantities
!.......................................................................

        call micxinit

 113  continue

      ieq_tot=0
      do ll=1,lrors
         ieq_tot=ieq_tot+inewjx_(ll) ! inewjx_() is defined in micxinit
         if (ll.eq.1) then
            ieq_(ll)=1
         else  ! Eqn no. at beginning of each flux surface:
            ieq_(ll)=ieq_(ll-1)+inewjx_(ll-1)
         endif
      enddo
      ieq_(lrors+1)=ieq_tot



!.......................................................................
!     adapt dyp5(iyh) and dy(iyh) such that sum(cynt2) over i is constant
!     on all l_ surfaces
!.......................................................................

      if (nchgdy .eq. 1) then
        do 114 ll=1,lrors
          call tdnflxs(ll)
          call wpchgdy
 114    continue
      endif

      ! YuP[Sept-2014] moved this part outside of ainitial:
      ! needed for proper work of FOW-logic.
      do ll=lrors,1,-1
        call tdnflxs(ll)
        call micxinil ! integration coefficients, mostly for cfpcoefn:
        ! gamma,alphan,cog,cosz,zmaxpsi (uses bbpsi,solrz from micxiniz)
        if (l_ .eq. lmdpln_) then
           ! calculate tau()==v*tauB_ZOW [also dtau(), vptb()]
           if(taunew.eq."enabled") then
              call baviorbt ! must be called after micxiniz
           else
              call baviorbto
           endif
        endif
      enddo

!....................................................................
!     Set some more radial mesh quantities.
      call tdrmshst !must be called AFTER: aingeom,tdtoarad,micxinit,micxinil
!....................................................................

!BH160530
!     Must be called after micxiniz and baviorbt/baviorbto:
!YuP      if (lossmode(1).eq.'simplban' .or. lossmode(1).eq.'simplbn1'
!YuP      +     .or.ndeltarho.ne.'disabled')
!YuP      +     call deltar
!YuP[07-2017] exclude simplban  - deltar is not needed.
      if (lossmode(1).eq.'simplbn1' &
           .or.ndeltarho.ne.'disabled') &
           call deltar

!.....................................................................
!l    1.2 Loop over spatial variable index for which the equations will be
!     solved. Compute initial distributions and sources etc.
!.....................................................................

      do 120 ll=lrors,1,-1

!......................................................................
!     determine local variables depending on variable index l_
!......................................................................

        call tdnflxs(ll)

!.....................................................................
!     transfer input data from CQL3d to CQL
!.....................................................................

        call tdstin

!.....................................................................
!     Call spatial variable "surface" initialization routine
!..................................................................

        call ainitial ! mrfn= is set here also (in vlf or vlh)

!..................................................................
!     Call a routine to initialize some urf module routines.
!..................................................................

!     Tried following, but gives problem of code stopping?? [BH041111].
!        if (cqlpmod.ne."enabled".and.urfmod.ne."disabled") call urfsetup
        if (cqlpmod.ne."enabled") call urfsetup ! mrfn= is set here also

!..................................................................
!     Copy some diagnostic quantities
!..................................................................

        call tdtoaray

 120  continue
!----------------------------------------------------------------------

!..................................................................
!     If running with the Brambilla code, create and dump and eqdsk
!     file for Brambilla ray tracing code to read in.
!..................................................................

      if (partner.eq."bramb") then
        call tdeqdsk
      endif

!....................................................................
!     Set some more radial mesh quantities.
!....................................................................

      call tdrmshst

!..................................................................
!     Call first-order radial orbit-shift module
!..................................................................

!BH160530      if (ndeltarho.ne."disabled") call deltar

!....................................................................
!     Call routines to set up radial integration coefficients, radial
!     diffusion and advection coefficients, certain flags and the radial
!     Chang-Cooper weights.
!....................................................................

      if (transp.eq."enabled") then
        if (cqlpmod .ne. "enabled") then
!     warning: these routines use lrors, not lrzmax
          call tdtrvint

!         Obtain radial diffusion d_rr, and optionally
!         pinch velocity d_r, from file.  In this case, the
!         d_rr and d_r can be multiplied by a time-dependent
!         scale factor at point of use, in trans.h and tdtravct.
!         The read in d_rr/d_r won't be changed.
          if (difus_io(1).eq."drrin") then
             call diffus_io(3) !Reads d_rr
             if (n_d_rr.ne.0) then
                !Getting time-dep scale factor at t=0
                do k=1,n_d_rr   !n_d_rr obtained in diffus_io
                   drrt(k)=difus_io_scale(k,1)
                enddo
             endif
          elseif(difus_io(1).eq."drrdrin") then
             pinch="disabled"   !Since getting d_r from file
             call diffus_io(4) !Reads d_rr and d_r
             if (n_d_rr.ne.0) then
                !Getting time-dep scale factor at t=0
                do k=1,n_d_rr   !n_d_rr obtained in diffus_io
                   drrt(k)=difus_io_scale(k,1)
                   drt(k)=difus_io_scale(k,2)
                enddo
             endif
          else
             call tdtrdfus
          endif

          call tdtrflg
          call tdtrwtl
        else
          call wpinitl
        endif
      endif
      call tddiag

!%OS
!     check some mesh values
      call wpmshchk
!%OS


!.....................................................................
!     Set up sigma-v calculation
!.....................................................................

      if (sigmamod.eq."enabled")  then
        icall="first"
        call sigv(icall)
      endif

!.....................................................................
!     Set up soft X-ray and NPA diagnostic.
!.....................................................................

      if (lrzmax .ge. 3) then
!**bh if (softxry.eq."enabled".and.kelecg.gt.0.and.eqmod.ne."enabled")
        if (softxry.ne."disabled".and.kelecg.gt.0) &
          then
          icall="first"
          if (noplots.ne."enabled1") then
             iplotsxr='yes'
          else
             iplotsxr='no'
          endif
          call tdsxray(icall,iplotsxr)
        endif
        if (npa_diag.ne."disabled".and.niong.ne.0) then
           icall="first"
!BH_to skip first call, per Vincent Tang:     call tdnpadiag(icall)
!BH100816:  No.      Also, actually for npa, icall has no effect,
!BH100816:  as calculation as setup repeated each call to tdnpadiag.
           call tdnpadiag(icall)
        endif
      endif

!..................................................................
!     Call the neutral beam source module
!..................................................................

      write(*,*)'tdinitl:call frnfreya, n=',n,'time=',timet
      call frnfreya(frmodp,fr_gyrop,beamplsep,beamponp,beampoffp, &
                    hibrzp,mfm1p,noplots)
       write(*,*) 'mfm1 (freya) = ',mfm1p
!       write(*,*) 'hibrzp(i,1,1),hibrzp(i,2,1),hibrzp(i,3,1)'
!       do i=1,mfm1p
!         write(*,'(i4,2x,0p9f9.4)') i, hibrzp(i,1,1),
!     >        hibrzp(i,2,1),hibrzp(i,3,1),
!     >        hibrzp(i,1,2),hibrzp(i,2,2),hibrzp(i,3,2),
!     >        hibrzp(i,1,3),hibrzp(i,2,3),hibrzp(i,3,3)
!       enddo
!      stop
!
!     Initialize ibeampon/ibeamponp (present and previous time step)
      if (frmodp.eq."enabled") then
         ibeampon=0   !Reset in tdchief, for pulse beam
         ibeamponp=1  !Indicates beam source calculated
      else
         ibeampon=0
         ibeamponp=0
      endif

!..................................................................
!     Plot out initial radial profiles
!..................................................................
      if (noplots.ne."enabled1") call tdpltmne

!..................................................................
!     Read RF diffusion coefficient data for multiple flux
!     surfaces if rdcmod="format1" or "aorsa" or "format2"
!..................................................................

      if (rdcmod.eq."format1" .or. rdcmod.eq."aorsa" &
                              .or. rdcmod.eq."format2") then
        call rdc_multi
      endif

!..................................................................
!     print initial output and plot equilibrium surfaces
!..................................................................
      call tdoutput(1)
      if (noplots.ne."enabled1" .and. eqmod.eq."enabled") then
         call tdplteq(0) !'0' means no rays (Here, no data on rays yet)
      endif
!$$$      if (eqmod .eq. "enabled") call tdplteq
!
      return
      end
end module tdinitl_mod