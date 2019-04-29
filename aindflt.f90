module aindflt_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine aindflt
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Set namelist input defaults for all namelist sections
!     except setup0/fsetup and frsetup.
!     Warning: should not set variables read in setup0/fsetup namelist
!     as aindflt is called AFTER the first read(2,setup0).
!     BH070305:  Some other constants derived from namelist
!                variables have been moved to new subroutine aindlft1.
!..................................................................


!     Set a few local constants, same as in subroutine ainsetpa.
      ep100=1.d+100
      zero=0.d0
      one=1.d0
!      pi=3.141592653589793d0
      pi=atan2(zero,-one)

      iy=200 ! default value; will be over-written by cqlinput value
      jx=300 ! default value; will be over-written by cqlinput value
      mx=3   ! default value; will be over-written by cqlinput value

      jfl=151  !!!min(201,jx)
      if (mod(jfl,2).eq.0) jfl=jfl-1
      ! jfl needed to be odd because of jpxyh=(jfl+1)/2 in pltprppr.f

      nmods=nmodsa ! YuP-101220: should be mrfn, but not known yet

      nso=0
      lz=lza

      ampfmod="disabled"
      nampfmax=2
      ampferr=1.d-3
      nonampf=0

      bootst="disabled"   ! analytic (Hinton and Haseltine)
                          ! bootstrap current
      bootcalc="disabled" !computational bootstrap current off.
      bootupdt="disabled" !updating 0th order distn for bs radial derv.
      bootsign=+1.0
      nonboot=2           !turn on computational bootstrap at n=nonboot.
      bremsrad="disabled"
      brfac=0.
      brfac1=0.
      brfacgm3=1.0
      isoucof=0
      faccof=1.e0
      bth=1000.
      btor=10000.
      chang="enabled"
      constr=1.d-3
      contrmin=1.d-12
      curr_edge=0.
      currerr=0.1 !0.1
      deltabdb=0.
      do k=1,ngena
         difus_type(k)="specify"
         difus_io(k)="disabled"
      enddo
      difus_io_file="drrin.nc"
      ndifus_io_t=0
      do ii=1,nbctimea
         difus_io_t(ii)=zero
         do k=1,ngena
            difus_io_drrscale(ii,k)=one
            difus_io_drscale(ii,k)=one
         enddo
      enddo
      difusr=1.d4
      difus_rshape(1)=1.0
      difus_rshape(2)=3.0
      difus_rshape(3)=3.0
      difus_rshape(4)=1.0
      difus_rshape(5)=-1.0
      difus_rshape(6)=0.0
      difus_rshape(7)=0.0
      difus_rshape(8)=0.0
      difus_vshape(1)=0.
      difus_vshape(2)=0.
      difus_vshape(3)=0.
      difus_vshape(4)=0.
      droptol=0.001d0
      dtr=5.d0
      dtr0=dtr
      do 3 i=1,ndtr1a
         dtr1(i)=0.
 3    continue
      esink=0.
      efflag='toroidal'
      efiter="enabled"
      efswtch="method1"
      efswtchn="disabled"
      efrelax=0.5
      efrelax1=0.5 !0.8
      elpar0=0.
      eoved=-.01
      enorm=200.
      enorme=enorm
      enormi=enorm
      epsthet=0.1
      eqmodel="power"
      eseswtch="disabled"
      f4d_out="disabled"
      gamaset=0.
      gamafac="disabled"
      gsla=270.
      ephicc=0.
      gslb=35.
      jhirsh=0
      kfrsou=0
      lbdry0="enabled"
      sbdry="bounded"
      scheck="enabled"
      iactst="disabled"
      idrop="No-Op"
      idskf="disabled"
      idskrf="disabled"
      ichkpnt="disabled"
      izeff="backgrnd"
      implct="enabled"
      ineg="disabled"
      isigtst=1
      do 6 i=1,6
        isigmas(i)=0
 6    continue
      isigsgv1=0
      isigsgv2=0
      kenorm=1
      lfil=30
      lmidpln=1
      lmidvel=0
      laddbnd=1
      colmodl=1
      locquas="disabled"
      machine="toroidal"
      manymat="disabled"
      meshy="free"
      nummods=1
!     if numixts= 1=>forw/back; -1=>back/forw for numindx=2
      numixts=1
      nchec=1
!     if nchgdy=1 adapt dy(ith)
      nchgdy=0
      ndeltarho="disabled"
      negyrg=0
      netcdfnm="disabled"
      netcdfshort="disabled"
      netcdfvecal="disabled"
      netcdfvecc="disabled"
      netcdfvece="disabled"
      netcdfvecrf="disabled"
      netcdfvecs="all"
      ncoef=1
      ncont=25
      nrstrt=1
!     if ngauss.ge.1 => analegco=disabled
!     good numbers are nlagran=4 and ngauss=4 or 6
!     max. nlagran allowed: 15
      nfpld=0
      ngauss=0
      nlagran=4
      nrf=0
      ngen=ngena
      nmax=nmaxa
      nonavgf=5
      nofavgf=10
      noncntrl=0
      nonel=0
      noffel=10000
      nonloss=0
      noffloss=10000
      nonvphi=10000
      noffvphi=10000
      nontran=0
      nofftran=10000
      nonelpr=10000
      noffelpr=0
!BH080305      do k=1,nmodsa
      do k=1,ngena
         nonrf(k)=0
         noffrf(k)=10000
      enddo
      nrskip=10
      numby=20
      do i=1,nplota
         nplot(i)=-10000
         nplt3d(i)=-10000
      enddo
      do i=1,nsavea
         nsave(i)=-10000
      enddo
      do 4 i=1,ndtr1a
         nondtr1(i)=-1
 4    continue
      npa_diag="disabled"
      atten_npa="enabled"
      nstop=5
      nstps=100
!     old way of integrating dens,cur in diaggnde
      oldiag="enabled"
      partner="disabled"
      profpsi="disabled"
      plt3d="enabled"
      pltd="enabled"
       !YuP[2018-02-07] New: pltd='color' and 'df_color'
       !for color contour plots
      pltdn="disabled"
      pltend="enabled"
      pltfvs="disabled"
      pltinput="enabled"
      pltlim="disabled"
      pltlimm=1.
      pltlos="disabled"
      pltso="disabled"
       !YuP[2018-02-07] New: pltso='color' and 'first_cl'
       !for color contour plots
      pltmag=1.
      pltsig="enabled"
      pltpowe="disabled"
      pltprpp="disabled"
      pltfofv="disabled"
      pltrst="enabled"
      pltstrm="disabled"
      pltflux="disabled"
      do 7 i=1,6
         pltflux1(i)=1.
 7    continue
      pltflux1(7)=0.
      pltvecal="disabled"
      pltvecc="disabled"
      pltvece="disabled"
      pltvecrf="disabled"
      plturfb="enabled" !YuP[2018-02-07] New: 'color' for color contour plots
      pltvflu="disabled"
      pltvs="rho"
      pltra="disabled"
      psimodel="axitorus"
      radcoord="sqtorflx"
      radmaj=100.
      eleccomp="enabled"
      radmin=50.
      relativ="enabled"
      relaxden=1.
      relaxtsp="disabled"
      rfacz=.7
      roveram=1.e-6
      rmirror=2.
      rzset="disabled"
      sigmamod="disabled"
      sigvcx=0.
      sigvi=0.
      simpbfac=1.d0
      symtrap="enabled"
      syncrad="disabled"
      softxry="disabled"
      soln_method="direct"
      tandem="disabled"
      taunew="disabled"
      tavg="disabled"
      do i=1,ntavga
         tavg1(i)=zero
         tavg2(i)=zero
      enddo
      tbnd(1)=.002
      do ll=2,lrorsa
         tbnd(ll)=0.
      enddo
      temp_den=0.0
      tfac=1.
      tfacz=1.
      thetd=0.0
      transp="disabled"
      pinch="disabled"
      advectr=1.0
      adimeth="disabled"
      nonadi= 5
      updown="symmetry"

      rdcmod="disabled"
      rdc_clipping="disabled"
      rdc_upar_sign=+1.
      nrdc=1
      rdcfile(1)="du0u0_input"
      nrdcspecies(1)=1
      rdcscale(1)=1.d0
      do i=2,nrdca
         rdcfile(i)="notset"
         nrdcspecies(i)=0
         rdcscale(i)=1.d0
      enddo

      urfmod="disabled"
      urfmult=1.0
      veclnth=1.0
      vdalp=.03
      vlfmod="disabled"
      vlfmodes=1.
      vlfnpvar="1/r"
      vlfbes="enabled"
      do k=1,nmodsa
         vlfharms(k)=1.
         dlndau(k)=1.
         vlfdnorm(k)=10.
         vlffreq(k)=.8e9
         vlfnp(k)=5.
         vlfdnp(k)=.2
         vlfddnp(k)=.1
         vlfnperp(k)=5.
         vlfharm1(k)=0.
         vlfeplus(k)=(0.,0.)
         vlfemin(k)=(0.,0.)
         vlfpol(k)=0.
         vlfdpol(k)=360.
         vlfddpol(k)=20.
         vlfparmn(k)=-ep100
         vlfparmx(k)=+ep100
         vlfprpmn(k)=0.0
         vlfprpmx(k)=+ep100
         vparmin(k)=-1.
         vparmax(k)=1.
         vprpmin(k)=0.
         vprpmax(k)=ep100
         vlh_karney=0.
         vlhpolmn(k)=0.
         vlhpolmx(k)=180.
         vlhprprp(k)="parallel"
      enddo
      vlhplse="disabled"
      vlhmod="disabled"
      vlhmodes=1.
      vlhpon=.1
      vlhpoff=.11
      vnorm=4.e10   !  Usually set through enorm
      vprprop="disabled"
      xfac=1.
      xlfac=1.
      xlpctlwr=.1
      xlpctmdl=.4
      xllwr=1./43.
      xlmdl=.25
      xpctlwr=.1
      xpctmdl=.4
      xlwr=1./43.
      xsinkm=1.
      xmdl=.25
      xprpmax=1.
      ylower=1.22
      yupper=1.275
      yreset="disabled"
      npwrzeff=1.
      mpwrzeff=1.
      npwrvphi=2.
      mpwrvphi=0.
      npwrelec=1.
      mpwrelec=1.
      npwrxj=1.
      mpwrxj=1.
      urfrstrt="disabled"
      urfwrray="disabled"
      xsink=0.

!.......................................................................
!     lrza arrays
!.......................................................................
      drya=1.d0/DBLE(lrza)
      do 100 ll=1,lrza
        rovera(ll)=.1
        rya(ll)=(ll-0.5)*drya
        elecfld(ll)=0.0
        zmax(ll)=1000.
 100  continue
      elecfld(0)=0.0
      zmax(0)=1000.
      rya(0)=0.
      rya(lrza+1)=1.

!.......................................................................
!     lrorsa arrays
!.......................................................................
      do 105 ll=1,lrorsa
!BH080122         irzplt(ll)=ll
         irzplt(ll)=0
 105  continue

!.......................................................................
!     nva arrays
!.......................................................................
      do nn=1,nva
         thet1(nn)=90.
         thet2(nn)=180.
         thet1_npa(nn)=90.
         thet2_npa(nn)=180.
         rd(nn)=zero
         thetd(nn)=zero
         x_sxr(nn)=zero
         z_sxr(nn)=zero
         rd_npa(nn)=100.d0
         thetd_npa(nn)=zero
         x_npa(nn)=zero
         z_npa(nn)=zero
      enddo
         rd(1)=100.d0


!DIR$ NEXTSCALAR
      mpwrsou(0)=1.
      npwrsou(0)=2.
      do 11 k=1,ngena
         torloss(k)="disabled"
         lbdry(k)="conserv"
         lossfile(k)="./prompt_loss.txt"
         lossmode(k)="disabled"
         enloss(k)=200.
         tauegy(k,0)=-1.
         negy(k)=0.0
         megy(k)=0.0
         regy(k)="disabled"
         gamegy(k)=0.0
         eparc(k,0)=-1.
         eperc(k,0)=-1.
         mpwrsou(k)=3.
         npwrsou(k)=2.
         ntorloss(k)=0.0
         mtorloss(k)=0.0
         paregy(k)=0.0
         peregy(k)=0.0
         pegy(k)=0.0
         tauloss(1,k)=0.
         tauloss(2,k)=0.
         tauloss(3,k)=0.
         do 8 i=1,6
            fpld(i,k)=0.
 8       continue
         fpld(7,k)=0.
         fpld(8,k)=1.e10
         fpld(9,k)=0.
         fpld(10,k)=pi
         do ll=1,lrza
            tauegy(k,ll)=-1.
            eparc(k,ll)=-1.
            eperc(k,ll)=-1.
         enddo
 11   continue
      do i1=1,negyrga
         do i2=1,2
            do k=1,ngena
               do ll=1,lrza
                  eegy(i1,i2,k,ll)=zero
                  jegy(i1,i2,k,ll)=0
               enddo
            enddo
         enddo
      enddo
      mpwr(0)=1.
      npwr(0)=2.
      do 10 k=1,ngen+nmax
        mpwr(k)=3.
        npwr(k)=2.
 10   continue
      do 9 k=1,ngen+nmax
        kpress(k)="enabled"
        kfield(k)="enabled"
 9    continue
      qsineut="disabled"
      trapmod="disabled"
      trapredc=0.
      scatmod="disabled"
      scatfrac=1.
      do 21 k=1,ngen+nmax
        reden(k,0)=1.
        denpar(k,0)=1.
        temp(k,0)=1.
        temppar(k,0)=1.
        reden(k,1)=1.
        denpar(k,1)=1.
        temp(k,1)=1.
        temppar(k,1)=1.
        bnumb(k)=1.
        fmass(k)=1.e-29
        nkconro(k)=0
        kspeci(1,k)=" "
        kspeci(2,k)=" "
        do ll=2,lrza
           reden(k,ll)=1.
           temp(k,ll)=1.
        enddo
        do ll=2,lza+1
          if (cqlpmod .ne. "enabled") then
            denpar(k,ll)=1.
            temppar(k,ll)=1.
          endif
        enddo
 21   continue
      nkconro(1)=1
      nkconro(2)=2
      nnspec=1
!..................................................................
!     Profile options are "parabola", "splines", and "asdex":
!..................................................................

      iprone="parabola"
      iprote="parabola"
      iproti="parabola"
      iprozeff="disabled"
      iprovphi="disabled"
      iproelec="parabola"
      ipronn="disabled"
      iprocur="parabola"
      tmdmeth="method1"
      njene=0
      njte=njene
      njti=njene
      enescal=1.
      tescal=1.
      tiscal=1.
      zeffscal=1.
      elecscal=1.
      vphiscal=1.
      bctimescal=1.d0

!..................................................................
!     acoef's specify ASDEX exponentail profiles
!     (ti profiles given by te profile, for "asdex" option).
!..................................................................

      acoefne(1)=-1.87
      acoefne(2)=-0.57
      acoefne(3)=-24.78
      acoefne(4)=-181.38
      acoefte(1)=7.51
      acoefte(2)=-13.45
      acoefte(3)=6.21
      acoefte(4)=-125.64

      zeffin(0)=1.0d0
      vphiplin(0)=0.0d0
      do 12  i=1,njenea
        ryain(i)=0.0d0
        elecin(i)=0.0d0
        tein(i)=0.0d0
        tiin(i)=0.0d0
        zeffin(i)=1.0d0
        vphiplin(i)=0.0d0
        difin(i)=0.0d0
 12   continue

      do k=1,npaproca
         npa_process(k)="notset"
         ennl(k)=5.0
         ennb(k)=1.e10
         do  i=1,njenea
            ennin(i,k)=0.0d0
         enddo
         ennscal(k)=1.0
      enddo
      npa_process(1)="cxh"

      do 13 k=1,ntotala
        do 14 i=1,njenea
          enein(i,k)=0.0d0
 14     continue
 13   continue

!..................................................................
!     Time-dependent profile quantities
!..................................................................

      nbctime=0
      do 30 i=1,nbctimea
         do 31 k=1,ntotala
            redenc(i,k)=0.d0
            redenb(i,k)=0.d0
            tempc(i,k)=0.d0
            tempb(i,k)=0.d0
 31      continue
         bctime(i)=DBLE(i-1)
         zeffc(i)=0.d0
         zeffb(i)=0.d0
         elecc(i)=0.d0
         elecb(i)=0.d0
         vphic(i)=0.d0
         vphib(i)=0.d0
         xjc(i)=0.d0
         xjb(i)=0.d0
         totcrt(i)=0.d0
 30   continue

      do 32 i=1,nbctimea
         do 33 l=1,njenea
            do 34 k=1,ntotala
               enein_t(l,k,i)=0.0d0
 34         continue
            tein_t(l,i)=0.0d0
            tiin_t(l,i)=0.0d0
            zeffin_t(l,i)=0.0d0
            elecin_t(l,i)=0.0d0
            xjin_t(l,i)=0.0d0
            vphiplin_t(l,i)=0.0d0
 33      continue
 32   continue

      ! For calc. of 1st-order orbit shift
      nr_delta=65
      nz_delta=65
      nt_delta=80 !Needs to be even
      ! For saving f4d== f(R,Z,u,theta) distribution:
      nr_f4d=20
      nz_f4d=21
      nv_f4d=20
      nt_f4d=20

      nen=nena
      nen_npa=nena
      mmsv=mx
      msxr=mx
      npaproc=1
      nv=1
      nv_npa=1
      enmin=5.
      enmax=50.
      fds=0.2
      enmin_npa=5.
      enmax_npa=50.
      fds_npa=0.2
      soucoord="disabled"
      nsou=1
      pltso="enabled"
      flemodel="fsa"
      knockon="disabled"
      komodel="mr"
      nkorfn=1
      nonko=10000
      noffko=10000
      soffvte=3.
      soffpr=0.5
      do 19 k=1,ngena
        do 20 m=1,nsoa
          nonso(k,m)=100000
          noffso(k,m)=100000
          do ll=1,lrza
             asor(k,m,ll)=0.
          enddo
          sellm1(k,m)=1.
          sellm2(k,m)=1.
          seppm1(k,m)=0.
          seppm2(k,m)=1.
          sem1(k,m)=0.
          sem2(k,m)=0.
          sthm1(k,m)=0.
          scm2(k,m)=0.
          szm1(k,m)=0.
          szm2(k,m)=1.
 20     continue
 19   continue

!     Some specific settings from cqlinput_help


!BH080125  DON'T reset this, as it conflicts with past
!BH080125  usage of asor.
!BH080125      do ll=1,lrza
!BH080125         asor(1,1,ll)=.25e+13
!BH080125         asor(1,2,ll)=3.25e+13
!BH080125      enddo

      scm2(1,1)=.001
      scm2(1,2)=10000.
      sellm1(1,1)=1.
      sellm1(1,2)=1.
      sellm2(1,1)=1.
      sellm2(1,2)=1.
      sem1(1,1)=1600.
      sem1(1,2)=0.
      sem2(1,1)=.5
      sem2(1,2)=25.
      seppm1(1,1)=1.
      seppm1(1,2)=1.
      seppm2(1,1)=1.
      seppm2(1,2)=1.
      sthm1(1,1)=5.
      sthm1(1,2)=0.
      szm1(1,1)=0.
      szm1(1,2)=0.
      szm2(1,1)=1.e+5
      szm2(1,2)=1.e+5

      do k=1,ngena
         do m=1,nsoa
            do ll=0,lrza
               sellm1z(k,m,ll)=sellm1(k,m)
               sellm2z(k,m,ll)=sellm2(k,m)
               seppm1z(k,m,ll)=seppm1(k,m)
               sem1z(k,m,ll)=sem1(k,m)
               sem2z(k,m,ll)=sem2(k,m)
               sthm1z(k,m,ll)=sthm1(k,m)
               scm2z(k,m,ll)=scm2(k,m)
               szm1z(k,m,ll)=szm1(k,m)
               seppm2z(k,m,ll)=seppm2(k,m)
               szm2z(k,m,ll)=szm2(k,m)
            enddo
            asorz(k,m,0)=0.
            do ll=1,lrza
               asorz(k,m,ll)=asor(k,m,ll)
            enddo
         enddo
      enddo



!.......................................................................
!l    4. Output option arrays
!.......................................................................

      do 400 i=1,noutpta
        nlotp1(i)=.false.
        nlotp2(i)=.false.
        nlotp3(i)=.false.
        nlotp4(i)=.false.
 400  continue
      nlotp1(4)=.true.
!BH070414      nlotp1(4)=.true.

!.......................................................................
!     5. Others, sometimes initialized later, but better do it before
!     reading namelist
!.......................................................................


!$$$c.......................................................................
!$$$c     6. Default values for finite orbit width (FOW) calculations
!$$$c.......................................................................
      fow="disabled" ! "disabled" is to use ZOW model
                     ! as the main model in CQL3D
!$$$      outorb="Not-detailed" ! outorb='detailed' or 'Not-detailed'
!$$$                ! (saving/not-saving data to a file for plotting)
!$$$      nmu =100  ! grid sizes for ad.ivariant mu
!$$$      npfi=100  ! and canonical momentum Pfi;
!$$$                ! to setup COM->R lookup table.
!$$$      nsteps_orb=50000 ! Max.number of time steps for orbit integration.
!$$$                ! Also used to trace Pfi=const levels for COM->R table
!$$$                ! in order to find intersections with mu=const levels.
!$$$      nptsorb=1 ! Number of points on a complete orbit
!$$$                ! (ityp=0 "main" orbit)
!$$$                ! from which ityp=1 "secondary" orbits are launched.
!$$$                ! ityp=1 orbit is stopped when it reaches the midplane.
!$$$                ! (Note: secondary orbits are not traced usually,
!$$$                ! see below, iorb2=0)
!$$$      i_orb_width=1 ! 1 -> Normal finite-orbit-width calculations.
!$$$                    ! 0 -> V_drift_perp is set to 0
!$$$                    ! (to mimic ZOW approximation)
!$$$      iorb2= 0 ! set to 1 to perform Runge-Kutta integration for tracing
!$$$               ! SECONDARY orbits to midplane; 0 - no RK tracing.
!$$$               ! This option (1) can be used for plotting orbits
!$$$               ! (together with outorb="detailed"),
!$$$               ! otherwise not needed.
!$$$      !=> MAIN orbits:
!$$$      ! u-grid:
!$$$      j0_ini= 20 !20 !2 !max(2, inc_j0/2)
!$$$      j0_end= 20 !20 !40
!$$$      inc_j0= 2  ! increment
!$$$      ! pitch-angle grid:
!$$$      i0_ini= 2  ! select as max(1, (inc_i0+1)/2)
!$$$      i0_end= iy-1
!$$$      inc_i0= 4  ! increment
!$$$      !=> SECONDARY orbits (launched from selected points on MAIN orbit)
!$$$      ! u-grid:
!$$$      j2_ini= 20 !20 !2  !max(2, inc_j2/2)
!$$$      j2_end= 20 !20 !40
!$$$      inc_j2= 2  !1 !2  ! increment
!$$$      ! pitch-angle grid:
!$$$      i2_ini= 1 !max(1, (inc_i2+1)/2) !1
!$$$      i2_end= iy-1
!$$$      inc_i2= 10 !4  ! increment
      return
      end
end module aindflt_mod
