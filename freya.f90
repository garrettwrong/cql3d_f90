!
!

!
!..................................................................
!     ONETWO DIVERGENCE
!     THE CALLING SEQUENCE TO FREYA IS ALTERED. CIVIC WILL ONLY ALLOW
!     61 ELEMENTS IN CALLING SEQUENCE. THEREFORE THE COMMUNICATION WILL
!     BE THROUGH FRCOMM
!..................................................................

      subroutine freya(ipts,mi,mj,codeid,rin,rmax,zax,zmin,zmax)
      use aminmx_mod, only : aminmx
      use bcast_mod, only : bcast
      use param_mod
      implicit none
      save
!---------------------------------------------------------------------
!
!     this subroutine calculates the particle and energy sources
!     due to neutral beam injection.
!
!---------------------------------------------------------------------
!     BH130329:  Adding capability to read a set of birth point files
!                from NUBEAM, to compare with results from freya
!                calculation.
!                New input variables in frsetup namelist are:
!                read_birth_pts="enabled", for use of this capability
!                birth_pts_files=input nbirth_pts_files filenames into
!                               the namelist, in the order which they
!                               will be used in cql3d. character*256.
!                               Max list length=24.  default='notset'
!                nbirth_pts_files=number of birth point files.
!                nbirth_pts= should be same as in all the NUBEAM files
!
!                A single beam is assumed, with up to three beam energy
!                components.
!
!                The shine through, and total power is read from the
!                files. nbirth_points is checked.  Power and shine
!                through numbers may vary from file to file.
!
!
!---------------------------------------------------------------------
!     the input quantities to freya are:
!
!     mb              Number of neutral beam injectors (.le.kb)
!     anglev(ib)      Vertical angle (degrees) between optical axis
!     and horizontal plane; a positive value indicates
!     particles move upward
!     angleh(ib)      Horizontal angle (degrees) between optical axis and
!     vertical plane passing through pivot point and
!     toroidal axis; a zero value denotes perpendicular
!     injection, while a positive value indicates par-
!     ticles move in the co-current direction
!     nsourc          Number of sources per beamline.
!     If 1, source is centered on beamline axis.
!     If nsourc=2, distinguish between the beamline
!     axis and the source centerline (optical axis).
!     The two sources are assumed to be mirror images
!     through the beamline axis.
!     In either case, the exit grid plane is perpendicula
!     to the beamline axis, and contains the source
!     exit grid center(s).
!     If nsourc=2, the alignment of the sources w.r.t.
!     the beamline axis is specified through bhofset,
!     bvofset, and bleni (described further below).
!     bvofset(ib)     Vertical offset from beamline axis to center
!     of each source (cm; used only for nsourc=2)
!     bhofset(ib)     Horizontal offset from beamline axis to center
!     of each source (cm; used only for nsourc=2)
!     bleni(ib)       Length along source centerline (source optical axis)
!     source to intersection point with the beamline axis
!     sfrac1(ib)      Fraction of source current per beamline coming
!     from upper source (used only for nsourc=2)
!     bcur(ib)        Total current (a) in ion beam (used only if bptor
!     is zero)
!     bptor(ib)       Total power (w) through aperture into torus; when
!     nonzero, bptor takes precedence over bcur
!     bshape(ib)      Beam shape
!     'circ' : circular
!     'rect' : rectangular
!     bheigh(ib)      Height of source (cm)
!     bwidth(ib)      Width of source (cm); diameter for
!     circular source.
!     bhfoc(ib)       Horizontal focal length of source (cm)
!     bvfoc(ib)       Vertical focal length of source (cm)
!     bhdiv(ib)       Horizontal divergence of source (degrees)
!     bvdiv(ib)       Vertical divergence of source (degrees)
!     ebkev(ib)       Maximum particle energy in source (keV)
!     fbcur(ie,ib)    Fraction of current at energy ebkeV/ie
!     iborb           Flag for modeling orbit effects on beam-
!     generated fast ions
!     1, model orbit effects
!     0, do not model orbit effects
!     npart           Number of particles followed into plasma
!     (suggest 10000)
!     npskip          Ratio of number of particles followed into plasma
!     to number of source particles (suggest 1)
!     naptr           Total number of aperatures encountered by a particle
!     as is moves from the source into the plasma chamber
!     Maximum is specified by parameter nap (=10).
!     First set of aperatures encountered by the particle
!     are assumed centered on the source axis, and subseq
!     aperatures are centered on the beamline axis; the
!     distinction is made through ashape.
!     ashape(iap,ib)  Aperture shape.
!     Prefix 's-' indicates source axis centered.
!     Prefix 'b-' indicates beamline axis centered.
!     's-circ'          'b-circ'
!     's-rect'          'b-rect'
!     's-vert'          'b-vert'
!     's-horiz'         'b-horiz'
!     'b-d3d'
!     (circ=circular aperature, rect=rectagular,
!     vert=limits vertical height of source particles,
!     horiz=limits horizontal height of source particles,
!     d3d= special DIII-D polygonal aperature)
!     'circ' : circular
!     'rect' : rectangular
!     aheigh(iap,ib)  Height of aperture (cm)
!     awidth(iap,ib)  Width of aperture (cm); diameter for circular
!     aperture
!     alen(iap,ib)    Length from source to aperature for 's-type' aperatur
!     and from exit grid plane along beamline axis for
!     'b-type' aperatures.
!     blenp(ib)       Distance along beamline axis from source exit
!     plane to the fiducial "pivot" point.
!     rpivot(ib)      Radial position of pivot (cm)
!     zpivot(ib)      Axial position of pivot (cm)
!     atw(k)          atomic mass of ion species k
!     codeid          flux surface geometry
!     "onedee"      : elliptical
!     anything else : nonelliptical
!     elong           elongation (height/width) of elliptical
!     cross-section plasma
!     ibion           ion index of beam species
!     mf              number of flux zones
!     mi              number of radial mesh points for nonelliptical
!     plasma
!     mj              number of axial mesh points for nonelliptical
!     plasma
!     nion            number of ion species
!     norb            flag controlling output of orbit information
!     >0, write out extensive orbit information
!     to unit norb
!     0, omit output of orbit information
!     <0, write out abbreviated orbit information
!     to unit iabs(norb)
!     pinsid(i)       poloidal magnetic flux (G-cm2) along horizontal
!     chord inside and through the magnetic axis vs.
!     a uniform mesh in major radius; pinsid is needed
!     only if iborb=1
!     potsid(i)       poloidal magnetic flux (G-cm2) along horizontal
!     chord outside and through the magnetic axis vs.
!     a uniform mesh in major radius; potsid(1) and
!     potsid(mf) are needed if codeid.ne.'onedee';
!     the entire potsid array is needed if iborb=1
!     psi(i,j)        poloidal magnetic flux (G-cm2) at mesh point i,j
!     (needed only if codid.ne.'onedee')
!     psivol(i)       volume (cm3) of flux zone i; depending upon
!     whether codeid.eq.'onedee' or not, psivol is
!     chosen such that either r or sqrt(psi-psiax)
!     varies a constant amount from one flux surface
!     to the next
!     rinsid(i)       major radius (cm) along horizontal chord inside
!     and through the magnetic axis vs. a uniform mesh
!     in sqrt(psi-psiax); rinsid is needed only if
!     iborb=1
!     rotsid(i)       major radius (cm) along horizontal chord outside
!     and through the magnetic axis vs. a uniform mesh
!     in sqrt(psi-psiax); rotsid(1) and rotsid(mf) are
!     needed if codeid.eq.'onedee'; the entire rotsid
!     array is needed if iborb=1
!     r(i)            maj radius mesh point i for nonelliptical plasma (cm)
!     rin             major radius of inside of vacuum vessel (cm)
!     rmax            maximum radial position of plasma (cm)
!     rpivot(ib)      radial position of pivot point (cm)
!     sfrac1(ib)      fraction of source current per beamline coming
!     from upper source (used only for nsourc=2)
!     sofset(ib)      vertical offset from optical axis to center
!     of each source (cm; used only for nsourc=2)
!     ONETWO DIVERGENCE: use w instead of z
!     w(j)            axial mesh point j for nonelliptical plasma (cm)
!     zax             axial position of magnetic axis (cm) for
!     elliptical plasma
!     zmin            minimum axial position of plasma (cm)
!     zmax            maximum axial position of plasma (cm)
!     zpivot(ib)      axial position of pivot point (cm)
!     zne(i)          electron density in zone i (cm-3)
!     zni(i,k)        density of ion species k in zone i (cm-3)
!     zte(i)          electron temperature in zone i (kev)
!     zzi(i,k)        charge number of ion species k in zone i
!     ONETWO DIVERGENCE
!     zeffctv(i)      z-effective (used with iexcit=-1)
!     iexcit       switch that controls the following
!     ONETWO DIVERGENCE
!     iexcit=-1 uses Stearn's average cross sections:otw. like =0
!     iexcit = 0   (default mode) normal freya run
!     iexcit = 1   use hexnb instead of crsecs to get cross sections
!     but neglect presence of excited sates in beam
!     iexcit = 2   use hexnb with excited beam states allowed.
!     inubpat     switch controlling 2-d beam deposition calculation
!     0, (default) bypass 2-d deposition
!     1, determine beam deposition on standard (r,w) grid
!     and output deposition of 3rd excited state
!     components to file 'beamdep'
!     2, same as inubpat=1, except bin deposition on
!     modified (r,w) grid whose dimensions are given by
!     'npat' (see below).
!     note:  if inubpat.gt.0, iexcit must be set equal to 2.
!     this will automatically be checked, iexcit rese
!     (if needed), and a message output to the crt.
!     npat        modified (r,w) dimension array for optional use if
!     inubpat.gt.0, containing:
!     npat(1)=number of new 'r' elements
!     npat(2)=number of new 'w' elements
!     default:  npat(1)=mi, npat(2)=mj
!
!     the output quantities are:
!
!     bion(ie,ib)     intensity of ion beam (particles/s)
!                     [BH: evidently out of the ion source, before
!                      neutralization]
!     bneut(ie,ib)    intensity of neutral beam (particles/s)
!                     [BH: Coming out of the neutralizer,
!                      a fraction of bion]
!     bpow(ie,ib)     beam power to aperture (w)
!     ebeam(ie,ib)    particle energy (kev)
!     fap(ie,ib)      fraction of beam stopped by aperture
!     fwall(ie,ib)    fraction of beam incident on wall (shinethrough)
!     forb(ie,ib)     fraction of beam lost on orbits
!     ftrapfi(i,ie,ib)fraction of trapped fast ions in each zone
!     fb11(ie,ib)     fraction of ions passing and axis-encircling
!     fb10(ie,ib)     fraction of ions passing and not encircling
!     fb01(ie,ib)     fraction of ions trapped and axis-encircling
!     fb00(ie,ib)     fraction of ions trapped and not encircling
!     fber(ie,ib)     fraction of ions trapped for which error was detec
!     hibrz(i,ie,ib)  normalized hot ion birth rate
!     hdepz(i,ie,ib)  normalized hot ion deposition rate
!     ipts            number of birth points to be plotted
!     xpts(ii)          x coordinate of birth point
!     ypts(ii)          y coordinate of birth point
!     zpts(ii)          w coordinate of birth point
!     vx(ii)          x component of birth velocity
!     vy(ii)          y component of birth velocity
!     vz(ii)          w component of birth velocity
!     wb11(ie,ib)     orbit width of fb11 ions (cm)
!     wb10(ie,ib)     orbit width of fb10 ions (cm)
!     wb01(ie,ib)     orbit width of fb01 ions (cm)
!     wb00(ie,ib)     orbit width of fb00 ions (cm)
!     zetaz(i,ie,ib)  average pitch angle cosine of deposited hot ions
!     angmpz(i,ie,ib) toroidal angular momentum density deposition rate.
!     (i.e. ang. momtm. deposited in shell per second /volume of she
!     angmpz contains only the toroidal component of ang momentum.
!---------------------------------------------------------------------
!     ONETWO DIVERGENCE
      include 'frcomm.h'

!     Automatic arrays for local dynamic mem, for case where NUBEAM
!     particle birth point list is used. nbirth_pts from frcomm.h
      integer, intent(out) :: ipts
      integer mi,mj ! in arg.list
      real(c_double) :: rin,rmax,zax,zmin,zmax ! in arg.list
      integer i,j,k,ib ! local

      real(c_double), allocatable :: x_nub(:),y_nub(:),z_nub(:)
      real(c_double), allocatable :: vx_nub(:),vy_nub(:),vz_nub(:)
      real(c_double), allocatable :: v_nub(:),en_nub(:)
      integer, allocatable :: ie_nub(:)  !energy cmpt,=1, 2, or 3
      character*128 filenm
      real(c_double), dimension(3) :: en_avg_cmpts_nub,pabs_cmpts_nub
      integer, dimension(3) :: nbirth_cmpts_nub
      real(c_double), dimension(3) :: birth_rate_cmpts_nub

      character*8 codeid
      real(c_double) :: psi(ki,kj),r(ki),w(kj)
      real(c_double) :: bpow(ke,kb),eb(ke,kb)
      equivalence(bpow,pbeam)
      equivalence(eb,ebeam)
      equivalence(r,xxx)
      equivalence(w,yyy)
      equivalence(psi,p)
      real(c_double) :: cangv(kb), cangh(kb), sangv(kb), sangh(kb)
      real(c_double) :: thetp(kb), thetpp(kb)
      real(c_double) :: costp(kb), sintp(kb), costpp(kb), sintpp(kb)
      real(c_double) :: vbeam(ke,kb)
      integer iatype(nap,kb)
!      dimension sgvxne(kz),sgxn(kz,ke,kb),sgxnmi(ke,kb),hxfrac(ke,kb)
      real(c_double) :: sgvxne(kz),sgxn(kcmp1,kz,kbe,ksge)
      real(c_double) :: sgxnmi(ke,kb),hxfrac(ke,kb)
!      real(c_double) :: sgxn1(kcmp1,kz,kbe,ksge)
      real(c_double) :: wt(kz), zeta_(kz),angmtp(kz)
      integer nmbrz(kz)
      real(c_double) :: rzpat(kix2,kjx2,ke,kb)
      real(c_double) pio180,dr,dz,dri,dzi,elong,pzone,rzone !local
      real(c_double) rmajor,drot,drin,drutp,deps,droti,drini !local
      real(c_double) drutpi,elongi,bntot,beff,ebev
      integer ncalls,ie,mb,norb,mim1,mjm1,kmin,kmax !local
      real(c_double) :: drpat,dzpat ! arg.in setrz()
      integer nrpat,nzpat  ! arg.in setrz()
      real(c_double) dnemin,dnemax,dtemin,dtemax
      real(c_double) dzemin,dzemax
      real(c_double) pinj_nub,pabs_nub ! arg. in read_nubeam_data()
      integer nshine_nub
      integer ipar,npar,nparx,newpar,ic,ii,iskip,isourc
      real(c_double) x0,y0,z0,vx0,vy0,vz0,xpos,ypos,zpos,rpos
      real(c_double) atwb  ! in arg. of inject1()
      integer myid ! in arg. of inject1()
      integer maxp,mlost,izone
      real(c_double) ebx,csgn,volume,tenter,smax,texit,zetai,vplane
      real(c_double) vtroid,vrad,angmtm
      real(c_double) wid,xnorm,xloss1,xloss2,bptorx
      integer ipass,iaxis,ier,izp
      
      data pio180 /0.017453293d0/
      data ncalls /0/
!
!.. new (JK)
!   kimp-> 2 for namei
!
      
      real(c_double) :: de_tk
      real(c_double) :: sgxnloc(kbe)
      real(c_double) :: zangrot(kz), znis(kz,kion), zzis(kz,kion)
!
      if (read_birth_pts.eq."enabled") then
      if (.NOT. ALLOCATED(x_nub)) then
      allocate(x_nub(nbirth_pts),y_nub(nbirth_pts),z_nub(nbirth_pts))
      allocate(vx_nub(nbirth_pts),vy_nub(nbirth_pts),vz_nub(nbirth_pts))
      allocate(v_nub(nbirth_pts),en_nub(nbirth_pts))
      allocate(ie_nub(nbirth_pts))
      endif
      endif

!----------------------------------------------------------------------
!     general freya initialization
      call bcast(rzpat,zero,kix2*kjx2*ke*kb)   !rzpat not used w cql3d
      mb=nbeams
      elong=0.d0
      norb=0
!
      do ib=1,ke
       do i=1,kb
         bneut(ib,i)=0.d0
         vbeam(ib,i)=0.d0
         sgxnmi(ib,i)=0.d0
         hxfrac(ib,i)=0.d0
       enddo
      enddo
      do i=1,kz
       do j=1,ke
         do k=1,kb
           hibrz(i,j,k)=0.d0
           hdepz(i,j,k)=0.d0
           ftrapfi(i,j,k)=0.d0
         enddo
       enddo
      enddo
      do i=1,kz
       do j=1,ke
         do k=1,kb
          sgxn(1,i,j,k)=0.d0
          sgxn(2,i,j,k)=0.d0
          sgxn(3,i,j,k)=0.d0
          sgxn(4,i,j,k)=0.d0
        enddo
       enddo
      enddo
      pzone=0.D0
      rzone=0.D0
!
!----------------------------------------------------------------------
!     set up some data for optional (r,w) deposition calculation
      inubpat=0  !Not used with cql3d
      if(inubpat.gt.0) then
        call setrz(npat,r(1),r(mi),w(1),w(mj),drpat,dzpat,nrpat,nzpat)
        iexcit=2
      endif

!----------------------------------------------------------------------
!     if nubeam list case, make sure only 1 beam
      if (read_birth_pts.eq."enabled") then
         nbeams=1
         mb=1
         if (codeid.eq.'onedee') then
            WRITE(*,*)'Freya: nubeam list not setup for this codeid'
            WRITE(*,*)'STOP'
            STOP
         endif
      endif

!
!     turn off orbit calculation in case of counter-injection; present
!     orbit model applies only for co-injection
      do 5 ib=1,mb
 5    if(angleh(ib).lt.5.) iborb = 0
!
!     initialize flux surface quantities
!     BH130327:
!          rotsid, rinsid equal 0, potsid(mf),potsid(1)set but other 0.
      mfm1 = mf-1
      rmajor = rotsid(1)
      drot   = (rotsid(mf)-rotsid(1))/mfm1
      drin   = (rinsid(1)-rinsid(mf))/mfm1
      drutp  = sqrt(potsid(mf)-potsid(1))/mfm1
!990131      drot   = amax1(drot,1.e-6)
!990131      drin   = amax1(drin,1.e-6)
!990131      drutp  = amax1(drutp,1.e-6)
!990131      elong  = amax1(elong,1.e-6)
      deps=1.d-6
      drot   = max(drot,deps)
      drin   = max(drin,deps)
      drutp  = max(drutp,deps)
      elong  = max(elong,deps)   !elong not properly set.
      droti  = one/drot
      drini  = one/drin
      drutpi = one/drutp
      elongi = one/elong
      if(codeid.ne.'onedee') then
        mim1  = mi-1
        mjm1  = mj-1
        dr    = (r(mi)-r(1))/mim1
        dz    = (w(mj)-w(1))/mjm1
        dri   = one/dr
        dzi   = one/dz
      endif
!
!     calculate sines and cosines of various angles
      do 15 ib=1,mb
        cangv(ib)=cos(anglev(ib)*pio180)
        cangh(ib)=cos(angleh(ib)*pio180)
        sangv(ib)=sin(anglev(ib)*pio180)
        sangh(ib)=sin(angleh(ib)*pio180)
 15   continue
!
!     calculate beam power; account for neutralizer efficiency
!
      atwb = atw(ibion)
      bntot = zero
      do 20 ib=1,mb
        do 21 ie=1,3
          ebeam(ie,ib) = ebkev(ib)/ie
          ebx = ebeam(ie,ib)/atwb
!     ONETWO DIVERGENCE   - delete reference to ionization efficiency
!     routines.
      call logint(ebx,beff)  ! JK  BH140501: Check if any diff.
!          beff=1.
          ebev = 1.e3*ebx
          vbeam(ie,ib) = 1.384e6*sqrt(ebev)  !adjstd for mass, cm/sec
          bion(ie,ib)  = 0.625e19*fbcur(ie,ib)*bcur(ib) ! molecular src
                                       !rate (as I understand, BH)
                                       !bcur need not be set
                                       !in nml, if bptor.ne.0.
                                       !In this case, bneut, bpow are
                                       !set for default bcur()=110.,
                                       !and renormalize for bptor below.
          bneut(ie,ib) = ie*beff*bion(ie,ib)  !ion rates at energy/ie.
          bpow(ie,ib)  = ebeam(ie,ib)*bneut(ie,ib)/0.625e16
!          write(*,*)'freya: ibion,atw(ibion),ie,ib,vbeam(ie,ib)',
!     &                      ibion,atw(ibion),ie,ib,vbeam(ie,ib)
 21     bntot = bntot + bneut(ie,ib)
 20   continue
!
!.......................................................................
!BH130329:  If read_birth_pts="enabled", then read in beam birth points
!           from a NUBEAM generated file, as described above.
!           There may be some unnecessary freya related calculation
!           since the NUBEAM data-read is shoe-horned in on top of freya
!
!           Else:  freya calculation.
!.......................................................................

      if (read_birth_pts.eq."enabled") then
!        read_nubeam_data reads data file and returns it
!        through the argument list.
!        NOTE:  The NUBEAM birth points are actually shifted from
!        the ion birth point, accounting for the distance to the
!        GC starting position of the particle.  That is, these
!        are birth points for particle guiding centers.
         ncalls=ncalls+1
         if (ncalls.gt.nbirth_pts_files) then
            write(*,*)'FREYA WARNING: Insufficient birth_pts_files'
            write(*,*)'FREYA WARNING: Stepping back to last file'
            ncalls=ncalls-1
         else
            write(*,*)'Freya: read_birth_pts case, ncalls=',ncalls
         endif
         filenm=trim(birth_pts_files(ncalls))
         call read_nubeam_data(filenm,nbirth_pts,atwb, &
                         nbirth_cmpts_nub,nshine_nub, &
                         pinj_nub,pabs_nub,x_nub,y_nub,z_nub, &
                         vx_nub,vy_nub,vz_nub,v_nub,en_nub,ie_nub, &
                         birth_rate_cmpts_nub,en_avg_cmpts_nub, &
                         pabs_cmpts_nub)

!        pinj_nub=Injected power (Watts) = freya bptor
         bptor(1)=pinj_nub
!        Set bcur,bion,bneut,bpow from Nubeam list
!        This is setup according to freya bptor, as above.
!        Here, bcur, here, refers to absorbed power, since
!        present Nubeam list doesn't contain a breakdown of the shine
!        through components.  [Fix later].
         bntot=0d0
         ib=1                   !Nubeam list only set up for 1 beam
         do ie=1,3
            bcur(ib)= bcur(ib)+birth_rate_cmpts_nub(ie)*1.6022e-19
            bneut(ie,ib)=birth_rate_cmpts_nub(ie)
            bion(ie,ib)=bneut(ie,ib)/ie
            bpow(ie,ib)=pabs_cmpts_nub(ie)
            bntot=bntot+bpow(ie,ib)  !Prob here, compared to freya
                                     !calc, since for freya should be
                                     !same units as bneut[see below use]
         enddo

      endif    ! On read_birth_pts.eq."enabled"

      if (read_birth_pts.ne."enabled") then  !Skip, if nubeam case

!     ONETWO DIVERGENCE
!..................................................................
!     Determine peak electron density, temp and zeff(lr_) for use with
!     iexcit=-1 (Stearn's formula)
!..................................................................

      call aminmx(zne,1,mfm1,1,dnemin,dnemax,kmin,kmax)
      call aminmx(zte,1,mfm1,1,dtemin,dtemax,kmin,kmax)
      call aminmx(zeffctv,1,mfm1,1,dzemin,dzemax,kmin,kmax)
      dtemax=dtemax/10.
!
!...calculate macroscopic cross sections
!     iexcit .le. 0 : Freeman and Jones
!                 5 : ADAS
!                 6 : Boley parameterization
!
!... new coding (JK)
!
      csgn=1.d0
      do i=1,kz
       zangrot(i)=0.d0
      enddo
!
!... In CQL3D, for nprim=1, nimp=1, main ions in zni(j,1), (j,2)
!    and zzi(j,1) and zzi(j,2) with general and maxwellian
!    distributions, respectively. The impurity ion is in
!    zni(j,3) and zzi(j,3). Also, Zeff has been copied into
!    zzi(j,nion+2=4). So, for this case we slip
!    the arrays such that (j,3) - > (j,2) for the impurity
!
      if(nprim.eq.1 .and. nimp.eq.1) then
       do j=1,kz
         znis(j,1)=zni(j,2)
         znis(j,2)=zni(j,3)
         zzis(j,1)=zzi(j,2)
         zzis(j,2)=zzi(j,3)
         zzis(j,3)=zzi(j,4)  ! Zeff profile
       enddo
      else
       do j=1,kz
         do i=1,kion
           znis(j,i)=zni(j,i)
           zzis(j,i)=zzi(j,i)
         enddo
       enddo
      endif
!
      write(*,*) 'calling nbsgxn ...'
      write(*,*) 'iexcit = ',iexcit
!       write(*,*) 'zte:'
!       write(*,'(5(2x,1pe14.5))') (zte(j),j=1,kz)
!       write(*,*) 'zti:'
!       write(*,'(5(2x,1pe14.5))') (zti(j),j=1,kz)
!       write(*,*) 'zne:'
!       write(*,'(5(2x,1pe14.5))') (zne(j),j=1,kz)
!       write(*,*) 'zni(j,1-3):'
!       write(*,'(5(2x,1pe14.5))') (znis(j,1),j=1,kz)
!       write(*,'(5(2x,1pe14.5))') (znis(j,2),j=1,kz)
!       write(*,'(5(2x,1pe14.5))') (znis(j,3),j=1,kz)
!       write(*,*) 'zzi(j,1-3):'
!       write(*,'(5(2x,1pe14.5))') (zzis(j,1),j=1,kz)
!       write(*,'(5(2x,1pe14.5))') (zzis(j,2),j=1,kz)
!       write(*,'(5(2x,1pe14.5))') (zzis(j,3),j=1,kz)
!       stop
!
        call nbsgxn (iexcit,namep,namei,mb,mfm1,ne_tk,nprim, &
                   nimp,nion,atwb,atw,ebkev,fe_tk,ibion,vbeam, &
                   zne,znis,zte,zti,zzis,de_tk,dtemax,dnemax,dzemax, &
                   hxfrac,sgxn,sgxnmi)
!      do ib=1,mb
!       do j=1,3
!         write(*,'(2i3,2x,1p1e12.4,a12)') j,ib,sgxnmi(j,ib),' sgxnmi'
!       enddo
!      enddo
!      stop
!
!... old coding (JK)
!
!      if(iexcit.le.0)then
!     ONETWO DIVERGENCE
!        call crsecs(atw,ebkev,ibion,ke,kz,mb,mfm1,nion,vbeam,
!     &    zne,zni,zte,zzi,sgvxne,sgxn,sgxnmi,iexcit,
!     &    dtemax,dnemax,dzemax)
!         do ib = 1, mb
!           do j = 1, 3
!             do i=1,mfm1
!        write(*,'(3i3,2x,1p1e10.4,a12)') i,j,ib,sgxn(i,j,ib), ' sgxn'
!             enddo
!           enddo
!         enddo
!         stop
!      else
!     ONETWO DIVERGENCE
!         write(*,*) 'freya: before frhexdrv'
!        call frhexdrv(mb,sgxn,sgxnmi,hxfrac)
!      endif
!... end old coding
!
      endif  !on read_birth_pts.ne."enabled"
!
!     calculate total plasma volume
      volume = 0.
      do i=1,mfm1
        write(*,*)'freya: i,psivol(i)=',i,psivol(i)
        volume = volume + psivol(i)
      enddo
      write(*,*)'freya: volume(sum of psivol)=',volume
      !Note: psivol is based on eqvol() which is calc-ed by subr.eqvolpsi.
      !Can be a little different from alternative definition
      ! of flux surface volume as setup in subr.tdrmshst.
!----------------------------------------------------------------------
!     begin loop over beams, ib
!----------------------------------------------------------------------

!     Following statment executed in subroutine frset.
!      if (read_birth_pts.eq."enabled") npart=nbirth_pts  !nubeam case

!BH110309      iskip =       1 + (npart-1)/1500
      maxp=1500000 !Formerly parameter giving max # of ions launched,
                  !according to code comment
      iskip =       1 + (npart-1)/maxp
      ic    = 0
      ipts  = 0

      do 200 ib=1,mb  ! 386 lines down to line 989

      if (read_birth_pts.ne."enabled") then  !Skip, if nubeam case
!
!     Determine aperature designators.
!     If nsourc=1, treat source axis centered aperatures equivalent
!     to beam axis centered.
!BH130915: Not sure of effect of this change. Reverting:do 40  i=1,naptr
         do 40  i=1,nap
            if (nsourc.eq.1) then
               if(ashape(i,ib).eq.'s-circ')  iatype(i,ib)=5
               if(ashape(i,ib).eq.'s-rect')  iatype(i,ib)=6
               if(ashape(i,ib).eq.'s-vert')  iatype(i,ib)=7
               if(ashape(i,ib).eq.'s-horiz') iatype(i,ib)=8
            elseif (nsourc.gt.1) then
               if(ashape(i,ib).eq.'s-circ')  iatype(i,ib)=1
               if(ashape(i,ib).eq.'s-rect')  iatype(i,ib)=2
               if(ashape(i,ib).eq.'s-vert')  iatype(i,ib)=3
               if(ashape(i,ib).eq.'s-horiz') iatype(i,ib)=4
            endif
            if(ashape(i,ib).eq.'b-circ')  iatype(i,ib)=5
            if(ashape(i,ib).eq.'b-rect')  iatype(i,ib)=6
            if(ashape(i,ib).eq.'b-vert')  iatype(i,ib)=7
            if(ashape(i,ib).eq.'b-horiz') iatype(i,ib)=8
            if(ashape(i,ib).eq.'b-d3d')   iatype(i,ib)=9
!        write(*,*) 'ashape(i,ib) = ',i,ib,ashape(i,ib)
!        write(*,*) 'iatype(i,ib) before rotate = ',iatype(i,ib)
 40      continue
!
!     Some angles for subroutine rotate
        thetp(ib)=atan2(bvofset(ib),sqrt(bleni(ib)**2-bvofset(ib)**2))
        costp(ib)=cos(thetp(ib))
        sintp(ib)=sin(thetp(ib))
        thetpp(ib)=atan2(bhofset(ib),sqrt(bleni(ib)**2-bvofset(ib)**2))
        costpp(ib)=cos(thetpp(ib))
        sintpp(ib)=sin(thetpp(ib))
!
        endif  !on read_birth_pts.ne."enabled"
!
!----------------------------------------------------------------------
!     begin loop over beam energy components, ie
!----------------------------------------------------------------------
        do 201 ie=1,3  !340 lines down
!         beam fractions
          fap(ie,ib)   = 0.
          fwall(ie,ib) = 0.
          forb(ie,ib)  = 0.
          fb11(ie,ib)  = 0.
          fb10(ie,ib)  = 0.
          fb01(ie,ib)  = 0.
          fb00(ie,ib)  = 0.
          fber(ie,ib)  = 0.
          do 110 i=1,mfm1
             ftrapfi(i,ie,ib)=0.
             hibrz(i,ie,ib) = 0.
             hdepz(i,ie,ib) = 0.
             angmpz(i,ie,ib) = 0.
             hicmz(i,ie,ib,1) = 0.
             hicmz(i,ie,ib,2) = 0.
             hicmz(i,ie,ib,3) = 0.
             zetaz(i,ie,ib) = 0.
 110      continue
          do 310 i=1,mfm1
             nmbrz(i)=0
 310      continue
!     --- nmbrz(i) counts ions born in zone i
!----------------------------------------------------------------------
!     begin loop over particles, for each beam(ib),energy(ie)
!----------------------------------------------------------------------
!
          if (read_birth_pts.ne."enabled") then
             npar = (bneut(ie,ib)/bntot)*npart   !Not nubeam list
          else
             npar=nbirth_cmpts_nub(ie)           !nubeam list
             ii=0
          endif
!
!BH130925          npar = (bneut(ie,ib)/bntot) * npart
          if(npar.eq.0) go to 201
          nparx = 0
          newpar = 0  ! JK
!
!.......................................................................
!     Loop over the particles, for each ie,ib
!.......................................................................
         do 180 ipar=1,npar  ! down 218 lines, to line 919
            if(mod(ipar-1,npskip).eq.0) newpar=1
! YuP[171103] Moved this line inside if()            if(newpar.eq.0) go to 120
!
          if (read_birth_pts.ne."enabled") then  !Skip if nubeam list
             ! FREYA normal calculations
             if(newpar.eq.0) go to 120 ! YuP[171103] moved from above
!
!... generate neutral particle at beam source
!
             call sorspt1(bshape,bheigh,bwidth,bhfoc,bvfoc,bhdiv, &
                  bvdiv,ib,ie,isourc,nsourc,sfrac1,vbeam,x0,y0,z0, &
                  vx0,vy0,vz0)
!
!... transform coordinates and advance particle to pivot point
!
            call rotate(naptr,iatype,aheigh,awidth,alen,bhofset, &
              bvofset,cangv,cangh,ib,isourc,costp,sintp,costpp, &
              sintpp,blenp,nsourc,sangv,sangh,rpivot,zpivot,mlost, &
              x0,y0,z0,vx0,vy0,vz0)
!      write(*,*)'freya,rotate:x0,y0,z0,vx0,vy0,vz0',x0,y0,z0,vx0,vy0,vz0
!
!     skip injection if particle is lost at aperture
!Check inj  write(*,*)'freya: ie,ipar,npar,mlost',ie,ipar,npar,mlost
! 120        if(mlost.ne.0) go to 160
  120       continue
            if (mlost .ne. 0)  go to 160
!
!... inject particle into plasma, i.e., follow particle from pivot
!     point into, through, or around plasma. Use old subroutine
!     inject for crsecs.
!
!       open(42,file='psi.dat',status='unknown')
!       write(42,*) 'psi data for DIII-D #122080 from cql3d'
!       write(42,*) '(nw,nh) = ',65,65
!       do i=1,65
!        do j=1,65
!          write(42,'(2i4,2x,1p1e15.7,a16)') i,j,psi(i,j),' psi'
!        enddo
!       enddo
!       close(42)
!       stop
!
            if (iexcit.le.0) then
             call inject_old(atw,codeid,drutpi,droti,dri,dzi, &
                   elongi,ib,ie,mfm1,mim1,mjm1,newpar,potsid(1), &
                   psi,r,rmajor,rin,rmax,sgxn,sgxnmi,x0,y0,z0, &
                   vx0,vy0,vz0,vbeam,w,zax,zmin,zmax,izone, &
                   pzone,rzone,rpos,xpos,ypos,zpos)
            else
             call inject1(atwb,codeid,de_tk,drutpi,droti,dri,ds_tk,dzi, &
                   elongi,ib,ie,kb,kbe,ksge,ke,kz,ki,mfm1,mim1, &
                   mjm1,ne_tk,newpar,nout,potsid(1),psi,r,rmajor,rin, &
                   rmax,sgxn,sgxnloc,sgxnmi,x0,y0,z0,vx0,vy0,vz0, &
                   vbeam,w,zangrot,zax,zmin,zmax,izone, &
                   pzone,rzone,rpos,xpos,ypos,zpos,myid,tenter, &
                   smax,texit)
            endif

!           Shift from particle birth point to guiding center point


          else    !On read_birth_pts.ne."enabled"  ! NUBEAM list

!           Nubeam list case.
!           Remember, for given ie (energy cmpt),
!           need to skip to next particle for the given ie.
            do
              ii=ii+1
              if (ie_nub(ii).eq.ie) exit
              if (ii.gt.npart) then
                WRITE(*,*)'Freya: Nubeam list inconsistency'
                STOP
              endif
            enddo
!
            xpos=x_nub(ii)
            ypos=y_nub(ii)
            zpos=z_nub(ii)
            vx0=vx_nub(ii)
            vy0=vy_nub(ii)
            vz0=vz_nub(ii)
            rpos=sqrt(xpos**2+zpos**2)
!           get pzone (i.e., psi-value) and radial izone
            call zone(drutpi,ki,mfm1, mim1, mjm1,dri,dzi,potsid(1),psi, &
                 r,w,xpos,ypos,zpos,pzone,izone)
!
            vbeam(ie,ib)=v_nub(ii)
!
          endif                  !On read_birth_pts.ne."enabled"
!
!... skip birth data if:  particle missed plasma
!
            if(izone.ge.mf) go to 170
!     For removing NBI source at all psi outside of psicutoff:
          !print*,'Rpos,Zpos, pzone,psicutoff=',rpos,zpos,pzone,psicutoff
          !pause
            if(psicutoff.ne.0.d0  .and. (-pzone).lt.psicutoff) go to 170
            ! pzone is increasing from center to edge.
            ! equilpsi is decreasing, and psicutoff is based on equilpsi
!
!... accumulate hot ion birth rate
!
            hibrz(izone,ie,ib) = hibrz(izone,ie,ib) + one
            hicmz(izone,ie,ib,1) = hicmz(izone,ie,ib,1) + sgxnloc(1)
            hicmz(izone,ie,ib,2) = hicmz(izone,ie,ib,2) + sgxnloc(2)
            hicmz(izone,ie,ib,3) = hicmz(izone,ie,ib,3) + sgxnloc(3)
!           write(*,*) 'hibrz(1,1,2) = ',hibrz(1,1,2)
!          write(*,*) 'hibrz = ',hibrz(izone,ie,ib)
!          write(*,*) 'hicmz-1 = ',hicmz(izone,ie,ib,1)
!
!... Calculate pitch angle cosine at birth point; accumulate average
!    pitch angle cosine. Calculate toroidal angular momentum deposited
!    in shell by each monte carlo ion.
!990131            zetai = amin1(zetai,1.)
!990131            zetai = amax1(zetai,-1.)
            zetai = (-ypos*vx0+xpos*vy0)/(rpos*vbeam(ie,ib))    ! Original
            vplane = SQRT (vx0**2+vy0**2)                       ! Onetwo
!BH130914            zetai  = csgn*(xpos*vy0-ypos*vx0)/(rpos*vplane)     ! Onetwo
            zetai = min(zetai,one)
            zetai = max(zetai,-one)
            vtroid=zetai*vbeam(ie,ib)                           ! Original
!BH130914            vtroid = csgn*zetai*vplane                          ! Onetwo
            vrad   = -SQRT (1.0 - zetai**2) * vplane            ! Onetwo
            angmtm=rpos*vtroid*atwb*1.673e-24
            angmpz(izone,ie,ib) = angmpz(izone,ie,ib) + angmtm  ! Onetwo
            if(iborb.ne.1) then
              zetaz(izone,ie,ib) = zetaz(izone,ie,ib) + zetai
!              angmpz(izone,ie,ib)=angmpz(izone,ie,ib) + angmtm
            endif
!
!... save occasional birth point for subsequent plotting
!
            ic = ic + 1
            if(mod(ic-1,iskip).eq.0) then
              ipts = ipts+1
              xpts(ipts) = xpos
              ypts(ipts) = ypos
              zpts(ipts) = zpos
              rpts(ipts)=sqrt(xpos**2+ypos**2)
              vx(ipts) = vx0
              vy(ipts) = vy0
              vz(ipts) = vz0
            endif
!
!... calculate (r,w) grid location
!
            if(inubpat.gt.0)then
              i=(rpos-r(1))/drpat + one
              j=(zpos-w(1))/dzpat + one
              rzpat(i,j,ie,ib) = rzpat(i,j,ie,ib) + one
            endif
!
! ----------------------------------------------------------------------
!... calculate orbit widths and orbit loss
!    [BH140501:  Might be good to replace this with the G.C. orbit-based
!                calcs.
!
            if(iborb.ne.0)then
              nparx  = nparx + 1
              call freyorb(atwb,codeid,drutpi,drini,droti,ic,iskip, &
                izone,mfm1, &
                norb,pinsid,potsid,pzone,rinsid,rotsid,rzone,rpos, &
                zetai,vbeam(ie,ib),zpos,ipass,iaxis,ier,izp,wid,wt, &
                zeta_,angmtp)
!
              if(ier.ne.0) then
                fber(ie,ib) = fber(ie,ib) + one
              elseif(izp.gt.mfm1) then
                go to 175
              else
                nmbrz(izp)=nmbrz(izp)+1
                if(ipass.eq.1 .and. iaxis.eq.1) then
                  fb11(ie,ib) = fb11(ie,ib) + one
                  wb11(ie,ib) = wb11(ie,ib) + wid
                elseif(ipass.eq.1 .and. iaxis.eq.0) then
                  fb10(ie,ib) = fb10(ie,ib) + one
                  wb10(ie,ib) = wb10(ie,ib) + wid
                elseif(ipass.eq.0 .and. iaxis.eq.1) then
                  fb01(ie,ib) = fb01(ie,ib) + one
                  wb01(ie,ib) = wb01(ie,ib) + wid
                  ftrapfi(izp,ie,ib)=ftrapfi(izp,ie,ib) + one
                elseif(ipass.eq.0 .and. iaxis.eq.0) then
                  fb00(ie,ib) = fb00(ie,ib) + 1.
                  wb00(ie,ib) = wb00(ie,ib) + wid
                  ftrapfi(izp,ie,ib)=ftrapfi(izp,ie,ib)+1.
                endif
              endif
!
!... accumulate hot ion deposition rate and average pitch angle cosine
!
              do 155 i=1,mfm1
                hdepz(i,ie,ib) = hdepz(i,ie,ib) + wt(i)
                angmpz(i,ie,ib)=angmpz(i,ie,ib) + wt(i)*angmtp(i)
 155          zetaz(i,ie,ib) = zetaz(i,ie,ib) + wt(i)*zeta_(i)
!
            endif    ! iborb=1
! ----------------------------------------------------------------------
!
            go to 180
!
!... accumulate particles that are not deposited in plasma
!
 160        fap(ie,ib) = fap(ie,ib) + one
            go to 180
 170        fwall(ie,ib) = fwall(ie,ib) + one
            go to 180
 175        forb(ie,ib) = forb(ie,ib) + one
!----------------------------------------------------------------------
!     end loop over particles
!----------------------------------------------------------------------
!
            newpar = 0
 180     continue

!
!... normalize average pitch angle cosine. normalize momentum and birth
!    mode in each shell to a single particle.
!
            write(*,*) '** end loop over particles **'
!          do i=1,mfm1
!        write(*,'(3i4,2x,1p1e10.4,a16)')i,ie,ib,hibrz(i,ie,ib),' hibrz'
!          enddo
!          stop
!
          do 182 i=1,mfm1
            xnorm = hibrz(i,ie,ib)
            if (xnorm.ne.0) then
              angmpz(i,ie,ib)  = angmpz(i,ie,ib)/xnorm   ! Onetwo
              hicmz(i,ie,ib,1) = hicmz(i,ie,ib,1)/xnorm  ! Onetwo
              hicmz(i,ie,ib,2) = hicmz(i,ie,ib,2)/xnorm  ! Onetwo
              hicmz(i,ie,ib,3) = hicmz(i,ie,ib,3)/xnorm  ! Onetwo
            endif
            if(iborb.ne.0) xnorm = hdepz(i,ie,ib)
            if(xnorm.ne.zero) xnorm = 1./xnorm
 182      zetaz(i,ie,ib) = xnorm*zetaz(i,ie,ib)
!
!... get the fraction of trapped ions in each zone. if orbit effects are
!    not turned on or itrapfi=0 then this effect is not included.
!
          do 300 i=1,mfm1
            nmbrz(i) = MAX0(nmbrz(i),1)
            ftrapfi(i,ie,ib) = ftrapfi(i,ie,ib)/nmbrz(i)
 300      continue
!
!... get true toroidal ang. mtm. density deposition rate by
!    multiplying by weight of each m.c. particle
!
          do 183 i=1,mfm1
            angmpz(i,ie,ib)=angmpz(i,ie,ib)*bneut(ie,ib)/ &
              (psivol(i)*npart)
 183      continue
!
!... normalize loss fractions and hot ion birth rate
!
          fap(ie,ib)   = fap(ie,ib) / npar
          fwall(ie,ib) = fwall(ie,ib) / npar
          forb(ie,ib)  = forb(ie,ib) / npar
          xloss1 = fap(ie,ib) + fwall(ie,ib)
          xloss2 = fap(ie,ib) + fwall(ie,ib) + forb(ie,ib)
          if(xloss1.ge.1.) go to 201
!
!         do i = 1, mfm1
!         write(*,'(i4,2x,1p2e12.4,a12)') i, psivol(i), hibrz(i,ie,ib),
!     .        ' hibrz-pre'
!         enddo
          do 190 i=1,mfm1
            hibrz(i,ie,ib) = hibrz(i,ie,ib)*volume &
              /((one-xloss1)*npar*psivol(i))
            if(iborb.eq.0) then
              hdepz(i,ie,ib) = hibrz(i,ie,ib)
            elseif(xloss2.lt.1.) then
              hdepz(i,ie,ib) = hdepz(i,ie,ib)*volume &
                /((one-xloss2)*npar*psivol(i))
            endif
 190      continue
!      if(ie.eq.1 .and. ib.eq.2) stop
!
!... normalize orbit widths and fractions
!
          if(iborb.eq.0) go to 201
          if(fb11(ie,ib).ne.zero) wb11(ie,ib) = wb11(ie,ib)/fb11(ie,ib)
          if(fb10(ie,ib).ne.zero) wb10(ie,ib) = wb10(ie,ib)/fb10(ie,ib)
          if(fb01(ie,ib).ne.zero) wb01(ie,ib) = wb01(ie,ib)/fb01(ie,ib)
          if(fb00(ie,ib).ne.zero) wb00(ie,ib) = wb00(ie,ib)/fb00(ie,ib)
          fb11(ie,ib) = fb11(ie,ib)/nparx
          fb10(ie,ib) = fb10(ie,ib)/nparx
          fb01(ie,ib) = fb01(ie,ib)/nparx
          fb00(ie,ib) = fb00(ie,ib)/nparx
          fber(ie,ib) = fber(ie,ib)/nparx
!
 201    continue  ! end loop over energy components
        write(*,*) '** end loop over energy components **'
 200  continue  ! end loop over beams
!
!----------------------------------------------------------------------
!
!     end loop over beams and components
!
!----------------------------------------------------------------------
      write(*,*) 'end loop over beams...'
      write(*,*) 'ie,ib,ipar = ',ie,ib,ipar
!
      if (read_birth_pts.ne."enabled") then
!
       write(*,*) 'hibrz(i,1,1),hibrz(i,2,1),hibrz(i,3,1)'
       do i=1,mfm1
         write(*,'(i4,2x,0p9f9.4)') i, hibrz(i,1,1), &
              hibrz(i,2,1),hibrz(i,3,1), &
              hibrz(i,1,2),hibrz(i,2,2),hibrz(i,3,2), &
              hibrz(i,1,3),hibrz(i,2,3),hibrz(i,3,3) 
       enddo
!       stop
!
!... renormalize currents and powers to bptor
!
      do 240 ib=1,mb
        if(bptor(ib).gt.0.) then
          bptorx = 0.
          do 210 ie=1,3
 210      bptorx = bptorx + (1.-fap(ie,ib))*bpow(ie,ib)
          if(bptorx.gt.0.) then
            xnorm = bptor(ib)/bptorx
            bcur(ib) = xnorm*bcur(ib)
            do 220 ie=1,3
              bion(ie,ib) = xnorm*bion(ie,ib)
              bneut(ie,ib) = xnorm*bneut(ie,ib)
 220        bpow(ie,ib) = xnorm*bpow(ie,ib)
            do 230 ie=1,3
              do 231 i=1,mfm1
 231          angmpz(i,ie,ib)=angmpz(i,ie,ib)*xnorm
 230        continue
          endif
        endif
 240  continue
!
      else
!
!     Set bcur,bion,bneut,bpow from Nubeam list [above]
!
      endif  !On read_birth_pts.ne."enabled"
!
!     write hot ion birth rate, and lost franction:
!     BH130407:  hibrz is not used in cql3d, except for
!       following printout.  For that purpose, it will be useful
!       (in the future) to adjust the calculation of izone, perhaps
!       to base it on the rho coordinate and rotsid,  and to
!       fill in the rotsid() array.
      write(*,*)
      write(*,*)'freya, hot ion birth rate vs rotsid, for ib=1, ie=1:3'
      write(*,*) (rotsid(i),i=1,mfm1)
      write(*,*) (hibrz(i,1,1),i=1,mfm1)
      write(*,*) (hibrz(i,2,1),i=1,mfm1)
      write(*,*) (hibrz(i,3,1),i=1,mfm1)
!
      write(*,*)
      write(*,*) 'fap(ie=1:3)',(fap(i,1),i=1,3)
      write(*,*) 'fwall(ie=1:3)',(fwall(i,1),i=1,3)
      write(*,*) 'forb(ie=1:3)',(forb(i,1),i=1,3)
      write(*,*)
!      write(*,*)
!      write(*,*)'freya:i,xpts(i),ypts(i),zpts(i),rpts(i),vx(i),vy(i),vz(i)'
!      do i=1,5
!         write(*,*) i,xpts(i),ypts(i),zpts(i),rpts(i),vx(i),vy(i),vz(i)
!         write(*,*) i,xpts(i),ypts(i),zpts(i),rpts(i),vx(i),vy(i),vz(i)
!      enddo
!
!     calculate neutral beam deposition density on (r,w) grid
      if (inubpat.eq.1 .and. codeid.ne.'onedee') then
        call frnbdep2(psi,mi,mj,r,w,potsid,mf,rzpat,nrpat,nzpat, &
        ke,mb,sgxn,vbeam,hxfrac,inubpat)
      endif
!
      return
      end subroutine freya


      subroutine read_nubeam_data(filenm,nbirth_pts,atwb, &
                         nbirth_cmpts_nub,nshine_nub, &
                         pinj_nub,pabs_nub,x_nub,y_nub,z_nub, &
                         vx_nub,vy_nub,vz_nub,v_nub,en_nub,ie_nub, &
                         birth_rate_cmpts_nub,en_avg_cmpts_nub, &
                         pabs_cmpts_nub) 
!
!     Read NUBEAM data file
!
      use iso_c_binding, only : c_double
      implicit none

!
!     input: filenm, nbirth_pts, atwb
!     output: nbirth_cmpts_nub,nshine_nub,pinj_nub,pabs_nub,
!             x_nub,y_nub,z_nub,vx_nub,vy_nub,vz_nub,
!             en_nub,ie_nub,birth_rate_cmpts_nub,en_avg_cmpts_nub,
!             pabs_cmpts_nub
!
!     Of the nbirth_pts (frset namelist, checked against nubeam list
!     value): nbirth_cmpts_nub(1:3) is breakdown versus energy cmpt
!     In addition, nshine_nub particles shine through (from nubeam list).
!     [Breakdown into energy components, presumably most full energy,
!      is not provided.]
!     pinj_nub is injected power into the torus(Watts)
!     pabs_nub is absorbed (Watts)
!     en_nub/ie_nub(1:nbirth_pts) is energy/energy-component for each
!        particle in the list
!     birth_rate_cmpts_nub,en_avg_cmpts_nub,pabs_cmpts_nub(1:3) are the
!        neutral particle birth rates, average energy (keV) and power
!        absorbed for each of the three energy components.
!

      character*128 filenm
      integer nbirth_pts ! input
      real(c_double) :: atwb
      real(c_double),dimension(nbirth_pts) :: x_nub,y_nub,z_nub
      real(c_double),dimension(nbirth_pts) :: vx_nub,vy_nub,vz_nub
      real(c_double),dimension(nbirth_pts) :: v_nub,en_nub
      integer, dimension(nbirth_pts) :: ie_nub
      real(c_double),dimension(3) :: en_avg_cmpts_nub,pabs_cmpts_nub
      integer, dimension(3) :: nbirth_cmpts_nub
      real(c_double),dimension(3) :: birth_rate_cmpts_nub
      real(c_double),dimension(3) :: birth_rate_cmpts_nub1
      
      integer iunit,i,kode ! local
      integer nbirth_pts_nub, nshine_nub ! to read in
      real(c_double) :: pinj_nub,pabs_nub,total_source_rate ! to read in
      real(c_double) :: ergtkev,enmin,enmax,end2,en_avg,p_avg_tot_nub ! local
      real(c_double) :: bcnst

      iunit=14
      open(unit=iunit,file=filenm,status='old',iostat=kode)
      if (kode.ne.0) then
         WRITE(*,*)'NUBEAM file', filenm,' cannot be opened'
         STOP
      endif

      read(iunit,100) pinj_nub
      read(iunit,101) pabs_nub
      read(iunit,102) nbirth_pts_nub,total_source_rate
      read(iunit,103) nshine_nub
      read(iunit,104)

 100  format(///////,37x,e13.6)   !Skip 7 lines plus 37 columns.
 101  format(37x,e13.6)
 102  format(/////////,3x,i8,15x,e13.6)
 103  format(19x,i8)
 104  format(///)
 105  format(6e14.7)

!     Check n_birth_pts_nub against namelist value
      if (nbirth_pts_nub .ne. nbirth_pts) then
         WRITE(*,*)'STOP: inconsistency with NUBEAM file n_birth_pts'
         WRITE(*,*)'Namelist n_birth_pts =',nbirth_pts
         WRITE(*,*)'NUBEAM   n_birth_pts =',nbirth_pts_nub
         STOP
      endif

!     Read the NUBEAM birth points
      ergtkev=1.6022d-09
      do i=1,nbirth_pts
         read(iunit,105)x_nub(i),y_nub(i),z_nub(i), &
                        vx_nub(i),vy_nub(i),vz_nub(i)
!        Calc beam velocity and energy
         v_nub(i)=vx_nub(i)**2+vy_nub(i)**2+vz_nub(i)**2
         en_nub(i)=0.5*atwb*1.6726e-24*v_nub(i)/ergtkev
         v_nub(i)=sqrt(v_nub(i))
      enddo

!     Close the file.
      close(iunit)

!     Find maximum and minimum particle energies, and use them to
!     calculate if full,half, or 1/3 energy cmpt (ie_nub=1,2,3, resp.)
!     Relative beam deposition rate at each component propto total
!     birth points at each energy.
!     [It is assumed that the list contains only particles born
!      within the plasma.  Add coding to check this.]
      enmin=minval(en_nub,nbirth_pts)
      enmax=maxval(en_nub,nbirth_pts)
      end2=0.5*enmax
      nbirth_cmpts_nub=0d0
      en_avg_cmpts_nub=0d0

      do i=1,nbirth_pts
         if (en_nub(i) .gt. 1.5*end2) then
            ie_nub(i)=1
            nbirth_cmpts_nub(1)=nbirth_cmpts_nub(1)+1
            en_avg_cmpts_nub(1)=en_avg_cmpts_nub(1)+en_nub(i)
         elseif (en_nub(i) .lt. 0.833*end2) then
            ie_nub(i)=3
            nbirth_cmpts_nub(3)=nbirth_cmpts_nub(3)+1
            en_avg_cmpts_nub(3)=en_avg_cmpts_nub(3)+en_nub(i)
         else
            ie_nub(i)=2
            nbirth_cmpts_nub(2)=nbirth_cmpts_nub(2)+1
            en_avg_cmpts_nub(2)=en_avg_cmpts_nub(2)+en_nub(i)
         endif
      enddo

!     Avg energy of each of deposited beam components.
      do i=1,3
         en_avg=en_avg_cmpts_nub(i)
         en_avg_cmpts_nub(i)=en_avg/nbirth_cmpts_nub(i)
      enddo

!     To obtain the normalization for birth_rate_cmpts_nub(), invert
!     pabs(Watt)=sum(i=1,3)[birth_rate_cmpts_nub(i)*en_avg_cmpts_nub(i)]
!                *bcnst
!     to obtain bcnst. Adjust keV energy to joules:*ergtkev*1.e-7(J/erg)
      p_avg_tot_nub=0d0
      do i=1,3
         p_avg_tot_nub=p_avg_tot_nub + nbirth_cmpts_nub(i)* &
                       en_avg_cmpts_nub(i)*ergtkev*1.d-7  !Watts
      enddo
      bcnst=pabs_nub/p_avg_tot_nub

!     Normalize for deposited power and deposition rate at each energy
      do i=1,3
         birth_rate_cmpts_nub(i)=bcnst*nbirth_cmpts_nub(i)
         pabs_cmpts_nub(i)=birth_rate_cmpts_nub(i)*en_avg_cmpts_nub(i)* &
                           ergtkev*1.e-7   !Watts
      enddo

!     Alternative calc of birth_rate_cmpts_nub(i)
      do i=1,3
         birth_rate_cmpts_nub1(i)=(total_source_rate/nbirth_pts_nub)* &
                                  nbirth_cmpts_nub(i)
      enddo

!     printout birth_rates, as check for consistency
      write(*,*)
      write(*,*)'birth rate(particles/sec) for each cmpt based on'
      write(*,*)' cmpt   total power    total birth rate'
      do i=1,3
        write(*,106) i,birth_rate_cmpts_nub(i),birth_rate_cmpts_nub1(i)
      enddo
 106  format(i3,5x,e13.6,5x,e13.6)

      return
      end subroutine read_nubeam_data


      subroutine zone(drutpi,ki,mfm1, mim1, mjm1,dri,dzi,psiax,psi,r,z, &
           xpos,ypos,zpos,pzone,izone)
      use iso_c_binding, only : c_double
      implicit none
      save

      real(c_double) :: psi(ki,*),r(*),z(*) ! in arg. list
      real(c_double) :: drutpi,xpos,ypos,zpos,pzone,dri,dzi,psiax !in arg. list
      integer ki,mfm1,mim1,mjm1,izone ! in arg. list
      integer i,j !local
      real(c_double) :: rpos,psix,ptest,area1,area2,area3,area4,dum !local


      rpos=sqrt(xpos**2+ypos**2)
 120  i=(rpos-r(1))*dri+1.
      j=(zpos-z(1))*dzi+1.
      if(i.gt.mim1) i=mim1
      if(j.gt.mjm1) j=mjm1

      psix  = min(psi(i,j), psi(i+1,j), psi(i,j+1), psi(i+1,j+1))
      ptest = (psix-psiax)*(drutpi/mfm1)**2
      if(ptest.lt.0.02) go to 124
      area1=(rpos-r(i))*(zpos-z(j))
      area2=(r(i+1)-rpos)*(zpos-z(j))
      area3=(r(i+1)-rpos)*(z(j+1)-zpos)
      area4=(rpos-r(i))*(z(j+1)-zpos)
      pzone=(area3*psi(i,j)+area4*psi(i+1,j)+area1*psi(i+1,j+1) &
        +area2*psi(i,j+1))*dri*dzi
      go to 126
 124  call pfit(psi(i-1,j-1), r(i-1), z(j-1), rpos, zpos, ki, pzone, &
        dum, dum)
 126  pzone = max(pzone,psiax)
      izone = sqrt(pzone-psiax)*drutpi + 1.

      return
      end subroutine zone
