!
!
      subroutine frinitl
      use param_mod
      implicit none
      include 'frname_decl.h'
      character*8 machinei
      
      integer i,j,k ! local
      real(c_double) :: ep100 ! local

      REAL RILIN

      ep100=1.d+100

!..................................................................
!     Jump past the "fr" namelist if "fr" is disabled.
!..................................................................

!     some preset to have clean namelist output
        do i=1,nap
          do j=1,kb
            ashape(i,j) = 'disabled'
          end do
        end do
        do j=1,kb
           bshape(j)='disabled'
        enddo

!..................................................................
!     This routine sets defaults and reads in data for the
!     neutral beam (NFREYA) code as it exists in ONETWO.
!..................................................................

!
!     NEUTRAL BEAM HEATING PARAMETERS
!--------------------------------------------------------------------

!******************************************************************
!     BELOW TO **** NOT USED IN CQL3D  **** ONETWO DIVERGENCE

!     timbplt    Times (up to 5) to produce data for Freya -like plots
!     of beam deposition. o/p is processed by nubplt.
!     Defaulted to off.
!     timbplt(1).le.time0.and.beamon.lt.time0 gives o/p
!     for initial time.
!     beamon      Time (sec) at which beam heating is turned on
!     btime       Time interval (sec) during which beam heating is on
!     ibcur       Flag for neutral beam driven current
!     1, include beam driven current
!     0, neglect beam driven current
!     ibcx        Flag for secondary charge exchange between fast ions
!     and thermal neutrals
!     1, include secondary charge exchange
!     0, neglect secondary charge exchange
!     ibslow      Flag for neutral beam slowing down
!     1, include neutral beam slowing down
!     0, neglect neutral beam slowing down
!     fionx       Allows testing for non-classical slowing down
!     (see subroutine slow1).
!     nameb       Name of neutral species in beam
!     'h', 'd', 't', 'dt'
!     Must be primary species #1 or #2.
!     relnub      Minimum relative change in ion density or electron
!     temperature between neutral beam calculations (suggest 0.10)
!      tfusbb      'thermal fusion balance for beams', fraction by which
!      the net energy gain from thermal fusion must exceed
!      the net energy loss for automatic beam turnoff.
!      If tfusbb = 0., automatic beam turnoff is not done.
!      iddcal     Flag controlling treatment of beam effects on fuscal,
!      the calculated fusion neutron rate:
!      0 = do not include knock-on or beam-d neutrons in fuscal
!      1 = include only knock-on neutrons in fuscal
!      2 = include only beam-d neutrons in fuscal
!      3 = include both knock-on and beam-d neutrons in fuscal
!     fdbeam     Fraction of deuterium atoms in neutral beam.

!     ABOVE NOT USED IN CQL3D** ONETWO DIVERGENCE**
!******************************************************************
!     ranseed    Starting seed for random number generator used in the
!     Freya determination of the beam deposition.
!     npart       Number of particles followed into plasma
!     (suggest 10000)
!     npskip     Ratio of number of particles followed into plasma
!     to number of source particles (suggest 1)
!     iborb      Flag for modeling orbit effects on beam-
!     generated fast ions
!     1, model orbit effects
!     0, do not model orbit effects
!     itrapfi    Flag for trapped fast ions.  If iborb=1 then
!     setting itrapfi=1 will modify the beam driven
!     current by the initial trapped ion fraction.
!     (Subsequent pitch angle diffusion is not!!! taken into
!     account).  itrapfi=0 neglects this effect.  If iborb=0,
!     itrapfi has no affect.  itrapfi=0 is default
!     iexcit    excited beam state option
!     0,  do not use excited state cross sections (default)
!     1,  use hexnb package but do not include excitations in
!     its calculation of cross sections
!     2,  use hexnb package, include excitations in calculation
!     of cross sections
!     inubpat   two-dimensional beam deposition option
!     0,  do not calculate beam deposition in (r,z) coordinates
!     1,  calculate beam deposition on eqdisk (r,z) coordinates,
!     write third excited state population to file 'beamdep'
!     2,  same as inubpat=1 except rescale (r,z) grid according
!     to npat (see below)
!     if inubpat.gt.0, iexcit is automatically set, iexcit=2
!     npat      modified (r,z) grid dimensions, used only for inubpat=2
!     npat(1)=number of elements in 'r' (<=2*nnra)
!     npat(2)=number of elements in 'z' (<=2*nnza)
!     defaults, npat(1)=nnra, npat(2)=nnza
!     mf not read in - ONETWO DIVERGENCE
!     mf         Number of flux zones plus 1 (max = 81)
!
!     In the following list the index ib designates the beam injector,
!     while ie refers to one of the three energy components.
!     iap refers to one of the aperatures.
!
!     nbeams          Number of neutral beam injectors (.le.kb)
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
!     naptr           Total number of aperatures encountered by a particle
!     as is moves from the source into the plasma chamber
!     Maximum is specified by parameter nap (=10).
!     First set of aperatures encountered by the particle
!     are assumed centered on the source axis, and subseq
!     aperatures are centered on the beamline axis; the
!     distinction is made through ashape.
!     anglev(ib)      Vertical angle (degrees) between optical axis
!     and horizontal plane; a positive value indicates
!     particles move upward
!     angleh(ib)      Horizontal angle (degrees) between optical axis and
!     vertical plane passing through pivot point and
!     toroidal axis; a zero value denotes perpendicular
!     injection, while a positive value indicates par-
!     ticles move in the co-current direction
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
!     'rect-lps':  rect. long pulse source (DIII-D only)
!     a choice of short or long pulse sources is
!     available by injector (not by source).  on
!     or both injectors may be long pulse by
!     setting one or both to 'rect-lps'
!     CAUTION:  DIII-D sources are defaulted to lps specs.
!     It is the user's responsibility to overide these for
!     sps configuration(s).
!     bheigh(ib)      Height of source (cm)
!     Default is bshape(1)=bshape(2)='rect-lps'
!     bwidth(ib)      Width of source (cm); diameter for
!     circular source.
!     bhfoc(ib)       Horizontal focal length of source (cm)
!     bvfoc(ib)       Vertical focal length of source (cm)
!     bhdiv(ib)       Horizontal divergence of source (degrees)
!     bvdiv(ib)       Vertical divergence of source (degrees)
!     ebkev(ib)       Maximum particle energy in source (keV)
!     fbcur(ie,ib)    Fraction of current at energy ebkeV/ie
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
!
!     hdepsmth        set this parm to a positive value (gt.0.0) to turn of
!     the smoothing of hibrz and hdepz in sub postnub.
!     If this option is used then enough zones must be spec
!     for adequate resolution (zones=number of radial grid
!     and enough injected neutrals must be followed to minimize
!     the statistical noise enough so that no greatly uneven
!     profiles result! this option was added because the smoothing
!     of the profiles by subroutine smooth can lead to unphysical
!     forms of the birth and deposition profiles.
!--------------------------------------------------------------------
!     ONETWO DIVERGENCE BELOW
      nprim=1  !number of primary species. Needs to be equal to ngen.
      nimp=0   !If impurities are present, should include at least =1.
      machinei="iter"
      if (machinei.eq."iter") go to 1660
      do 1650 i=1,5
        timbplt(i)=1000.
 1650 continue
      beamon = 1000.  ! Not used for anything
      btime  = 0.     ! Not used for anything
      ibcur = 1
      itrapfi=0
      ibcx = 1
      ibslow = 1
      fionx = 0.0
      nameb = 'h'
      relnub=.10
!NLremoved      tfusbb = 0.
!NLremoved      iddcal = 3
!NLremoved      fdbeam = 0.150e-3

 1660 continue
!BH091020:  Adding above defaults which were skipped for some unknown
!BH091020:  reason, for machinei='iter'
!BH091020:  Simply zeroing.
      beamon = 0.  ! Not used for anything
      btime  = 0.  ! Not used for anything
      ibcur = 0
      itrapfi=0
      ibcx = 0
      ibslow = 0
      fionx = 0.0
!      nameb = 'h'
      relnub=0.

      nameb = 'd'
      ranseed=7**7
!     ONETWO DIVERGENCE
      npart= 150000 ! YuP-101211: default value. Was: npart=nap
      hdepsmth=-1.0
!---smoothing is normally on by above line
      npskip = 1
!     ONETWO DIVERGENCE BELOW
      iborb = 0
      inubpat=0
      npat(1)=nnra
      npat(2)=nnza
      mf = 41
      nbeams = 1
!
!     ONETWO DIVERGENCE BELOW
      nsourc=1
!
      if (ke.lt.3) then
         WRITE(*,*)'frinitl:  parameter ke needs to be .ge.3'
         stop
      endif
!
!     DIII beam input
      if(machinei.ne.'doub-iii')go to 9020
      naptr=2
      do 1700 i=1,kb
        anglev(i)  = 0.
        angleh(i)  = 13.5
        bvofset(i)=39.75
        bhofset(i)=0.0
        bleni(i)=553.88
        bcur(i)    = 110.
        bptor(i)   = 3.5e6
        bshape(i)  = 'rect'
        bheigh(i)  = 10.
        bwidth(i)  = 40.
        bhfoc(i)   = 480.
        bvfoc(i)   = 550.
        bhdiv(i)   = 1.4
        bvdiv(i)   = .45
        ebkev(i)   = 80.0
        if (ke.lt.3) then
           WRITE(*,*)'frinitl:  parameter ke needs to be .ge.3'
           stop
        else
           fbcur(1,i) = .6
           fbcur(2,i) = .3
           fbcur(3,i) = .1
        endif
        if (ke.ge.4) then
           do k=4,ke
              fbcur(k,i)=0.
           enddo
        endif
        sfrac1(i)  = 0.5
        blenp(i) = 486.13
        rpivot(i)  = 270.
        zpivot(i)  = 89.
        if (nap.lt.2) then
           WRITE(*,*)'frinitl:  parameter nap needs to be .ge.2'
           stop
        else
           ashape(1,i)='s-vert'
           ashape(2,i)='b-horiz'
           aheigh(1,i)=10.1
           aheigh(2,1)=0.0
           awidth(1,i)=0.0
           awidth(2,i)=32.
           alen(1,i)=442.*1.0026
           alen(2,i)=456.*1.0026
        endif
        if (nap.ge.3) then
           do k=3,nap
              ashape(k,i)='disabled'
              aheigh(k,i)=0.0
              awidth(k,i)=0.0
              alen(k,i)=0.0
           enddo
        endif
 1700 continue
      go to 9030
!
!     ITER BEAM INPUT.
 9020 naptr=4
      do 1705  i=1,kb
        anglev(i)=0.0
        angleh(i)=19.5
        bvofset(i)=0.0
        bhofset(i)=42.074
        bleni(i)=556.808
        bcur(i)=110.
        bptor(i)=10.e6
        bshape(i)='rect-lps'
        bheigh(i)=48.
        bwidth(i)=12.
        bhdiv(i)=.50
        bvdiv(i)=1.3
        if (ke.lt.3) then
           WRITE(*,*)'frinitl:  parameter ke needs to be .ge.3'
           stop
        else
           fbcur(1,i)=0.7
           fbcur(2,i)=0.2
           fbcur(3,i)=0.1
        endif
        if (ke.ge.4) then
           do k=4,ke
              fbcur(k,i)=0.0
           enddo
        endif
        bhfoc(i)=ep100
        bvfoc(i)=1.e3
        ebkev(i)=75.
        sfrac1(i)=0.5d0
        blenp(i)=539.
        rpivot(i)=286.6
        zpivot(i)=0.0d0
        if (nap.lt.4) then
           WRITE(*,*)'frinitl:  parameter nap needs to be .ge.4'
           stop
        else
           ashape(1,i)='s-rect'
           ashape(2,i)='s-rect'
           ashape(3,i)='b-d3d'
           ashape(4,i)='b-circ'
           aheigh(1,i)=47.8
           aheigh(2,i)=48.0
           aheigh(3,i)=0.0
           aheigh(4,i)=0.0
           awidth(1,i)=13.8
           awidth(2,i)=17.7
           awidth(3,i)=0.0
           awidth(4,i)=50.9
           alen(1,i)=186.1
           alen(2,i)=346.0
           alen(3,i)=449.0
           alen(4,i)=500.
        endif
        if (nap.ge.5) then
           do k=5,nap
              ashape(k,i)='disabled'
              aheigh(k,i)=0.0
              awidth(k,i)=0.0
              alen(k,i)=0.0
           enddo
        endif

 1705 continue
 9030 continue
!---
!---following parms are used in subroutine hexnb
!---
      kdene=1
      kdeni=1
      kdenz=1
      ksvi=0
      ksvz=0
      ksve=0
      krad=1
      ngh=10
      ngl=10
      iexcit=0
      ilorent=0
      mstate=4

!     note izstrp=0 implies coronal equilibrium for impurity j in hexnb
      do 1720 j=1,kimp
        izstrp(j)=0
 1720 continue


      frmod="disabled"
      !gyro-radius correction for NBI deposition
      fr_gyro="disabled"
      frplt="enabled"
      nfrplt=300
      bmsprd=.03
      multiply="disabled"
      multiplyn=0

      beamplse="disabled"
      beampon=0.d0
      beampoff=0.d0

!     defaults for reading NUBEAM particle birth pt list
      read_birth_pts="disabled"
      nbirth_pts_files=1
      nbirth_pts=250000
      do i=1,24
         birth_pts_files(i)="notset"
      enddo

!     For removing NBI source at all psi outside of psicutoff:
      psicutoff=0.d0 ! if 0.0, no removal is done
      ! The value of psicutoff can be determined from
      ! screen output of rho and psi values.
      ! Select one of psi(lr) values, set it in cqlinput.


      return
      end subroutine frinitl
