module ainsetpa_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine ainsetpa
      use param_mod
      use cqlcomm_mod
      use cqlconf_mod, only : setup0
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine should be called after (first) namelist setup0 is read.
!     Determines values of setup0%lrz, lrors, and setup0%ls according to model
!     if setup0%lrzdiff.eq."enabled":
!     setup0%lrindx: index of radial surface on which equations are solved
!     indxlr: index of array setup0%lrindx giving lr_=> indxlr(setup0%lrindx(l=1,setup0%lrz))=l.
!     indxlr used to save space (arrays source, gone, egylosa,..)
!     (indxlr labels the radial surfaces from 1 to setup0%lrz)
!     if setup0%lsdiff.eq."enabled" (operative for setup0%cqlpmod="enabled"):
!     setup0%lsindx: index of orbit position s on which equ. are solved
!     indxls: index of array setup0%lsindx giving ls_=> indxls(setup0%lsindx(l=1,setup0%ls))=l.
!     (indxls labels the orbit position from 1 to setup0%ls).
!     setup0%lrindx(1)=namelist input
!     setup0%lrindx(2:setup0%ls)=setup0%lrindx(1), indxlr(setup0%lrindx(1))=1
!
!     Also sets values of some physical constants,
!     and some real(c_double) numbers used in subroutine arguments.
!..................................................................


!MPIINSERT_INCLUDE

!.......................................................................
!l    1. General specification
!.......................................................................

      if (setup0%lrz.eq.0) setup0%lrz=lrza
!000720      if (lrza.lt.setup0%lrz) setup0%lrz=lrza
      if (lrza.lt.setup0%lrz) stop 'ainsetpa: setup0%lrz.lt.lrza'
!000720      if (setup0%lrzmax .gt. lrza) setup0%lrzmax=lrza
      if (setup0%lrzmax .gt. lrza) stop 'ainsetpa: setup0%lrzmax.gt.lrza'
      if (setup0%ls.eq.0) setup0%ls=lsa
!000720      if (lsa.lt.setup0%ls) setup0%ls=lsa
      if (lsa.lt.setup0%ls) stop 'ainsetpa: lsa.lt.setup0%ls'
!000720      if (setup0%lsmax .gt. lsa) setup0%lsmax=lsa
      if (setup0%lsmax .gt. lsa) stop  'ainsetpa: setup0%lsmax .gt. lsa'

      if (setup0%lrzdiff .ne. "enabled") setup0%lrzmax=setup0%lrz
      if (setup0%lrzmax .le. setup0%lrz) then
        setup0%lrzdiff="disabled"
        setup0%lrzmax=setup0%lrz
        do 100 l=1,setup0%lrz
          setup0%lrindx(l)=l
 100    continue
      endif

      if (setup0%lsdiff .ne. "enabled") setup0%lsmax=setup0%ls
      if (setup0%lsmax .le. setup0%ls) then
        setup0%lsdiff="disabled"
        setup0%lsmax=setup0%ls
        do 101 l=1,setup0%ls
          setup0%lsindx(l)=l
 101    continue
      endif

      if (setup0%lrzdiff .eq. "enabled") then

!     shift setup0%lrindx if setup0%lrindx(0).ne.0
!     (assume error in namelist:setup0%lrindx=... instead of setup0%lrindx(1)=...)
        if (setup0%lrindx(0) .ne. 0) then
!MPIINSERT_IF_RANK_EQ_0
          WRITE(*,'(/"  WARNING: setup0%lrindx has been shifted by one index" &
            ," as setup0%lrindx(0).ne.0",/)')
!MPIINSERT_ENDIF_RANK
!dir$ novector
          do 102 ll=setup0%lrz,1,-1
            setup0%lrindx(ll) = setup0%lrindx(ll-1)
 102      continue
          setup0%lrindx(0)=0
        endif
      endif
      do 110 l=1,setup0%lrz
        indxlr(setup0%lrindx(l))=l
 110  continue

      if (setup0%cqlpmod .ne. "enabled") then

!.......................................................................
!l    2. CQL3D: CQL or CQL3D (radial coordinate as variable)
!.......................................................................

        lrors=setup0%lrz
        setup0%ls=1
        setup0%lsmax=1
        setup0%lsdiff="disabled"
        setup0%lsindx(1)=1

      else

!.......................................................................
!     3. CQLP (parallel): coordinate s along B is variable
!        instead of radius
!.......................................................................
 300    continue

        lrors=setup0%ls
        if (setup0%lsdiff .ne. "enabled") then
          do 310 ll=1,setup0%ls
            setup0%lsindx(ll)=ll
 310      continue
        endif
        do 311 l=1,setup0%ls
          indxls(setup0%lsindx(l))=l
 311    continue

!     so far, CQLP can only treat one flux surface
!     (need 4D storage arrays)
        setup0%lrz=1
        if (setup0%lrzdiff .eq. "disabled") then
          setup0%lrzmax=setup0%lrz
          setup0%lrindx(1)=1
          indxlr(1)=1
        endif
        do 315 l=1,setup0%ls
          setup0%lrindx(l)=setup0%lrindx(1)
 315    continue

      endif


!.......................................................................
!l    Set some physical, numerical and normalization constants.
!l    Numerical constants set to help maintain real(c_double) arithmetic.
!.......................................................................
      zero=0.d0
      one=1.d0
      two=2.d0
      three=3.d0
      four=4.d0
      sevenhun=700.d0
      one_=1.d0
      em6=1.d-6
      em8=1.d-8
      em10=1.d-10
      em12=1.d-12
      em14=1.d-14
      em37=1.d-37
      em40=1.d-40
      em90=1.d-90
      em100=1.d-100
      em120=1.d-120
      em300=1.d-300
      ep37=1.e+37
      ep90=1.d+90
      ep100=1.d+100
      half=one/two
      third=one/three

!      pi=3.141592653589793d0
      pi=atan2(zero,-one)
      pio180=pi/180.d0
      twopi=two*pi
      fourpi=4.0*pi
      stopm=-sqrt(two/pi)
      charge=4.8032d-10
      alp=7.2974d-3
      r0=2.817940285d-13
      pio2=pi*half
      ergtkev=1.6022d-09  !energy(ergs) associated with 1 keV
      clight=2.99792458d10
      clite2=clight**2
      proton=1.67262158d-24
      restmkev=510.998902d0

!.......................................................................
!     Additional variable initialization
!.......................................................................

      do 111 ll=0,lrza
        lmdpln(ll)=ll
 111  continue

      nonch=nstop+1 !!! Before 101220: nonch=noncha (==2000)

      return
      end subroutine ainsetpa

end module ainsetpa_mod
