module ainsetpa_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine ainsetpa
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!..................................................................
!     This routine should be called after (first) namelist setup0 is read.
!     Determines values of lrz, lrors, and ls according to model
!     if lrzdiff.eq."enabled":
!     lrindx: index of radial surface on which equations are solved
!     indxlr: index of array lrindx giving lr_=> indxlr(lrindx(l=1,lrz))=l.
!     indxlr used to save space (arrays source, gone, egylosa,..)
!     (indxlr labels the radial surfaces from 1 to lrz)
!     if lsdiff.eq."enabled" (operative for cqlpmod="enabled"):
!     lsindx: index of orbit position s on which equ. are solved
!     indxls: index of array lsindx giving ls_=> indxls(lsindx(l=1,ls))=l.
!     (indxls labels the orbit position from 1 to ls).
!     lrindx(1)=namelist input
!     lrindx(2:ls)=lrindx(1), indxlr(lrindx(1))=1
!
!     Also sets values of some physical constants,
!     and some real*8 numbers used in subroutine arguments.
!..................................................................


!MPIINSERT_INCLUDE

!.......................................................................
!l    1. General specification
!.......................................................................

      if (lrz.eq.0) lrz=lrza
!000720      if (lrza.lt.lrz) lrz=lrza
      if (lrza.lt.lrz) stop 'ainsetpa: lrz.lt.lrza'
!000720      if (lrzmax .gt. lrza) lrzmax=lrza
      if (lrzmax .gt. lrza) stop 'ainsetpa: lrzmax.gt.lrza'
      if (ls.eq.0) ls=lsa
!000720      if (lsa.lt.ls) ls=lsa
      if (lsa.lt.ls) stop 'ainsetpa: lsa.lt.ls'
!000720      if (lsmax .gt. lsa) lsmax=lsa
      if (lsmax .gt. lsa) stop  'ainsetpa: lsmax .gt. lsa'

      if (lrzdiff .ne. "enabled") lrzmax=lrz
      if (lrzmax .le. lrz) then
        lrzdiff="disabled"
        lrzmax=lrz
        do 100 l=1,lrz
          lrindx(l)=l
 100    continue
      endif

      if (lsdiff .ne. "enabled") lsmax=ls
      if (lsmax .le. ls) then
        lsdiff="disabled"
        lsmax=ls
        do 101 l=1,ls
          lsindx(l)=l
 101    continue
      endif

      if (lrzdiff .eq. "enabled") then

!     shift lrindx if lrindx(0).ne.0
!     (assume error in namelist:lrindx=... instead of lrindx(1)=...)
        if (lrindx(0) .ne. 0) then
!MPIINSERT_IF_RANK_EQ_0
          WRITE(*,'(/"  WARNING: lrindx has been shifted by one index" &
            ," as lrindx(0).ne.0",/)')
!MPIINSERT_ENDIF_RANK
!dir$ novector
          do 102 ll=lrz,1,-1
            lrindx(ll) = lrindx(ll-1)
 102      continue
          lrindx(0)=0
        endif
      endif
      do 110 l=1,lrz
        indxlr(lrindx(l))=l
 110  continue

      if (cqlpmod .ne. "enabled") then

!.......................................................................
!l    2. CQL3D: CQL or CQL3D (radial coordinate as variable)
!.......................................................................

        lrors=lrz
        ls=1
        lsmax=1
        lsdiff="disabled"
        lsindx(1)=1

      else

!.......................................................................
!     3. CQLP (parallel): coordinate s along B is variable
!        instead of radius
!.......................................................................
 300    continue

        lrors=ls
        if (lsdiff .ne. "enabled") then
          do 310 ll=1,ls
            lsindx(ll)=ll
 310      continue
        endif
        do 311 l=1,ls
          indxls(lsindx(l))=l
 311    continue

!     so far, CQLP can only treat one flux surface
!     (need 4D storage arrays)
        lrz=1
        if (lrzdiff .eq. "disabled") then
          lrzmax=lrz
          lrindx(1)=1
          indxlr(1)=1
        endif
        do 315 l=1,ls
          lrindx(l)=lrindx(1)
 315    continue

      endif


!.......................................................................
!l    Set some physical, numerical and normalization constants.
!l    Numerical constants set to help maintain real*8 arithmetic.
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
      end
end module ainsetpa_mod
