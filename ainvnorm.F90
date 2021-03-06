! Copyright 2019 Garrett Wright, Princeton Plasma Physics Laboratory,
!    contracted by the U.S. Department of Energy (DE-AC02-09CH11466).
!
! This file is part of cql3d_f90. See LICENSE.
!
! cql3d_f90 is free software: you can redistribute it and/or modify it
! under the terms of the GNU Affero General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! cql3d_f90 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with cql3d_f90.  If not, see <https://www.gnu.org/licenses/>.

module ainvnorm_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use diagwrng_mod, only : diagwrng

  !---END USE

!
!

contains

      subroutine ainvnorm
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!.....................................................................
!     Determine mesh normalization constant vnorm.
!     vnorm is the maximum velocity for a non-
!     relativistic mesh. enorm is the maximum energy.
!     For a relativistic mesh vnorm goes over to maximum
!     momentum/unit rest mass.
!     Running electrons and ions in tandem as general species
!     requires special treatment.
!.....................................................................

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif


      if (tandem.eq."enabled") then
        xfac=-1.
        icase=1
        enorm=enorme
        if (relativ.eq."disabled") then
          vnorm=sqrt(enorm*ergtkev*2./fmass(kelecg))
        else
          vnorm=sqrt((enorm*ergtkev/(fmass(kelecg)*clite2)+1.)**2-1.)* &
            clight
        endif
        xlwr=sqrt(enormi*fmass(kelecg)/(enorme*fmass(kionn)))
        xpctlwr=.65
        xmdl=1.
        xpctmdl=.35
#ifdef __MPI
      if(mpirank.eq.0) then
#endif
        if(setup0%verbose>0) WRITE(*,*)' WARNING/ainvnorm: For tandem=enabled, &
         xfac,xlwr,xpctlwr,xmdl,xpctmdl are reset'
        if(setup0%verbose>0) WRITE(*,*)' ainvnorm: xfac,xlwr,xpctlwr,xmdl,xpctmdl=', &
                              xfac,xlwr,xpctlwr,xmdl,xpctmdl
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

      elseif((kenorm.ge.1).and.(kenorm.le.ngen).and.(enorm.gt.0.)) then
        icase=2
        if (relativ.eq."disabled") then
          vnorm=sqrt(enorm*ergtkev*2./fmass(kenorm))
        else
          vnorm=sqrt((enorm*ergtkev/(fmass(kenorm)*clite2)+1.)**2-1.)* &
            clight
        endif
      else if (vnorm.gt.0.) then
        icase=3
        if (relativ.eq."disabled") then
          enorm=fmass(1)*vnorm*vnorm/ergtkev*0.5
        else
          enorm=fmass(1)*clite2*(sqrt(1.+vnorm**2/clite2)-1.)/ergtkev
        endif
      else
        call diagwrng(1)
      endif
      vnorm2=vnorm*vnorm
      vnorm3=vnorm2*vnorm
      vnorm4=vnorm2*vnorm2
      cnorm=clight/vnorm
      cnorm2=cnorm**2
      cnorm3=cnorm**3
      if (relativ.eq."disabled") then
        cnorm2i=0. !YuP[07-2016] Bad way to control relativ.eq."disabled"
        cnormi=0.  !YuP[07-2016] Bad way to control relativ.eq."disabled"
      else
        cnorm2i=1./cnorm2
        cnormi=1./cnorm
      endif

!.....................................................................
!     energy conversion factor...
!.....................................................................

!DIR$ NEXTSCALAR
      do 12 k=1,ntotal
        fions(k)=.5*fmass(k)*vnorm**2/ergtkev
 12   continue

!..................................................................
!     Determine some mass ratios and normalization constants.
!..................................................................

      alp=7.3e-3
      r0=2.82e-13
      gacon2=2.*alp/(r0*fmass(kelec)*clight)
      do 10 i=1,ntotal
!DIR$ NEXTSCALAR
        do 11 k=1,ntotal
          gamt(i,k)=fmass(i)*fmass(k)*gacon2/(fmass(i)+fmass(k))
 11     continue
 10   continue
!DIR$ NEXTSCALAR
      do 40 k=1,ngen
        gam1=4.*pi*(charge*bnumb(k))**4/fmass(k)**2
        tnorm(k)=vnorm**3/(gam1*one_)
 40   continue
      return
      end subroutine ainvnorm


end module ainvnorm_mod
