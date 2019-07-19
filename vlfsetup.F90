! Copyright 2019 Garrett Wright, Princeton Plasma Physics Laboratory,
!    contracted by the U.S. Department of Energy (putnumberhere).
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

module vlfsetup_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use vlfalloc_mod, only : vlfalloc
  use zcunix_mod, only : zzbeslri

  !---END USE

!
!

contains

      subroutine vlfsetup
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine sets up arrays for vlf urf modules
!..................................................................


!...................................................................
!     mrf is the number of RF modes  being calculated.
!     nharms(1:mrf)) is the number of cyclotron harmonics calculated,
!     starting at nharm1.
!     Both nharms().gt.1 and mrf.gt.1 are now permitted by  the
!     storage scheme [BH, 060314].
!...................................................................

      pointer besl
      dimension besl(:) ! Now it's a local working array

      mrf=vlfmodes   !i.e., number of wave types

      if (mrf.gt.nmodsa) then
         write(*,*)'STOP:    vlfmodes.gt.nmodsa'
      endif

      deps=1.e-10
      do k=1,mrf
         nharms(k)=vlfharms(k)+deps
         nharm1(k)=vlfharm1(k)+deps
      enddo

!BH060314      one=1.
!BH060314      nharms=max(one,vlfharms)

!     mrfn is the number of wave "modes", that is, the sum over
!     wave types of the number of harmonics for each wave type
      mrfn=0
      do k=1,mrf
         if (nharms(k).eq.0) then
            mrfn=mrfn+1
         else
            mrfn=mrfn+nharms(k)
         endif
      enddo

      if (mrfn.gt.nmodsa) then
         write(*,*)'vlfsetup: mrfn>nmodsa.  mrfn,nmodsa=',mrfn,nmodsa
         STOP 'Increase nmodsa.'
      endif

!BH020913, added following line.
!BH060314      nharm1=vlfharm1(1)

!BH060314      mrfn=max0(mrf,nharms)

!BH060314c     Test for non-implemented options:
!BH060314      if (vlfmodes.gt.1.and.vlfharms.gt.1) then
!BH060314        print 102
!BH060314 102    format(' vlfmodes.gt.1.and.vlfharms.gt.1')
!BH060314        stop ' in vlfsetup'
!BH060314      endif
!BH060314      if (nharms.gt.nmodsa)  stop 'nharms.gt.nmodsa'

!.......................................................................
!
!     Set up table krfn(1:mrfn) pointing to wave type index.
!     Set up table irfn(1:mrf) pointing to the wave mode index of
!       the lowest harmonic for each wave type.
!     Set up nharm(1:mrfn).
!     NOTE:  Cannot set these rrays earlier in code because
!            nharm is possibly read above from ray data files,
!            giving value for nharm1.
!
!.......................................................................

      k=0
      do krf=1,mrf
         do kk=1,nharms(krf)
            k=k+1
            krfn(k)=krf
            if (kk.eq.1) irfn(krf)=k
            nharm(k)=nharm1(krf)+(kk-1)
         enddo
      enddo

      write(*,*)
      write(*,*)'vlfsetup: mrf = ',mrf
      write(*,*)'vlfsetup: mrfn = ',mrfn
      write(*,*)'vlfsetup: irfn = ',irfn
      write(*,*)'vlfsetup: krfn = ',krfn
      write(*,*)'vlfsetup: nharm1 = ',nharm1
      write(*,*)'vlfsetup: nharms = ',nharms
      write(*,*)'vlfsetup: nharm = ',nharm


!     Duplicate ray data in sets nharms(krf) long using the first ray data
!     for each ray type krf, if there is more than one harmonic:
!     Need to do this duplication in reverse order, k=mrfn,1,-1
!     in order not to overwrite data initially stored as krf=1,mrf.

      do k=mrfn,1,-1
         freqcy(k)=vlffreq(krfn(k))
         omega(k)=twopi*vlffreq(krfn(k))
         vlfnperp(k)=vlfnperp(krfn(k))
         vlfnp(k)=vlfnp(krfn(k))
         vlfdnp(k)=vlfdnp(krfn(k))
         vlfddnp(k)=vlfddnp(krfn(k))
         vlfeplus(k)=vlfeplus(krfn(k))
         vlfemin(k)=vlfemin(krfn(k))
         vlfpol(k)=vlfpol(krfn(k))*pi/180.
         vlfdpol(k)=vlfdpol(krfn(k))*pi/180.
         vlfddpol(k)=vlfddpol(krfn(k))*pi/180.
         vlfdnorm(k)=vlfdnorm(krfn(k))
      enddo

      write(*,*)'vlfsetup: vlfdnorm(k),k=1,mrfn ', &
                          (vlfdnorm(k),k=1,mrfn)

      vlfpol_inrange=0.d0
      !YuP[03-2016] Scan all poloidal angles along flux surface:
      ! if vlfpol is outside of range of all pol(), print warning.
      !Note: in a mirror machine the range of pol() is limited,
      !      usually less than [-90;+90] degrees.
      ! In a tokamak (closed surfaces), it is [0;180] or [0;360] degrees
      ! depending on eqsym value.
      do k=1,mrfn
      do l=1,lz
      if( abs(pol(l,lr_)-vlfpol(k)) .lt. 0.5*vlfdpol(k) )then
        ! ok, this pol() is near vlfpol, within the +/- vlfdpol/2 range
        vlfpol_inrange=1.0 ! value is changed
      endif
      enddo
      enddo
      if(vlfpol_inrange.eq.0.d0)then
        ! the above condition was not met for any vlfpol(k)
        WRITE(*,*)'vlfsetup: vlfpol is outside of pol.angle range.'
        WRITE(*,*)'vlfsetup: pol(1,setup0%lrz)=',pol(1,setup0%lrz), &
                          '  pol(lz,setup0%lrz)=',pol(lz,setup0%lrz)
        WRITE(*,*)'vlfsetup: vlfpol(k)=',vlfpol ! print for all k
        stop
      endif

!..................................................................
!     Allocate arrays
!..................................................................
      jjx=((jx-1)/ibytes)*ibytes+ibytes

      call vlfalloc

      do 15 j=1,jx
        sx(j)=x(j)/gamma(j)
 15   continue

!..................................................................
!     "l" refers to poloidal position z(l,lr_)
!     lr_ is the flux surface label.
!..................................................................

      do 20 l=1,lz
        do 21 i=1,iyh
          cosmz(i,l,lr_)=-.5*(cosz(i+1,l,lr_)+cosz(i,l,lr_))
          cosmz(iy+1-i,l,lr_)=-cosmz(i,l,lr_)
 21     continue
        cosmz(iy,l,lr_)=-cosz(iy,l,lr_)


 20   continue

!..................................................................
!     Fill in bessel function table.
!..................................................................
      nharmx = 1
      do k=1,mrfn
         nharmx=MAX(nharmx,nharm(k))
      enddo
      write(*,*)'vlfsetup:  Before allocate'
      allocate(besl(nharmx+2),STAT=istat) !-YuP->added
      call bcast(besl,zero,SIZE(besl))     !-YuP->added
      allocate(jbm1(nbssltbl,mrfn),STAT=istat) !-YuP-> modified to real
      call bcast(jbm1,zero,SIZE(jbm1))          !-YuP-> modified to real
      allocate(jb0(nbssltbl,mrfn),STAT=istat)  !-YuP-> modified to real
      call bcast(jb0,zero,SIZE(jbm1))           !-YuP-> modified to real
      allocate(jbp1(nbssltbl,mrfn),STAT=istat) !-YuP-> modified to real
      call bcast(jbp1,zero,SIZE(jbm1))          !-YuP-> modified to real
      write(*,*)'vlfsetup:  After allocate'
!..................................................................
!     Loop over excitation modes
!..................................................................

      do 200 k=1,mrfn

!..................................................................
!     Find maximum k-perp, accounting for aspect ratio
!..................................................................

!BH000416        xkp=vlfnperp(k)*freqcy(k)*twopi/clight*solr(1,lr_)/
!BH000416     +        solr(lorbit(lr_),lr_)
        xkp=vlfnperp(k)*freqcy(k)*twopi/clight*solrz(1,lr_)/ &
              solrz(lz,lr_)

!..................................................................
!     Guess at minimum cyclotron frequency - leave bpol out of estimate
!     this will decrease guess.
!..................................................................

        bvalmin=fpsi(setup0%lrindx(lrors))/solr(lorbit(lr_),lr_)
        !YuP[03-2016] But in a mirror machine, Btor=0 (fpsi=0).
        !write(*,*)'vlfsetup: krf,bvalmin=',k,bvalmin
        bvalmin=bmidplne(lr_) !YuP[03-2016]

        wcemin=abs(bnumb(1))*charge*bvalmin/(fmass(1)*clight)

!..................................................................
!     Maximum argument in table will be argmax - give it a 40% boost
!..................................................................

        argmax=xkp*vnorm/wcemin*1.4
        write(*,'(a,3i6,2e12.3)')'vlfsetup: n,lr_,krf,bvalmin,argmax=', &
          n,lr_,k,bvalmin,argmax
        bsslstp(k)=argmax/(nbssltbl-1)
        arg=0.
        do 30 i=1,nbssltbl
          call zzbeslri(arg,nharm(k)+2,0,besl,ncalc)
          if (ncalc.ne.nharm(k)+2) stop ' in vlfsetup'
          if(nharm(k).eq.0)  then
            jbm1(i,k)=-besl(2)
            jb0(i,k)=besl(1)
            jbp1(i,k)=besl(2)
          elseif (nharm(k).gt.0)  then
            jbm1(i,k)=besl(nharm(k))
            jb0(i,k)=besl(nharm(k)+1)
            jbp1(i,k)=besl(nharm(k)+2)
          endif
          arg=arg+bsslstp(k)
 30     continue
        argmax=arg-bsslstp(k)

 200  continue ! k=1,mrfn

      return
      end subroutine vlfsetup


end module vlfsetup_mod
