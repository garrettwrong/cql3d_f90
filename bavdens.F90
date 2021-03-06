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

module bavdens_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use exlin_mod, only : exlin
  use r8subs_mod, only : dcopy
  use psif_mod, only : psif

  !---END USE

!

contains

      subroutine bavdens(k)
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..........................................................
!     This routine computes the bounce average of certain
!     density related quantities.
!.............................................................


      do 1 i=1,iy
        bavdn(i,lr_)=reden(k,lr_)
        bavpd(i,lr_)=batot(i,lr_)*sinn(i,lmdpln_)*reden(k,lr_)
 1    continue

!..................................................................
!     Exit if locquas.ne. "enabled"
!..................................................................

      if (kelecg.gt.0.or.locquas.eq."disabled".or.kelecm.ne.k) return
      call bcast(temc3,zero,iy)
      call bcast(temc4,zero,iy)


      lrange=lz
      if (numclas .eq. 1) lrange=lz/2+1


      do 30 l=1,lrange
        do 31 i=1,imax(l,lr_)

          ! 1. All orbits that reach/pass a given poloidal node l:
          ! passing (i<itl), or trapped that could reach/pass this l;
          ! also includes Last Trapped, with tip at l=lz_bmax
          !(for such particle, lmax(itl)=lz_bmax; see micxinil)
            ax=dtau(i,l,lr_)/tau(i,lr_)
            y1=densz(l,ngen+1,negyrg,lr_)
            temc3(i)=temc3(i)+y1*ax
            temc4(i)=temc4(i)+(y1/bbpsi(l,lr_))*ax

          ! 2a. Trapped, with tips between l and l+1 (ABOVE midplane):
          if (l.eq.lmax(i,lr_) .and. l.lt.lz_bmax(lr_)) then
            ! Add contribution from orbit tip:
            ax=dtau(i,l+1,lr_)/tau(i,lr_)
            y1=densz(l,  ngen+1,negyrg,lr_)
            y2=densz(l+1,ngen+1,negyrg,lr_)
            xx=zboun(i,lr_)
            qq=exlin(y1,y2,z(l,lr_),z(l+1,lr_),xx)
            temc3(i)=temc3(i)+qq*ax
            temc4(i)=temc4(i)+(qq/psif(zboun(i,lr_)))*ax
          endif

          ! 2b. Trapped, with tips between l and l-1 (BELOW midplane):
          if (l.eq.lmax(i+iyh,lr_) .and. l.gt.lz_bmax(lr_)) then
            ! Add contribution from orbit tip:
            ax=dtau(i,l-1,lr_)/tau(i,lr_)
            y1=densz(l,  ngen+1,negyrg,lr_)
            y2=densz(l-1,ngen+1,negyrg,lr_)
            xx=zboun(i,lr_)
            qq=exlin(y1,y2,z(l,lr_),z(l-1,lr_),xx)
            temc3(i)=temc3(i)+qq*ax
            temc4(i)=temc4(i)+(qq/psif(zboun(i,lr_)))*ax
          endif

 31     continue
 30   continue


      ! Symmetrize around pitch=pi/2
      do 40 i=1,iyh
        iii=iy+1-i
        temc3(iii)=temc3(i)
        temc4(iii)=temc4(i)
 40   continue
      call dcopy(iy,temc3,1,bavdn,1)


      ! Symmetrize around pitch=pi/2
      do 50 i=2,iyh
        iii=iy+1-i
        bavpd(i,lr_)=(temc4(i)/sinn(i,lmdpln_)**2 &
          -temc3(i))*tann(i,lmdpln_)**2*sinn(i,lmdpln_)
        bavpd(iii,lr_)=+bavpd(i,lr_)
 50   continue
      bavpd(1,lr_)=0.
      bavpd(iy,lr_)=0.


      return
      end subroutine bavdens


end module bavdens_mod
