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

module tdtrvint_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine tdtrvint
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     This routine creates some integration coefficients on the
!     transport mesh.
!..............................................................


      if (transp.eq."disabled") return

      call bcast(vpint_,zero,iy*setup0%lrzmax)
      call bcast(cynt2_,zero,iy*setup0%lrz)
      call bcast(vptb_,zero,iy*(setup0%lrzmax+1))
      call bcast(bovcos,zero,(iy+1)*(setup0%lrz+1))
      call bcast(cosovb,zero,(iy+1)*(setup0%lrz+1))

      do 10 l=1,lrors
        ilr = setup0%lrindx(l)
        itl=itl_(l)
        itu=itu_(l)
        iyh=iyh_(l)
        do 20 i=1,iy_(l)
          vpint_(i,ilr)=vpint(i,ilr)
          cynt2_(i,l)=cynt2(i,l)
          vptb_(i,ilr)=vptb(i,ilr)
 20     continue

!BH070419:   removing special itl,itu treatment for ipacktp=0
        if (ipacktp.eq.3) then

        vpint_(itl-2,ilr)=vpint(itl-2,ilr)+vpint(itl-1,ilr)+ &
          .5*vpint(itl,ilr)
        vpint_(itl-1,ilr)=0.
        vpint_(itl,ilr)=0.
        vpint_(itl+1,ilr)=0.
        vpint_(itu+1,ilr)=0.
        vpint_(itu,ilr)=0.
        vpint_(itu-1,ilr)=0.
        vpint_(itl+2,ilr)=vpint(itl+2,ilr)+vpint(itl+1,ilr)+ &
          .5*vpint(itl,ilr)
        vpint_(itu-2,ilr)=vpint_(itl+2,ilr)
        vpint_(itu+2,ilr)=vpint_(itl-2,ilr)
        cynt2_(itl-2,l)=cynt2(itl-2,l)+cynt2(itl-1,l)+.5*cynt2(itl,l)
        cynt2_(itl-1,l)=0.
        cynt2_(itl,l)=0.
        cynt2_(itl+1,l)=0.
        cynt2_(itu+1,l)=0.
        cynt2_(itu,l)=0.
        cynt2_(itu-1,l)=0.
        cynt2_(itl+2,l)=cynt2(itl+2,l)+cynt2(itl+1,l)+cynt2(itl,l)*.5
        cynt2_(itu-2,l)=cynt2_(itl+2,l)
        cynt2_(itu+2,l)=cynt2_(itl-2,l)
        vptb_(itl-2,ilr)=vpint_(itl-2,ilr)/cynt2_(itl-2,l)
        vptb_(itl+2,ilr)=vpint_(itl+2,ilr)/cynt2_(itl+2,l)
        vptb_(itu-2,ilr)=vptb_(itl+2,ilr)
        vptb_(itu+2,ilr)=vptb_(itl-2,ilr)
        vptb_(itl-1,ilr)=0.
        vptb_(itl,ilr)=0.
        vptb_(itl+1,ilr)=0.
        vptb_(itu+1,ilr)=0.
        vptb_(itu,ilr)=0.
        vptb_(itu-1,ilr)=0.

        elseif (ipacktp.ne.0) then
           write(*,*)'STOP in tdtrvint:  Check ipacktp'
           stop
        endif

 10     continue

!..............................................................
!     Define the Jacobian of the transformation between (E,mu) and
!     and (theta,v) coordinates here. The procedure utilized below
!     is done to assure conservation in the transport difference
!     scheme. The evaluation coss(i,l)/bmidplne(l) is to be
!     avoided. In any event, we exclude the factor m**2/gamma
!     which does not depend on pitch angle.
!..............................................................

      do 80 l=1,lrors
        if (meshy.eq."fixed_mu") then
          do 30 i=1,iytr(l)/2
            m=iytr(lrors)+1-i
            cosovb(idx(i,l),l)=cynt2_(idx(i,lrors),lrors)/ &
              cynt2_(idx(i,l),l)*coss(idx(i,lrors),lrors)/ &
              bmidplne(setup0%lrindx(lrors))
            cosovb(idx(m,l),l)=cynt2_(idx(m,lrors),lrors)/ &
              cynt2_(idx(m,l),l)*coss(idx(m,lrors),lrors)/ &
              bmidplne(setup0%lrindx(lrors))
 30       continue
        else if (meshy.eq."fixed_y") then
          do 40 i=1,iytr(lrors)
            il=idx(i,l)
            cosovb(il,l)=1.
            bovcos(il,l)=1.
 40       continue
        endif
 80   continue
      if (meshy.ne."fixed_mu") go to 500
      do 60 l=1,lrors-1
        do 50 i=1,iytr(l)/2
          ii=iytr(lrors)+1-i
          i2=idx(ii,l+1)
          i1=idx(ii,l)

          ill=idx(i,l+1)
          il=idx(i,l)
          if(meshy.eq."fixed_mu") then
            bovcos(il,l)=.5*(1./cosovb(il,l)+ 1./cosovb(ill,l+1))
            bovcos(i1,l)=.5*(1./cosovb(i1,l)+ 1./cosovb(i2,l+1))
          else
            bovcos(il,l)=one
            bovcos(i1,l)=one
          endif
 50     continue
 60   continue
 500  continue
      do 70 i=1,iymax
        bovcos(i,0)=0.
        cosovb(i,0)=0.
 70   continue
      return
      end subroutine tdtrvint

end module tdtrvint_mod
