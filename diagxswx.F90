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

module diagxswx_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use coefmidt_mod, only : coefmidt
  use coefmidv_mod, only : coefmidv
  use coefstup_mod, only : coefstup
  use diagdens_mod, only : diagdens
  use diagentr_mod, only : diagentr
  use diagentr_mod, only : gfu
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine diagxswx(k)
      use advnce_mod !here: in diagimpd(). To get gfu(), gfi()
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     routine computes density gain due to partlcle sources and losses
!     during the during velocity sweep in exsweep.
!     gains are positive - losses are negative.
!     1- density gained due to setting negative values of
!     distribution function to 0 during velocity sweep.
!     3- particle source contribution.
!     4- particle loss due to flux at high velocity end.
!     5- bad orbits (part of krook operator)
!     6- toroidal losses(part of krook operator)
!     Further - the normalized velocity flux at mesh point x(j+.5) is pu
!     into vflux(j,k,l_).
!..................................................................


!..................................................................
!     sgain(i,k) will contain the local gain in density due to
!     process i for species k during this time step.
!..................................................................

      sgain(1,k)=0.
      sgain(2,k)=0.
      sgain(3,k)=0.
      sgain(4,k)=0.
      sgain(5,k)=0.
      sgain(6,k)=0.
      sgain(7,k)=0.
      sgain(8,k)=0.
      call bcast(tam5,zero,jx)
      call bcast(tam6,zero,jx)
      call bcast(tam7,zero,jx)
      call bcast(vflux(1:jx,k,l_),zero,jx)

!..................................................................
!     Collect contributions from various pieces of the Krook operator.
!..................................................................

      do 220 i=1,iy
        call bcast(tam8,zero,jx)
        do 210 j=1,jx
          tam5(j)=tam5(j)+vptb(i,lr_)*dtr*gon(i,j)* &
            temp2(i,j)*cynt2(i,l_)
          tam6(j)=tam6(j)-vptb(i,lr_)*dtr/taulos(i,j,indxlr_) &
            *temp2(i,j)*cynt2(i,l_)
          tam7(j)=tam7(j)-vptb(i,lr_)*dtr*tam8(j) &
            *temp2(i,j)*cynt2(i,l_)
 210    continue
 220  continue

!..................................................................
!     integrate over velocity.
!..................................................................

      do 230 j=1,jx
        sgain(5,k)=sgain(5,k)+tam5(j)*cint2(j)*one_
        sgain(6,k)=sgain(6,k)+tam6(j)*cint2(j)*one_
        sgain(7,k)=sgain(7,k)+tam7(j)*cint2(j)*one_
 230  continue

!..................................................................
!     Determine the fraction of particles leaving a sphere of
!     radius x(j) in a tauee(lr_) time. Also determine the density loss
!     due to particles leaving the domain at the high velocity end.
!..................................................................

      do 310 i=1,iy
        sgain(4,k)=sgain(4,k)+one_*gfu(i,jx,k)*cynt2(i,l_)*dtr
        do 311 j=1,jx
          vflux(j,k,l_)=vflux(j,k,l_)+one_*gfu(i,j,k)*cynt2(i,l_)
 311    continue
 310  continue
      do 320 j=1,jx
        vflux(j,k,l_)=vflux(j,k,l_)/xlndn(k,lr_)*tauee(lr_)
 320  continue

!..................................................................
!     Particle source term.
!..................................................................

      sgain(3,k)=xlncur(k,lr_)*dtr*0.5

!..................................................................
!     if ineg  .eq. "enabled" set negative values of f to zero.
!..................................................................

      if (n .eq. 0 .or. n/nchec*nchec .ne. n)  go to 90

!..................................................................
!     determine energy transfer diagnostics
!..................................................................

      entr(k,4,l_)=0.
      do 120 lefct=-1,8
        call diagentr(lefct,k)
 120  continue
      call diagentr(11,k)
      call diagentr(12,k)
      call coefstup(k)
      call coefmidv(da,1)
      call coefmidv(db,2)
      call coefmidv(dc,3)
      call coefmidt(dd,1)
      call coefmidt(de,2)
      call coefmidt(df,3)
 90   continue
      call dcopy(iyjx2,temp2(0:iy+1,0:jx+1),1,temp1(0:iy+1,0:jx+1),1)
      if (ineg .eq. "disabled") go to 400

      do 410 j=1,jx
        do 411 i=1,iy
          if(temp2(i,j) .lt. 0.) then
            temp1(i,j)=zero
            temp4(i,j)=-temp2(i,j)
          else ! temp2(i,j) .ge. 0.
            temp1(i,j)=temp2(i,j)
            temp4(i,j)=zero
          endif
 411    continue
 410  continue

!..................................................................
!     routine diagdens will compute density gained by setting negative
!     values of distribution to 0.
!..................................................................

      call diagdens(xline,xmidp,eline)
      engain(k)=eline*one_
      sgain(1,k)=xline*one_
 400  continue

!..................................................................
!     if debugging,compute density at end of velocity split
!..................................................................

      if (iactst .eq. "disabled") go to 500
      call dcopy(iyjx2,temp1(0:iy+1,0:jx+1),1,temp4(0:iy+1,0:jx+1),1)
      call diagdens(yline,ymidd,eline)
      yline=yline*one_
 500  continue
      return
      end subroutine diagxswx


end module diagxswx_mod
