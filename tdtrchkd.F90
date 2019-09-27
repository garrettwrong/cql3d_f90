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

module tdtrchkd_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine tdtrchkd(f1,vp,densty)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.............................................................
!     This routine computes the density of species f1, given
!     that vpint=vp(setup0%lrindx). This will work for arrays defined on
!     radial as well as the velocity mesh.
!.............................................................

      dimension f1(0:iyp1,0:jxp1,ngen,0:*),vp(iy,lrza)
      dimension densty(ngen,lrza)
      if (iactst.eq."disabled") return

      do 10 k=1,ngen
        do 15 l=1,lrors
          call bcast(tam1,zero,jx)
          do 20 i=1,iy_(l)
            do 30 j=1,jx
              tam1(j)=tam1(j)+f1(i,j,k,l)*vp(i,setup0%lrindx(l))
 30         continue
 20       continue
          densty(k,l)=0.
          do 40 j=1,jx
            densty(k,l)=densty(k,l)+tam1(j)*cint2(j)
 40       continue
 15     continue
 10   continue
      return
      end subroutine tdtrchkd

end module tdtrchkd_mod
