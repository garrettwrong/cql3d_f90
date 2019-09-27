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

module tdtrrtov_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdtrchkd_mod, only : tdtrchkd

  !---END USE

!
!

contains

      subroutine tdtrrtov(f1)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     This routine interpolates onto the full velocity mesh
!     from the transport velocity mesh.
!     f1 ======> f1,    THAT IS, just fix up itl+/-1,itu+/-1
!..............................................................

      dimension f1(0:iyp1,0:jxp1,ngen,0:*),denrad(ngen,lrorsa)
      dimension denvel(ngen,lrorsa)

      call tdtrchkd(f1,vpint_,denrad)

!BH070419:   removing special itl,itu treatment for ipacktp=0
        if (ipacktp.eq.3) then

      do 10 k=1,ngen
        do 11 l=1,lrors
          ilr=setup0%lrindx(l)
          itl=itl_(l)
          itu=itu_(l)
          do 12 j=1,jx
            fact1=(vpint_(itl-2,ilr)-vpint(itl-2,ilr))*f1(itl-2,j,k,l) &
              +2.*(vpint_(itl+2,ilr)-vpint(itl+2,ilr))*f1(itl+2,j,k,l) &
              +(vpint_(itu+2,ilr)-vpint(itu+2,ilr))*f1(itu+2,j,k,l)
            fact2=vpint(itl-1,ilr)*f_lm(j,k,l)+ &
              2.*vpint(itl,ilr)+2.*vpint(itl+1,ilr)* &
              f_lp(j,k,l)+vpint(itu+1,ilr)*f_up(j,k,l)
            f1(itl,j,k,l)=fact1/fact2
            f1(itu,j,k,l)=f1(itl,j,k,l)
            f1(itl-1,j,k,l)=f_lm(j,k,l)*f1(itl,j,k,l)
            f1(itu+1,j,k,l)=f_up(j,k,l)*f1(itl,j,k,l)
            f1(itl+1,j,k,l)=f_lp(j,k,l)*f1(itl,j,k,l)
            f1(itu-1,j,k,l)=f1(itl+1,j,k,l)
 12       continue
 11     continue
 10   continue

        elseif (ipacktp.ne.0) then
           write(*,*)'STOP in tdtrrtov:  Check ipacktp'
           stop
        endif

      call tdtrchkd(f1,vpint,denvel)
      return
      end subroutine tdtrrtov

end module tdtrrtov_mod
