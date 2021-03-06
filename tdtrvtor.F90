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

module tdtrvtor_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use r8subs_mod, only : dcopy
  use tdtrchkd_mod, only : tdtrchkd

  !---END USE

!
!

contains

      subroutine tdtrvtor(f1,f2)
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     This routine takes the distribution f1 defined on the full
!     velocity mesh and interpolates it onto the velocity mesh
!     on which the transport is done. The new distribution is f2
!     and it has the same density as f1 (particle conserving).
!..............................................................


      dimension f1(0:iyp1,0:jxp1,ngen,0:*)
      dimension f2(0:iyp1,0:jxp1,ngen,0:*)
      dimension denrad(ngen,lrorsa),denvel(ngen,lrorsa)

      call dcopy(iyjx2*ngen*(lrors+1),f1(0:iy+1,0:jx+1,1:ngen,0:lrors),1, &
                                      f2(0:iy+1,0:jx+1,1:ngen,0:lrors),1)

      call tdtrchkd(f1,vpint,denvel)

!BH070419:   removing special itl,itu treatment for ipacktp=0
!      write(*,*)'tdtrvtor: ipacktp', ipacktp
      if (ipacktp.eq.3) then

      do 10 k=1,ngen
        do 20 l=1,lrors
          ilr=setup0%lrindx(l)
          itl=itl_(l)
          itu=itu_(l)
          iyy=iy_(l) !-YuP-101215: Don't use iy=; it's in common /params/
                     ! Don't let overwrite the cqlinput value!
          do 30 j=1,jx
            do 40 i=itl-1,itl+1
              i_=iyy+1-i
              f2(i,j,k,l)=0.
              f2(i_,j,k,l)=0.
 40         continue

            if (f1(itl,j,k,l).ne.zero) then
              f_lm(j,k,l)=f1(itl-1,j,k,l)/f1(itl,j,k,l)
              f_lp(j,k,l)=f1(itl+1,j,k,l)/f1(itl,j,k,l)
              f_up(j,k,l)=f1(itu+1,j,k,l)/f1(itl,j,k,l)
            else
              f_lm(j,k,l)=1.0
              f_lp(j,k,l)=1.0
              f_up(j,k,l)=1.0
            endif

            f2(itl-2,j,k,l)=(f1(itl-2,j,k,l)*vpint(itl-2,ilr)+ &
              f1(itl-1,j,k,l)*vpint(itl-1,ilr)+ &
              f1(itl,j,k,l)*vpint(itl,ilr)*.5)/vpint_(itl-2,ilr)
            f2(itl+2,j,k,l)=(f1(itl+2,j,k,l)*vpint(itl+2,ilr)+ &
              f1(itl+1,j,k,l)*vpint(itl+1,ilr)+ &
              f1(itl,j,k,l)*vpint(itl,ilr)*.5)/vpint_(itl+2,ilr)
            f2(itu+2,j,k,l)=(f1(itu+2,j,k,l)*vpint(itu+2,ilr)+ &
              f1(itu+1,j,k,l)*vpint(itu+1,ilr)+ &
              f1(itu,j,k,l)*vpint(itu,ilr)*.5)/vpint_(itu+2,ilr)
            f2(itu-2,j,k,l)=f2(itl+2,j,k,l)
 30       continue
 20     continue
 10   continue

      elseif (ipacktp.ne.0) then
         if(setup0%verbose>0) write(*,*)'STOP in tdtrvtor:  Check ipacktp'
         stop
      endif

      call tdtrchkd(f2,vpint_,denrad)
      return
      end subroutine tdtrvtor

end module tdtrvtor_mod
