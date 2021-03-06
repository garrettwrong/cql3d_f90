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

module coefwtj_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use coefmidv_mod, only : coefmidv

  !---END USE

!
!

contains

      subroutine coefwtj(k)
      use param_mod
      use cqlcomm_mod
      use advnce_mod !here: in coefwtj(). Maybe not needed?
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine computes w = dy(i,l_)*dd(i,j)/df(i,j) and then sets
!     g(w)=di(i,j,k,l_)=1./w - (1./(exp(w)-1.))
!     Note g(minus infinity)=1.; g(0)=.5; g(plus infinity)=0.
!     This is the Chang-Cooper algorithm (2-d form ala Karney)
!     This routine actually has modified the above by subtracting
!     off the RF contribution to the diffusion (see the coding).
!     This procedure, while not yet
!     justified theoretically, has kept the distribution from
!     going negative for strongly driven problems.
!..................................................................



!-YuP      call bcast(dj(1,0,k,l_),half,iyjxp1) ! could it be error?
      call bcast(dj(0:iy+1,0:jx,k,l_),half,(iy+2)*(jx+1))
         !Note:  dj(0:iy+1,0:jx,1:ngen,lrors)

      if (chang.ne."disabled") then
!..................................................................
!     RF diffusion alone is in dbb, total in db.
!..................................................................

        op=one+em12
        call coefmidv(dbb,2)


        do 10 j=1,jx

          do 20 i=1,iy
            temc1(i)=dx(j)*da(i,j)*op*db(i,j)/(op*db(i,j) &
              -dbb(i,j))**2
 20       continue

!.................................................................
!     Limit magnitudes so exp or 1/temc1 do not blow up
!.................................................................

          do 81 i=1,iy
            if(abs(temc1(i)).gt.em6) then
              temc3(i)=temc1(i)
            else
              temc3(i)=em6
            endif
 81       continue

          do 82 i=1,iy
            if(temc3(i).lt.sevenhun) then
              temc2(i)=temc3(i)
            else
              temc2(i)=sevenhun
            endif
 82       continue

          do 83 i=1,iy
            if(temc2(i).le.-sevenhun) then
              temc2(i)=-sevenhun
            endif
 83       continue

!.................................................................
!     Evaluate the Chang-Cooper weight
!.................................................................

          do 84 i=1,iy
            dj(i,j,k,l_)=1.d0/temc2(i)-(1.d0/(exp(temc2(i))-1.d0))
 84       continue

!...............................................................
!     Limit for small temc1 follows...
!...............................................................

          do 21 i=1,iy
            wsub=(3.d0+temc1(i))/(2.d0+temc1(i))/3.d0
            if(temc2(i).eq.em6) then
              dj(i,j,k,l_)=wsub
            endif
 21       continue

!..............................................................
!     Limit for large positive or negative temc1 follows
!..............................................................

          do 22 i=1,iy
            if (temc2(i).eq.sevenhun) then
              dj(i,j,k,l_)=1.d0/temc1(i)
            elseif (temc2(i).eq.-sevenhun) then
              dj(i,j,k,l_)=1.d0+1.d0/temc1(i)
            endif
 22       continue

 10     continue ! j=1,jx

      endif ! if (chang.ne."disabled")


!.......................................................................
!     Ensures correct differentiation at end of intervals
!.......................................................................

      do 30 i=1,iy
        dj(i,jx,k,l_)=1.d0
        dj(i,0,k,l_)=0.d0
 30   continue


      return
      end subroutine coefwtj


end module coefwtj_mod
