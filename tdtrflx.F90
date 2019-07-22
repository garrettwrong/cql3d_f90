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

module tdtrflx_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE


!
!

contains

      subroutine tdtrflx
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : cvmgt

      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     This routine checks the accuracy of the radial time advancement
!     by integrating the flux at the outside. This should equal
!     the change in the number of particles in the plasma over
!     the last time step.
!..............................................................


      include 'trans.h'

      call bcast(fxsp,zero,iyjx2*ngen*lrors)  !BH080428, why?
      flxout=0.
      if (soln_method.ne.'it3drv') then
      do 10 k=1,ngen
        do 20 j=1,jx
          do 30 i=1,iytr(lrors)
            l=lrors-1
            id=idx(i,l)
            if (l.ne.lpt(i)) then
              flxout=cosovb(id,l)*(h_r(l)*bovcos(id,l) &
                *sfu(i,j,k,l))*cynt2_(id,l)*cint2(j)*4.*pi**2*radmaj &
                +flxout
            else if (l.eq.lpt(i)) then
              id=idx(i,l)
              i_=iytr(lrors)+1-i
              i_d=idx(i_,l)
              flxout=flxout+(.5*cosovb(id,l)*h_r(l)*bovcos(id,l) &
                *sfu(i,j,k,l)+ &
                .5*cosovb(i_d,l)*h_r(l)*bovcos(i_d,l)*sfu(i_,j,k,l)) &
                *cynt2_(i_d,l)*cint2(j)*4.*pi**2 &
                *radmaj
            endif
 30       continue
 20     continue
 10   continue

      else   !Using sfup rather than sfu [BH:probably can make coding
             !more efficient].

      do k=1,ngen
        do j=1,jx
          do i=1,iytr(lrors)
            l=lrors-1
            id=idx(i,l)
            if (l.ne.lpt(i)) then
              flxout=cosovb(id,l)*(h_r(l)*bovcos(id,l) &
                *sfup(i,j,k,l))*cynt2_(id,l)*cint2(j)*4.*pi**2*radmaj &
                +flxout
            else if (l.eq.lpt(i)) then
              id=idx(i,l)
              i_=iytr(lrors)+1-i
              i_d=idx(i_,l)
              flxout=flxout+(.5*cosovb(id,l)*h_r(l)*bovcos(id,l) &
                *sfup(i,j,k,l)+ &
                .5*cosovb(i_d,l)*h_r(l)*bovcos(i_d,l)*sfup(i_,j,k,l)) &
                *cynt2_(i_d,l)*cint2(j)*4.*pi**2 &
                *radmaj
            endif
          enddo
        enddo
      enddo
      endif

      flxout=flxout*dttr
      return
      end subroutine tdtrflx


end module tdtrflx_mod
