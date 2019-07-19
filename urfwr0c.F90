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

module urfwr0c_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfwr0c(x,nrayelt,nray,nrayelts)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      complex*16 x,c16i
      dimension x(nrayelts,*),nrayelt(*)
!....................................................................
!     For formatted o/p purposes, values of abs(x) .lt. 1.e-300
!     are set equal to 0. (Actually, .gt.1.e-999 OK for format,
!     but some fortrans restricted to .lt.1.e-37).
!....................................................................

      c16i=(0.d0,1.d0)
      deps=1.d-300
      zero=0.d0
      do 10  iray=1,nray
        do 11  is=1,nrayelt(iray)

           z16r=0.5*(x(is,iray)+conjg(x(is,iray)))
           z16i=-c16i*0.5*(x(is,iray)-conjg(x(is,iray)))
           if (abs(z16r).lt.deps) &
                x(is,iray)=cmplx(zero, z16i)

           if (abs(z16i).lt.deps) &
                x(is,iray)=cmplx(z16r, zero)

           if (abs(z16r).lt.deps .and. abs(z16i).lt.deps) &
                x(is,iray)=cmplx(zero, zero)

 11     continue
 10   continue

      return
      end subroutine urfwr0c

end module urfwr0c_mod
