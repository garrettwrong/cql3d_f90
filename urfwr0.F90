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

module urfwr0_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfwr0(x,nrayelt,nray,nrayelts)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      dimension x(nrayelts,*),nrayelt(*)
      save em300
      data em300 /1.d-300/
!....................................................................
!     For formatted o/p purposes, values of abs(x) .lt. 1.e-300
!     are set equal to 0. (Formats ok to 1.e-999, but real(c_double)
!     restricted to .gt. 1.e-327 in some fortrans).
!....................................................................

      do 10  iray=1,nray
        do 11  is=1,nrayelt(iray)
          if (abs(x(is,iray)).lt.em300)  x(is,iray)=0.0
 11     continue
 10   continue

      return
      end subroutine urfwr0

end module urfwr0_mod
