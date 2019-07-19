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

module aminmx_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE


contains

      subroutine aminmx(array,ifirst,ilast,istride,amin,amax, &
        indmin,indmax)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!     compute max and min with indices

      dimension array(ilast)
!
      amin = array(ifirst)
      amax = array(ifirst)
      indmin = ifirst
      indmax = ifirst
      do i=ifirst+istride,ilast,istride
        if (array(i) .lt. amin) then
          amin = array(i)
          indmin = i
        end if
        if (array(i) .gt. amax) then
          amax = array(i)
          indmax = i
        end if
      end do
!
      return
      end subroutine aminmx

end module aminmx_mod
