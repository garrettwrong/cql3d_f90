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

module tdpro_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      real(c_double) function tdpro(psi,rplasm,acoef)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!..................................................................
!     Calculate ASDEX YAG1 type profiles.
!     acoef(i),i=1,4  must be specified.
!..................................................................

      dimension acoef(4)

      x = rplasm*psi
      x2=x*x
      arg = acoef(4)
      do 10  i=1,3
        arg = arg*x2 + acoef(4-i)
 10   continue
      tdpro = exp(arg)
      return
      end function tdpro

end module tdpro_mod
