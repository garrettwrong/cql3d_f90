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

module cfpsymt_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine cfpsymt
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Symmetrize about pi/2  (f.p. coefficients)
!..................................................................

      save


      do 1 k=1,ngen
        do 2 i=1,iyh
          ii=iy+1-i
          do 3 j=1,jx
            cal(ii,j,k,l_)=cal(i,j,k,l_)
            cbl(ii,j,k,l_)=cbl(i,j,k,l_)
            cfl(ii,j,k,l_)=cfl(i,j,k,l_)
            ccl(ii,j,k,l_)=-ccl(i,j,k,l_)
            cdl(ii,j,k,l_)=-cdl(i,j,k,l_)
            cel(ii,j,k,l_)=-cel(i,j,k,l_)
 3        continue
 2      continue
 1    continue
      return
      end subroutine cfpsymt

end module cfpsymt_mod
