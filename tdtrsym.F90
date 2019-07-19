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

module tdtrsym_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine tdtrsym
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      do 5 k=1,ngen
        do 10 l=1,lrors-1
          do 20 i=itl_(l),iyh_(l)
            i_=iy_(l)+1-i
            do 30 j=1,jx
              frn(i,j,k,l)=(frn(i,j,k,l)+frn(i_,j,k,l))*.5
              frn(i_,j,k,l)=frn(i,j,k,l)
 30         continue
 20       continue
 10     continue
 5    continue
      return
      end subroutine tdtrsym

end module tdtrsym_mod
