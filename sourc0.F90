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

module sourc0_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine sourc0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     define source uniquely at v=0
!..................................................................


      do 10 k=1,ngen
        s=0.
        u=0.
        do 11 i=1,iy
          u=u+source(i,1,k,indxlr_)*cynt2(i,l_)*vptb(i,lr_)
          s=s+cynt2(i,l_)*vptb(i,lr_)
 11     continue
        do 12 i=1,iy
          source(i,1,k,indxlr_)=u/s
 12     continue
 10   continue
      return
      end subroutine sourc0

end module sourc0_mod
