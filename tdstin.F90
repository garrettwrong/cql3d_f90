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

module tdstin_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine tdstin
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!


!..................................................................
!     This routine initializes flux surface dependent data so the
!     2-d code initialization routine ainitial can be called.
!..................................................................

      do 2 k=1,ngen
        do 3 m=1,nso
          sellm1(k,m)=sellm1z(k,m,lr_)
          seppm1(k,m)=seppm1z(k,m,lr_)
          sellm2(k,m)=sellm2z(k,m,lr_)
          seppm2(k,m)=seppm2z(k,m,lr_)
          sem1(k,m)=sem1z(k,m,lr_)
          sem2(k,m)=sem2z(k,m,lr_)
          sthm1(k,m)=sthm1z(k,m,lr_)
          scm2(k,m)=scm2z(k,m,lr_)
          szm1(k,m)=szm1z(k,m,lr_)
          szm2(k,m)=szm2z(k,m,lr_)
          asor(k,m,lr_)=asorz(k,m,lr_)
 3      continue
 2    continue

      return
      end subroutine tdstin

end module tdstin_mod
