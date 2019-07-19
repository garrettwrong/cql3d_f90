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

module tdxin23d_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use profaxis_mod, only : profaxis

  !---END USE

!
!

contains

      subroutine tdxin23d(a,rya,klrz,ngn,nso,k,kk,expn1,expm1)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Fill in input arrays between center and edge parabolically.
!..................................................................


      dimension rya(0:klrz), a(ngn,nso,0:klrz)
      data em90 /1.d-90/
      if (abs(a(k,kk,0)) .le. em90) a(k,kk,0)=em90
      dratio=a(k,kk,1)/a(k,kk,0)
      do 1 ll=1,klrz
        call profaxis(rn,expn1,expm1,dratio,rya(ll))
        a(k,kk,ll)=a(k,kk,0)*rn
 1    continue
      return
      end subroutine tdxin23d


end module tdxin23d_mod
