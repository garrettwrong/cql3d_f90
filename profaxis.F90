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

module profaxis_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine profaxis(rn,expn1,expm1,dratio,rova)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!---------------------------------------------------------------------
!     Expands "parabolic" profiles by computing the ratio rn
!     of the local parameter value at normalized radius rova to the central
!     value, given exponents expn1 and expm1 of the parabola, and
!     the ratio dratio of edge to central parameter value.
!---------------------------------------------------------------------

      if(abs(dratio-1.).lt.1.e-12) then
        rn=1.0
      else
        rn=dratio+(1.0-dratio)*(1.-rova**expn1)**expm1
      endif
      return
      end subroutine profaxis

end module profaxis_mod
