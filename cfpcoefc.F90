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

module cfpcoefc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use cfpcoefn_mod, only : cfpcoefn
  use cfpcoefr_mod, only : cfpcoefr

  !---END USE

!
!

contains

      subroutine cfpcoefc
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!..................................................................
!     Subroutine to handle coefc subroutines
!..................................................................
!
!
      if(relativ.eq."fully") then
        call cfpcoefr
      else
        call cfpcoefn
      endif
!
      return
      end subroutine cfpcoefc


end module cfpcoefc_mod
