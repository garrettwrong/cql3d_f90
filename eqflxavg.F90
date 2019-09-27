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

module eqflxavg_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use eqorbit_mod, only : eqorbit

  !---END USE

!
!

contains

      subroutine eqflxavg(epsicon_,a,flxavg_,flxavgd_)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      dimension a(*)
!dir$ nobounds
!..................................................................
!     This routine returns the flux surface average of a.
!     It works for both updown and non-up-down symmetric cases.
!..................................................................

      if (eqorb.eq."enabled") then
        call eqorbit(epsicon_)
      endif
      sum1=0.
      sum2=0.

!..................................................................
!     eqdell=dl; eqbpol=B-poloidal; both defined on constant phi
!     flux surface.
!..................................................................
      do 10 l=2,lorbit_
        sum1=sum1+eqdell_(l)/eqbpol_(l)
        sum2=sum2+(a(l)+a(l-1))*.5*eqdell_(l)/eqbpol_(l)
 10   continue
      flxavg_=sum2/sum1
      flxavgd_=sum1
      return
      end subroutine eqflxavg


end module eqflxavg_mod
