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

module tdtrcon_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdtrflx_mod, only : tdtrflx

  !---END USE

!
!

contains

      subroutine tdtrcon
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!..............................................................
!     Compute conservation constant.
!..............................................................
      data sgainr /0.0/
! jakub urban 110708: commented out for g95 compiler
! the above blanket save statement should do the job
!      save sgainr

!..............................................................
!     Compute original number of particles in tokamak.
!..............................................................
      if (n.eq.1) then
        total0=0.
        do 10 l=1,setup0%lrz
          ilr=setup0%lrindx(l)
          total0=xlndn0(ilr)/zmaxpsi(ilr)*dvol(ilr)+total0
 10     continue
      endif

!..............................................................
!     Compute total number of particles in device now.
!...............................................................

      total=0.
      do 30 l=1,setup0%lrz
        ilr=setup0%lrindx(l)
        do 40 k=1,ngen
          total=xlndn(k,ilr)/zmaxpsi(ilr)*dvol(ilr)+total
 40     continue
 30   continue

!................................................................
!     Call routine to compute number of particles lost at limiter
!     this time step.
!................................................................

      call tdtrflx
      sgainr=sgainr+flxout
      conserv=(total-total0-sgainr)/(total*.5+total0*.5)
      return
      end subroutine tdtrcon


end module tdtrcon_mod
