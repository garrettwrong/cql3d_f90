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

module urfwrong_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfwrong(kerr)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)


!..................................................................
!     Flags errors in urf model - and terminates execution.
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

! print error messages - on mpirank.eq.0 only
#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      if (kerr.eq.1) then
        WRITE(*,10010)
      elseif (kerr.eq.2) then
        WRITE(*,10020)
      elseif (kerr.eq.3) then
        WRITE(*,10030)
      elseif (kerr.eq.4) then
        WRITE(*,10040)
      elseif (kerr.eq.5) then
        WRITE(*,10050)
      elseif (kerr.eq.6) then
        WRITE(*,10060)
      elseif(kerr.eq.7) then
        WRITE(*,10070)
      else
        WRITE(*,10990)
      endif
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

      stop 'urfwrong:' ! stop at all MPI cores

10010 format("subroutine urfbes: insufficient elements -  Bessel table")
10020 format("harmonic number...nharm(k) must be .le. nharma")
10030 format("subroutine urfb0 - unable to determine poloidal angle.")
10040 format("subroutine urfb0 - cannot locate the resonance region", &
        "for v-perp=0")
10050 format(" ")
10060 format("subroutine urfread: too many ray elements provided by", &
        /, "ray tracing code for some ray. Not enough RAM memory")
10070 format("failure in urfpack, can't find resonance region..")
10990 format("unspecified error originating from urf module.")
      return
      end subroutine urfwrong

end module urfwrong_mod
