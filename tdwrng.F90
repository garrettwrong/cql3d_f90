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

module tdwrng_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE



!-----------------------------------------------------------------

contains

      subroutine tdwrng(kerr)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!...
!
!mnt  diagnostic error messages
!
!...

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
      elseif (kerr.eq.7) then
        WRITE(*,10070)
      elseif (kerr.eq.8) then
        WRITE(*,10080)
      elseif (kerr.eq.9) then
        WRITE(*,10090)
      elseif (kerr.eq.10) then
        WRITE(*,10100)
      else
        WRITE(*,10990)
      endif
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

      stop ! stop at all MPI cores

!      if (kerr.eq.0) return
10010 format("machine = toroidal only  for 3-d calc. in cqlinput")
10020 format("wdmodel =ech1 only for 3-d calc. in cqlinput")
10030 format("rya(ll).lt.rya(ll-1) for an ll; rzset error (cqlinput)")
10040 format ("subroutine tdsxr - thpol out of bounds")
10050 format("subroutine tdsxr - field line does not intersect viewing", &
        " angle")
10060 format("subroutine tdtodskr - hpdealc returns error message")
10070 format("subroutine ???? - iostatus returns error message")
10080 format("transport model- meshy=fixed_mu.  Insufficient resolution" &
        "in the trapped region for some l. Set tfac larger")
10090 format("transport model - transport mesh falls EXACTLTY on p/t",/, &
        "boundary. Perturb tfac")
10100 format("STOP in tdtry: iactst=abort forces stop. error>1.e-8")
10990 format("unspecified")
      return
      end subroutine tdwrng

end module tdwrng_mod
