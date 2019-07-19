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

module wpwrng_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine wpwrng(kerr)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!...
!
!mnt  diagnostic error messages  + chkpnt if ichkpnt.ne."disabled"
!
!...

#ifdef __MPI
      include 'mpilib.h'
#endif

      character cptline*80
      integer getpid

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
      elseif (kerr.eq.11) then
        WRITE(*,10110)
      elseif (kerr.eq.12) then
        WRITE(*,10120)
      elseif (kerr.eq.13) then
        WRITE(*,10130)
      elseif (kerr.eq.14) then
        WRITE(*,10140)
      elseif (kerr.eq.15) then
        WRITE(*,10150)
      elseif (kerr.eq.16) then
        WRITE(*,10160)
      elseif (kerr.eq.17) then
        WRITE(*,10170)
      elseif (kerr.eq.99) then
        WRITE(*,10990)
      else
        WRITE(*,10990)
      endif

!..................................................................
!     Check-point this job, if ichkpnt.ne."disabled" into file ichkpnt.
!..................................................................
      if (ichkpnt.ne."disabled") then
!BH080118        i=getpid(0)
        WRITE(cptline,100) i,ichkpnt
 100    format("chkpnt -p ",i5," -f ",a8)
        PRINT *,cptline
!        i=dropfile(ichkpnt)
      endif
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

      stop 'wpwrng:' ! stop at all MPI cores

10010 format(" setup0%cqlpmod and colmodl=4 not allowed")
10020 format(" Cannot have periodic mesh along B and non full l mesh " &
        ,"for some pitch-angle")
10030 format(" lmidpln .ne. 1 not yet valid with CQLP")
10040 format(" Cannot run parallel transport with setup0%lsdiff=enabled")
10050 format(" error in micxinil in calculating imax")
10060 format(" error in micxinil in calculating ill and cosz")
10070 format(" updown .ne. symmetry with 9<numods<20")
10080 format(" sbdry .ne. periodic with 9<numods<20")
10090 format(" setup0%ls odd with 9<numods<20")
10100 format(" lz odd with 9<numods<20")
10110 format(" iymax .ne. iy_(1)")
10120 format(" l_lower(i).ne.1")
10130 format(" l_upper(i).ne.setup0%ls for i<itl_(1) or i>itu_(1)")
10140 format(" i.ne.iyh_(l_upper(i))")
10150 format(" wptrafx called with a non-constant y mesh")
10160 format(" numixts should equal +1 or -1 with numindx=2 option")
10170 format(" lmidvel should equal 0 with numindx=4 (centered scheme)")
10990 format(" unspecified in wpwrng")

      end subroutine wpwrng

end module wpwrng_mod
