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

module tdnpadiag_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdnpa0_mod, only : tdnpa0

  !---END USE

!
!

contains

      subroutine tdnpadiag(icall)
        use cqlconf_mod, only : setup0
        use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      character*8 icall

!..................................................................
!     sets up call to NPA diagnostic, Version 1.0
!..................................................................

      character*8 iplotnbi

      write(*,*)
      write(*,*)'tdnpadiag, time step ',n


!     Call npa routines to calc and plot output....

      if (setup0%noplots.eq."enabled1") then
         iplotnbi='no'
      else
         iplotnbi='yes'
      endif

      do 1 l=1,setup0%lrzmax
        tr1(l)=reden(kelec,l)
 1    continue
      call tdnpa0(rrz,tr1(1),icall,iplotnbi)

      return
      end subroutine tdnpadiag


end module tdnpadiag_mod
