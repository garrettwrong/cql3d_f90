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

module eqalloc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine eqalloc
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!...................................................................
!     Allocate allocatable arrays
!...................................................................

!dir$ nobounds

!..................................................................
!     A check on allocations is sucessful entering then exiting
!     the subroutine.
!..................................................................
      write(*,*)'eqalloc:  Entering eqalloc'

      lnlfield=lfielda*setup0%lrzmax
      lndumeq=4*lnlfield
      allocate(drpmconz(setup0%lrzmax),STAT=istat)
      call bcast(drpmconz,zero,SIZE(drpmconz))
      allocate(eqdell(lfielda,setup0%lrzmax),STAT=istat)
      call bcast(eqdell,zero,SIZE(eqdell))
      allocate(eqbpol(lfielda,setup0%lrzmax),STAT=istat)
      call bcast(eqbpol,zero,SIZE(eqbpol))
      allocate(solr(lfielda,setup0%lrzmax),STAT=istat)
      call bcast(solr,zero,SIZE(solr))
      allocate(solz(lfielda,setup0%lrzmax),STAT=istat)
      call bcast(solz,zero,SIZE(solz))

      write(*,*)'eqalloc:  Leaving eqalloc'

      return
      end subroutine eqalloc


end module eqalloc_mod
