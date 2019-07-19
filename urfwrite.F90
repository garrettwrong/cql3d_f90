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

module urfwrite_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use urfwrite__mod, only : urfwrite_
  !---END USE

!
!

contains

      subroutine urfwrite
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine performs the converse operation on urfread:
!     i.e., a disk file "rayop" in written or updated with ray tracing
!     data.
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

#ifdef __MPI
      if(mpirank.ne.0) return
#endif

      krf=0
      if (lh.eq."enabled") then
        krf=krf+1
        open(unit=20,file='raylh',status='old')
        call urfwrite_(krf,20)
      endif
      if (ech.eq."enabled") then
        krf=krf+1
        open(unit=23,file='rayech',status='old')
        call urfwrite_(krf,23)
      endif
      if (fw.eq."enabled") then
        krf=krf+1
        open(unit=24,file='rayfw',status='old')
        call urfwrite_(krf,24)
      endif
      return
      end subroutine urfwrite


end module urfwrite_mod
