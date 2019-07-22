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

module bcast_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

! XXX every thing in this file and related to it should be deleted
! YuP: basically I agree. I guess the intent was to use a fast procedure with unrolled loops,
! similar to dcopy().
! Never accomplished.
! But now there are too many lines with "bcast" in the source.

!
!

contains

      subroutine bcast(a,val,n)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Temporary bcast routine until I can find UNICOS equivalent
!..................................................................

      dimension a(n)
      do 100 i=1,n
        a(i)=val
 100  continue
      return
      end subroutine bcast
!
!
      subroutine ibcast(ia,ival,n)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Temporary bcast routine until I can find UNICOS equivalent
!..................................................................

      dimension ia(n)
      do 100 i=1,n
        ia(i)=ival
 100  continue
      return
      end subroutine ibcast

! NME bcast routine for complex arrays
      subroutine ccast(c,cval,n)
      implicit integer (i-n), complex*16 (c)
      dimension c(n)
      do 100 i=1,n
         c(i)=cval
 100  continue
      return
      end subroutine ccast

end module bcast_mod
