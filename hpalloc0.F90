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

module hpalloc0_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine hpalloc0(iptr,ln)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!     NO_OP:  Not needed after 081216, after changing from
!             CRAY pointers to f90 pointers.

!..................................................................
!     Allocates memory AND zeroes out the array.
!..................................................................
!990131      integer err
!$$$CPTR>>>DELETE PTR-HPALLOC
!$$$      pointer(iptr,a(1))
!$$$
!$$$c_cray      call hpalloc(iptr,ln,err,1)
!$$$c_pc Unix libU77.a library.  malloc allocates bytes.
!$$$c_pc      iptr=malloc(8*ln)
!$$$c_pc      if (iptr.eq.0) stop 'Problem with malloc in hpalloc0'
!$$$      iptr=malloc(8*ln)
!$$$      if (iptr.eq.0) then
!$$$         write(*,*)'hpalloc0: ln=',ln
!$$$         stop 'Problem with malloc in hpalloc0'
!$$$      endif
!$$$
!$$$      zero=0.d0
!$$$      call bcast(a,zero,ln)
!$$$CPTR<<<END PTR-HPALLOC
      return
      end
end module hpalloc0_mod
