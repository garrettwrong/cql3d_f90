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

module tdsxray_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdsxr0_mod, only : tdsxr0

  !---END USE

!
!

contains

      subroutine tdsxray(icall,iplotsxr)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     sets up call to soft-x-ray analyzer
!..................................................................

      character*8 icall,iplotsxr

      do 1 l=1,setup0%lrzmax
        tr1(l)=reden(kelec,l)
 1    continue
!BH081106:  In some radcoord cases, rrz contains normalized radial
!BH081106:  coord data, and is not suitable for eqmod.ne."enabled"
!BH081106:  circ plasma model, or the eqmod.eq."enabled" eqdsk
!BH081106:  equilibria.
      if (eqmod.ne.'enabled') then
         if (radcoord.ne.'sqtorflx') &
            write(*,*)'tdsxray: WARNING, check our radial coord further'
         call tdsxr0(rrz,tr1(1),icall,iplotsxr)
      else
         call tdsxr0(rpmconz,tr1(1),icall,iplotsxr)
      endif
      return
      end subroutine tdsxray


end module tdsxray_mod
