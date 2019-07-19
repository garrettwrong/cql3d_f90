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

module wpcthta_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine wpcthta
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine calculates the coefficient cthta(i,j) encounting for
!     the mu*grad_parallel B force, when setup0%cqlpmod=enabled
!..................................................................

!.......................................................................

      do 100 j=1,jx
        if (mod(nummods,10).le.4 .or. lmidvel.eq.0) then
          ztra=-0.5*x(j)*vnorm*psisp(l_)/psis(l_)
        else
          ztra=-0.5*x(j)*vnorm*(psisp(l_)+psisp(l_+1)) &
            /(psis(l_)+psis(l_+1))
        endif
        do 110 i=1,iy
          cthta(i,j)=ztra*sinn(i,l_)*dyi(i,l_)
 110    continue
 100  continue

      return
      end subroutine wpcthta

end module wpcthta_mod
