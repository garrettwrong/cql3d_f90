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

module coefefad_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine coefefad(k)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine adds in the contribution of the D.C. electric
!     field to the coefficients used for time advancement.
!..................................................................


!.......................................................................
!l    Include electrostatic field
!.......................................................................

      do 20 i=1,iy
        ztrda=-elparnw(l_)*charge/vnorm*coss(i,l_)
        ztrdd=elparnw(l_)*charge/vnorm*sinn(i,l_)**2
        do 21 j=1,jx
          da(i,j)=da(i,j)+bnumb(k)/fmass(k)*(elecfld(lr_)*cex(i,j,l_)+ &
            ztrda*xsq(j))
          dd(i,j)=dd(i,j)+bnumb(k)/fmass(k)*(elecfld(lr_)*cet(i,j,l_)+ &
            ztrdd*x(j))
 21     continue
 20   continue

      return
      end subroutine coefefad

end module coefefad_mod
