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

module pltcycl_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine pltcycl(iymn,iymx,ymn,ymx)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!-------------------------
!  find the minimum and maximum log-cycle for ymn,ymx
!-----------------------------
      iymn=0
      iymx=0
      imn=0
      imx=0
10    imn1=imn+1
      imn2=imn-1
      if (ymn.lt.10.**imn) then
        iymn=imn2
        imn=imn2
        if (ymn.ge.10.**imn2) then
          go to 20
        else
          go to 10
        endif
      endif
      if (ymn.ge.10.**imn) then
        if (ymn.ge.10.**imn1) then
          iymn=imn1
          imn=imn1
          go to 10
        else
          go to 20
        endif
      endif
20    imx1=imx+1
      imx2=imx-1
      if (ymx.gt.10.**imx) then
        iymx=imx1
        imx=imx1
        if (ymx.le.10.**imx1) then
          go to 30
        else
          go to 20
        endif
      endif
      if (ymx.le.10.**imx) then
        if (ymx.le.10.**imx2) then
          iymx=imx2
          imx=imx2
          go to 20
        else
          go to 30
        endif
      endif
30    continue
      return
      end subroutine pltcycl

end module pltcycl_mod
