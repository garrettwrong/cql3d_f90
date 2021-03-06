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

module urfedge_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfedge(i1,i2,c1,c2,lp,lrad,jval)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..........................................................
!     Computes weights for diffusion coefficients if spectral
!     width of ray element is small.
!..........................................................


      if (c2 .lt. -1.) c2=-1.
      if (c2 .gt. 1.) c2=1.
      if (c1 .lt. -1.) c1=-1.
      if (c1 .gt. 1.)  c1=1.

      if (i1.ne.i2) then

        if (i1.ne.1) then
          fact1=abs(cosmz(i1-1,lp,lrad)-cosmz(i1,lp,lrad))
          del1=abs(c1-cosmz(i1-1,lp,lrad))
        else           !Not possible
          fact1=abs(-1.-cosmz(i1,lp,lrad))
          del1=abs(-1.-c1)
        endif
        if (i2.ne.1) then
          del2=abs(c2-cosmz(i2,lp,lrad))
          fact2=abs(cosmz(i2,lp,lrad)-cosmz(i2-1,lp,lrad))
        else
!*bh*940227del2=abs(c2-cosmz(1,lp,lrad))
!*bh*940227fact2=abs(-1.-cosmz(1,lp,lrad))!closetozero.
          del2=1.0
          fact2=1.0
        endif
        scal1=del1/fact1
        scal2=del2/fact2
      else
        del1=abs(c2-c1)
        if (i1.ne.1) then
          fact1=abs(cosmz(i1-1,lp,lrad)-cosmz(i1,lp,lrad))
        else             !Not possible
          fact1=abs(-1.-cosmz(i1,lp,lrad))
        endif
        scal1=del1/fact1
        scal2=scal1
      endif
      ifct1(jval)=65535*scal1+.5
      ifct2(jval)=65535*scal2+.5
      return
      end subroutine urfedge

end module urfedge_mod
