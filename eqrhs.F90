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

module eqrhs_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use zcunix_mod, only : terp1
  use zcunix_mod, only : terp2

  !---END USE

!
!

contains

      subroutine eqrhs(neq,t,yv,ydot)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      dimension yv(2), ydot(2)
!..................................................................
!     This routine is used with O.D.E. solver LSODE. It provides
!     the r.h.s. of the set of coupled O.D.E.'s .
!..................................................................

      dpsidr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,1,0)
      dpsidz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,0,1)
      rbval=sqrt(dpsidr**2+dpsidz**2+(fpsi_)**2)
!
!     Calc Br and Bz for integration along the  poloidal flux surface.
!     But, omit cursign factor, so start in pos-Z dirn.
!
      ydot(1)=dpsidz/(rbval)
      ydot(2)=-dpsidr/(rbval)
      return
      end subroutine eqrhs

end module eqrhs_mod
