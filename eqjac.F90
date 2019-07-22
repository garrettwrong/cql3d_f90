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

module eqjac_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use zcunix_mod, only : terp1
  use zcunix_mod, only : terp2

  !---END USE

!
!

contains

      subroutine eqjac(neq,t,yv,ml,mu,pd,nrowpd)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      dimension yv(neq),pd(nrowpd,neq)
!..................................................................
!     This routine is used with the O.D.E. solver LSODE. It
!     provides the Jacobian of the system of equations.
!..................................................................

!..................................................................
!     Provide the various derivatives of epsi (first and second
!     order).
!..................................................................

      dpsidr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,1,0)
      dpsidz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,0,1)
      d2psidrr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,2,0)
      d2psidrz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,1,1)
      d2psidzz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,0,2)

!..................................................................
!     Compute B and it's first order derivatives.
!..................................................................

      alpha=(dpsidr**2+dpsidz**2+(fpsi_)**2)
      sqrtalp=sqrt(alpha)
      dalphadr=2.*(dpsidz*d2psidrz+dpsidr*d2psidrr)
      dalphadz=2.*(dpsidz*d2psidzz+dpsidr*d2psidrz)

!..................................................................
!     Compute the 4 elements of the Jacobian.
!..................................................................

      pd(1,1)=(d2psidrz-dpsidz*dalphadr/alpha)/alpha
      pd(1,2)=(d2psidzz-dpsidz/alpha*dalphadz)/alpha
      pd(2,1)=(-d2psidrr+dpsidr/alpha*dalphadr)/alpha
      pd(2,2)=(-d2psidrz+dpsidr/alpha*dalphadz)/alpha
      return
      end subroutine eqjac

end module eqjac_mod
