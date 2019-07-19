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

module firstdrv_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine firstdrv(x,f,fprim,n)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      dimension x(n),f(n),fprim(n)

!..................................................................
!     Compute the first derivative of the function f(x,l_) at the x(l)
!     mesh points and put the solution into the array fprim(l), for
!     l=1,...,n     Taylor expansion used is:
!     f(x+h(i,l_))=f(x,l_) + h(i)*f'(x) + h(i)**2/2*f"(x) for i=1,2
!     the x mesh need not be evenly spaced.
!..................................................................

      do 10 l=2,n-1
        h1=x(l-1)-x(l)
        h2=x(l+1)-x(l)
        a=h2**2/h1**2
        fprim(l)=(f(l-1)*a - f(l+1) - f(l)*(a-1))/(a*h1-h2)
 10   continue

!..................................................................
!     End points extrapolate, but expansion is same as above.
!..................................................................

      h1=x(2)-x(1)
      h2=x(3)-x(1)
      a=h2**2/h1**2
      fprim(1)=(f(2)*a - f(3) - f(1)*(a-1))/(a*h1-h2)
      h1=x(n-1)-x(n)
      h2=x(n-2)-x(n)
      a=h2**2/h1**2
      fprim(n)=(f(n-1)*a - f(n-2) - f(n)*(a-1))/(a*h1-h2)
      return
      end subroutine firstdrv

end module firstdrv_mod
