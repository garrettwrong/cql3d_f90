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

module impnorm_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use iso_c_binding, only : c_double

  !---END USE

contains

  subroutine impnorm(xnorm,a,rhs,nn)
    !implicit integer (i-n), real(c_double) (a-h,o-z)
    implicit none

    !..................................................................
    !     This routine normalizes the matrix a so that the maximum
    !     coefficient for each equation is of order 1.
    !..................................................................
    real(c_double) :: a(nn)
    real(c_double) :: xnorm
    real(c_double) :: rhs
    integer :: i
    integer :: nn


    xnorm=0.d0
    do i=1,nn
       xnorm=xnorm+dabs(a(i))
    end do
    if (xnorm.gt.0.d0) then
       rhs=rhs/xnorm
       do i=1,nn
          a(i)=a(i)/xnorm
       end do
    endif

    return
  end subroutine impnorm

end module impnorm_mod
