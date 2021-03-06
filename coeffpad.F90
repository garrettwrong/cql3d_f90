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

module coeffpad_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine coeffpad(k)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!      real(c_double),dimension(iy):: prnt1,prnt2,prnt3
!      real(c_double),dimension(iy):: prnt4,prnt5,prnt6

!..................................................................
!     This routine adds in the collisional contribution to the
!     coefficients employed in time advancement..
!     scatfrac=0. disables pitch angle scattering (along with
!       mx=0, see cqlinput_help. scatfrac=1. by default).
!..................................................................


      do 10 j=1,jx
        do 11 i=1,iy
          da(i,j)=da(i,j)+cal(i,j,k,l_)
          db(i,j)=db(i,j)+cbl(i,j,k,l_)
          dc(i,j)=dc(i,j)+scatfrac*ccl(i,j,k,l_)
          dd(i,j)=dd(i,j)+cdl(i,j,k,l_)
          de(i,j)=de(i,j)+scatfrac*cel(i,j,k,l_)
          df(i,j)=df(i,j)+scatfrac*cfl(i,j,k,l_)
 11     continue
 10   continue
!      do i=1,iy
!         prnt1(i)=da(i,2)
!         prnt2(i)=db(i,2)
!         prnt3(i)=dc(i,2)
!         prnt4(i)=dd(i,2)
!         prnt5(i)=de(i,2)
!         prnt6(i)=df(i,2)
!      enddo


      return
      end subroutine coeffpad

end module coeffpad_mod
