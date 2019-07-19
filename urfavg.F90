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

module urfavg_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfavg
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      alpha=1.
      if (n.gt.2) then
        alpha3=.35
        alpha=(1.-alpha3)/(nstop)*(n-3)+alpha3
      endif

      if (n.eq.3) then
!cc        call dcopy(iyjx2*ngen*lrors,f,1,g_,1)
      do l=1,lrors
         do k=1,ngen
            do j=0,jxp1
               do i=0,iyp1
                  g_(i,j,k,l)= f(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      endif

      do k=1,ngen
         do 10 l=1,lrors
            do 20 j=1,jx
               do 30 i=1,iy
                  g_(i,j,k,l)=alpha*f(i,j,k,l)+(1.-alpha)*g_(i,j,k,l)
 30            continue
 20         continue
 10      continue
      enddo

      return
      end subroutine urfavg

end module urfavg_mod
