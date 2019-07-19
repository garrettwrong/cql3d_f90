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

module diagxswt_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use diagdens_mod, only : diagdens
  use diagentr_mod, only : diagentr
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine diagxswt(k)
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Routine computes density gain due to particle sources and losses
!     during the during theta sweep in exsweep.
!     gains are positive - losses are negative.
!     2- density gained due to setting negative values of
!     distribution function to 0 during theta sweep.
!     3- particle source contribution.
!..................................................................



!..................................................................
!     Add in the source contribution.
!..................................................................

      sgain(3,k)=xlncur(k,lr_)*.5*dtr+sgain(3,k)
      call dcopy(iyjx2,temp2(0:iy+1,0:jx+1),1,temp1(0:iy+1,0:jx+1),1)
      s=0.
      if (ineg .eq. "disabled") go to 350

!..................................................................
!     if desired set negative values of distribution function = 0.
!..................................................................

      do 300 j=1,jx
        do 301 i=1,iy
          if(temp2(i,j) .lt. 0.) then
            temp1(i,j)=zero
            temp4(i,j)=-temp2(i,j)
          else ! temp2(i,j) .ge. 0.
            temp1(i,j)=temp2(i,j)
            temp4(i,j)=zero
          endif
 301    continue
 300  continue
      call diagdens(xline,xmidp,eline)
      engain(k)=engain(k)+eline*one_
      sgain(2,k)=xline*one_
 350  continue
      call dcopy(iyjx2,temp1(0:iy+1,0:jx+1),1,temp4(0:iy+1,0:jx+1),1)

!..................................................................
!     Compute power from df/dt and from setting neg f to zero.
!..................................................................

      if (n .gt. 0 .and. n/nchec*nchec .eq. n) then
        call diagentr(9,k)
        call diagentr(10,k)
        call dcopy(iyjx2,temp4(0:iy+1,0:jx+1),1,temp1(0:iy+1,0:jx+1),1)
      endif
      if (iactst.eq."disabled") go to 500

!..................................................................
!     if debugging, compute density at end of theta split
!..................................................................

      call diagdens(yline,ymidd,eline)
      yline=yline*one_
 500  continue
      return
      end subroutine diagxswt


end module diagxswt_mod
