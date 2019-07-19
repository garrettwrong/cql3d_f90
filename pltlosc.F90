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

module pltlosc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use pltdf_mod, only : pltcont

  !---END USE

!
!

contains

  subroutine pltlosc
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!     Plot contours of the loss region..
!
!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!

      if (setup0%noplots.eq."enabled1") return

      do 100 k=1,ngen
        suu=0.
        do 92001 i=1,iy
          do 92002 j=1,jx
            if(gone(i,j,k,indxlr_).lt.-.9) then
              temp1(i,j)=vnorm*x(j)*f(i,j,k,l_)/tau(i,lr_)
            elseif (gone(i,j,k,indxlr_).gt.0.) then
              temp1(i,j)=0.
            elseif (gone(i,j,k,indxlr_) .le. 0.) then
              temp1(i,j)=vnorm*x(j)*f(i,j,k,l_) &
                *(-gone(i,j,k,indxlr_))/tau(i,lr_)
            else
              temp1(i,j)=gone(i,j,k,indxlr_)*f(i,j,k,l_)
            endif
            if (temp1(i,j).ne.zero) suu=temp1(i,j)
92002     continue
92001   continue
        if (suu.eq.0.) go to 92003

        write(t_,588) k
 588    format("Loss due to lossmode(k) and torloss(k), k=",i5)
#ifndef NOPGPLOT
        CALL PGPAGE
#endif
        call pltcont(k,1,t_,8) ! itype=8 for pltlosc
        !call GSCPVS(.5,.4)
!$$$        call gxglfr(0)
!$$$        call gscpvs(.15,.85)
!$$$        write(t_,560)
!$$$        call gptx2d(t_)
!$$$ 560    format("Contour values:")
!$$$        write(t_,570) (temp2(jc,1),jc=1,(ncont/2)*2)
!$$$        call gptx2d(t_)
!$$$ 570    format((1x,e16.6,5x,e16.6),"$")
92003   continue
 100  continue
      return
      end subroutine pltlosc


end module pltlosc_mod
