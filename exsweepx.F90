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

module exsweepx_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use diagentr_mod, only : gfu
  use diagwrng_mod, only : diagwrng

  !---END USE

!
!

contains

      subroutine exsweepx(k)
      use param_mod
      use cqlcomm_mod
      use advnce_mod !here: in exsweepx(). To get gfu()
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine checks the solution obtained after the theta split
!     by differentiation of the solution in time and comparison with
!     the right hand side. temp4(i,j) will contain the l.h.s.and will b
!     multiplied by cint2(j)*cynt2(i,l_) so that it will contain the
!     the change in the local particle number density. temp4(i,j) will
!     be the r.h.s.  and the two arrays should (ideally) be identical.
!     Each array is summed for purposes of comparison.
!..................................................................


      sumleft=0.
      sumright=0.
      do 90 i=1,iy
        do 100 j=1,jx

!..................................................................
!     Differentiate the flux gfu
!..................................................................

          temp5(i,j)=gfu(i,j,k)
          temp3(i,j)=(gfu(i,j,k)-gfu(i,j-1,k))*cynt2(i,l_) &
!..................................................................
!     Krook operator..
!..................................................................
            +cah(i,j)*temp2(i,j)*vptb(i,lr_)*cint2(j)*cynt2(i,l_) &
!..................................................................
!     Particle source...
!..................................................................
            +.5*so(i,j)*vptb(i,lr_)*cint2(j)*cynt2(i,l_)

          temp3(i,j)=temp3(i,j)*dtr*one_

!..................................................................
!     Left hand side..
!..................................................................

          temp4(i,j)=(temp2(i,j)-temp1(i,j)) &
            *vptb(i,lr_)*cint2(j)*cynt2(i,l_) &
            *one_
          sumleft=sumleft+temp4(i,j)
          sumright=sumright+temp3(i,j)
 100    continue
 90   continue
      error=(sumleft-sumright)/xlndn(k,lr_)
      if (iactst.eq."abort") then
        if (error.gt.1.e-8) call diagwrng(15)
      endif
      return
      end subroutine exsweepx


end module exsweepx_mod
