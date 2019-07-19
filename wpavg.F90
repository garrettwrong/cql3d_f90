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

module wpavg_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine wpavg
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     Compute flux surface average of various quantities for
!     CQLP case
!..............................................................


!.......................................................................
!l    1. Surface average of parallel current and parallel electric field
!     <j_par/R>, <E_par/R> with <a>=int(ds a/B) / int(ds/B)
!.......................................................................

      zcuravg=0.0
      zcuravg2=0.0
      zeleavg=0.0
      zeleavg2=0.0
      z1oravg=0.0
      zflxavg=0.0
      zflxavg2=0.0
      ilr=setup0%lrindx(1)
      zelcof=elecfld(ilr)/300.*rmag*fpsi(ilr)/bmidplne(ilr)**2/3.e+09
      zelcof2=elecfld(ilr)/300.*rmag*fpsi(ilr)/bmidplne(ilr)/3.e+09
      do 100 l=1,setup0%ls
        zcuravg=zcuravg+dsz(l)*currmtpz(l)/psis(l)/bmidplne(ilr)/ &
          solrs(l)
        zcuravg2=zcuravg2+dsz(l)*currmtpz(l)
        zeleavg=zeleavg+dsz(l)*zelcof/(solrs(l)*psis(l))**2/solrs(l)
        zeleavg2=zeleavg2+dsz(l)*zelcof2/solrs(l)**2/psis(l)
        z1oravg=z1oravg+dsz(l)/psis(l)/bmidplne(ilr)/solrs(l)
        zflxavg=zflxavg+dsz(l)/psis(l)/bmidplne(ilr)
        zflxavg2=zflxavg2+dsz(l)
 100  continue
      zcuravg=zcuravg/z1oravg
      zeleavg=zeleavg/z1oravg

      write(6,'(/" surface averages:  <j_par/R>/<1/R>= ",1pe13.4,/ &
        "                    <E_par/R>/<1/R>= ",1pe13.4,/ &
        "                    <j_par*B>      = ",1pe13.4,/ &
        "                    <E_par/B>      = ",1pe13.4,/ &
        "        <E_par/R>/<j_par/R>/sptz(1)= ",1pe13.4,/ &
        "        <E_par*B>/<j_par*B>/sptz(1)= ",1pe13.4,/ &
        "  <1/R>= ",1pe13.4,"   flxavg= ",1pe13.4, &
        " flxavg2= ",1pe13.4,"  n=",i4)') &
        zcuravg,zeleavg,zcuravg2,zeleavg2,zeleavg/zcuravg/sptzr(1), &
        zeleavg2/zcuravg2/sptzr(1),z1oravg,zflxavg,zflxavg2,n

      return
      end subroutine wpavg

end module wpavg_mod
