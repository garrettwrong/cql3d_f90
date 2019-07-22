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

module soucrit_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use cfpgamma_mod, only : cfpgamma

  !---END USE

!
!

contains

      subroutine soucrit
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : luf
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!   scchiu, 9609..
!  calculate the critical momentum-per-mass.
!  Since jxcrit is used for calculation of relativistic runaways,
!   it is chosen to be no larger than for particles at 3.*clight.
!   This is to reduce a problem of definition which occurs
!   transiently for abs(E) < abs(E_crit_for_runaway).
!
      call cfpgamma
      do 10 k=1,ngen
        fack=abs(elecfld(lr_))/reden(k,lr_)*1.e16/0.0918 &
                    *18./gama(k,k)
        eoe0(k,lr_)=fack

            if((fack-1.).eq.zero)then
             write(*,*)'soucrit: (fack-1.)=0=',(fack-1.)
             !pause
            endif

        fack1=1./(fack-1.)
        if (fack1.le.0.d0) then
          ucrit(k,lr_)=1.e20
        else
          ucrit(k,lr_)=clight*sqrt(fack1)/vnorm
        endif
!  Take runaway electrons to have momentum per mass beyond the
!  minimum of 3.*clight or ucrit:
!990131        xcrit=amin1(3.*clight/vnorm,ucrit(k,lr_))
        xcrit=min(3.*clight/vnorm,ucrit(k,lr_))
        jxcrit(k,lr_)=luf(xcrit,x,jx)
10    continue
      return
      end subroutine soucrit


end module soucrit_mod
