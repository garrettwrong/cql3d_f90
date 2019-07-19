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

module eqfpsi_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use eqwrng_mod, only : eqwrng
  use zcunix_mod, only : terp1

  !---END USE

!
!

contains

      subroutine eqfpsi(psval,fpsi__,fppsi__)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine provides f(psi) to model the toroidal
!     magnetic field. For cases that eqsource="ellipse"
!     the f is ad-hoc and is determined through the namelist
!     model, fpsimodel. In the case that eqsource="filename", then
!     a file exists on disk which provides f and the equilibrium
!     psi. As of 9/21/88 filename=eqdsk or topeol.
!     Also provided is the derivative df/dpsi, fppsi.
!..................................................................

      if (eqsource.eq."ellipse") then
        if (fpsimodl.eq."constant") then
          fpsi__=btor*radmaj
          fppsi__=0.
        else
          call eqwrng(7)
        endif
      else
        itab(1)=1
        itab(2)=1
        itab(3)=0
        call terp1(nfp,psiar,fpsiar,d2fpsiar,psval,1,tab,itab)
        fpsi__=tab(1)
        fppsi__=tab(2)
      endif
      return
      end subroutine eqfpsi


!
!
      subroutine eqppsi(psval,ppsi__,pppsi__)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine provides p(psi) to model the plasma
!     pressure. For cases that eqsource="ellipse"
!     the p i zero. In the case that eqsource="filename", then
!     a file exists on disk which provides f and the equilibrium
!     psi. As of 9/21/88 filename=eqdsk or topeol.
!     Also provided is the derivative dp/dpsi, pppsi.
!..................................................................

      if (eqsource.eq."ellipse") then
         ppsi__=0.
         pppsi__=0.
      else
        itab(1)=1
        itab(2)=1
        itab(3)=0
        call terp1(nfp,psiar,prar,d2prar,psval,1,tab,itab)
        ppsi__=tab(1)
        pppsi__=tab(2)
      endif
      return
      end subroutine eqppsi

end module eqfpsi_mod
