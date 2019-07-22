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

module tdtoarad_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine tdtoarad
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!     Storing values on whole mesh setup0%lrzmax and equilibrium quantities
!     recall:
!     lr_     : radial index
!     l_      : spatial variable index
!     lmdpln_ : midplane at l_=lmdpln_
!


!..................................................................
!     BH091022:  Looks like a somewhat orphaned idea.  The "-z"
!                variable here a not used much elsewhere, and
!                original (non-z) variables could be used.
!                Seems to create confusion of variable names.
!                Clean it up sometime.....
!                No references to l_ or lmdpln_?
!..................................................................

!..................................................................
!     Copies flux surface dependent quantities into radial arrays
!..................................................................

!..................................................................
!     rfpwrz(k,lr_) is the local rf power deposited in species k
!     at flux surface lr_ in Watts/cc (averaged over the orbit)
!     currz(k,lr_) is the current of species k at lr_ in Amps (averaged
!     currtz(lr_) is the sum of currz(k,lr_) over ionic general species.
!     currmtz(lmdpln_) is the sum of the midplane currents over ionic species
!     curr(m)tpz(lr_) is curr(m)tz plus the electron contribution to
!     the current.
!..................................................................


      do 21 k=1,ngen
        if (nso.gt.0) then
          numsrce=nso
        else
          numsrce=1
        endif
        do 23 m=1,numsrce
          asorz(k,m,lr_)=asor(k,m,lr_)
 23     continue
 21   continue

      if (eqmod.eq."enabled") then
        area(lr_)=areacon(lr_)
        vol(lr_)=volcon(lr_)
        equilpsi(lr_)=epsicon(lr_)
        bmdplne(lr_)=bmidplne(lr_)
        onovrpz(lr_,1)=onovrp(1,lr_)
        onovrpz(lr_,2)=onovrp(2,lr_)
        aspin(lr_)=(rpcon(lr_)-rmcon(lr_))/(rpcon(lr_)+rmcon(lr_))
        rmconz(lr_)=rmcon(lr_)
        rpconz(lr_)=rpcon(lr_)
        fpsiz(lr_)=fpsi(lr_)
        bpolsqaz(lr_)=bpolsqa(lr_)

!..................................................................
!     Compute the plasma pressure (cgs)
!     if (kpress(k).ne. "enabled") this species is not used
!     in calculation of total pressure.
!..................................................................

        phot=0.
        do 10 k=1,ngen
          if (kpress(k).eq."disabled") go to 10
          phot=phot+reden(k,lr_)*energy(k,lr_)*2./3.*ergtkev
 10     continue
        prest(lr_)=phot
        do 20 k=ngen+1,ntotal
          if (kpress(k).eq."disabled") go to 20
          prest(lr_)=prest(lr_)+reden(k,lr_)*energy(k,lr_)*2./3.* &
            ergtkev
 20     continue
 30     continue
      endif
      return
      end subroutine tdtoarad

end module tdtoarad_mod
