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

module eqcoord_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use eqfndpsi_mod, only : eqfndpsi
  use eqrhopsi_mod, only : eqrhopsi

  !---END USE



contains

      subroutine eqcoord
      use param_mod
      use cqlcomm_mod
      use eqrhopsi_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      character*8 generate


!..................................................................
!     This routine controls the translation of the equilibrium psi
!     and f =(R*Btor) data into information specific to a flux surface
!     characterized by radial coordinate rovera(lr_). We have
!     erhocon(lr_)=rovera(lr_)*rhomax where rhomax=eqrho(nconteqn).
!     For a circular flux surface rhomax would correspond to radmin.
!     !!!! for small aspect ratio only  !!!!
!..................................................................

!..................................................................
!     Generate an array eqrho(j), j=1,..., nconteqn. This array
!     contains the toroidal radial coordinate, rho, corresponding
!     to the psi array eqpsi(j). This will be used later as a
!     basis for spline interpolation.
!..................................................................
      generate="disabled"
      if (lr_.eq.setup0%lrzmax) then
        generate="enabled"
        call eqrhopsi(generate)
      endif

!..................................................................
!     Next determine the radial coordinate of interest.
!..................................................................

      if (rovera(lr_).ge.0) then
        erhocon(lr_)=rovera(lr_)*rhomax
      endif

!..................................................................
!     Determine the psi value (epsicon(lr_)) corresponding to
!     erhocon(lr_) (erhocon(lr_) is passed in a common block).
!..................................................................

      call eqfndpsi(psides,areades,volum)
      epsicon(lr_)=psides
      areacon(lr_)=areades
      volcon(lr_)=volum
!-----------------------------------------------------------------------
!     Determine some geometrical factors needed to calculate the aspect
!     ratio and the elongation
!-----------------------------------------------------------------------
      rgeom(lr_) =0.5 * (rpcon(lr_) - rmcon(lr_))
      r0geom(lr_)=0.5 * (rpcon(lr_) + rmcon(lr_))
      call aminmx(solz_,1,lorbit(lr_),1,zmincon1,zmaxcon1,kmin,kmax)
      zgeom(lr_)=zmaxcon1

!..................................................................
!     If the density and temperature profiles are to have a strict
!     psi dependence, do that here. This will overwrite the default
!     rho profile determined in subroutine tdxinit
!..................................................................

      if (profpsi.eq."enabled".and. setup0%lrzmax.gt.1) then
        do 1000 k=1,ntotal
          reden(k,lr_)=reden(k,0)*(epsicon(lr_)/psimag)**.5
          temp(k,lr_)=(epsicon(lr_)/psimag)*temp(k,0)
 1000   continue
      endif
      return
      end subroutine eqcoord


end module eqcoord_mod
