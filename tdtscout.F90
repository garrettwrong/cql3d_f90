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

module tdtscout_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdinterp_mod, only : tdinterp

  !---END USE

!
!

contains

      subroutine tdtscout
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine creates a file 'tscout'; data
!     is read in from the file created by the TSC code. The file
!     created here has the power and the current as a function of
!     rho_, the radial coordinate given by TSC.
!..................................................................

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

      dimension powtsc(nrada),currtsc(nrada)

#ifdef __MPI
      if(mpirank.ne.0) return
#endif

      if (eqsource.ne."tsc") return

!..................................................................
!     Create output file.
!..................................................................

      open(unit=12,file='tscout',delim='apostrophe',status='unknown')

!..................................................................
!     interpolate to TSC mesh
!..................................................................

      do 10 l=1,setup0%lrzmax
        tr(l)=sorpw_rf(kelecg,l)
 10   continue
      call tdinterp("free","free",rya(1:ubound(rya,1)),tr(1),setup0%lrzmax,rho_,powtsc, &
        npsitm)
      call tdinterp("free","free",rya(1:ubound(rya,1)),currtpz(1),setup0%lrzmax,rho_,currtsc &
        ,npsitm)

!..................................................................
!     Now write to disk.
!     powtsc is the power in Watts/cc
!     currtsc is the current in Amps/cm**2
!..................................................................

      write(12,100) (powtsc(l),l=1,npsitm)
      write(12,100) (currtsc(l),l=1,npsitm)
 100  format(5e16.6)
      close(unit=12)
      return
      end subroutine tdtscout

end module tdtscout_mod
