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

module wpchgdy_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine wpchgdy
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     Adapt dyp5(iyh) and dy(iyh) and related quantities, such that
!
!     sum(i=1,iyh) [cynt2(i,l_)] = sum(i=1,iyh) [cynt2(i,ll)]
!     where ll=1 for CQLP and setup0%lrz for CQL3D. This may be useful for
!     meshy=fixed_mu option
!.......................................................................

!.......................................................................
!l    1. Re-define dyp5(iyh,l_)
!.......................................................................

      ilref=1
      if (setup0%cqlpmod .ne. "enabled") ilref=setup0%lrz
      if (l_ .eq. ilref) go to 999

      zsum0=0.0
      do 100 i=1,iyh_(ilref)
        zsum0=zsum0+cynt2(i,ilref)/twopi
 100  continue
      zsuml=0.0
      do 110 i=1,iyh-1
        zsuml=zsuml+cynt2(i,l_)/twopi
 110  continue
      dyp5(iyh,l_)=2.*(zsum0-zsuml)/sinn(iyh,l_) - dym5(iyh,l_)

!.......................................................................
!l    2. Correct related quantities
!.......................................................................

      dy(iyh,l_)=0.5*(dym5(iyh,l_)+dyp5(iyh,l_))
      dy(iyh+1,l_)=dy(iyh,l_)
      eyp5(iyh,l_)=1./dyp5(iyh,l_)
      dym5(iyh+1,l_)=dyp5(iyh,l_)
      eym5(iyh+1,l_)=eyp5(iyh,l_)
      cynt2(iyh,l_)=dy(iyh,l_)*twopi*sinn(iyh,l_)
      cynt2(iyh+1,l_)=cynt2(iyh,l_)
      twoint(l_)=0.
      do 200 i=1,iy
        twoint(l_)=twoint(l_)+cynt2(i,l_)
 200  continue

 999  return
      end subroutine wpchgdy

end module wpchgdy_mod
