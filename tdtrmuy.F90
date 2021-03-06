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

module tdtrmuy_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use micgetr_mod, only : micgetr

  !---END USE

!
!

contains

      subroutine tdtrmuy
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     This routine sets up the template pitch angle mesh to be
!     emulated by the actual pitch angle meshes on all the
!     flux surfaces.
!..............................................................

      ipacktp=0   ! y() mesh pts around t-p bndry not removed and
                  ! the v-space and transport y-meshes are the same.
                  ! Other option is =3 giving number of points
                  ! removed around each tp bndry.  This was the
                  ! only option pre-070419.
      if (meshy.eq."free") return
!BH070419      if (setup0%cqlpmod .ne. "enabled") iytr(lrors)=iy-6
      if (setup0%cqlpmod .ne. "enabled" .and. tfac.ge.0.) then
         ipacktp=3
         iytr(lrors)=iy-2*ipacktp  !Setting up transport mesh with bins
                                  !at and neighboring t-p bndry removed.
                                  !Used in splitting algorithm for
                                  !transp='enabled'.
      else  ! i.e., setup0%cqlpmod.eq."enabled" .and./.or. tfac.lt.0.
         ipacktp=0
         iytr(lrors)=iy
      endif

      iy_(lrors)=iy
      iymax=iy_(lrors)
      iyh_(lrors)=iy/2
      iyh=iyh_(lrors)
      yreset="disabled"

!..............................................................
!     Determine a normalized mu mesh if radial derivatives are to
!     be at constant perpendicular adiabatic invariant (mu).
!     Otherwise determine a theta mesh which will be the same
!     on all flux surfaces (except for the p/t region to be added in
!     subroutine tdtry).
!..............................................................

      top=pi*half
      iyhtr=iytr(lrors)/2
      iyhtrp=iyhtr+1

!     Note: mun might be redefined in subroutine wptrmuy

      mun(1)=0.
      hmu=top/iyhtr
      hmufg=abs(tfac)*hmu   ! ?? or abs(tfac) ??
      mun(2)=hmufg
      call micgetr(iyhtrp,top,hmufg,ram,ksingul)
      do 10 i=3,iyhtr
        mun(i)=mun(i-1)+ram*(mun(i-1)-mun(i-2))
 10   continue

      if (meshy.eq."fixed_mu") then
        do 40 i=2,iyhtr
          mun(i)=sin(mun(i))**2
 40     continue
      endif

      return
      end subroutine tdtrmuy

end module tdtrmuy_mod
