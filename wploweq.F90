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

module wploweq_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

  subroutine wploweq
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     Determine equilibrium quatities in the lower half cross-section
!     assuming an up/down symmetric flux surface. The total mesh
!     points is setup0%lsmax, thus the mesh point at theta_pol=pi is setup0%lsmax/2+1
!.......................................................................


!.......................................................................
!l    1. Copy values of top half flux-surface onto lower half
!     Assumes there is only one flux surface lr_ when setup0%cqlpmod=enabled
!.......................................................................

      lr_=setup0%lrindx(1)
      zmaxtot=2.*zmax(lr_)

!     lz meshes
      do 100 l=lz/2+2,lz
        z(l,lr_)=zmaxtot-z(lz-l+2,lr_)
        pol(l,lr_)=2.*pi-pol(lz-l+2,lr_)
!     symmetric quantities
        bbpsi(l,lr_)=bbpsi(lz-l+2,lr_)
 100  continue

!     setup0%ls meshes
      do 110 l=setup0%ls/2+2,setup0%ls
        sz(l)=zmaxtot-sz(setup0%ls-l+2)
!     symmetric quantities
        psis(l)=psis(setup0%ls-l+2)
        psipols(l)=psipols(setup0%ls-l+2)
        solrs(l)=solrs(setup0%ls-l+2)
        dsz(l)=dsz(setup0%ls-l+2)
        dszp5(l)=dszm5(setup0%ls-l+2)
        dszm5(l)=dszp5(setup0%ls-l+2)
        eszp5(l)=eszm5(setup0%ls-l+2)
        eszm5(l)=eszp5(setup0%ls-l+2)
!     anti-symmetric quantities
        psisp(l)=-psisp(setup0%ls-l+2)
        solzs(l)=-solzs(setup0%ls-l+2)
 110  continue
!     fix some quantities (periodic conditions assumed)
      z(lz+1,lr_)=zmaxtot
      pol(lz+1,lr_)=2.*pi
      bbpsi(lz+1,lr_)=bbpsi(1,lr_)
      sz(setup0%ls+1)=zmaxtot
      psis(0)=psis(setup0%ls)
      psis(setup0%ls+1)=psis(1)
      psisp(0)=psisp(setup0%ls)
      psisp(setup0%ls+1)=psisp(1)
      psipols(0)=psipols(setup0%ls)
      psipols(setup0%ls+1)=psipols(1)
      solrs(0)=solrs(setup0%ls)
      solrs(setup0%ls+1)=solrs(1)
      solzs(0)=solzs(setup0%ls)
      solzs(setup0%ls+1)=solzs(1)
      dszp5(0)=dszp5(setup0%ls)
      eszp5(0)=eszp5(setup0%ls)
      dszp5(setup0%ls/2+1)=sz(setup0%ls/2+2)-sz(setup0%ls/2+1)
      eszp5(setup0%ls/2+1)=1./dszp5(setup0%ls/2+1)
      dszm5(1)=dszp5(0)
      eszm5(1)=1./dszm5(1)
      dszm5(setup0%ls+1)=dszm5(1)
      dszm5(0)=dszm5(setup0%ls)
      eszm5(setup0%ls+1)=eszm5(1)
      eszm5(0)=eszm5(setup0%ls)
      dsz(0)=0.5*(dszm5(0)+dszp5(0))
      dsz(1)=0.5*(dszm5(1)+dszp5(1))
      dsz(setup0%ls)=dsz(0)
      dsz(setup0%ls+1)=dsz(1)
      dsz(setup0%ls/2+1)=0.5*(dszm5(setup0%ls/2+1)+dszp5(setup0%ls/2+1))

      do 111 k=1,ntotal
        do 112 l=setup0%ls/2+2,setup0%ls
          denpar(k,l)=denpar(k,setup0%ls-l+2)
          temppar(k,l)=temppar(k,setup0%ls-l+2)
 112    continue
 111  continue

      return
      end subroutine wploweq

end module wploweq_mod
