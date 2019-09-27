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

module soup_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use psif_mod, only : psif
  use r8subs_mod, only : dcopy
  use r8subs_mod, only : dscal
  use soup0_mod, only : soup0

  !---END USE

!
!

contains

      subroutine soup(cosi,l,kk,m)
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dscal, dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Given a point z(l,lr_) along the field line z, and given
!     the local cosine of the pitch angle of the particle,cosi,
!     compute a velocity source profile soupp(j,lr_),j=1,jx
!     for source number m of species kk.
!..................................................................


      dimension ifirst(lrza)
      data ifirst/lrza*0/

!..................................................................
!     For the case that a Gaussian profile in polar coordinates
!     is to be computed determine some exponentials depending on
!     speed alone and store them in sovt(j,kk,m,lr_)
!..................................................................

      if (ifirst(lr_) .eq. 0) call soup0
      ifirst(lr_)=1
      zl=z(l,lr_)

!..................................................................
!     ctl is cos(theta(z(l,lr_))) for the zero banana width pinch orbit.
!     In other words the trapped/passing boundary maps to acos(ctl).
!..................................................................

      ctl=sqrt(abs(1.-psif(zl)*(1.-coss(itl,lmdpln_)**2)))
      thtl=acos(ctl)
      if (soucoord.eq."cart") then
        sini=sqrt(abs(1.-cosi*cosi))

!..................................................................
!     Determine the Gaussian profile next
!     tam7,10 contain contribution from cosi and tam8,9 contribution fro
!     -cosi (reflected about pi/2).
!..................................................................

        do 10 j=1,jx
          tam13(j)=x(j)*cosi
          tam12(j)=-tam13(j)
          tam11(j)=x(j)*sini
          tam10(j)=-(tam13(j)-sxllm1(kk,m,lr_))**2/sxllm2(kk,m,lr_)
          tam9(j)=-(tam12(j)-sxllm1(kk,m,lr_))**2/sxllm2(kk,m,lr_)
          tam6(j)=-(tam11(j)-sxppm1(kk,m,lr_))**2/sxppm2(kk,m,lr_)
          tam7(j)=exp(tam10(j)+tam6(j))
          tam8(j)=exp(tam9(j)+tam6(j))
 10     continue

!..................................................................
!     Now for polar coordinate case.
!..................................................................

      else
        facc=-(cosi-cosm1(kk,m,lr_))**2/cosm2(kk,m,lr_)
        faccr=-(-cosi-cosm1(kk,m,lr_))**2/cosm2(kk,m,lr_)
        q1=exp(facc)
        q2=exp(faccr)
        do 11 j=1,jx
          tam7(j)=q1*sovt(j,kk,m,lr_)
          tam8(j)=q2*sovt(j,kk,m,lr_)
 11     continue
      endif

!..................................................................
!     Symmetrize in the trapped region.
!..................................................................

      if (abs(cosi).le.ctl) then
        do 12 j=1,jx
          soupp(j,lr_)=(tam7(j)+tam8(j))*.5
 12     continue
      else
        call dcopy(jx,tam7,1,soupp(1:jx,lr_),1)
      endif

!..................................................................
!     Give the desired z(l,lr_) dependence to the source current
!..................................................................

      if (isounor .ne. 1) then
        facz=exp(-(zl-zm1(kk,m,lr_))**2/zm2(kk,m,lr_))
        call dscal(jx,facz*sounor(kk,m,l,lr_)*asor(kk,m,lr_),soupp(1:jx,lr_),1)
      endif
      return
      end subroutine soup


end module soup_mod
