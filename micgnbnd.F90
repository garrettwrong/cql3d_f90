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

module micgnbnd_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use micgmbnd_mod, only : micgmbnd

  !---END USE

!
!

contains

      subroutine micgnbnd
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................

!     This routine computes the "average" bounce time
!     in the region of velocity space (a sliver) about the
!     pass/trapped boundary. Note that the bounce time is infinite
!     on the p/t boundary (y(itl,l_)) but the singularity is
!     integrable.
!     See the user manual for mathematical details.
!..................................................................


      sp=(y(itl,l_)+y(itl-1,l_))*.5
      sm=(y(itl+1,l_)+y(itl,l_))*.5
      dp=cos(sp)/coss(itl,l_)-1.
      dm=-(cos(sm)/coss(itl,l_)-1.)

      ilzhfs=lz
      if (numclas .eq. 1) ilzhfs=lz/2+1

!..................................................................
!     For standard circular cross section model..
!..................................................................

      if (psimodel.eq."axitorus") then
        betta(0,lr_)=0.
        betta(1,lr_)=1.
        alm(0,lr_)=1.
        alm(1,lr_)=-.5
        do 10 m=2,mbet
          betta(m,lr_)=-((2*m-3)*betta(m-2,lr_)-(4.*m-4.) &
            *betta(m-1,lr_))/(2*m-1)
          alm(m,lr_)=(2*m-3)/(2.*m)*alm(m-1,lr_)
 10     continue
        ak0=1.38629
        bk0=.5
 20     continue
        s=0.
        do 30 m=0,mbet
          s=s-alm(m,lr_)*coss(itl,l_)**(2*m+1)*betta(m,lr_)*(dm+dp)
 30     continue
        s=s+coss(itl,l_)*sinn(itl,l_) &
!990131     1    *(dm*(ak0-bk0*(alog(2.*dm)-1.))+dp*(ak0-bk0*(alog(2.*dp)-1.))) &
          *(dm*(ak0-bk0*(log(2.*dm)-1.))+dp*(ak0-bk0*(log(2.*dp)-1.)))
        xlbnd(lr_)=s*8.*zmax(lr_)

!..................................................................
!     For non-circular geometries use alternate algorithm.
!............................................................

      else
        xlbnd1=4.*pi*coss(itl,l_)*coss(itl,l_)*tau(itl,lr_)
        call micgmbnd(xlbmd2)
        xlbnd(lr_)=(dp+dm)*(xlbnd1)+xlbmd2
      endif
      xu=xlbnd(lr_)/(twopi*sinn(itl,l_)*dy(itl,l_)*coss(itl,l_)*2.)
      dtau(itl,ilzhfs,lr_)=xu-tau(itl,lr_)
      dtau(itu,ilzhfs,lr_)=dtau(itl,ilzhfs,lr_)
      if (lmax(itl,lr_)+1 .eq. ilzhfs) go to 60
      if (lmax(itl,lr_)+1 .gt.lz) stop 'in micgnbnd: need to check this'
      dtau(itl,lmax(itl,lr_)+1,lr_)=0.
      dtau(itu,lmax(itl,lr_)+1,lr_)=0.
 60   continue
      lmax(itl,lr_)=ilzhfs
      lmax(itu,lr_)=ilzhfs
      tau(itl,lr_)=xu
      tau(itu,lr_)=xu
      do 40 i=1,iy
        vptb(i,lr_)=abs(coss(i,l_))*tau(i,lr_)
 40   continue
      return
      end subroutine micgnbnd


end module micgnbnd_mod
