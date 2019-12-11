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

module vlhd_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      real(c_double) function vlhd(vll,vprp,thp,nmod)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!     Determine the local phenomenological lower hybrid
!     diffusion coefficient as a function of local (in z)
!     v-parallel, v-perpendicular and poloidal angle.
!
!
!


      if (vprprop .eq. "enabled") then
        xvpr=1.
      else
        xvpr=0.
      endif
      if (vdalp .lt. 1.e-8) then
         if(setup0%verbose>0) write(*,*)'vlhd: ***WARNING***  vdalp.lt.1.e-8, reset to 0.03'
         vdalp=.03
      endif
!
!     determine the vparallel interval at R=R0
!
      vpmax=vparmax(nmod)*clight
      vpmin=vparmin(nmod)*clight
!
!
!     determine the vparallel range of nonzero D at this
!     orbit point.
!
      rovr0=(radmaj+xvpr*radmin*rovera(lr_)*cos(thp))/radmaj
      vmin=vpmin*rovr0
      vmax=vpmax*rovr0
!
!     determine the vperp interval at R=R0
!
      vprmax=vprpmax(nmod)*clight
      vprmin=vprpmin(nmod)*clight
!
!
!     determine the vperp range of nonzero D at this
!     orbit point.
!
      vprmin=vprmin*rovr0
      vprmax=vprmax*rovr0

!     Determine the local diffusion coefficient.
!
      delv=(vmax-vmin)*vdalp

      if (vmin .gt. vll .or. vmax .lt. vll) then
        vlhd=0.

      elseif (vprmin.gt.vprp .or. vprmax.lt.vprp) then
        vlhd=0.

      elseif (vlhpolmn(nmod)/57.29577.gt.thp &
              .or. vlhpolmx(nmod)/57.29577.lt.thp) then
        vlhd=0.

      else
!       the nonzero diffusion coefficient..
!
        dene=reden(kelec,lr_)
        te=temp(kelec,lr_)*1.e3
!990131        xlog=24.-alog(sqrt(dene)/te)
        xlog=24.-log(sqrt(dene)/te)
        if (kelecg.eq.1) then
          taucoll=3.44074e5*te**1.5/(zeff(lr_)*dene*xlog)
        elseif (kiong(1).eq.1) then
          ti=temp(kiong(1),lr_)*1.e3
          taucoll=2.08507e7*ti**1.5 &
              /(zeff(lr_)*dene*bnumb(kiong(1))**2*xlog)
        else
          stop 'in vlhd'
        endif
        difus=dlndau(nmod)/taucoll*(vth(1,lr_)/vnorm)**2

!BH021118:  Karney uses following with vlh_karney=1.,
!           see, Karney and Fisch, Phys. Fl. 28, p. 116 (1985).
        if (vlh_karney.ne.0.) then
           clight2=clight*clight
           uu=vll**2+vprp**2
           if (relativ.eq.'disabled') then
              uu=sqrt(uu)
           else
              uu=sqrt(uu/(1.-uu/clight2))
           endif
           factor=(1./(1.+uu/vth(1,lr_)))**vlh_karney
           difus=factor*difus
!           write(*,*)'vlhd: factor=',factor
        endif

        if (vmin.le.vll .and. vll.le. vmin+delv) then
          vlhd=difus*(1.-cos((vll-vmin)/delv*pi))*.5
        elseif (vmax-delv .le. vll .and. vll .le. vmax) then
          vlhd=difus*(1.-cos((vmax-vll)/delv*pi))*.5
        else
          vlhd=difus
        endif
      endif

      return
      end function vlhd

end module vlhd_mod
