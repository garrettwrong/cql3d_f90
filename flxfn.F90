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

module flxfn_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use diagwrng_mod, only : diagwrng

  !---END USE

!
!

contains

      real(c_double) function flxfn(ra)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!mnt  This routine calculates the poloidal magnetic field and the
!mnt  flux function at (normalized) radius ra (r/radmin) based on a
!mnt  gaussian current profile of specified width.


      data idata/0/
!BH151202 Attempting to get code working with machine.eq.'mirror'
!BH151202      if (machine.ne."toroidal")call diagwrng(3)
      if (machine.ne."toroidal".and.machine.ne."mirror")call diagwrng(3)
      if (idata.eq.1) go to 10
      currwth=.5
      currnorm=1.-exp(-(1./currwth)**2)
      psi0=-radmin*radmaj*bth/currnorm*0.5
      wth2=1./currwth/currwth
      idata=1
 10   continue
      xra=ra/currwth
      xrbtemp=xra*xra
      if (xrbtemp.lt.1.e-13) then
        prf=xrbtemp/currnorm
      else
        expfact=exp(-xrbtemp)
        prf=(1.-expfact)/currnorm
      endif
      bthr(lr_)=0.
      atemp=sqrt(radmaj**2-(ra*radmin)**2)
      btemp=atemp/(radmaj**2-radmin**2)**.5
      if (ra.gt.0.) bthr(lr_)=btemp*bth*prf/ra
      robtemp=ra*ra
      wth2n=1.
      robtempn=1.
      denn=1.
      next=1
      last=1
      sum=0.
 20   continue
      wth2n=wth2n*wth2
      robtempn=robtempn*robtemp
      denn=-denn*next*next/last
      term=wth2n*(1.-robtempn)/denn
      sum=sum+term
      last=next
      next=next+1
      if (next.lt.100) go to 20
      flxfn=psi0*sum
      return
      end function flxfn


end module flxfn_mod
