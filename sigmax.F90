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

module sigmax_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine sigmax(kk)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!mnt  generate a maxwellian with same density and temp as species
!mnt  "kk".


      if(kk.eq.0) then
         if(setup0%verbose>0) WRITE(*,*)'sigmax: invalid species number kk=0'
         stop
      endif

      fxllm2=2./3.*energy(kk,lr_)/fions(kk)
      fxppm2=fxllm2
      do 10 i=1,iy
      do 10 j=1,jx
      xll=x(j)*coss(i,l_)
      xpp=x(j)*sinn(i,l_)
      facx=-xll**2/fxllm2
      facy=-xpp**2/fxppm2
      temp3(i,j)=exp(facx+facy)
10    continue
      call bcast(tam1,zero,jx)
      do 20 i=1,iy
      do 21 j=1,jx
      tam1(j)=tam1(j)+temp3(i,j)*cynt2(i,l_)
21    continue
20    continue
      gn=0.
      do 22 j=1,jx
      gn=gn+tam1(j)*cint2(j)
22    continue
!BH120314      Following sn gives int-d**3x{temp3}=1
!BH120314      sn=xlndn(kk,lr_)/zmaxpsi(lr_)/(reden(kk,lr_)*gn)
      sn=reden(kk,lr_)/gn
30    continue
!BH120314      do 23 i=0,iyjx2-1
!BH120314      temp3(i,0)=temp3(i,0)*sn
!BH12031423    continue
      do j=0,jx+1
         do i=0,iy+1
            temp3(i,j)=temp3(i,j)*sn
         enddo
      enddo

      return
      end subroutine sigmax

end module sigmax_mod
