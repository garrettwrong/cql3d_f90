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

module impchk_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use diagentr_mod, only : gfi
  use diagwrng_mod, only : diagwrng

  !---END USE

!
!

contains

      subroutine impchk(k)
      use param_mod
      use cqlcomm_mod
      use advnce_mod !here: in impchk(). To get gfi(), qz(),ry(), hfi(), etc.
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine computes the density gain (from dF/dt) and compares i
!     locally to the density gain from the r.h.s. of the F.P. equation.
!     sumleft and sumright are the integrated totals. error is the
!     normalized difference and should be as close to roundoff as
!     possible.
!..................................................................

!%OS
      dimension zermx(51),imx(51),jmx(51),zsumj3(iy),zsumj4(iy)
      dimension zdiffj(iy),zdiffi(jx)
      dimension zdns(lrorsa)
!%OS
      fpithta(i,j)=f(i+1,j,k,l_)*(1.-dithta(i,j,l_)) + &
                   f(i  ,j,k,l_)*dithta(i,j,l_)

!.......................................................................

      sumleft=0.
      sumright=0.
      call bcast(tam1,zero,jx)
      zdns(l_)=0.0

!..................................................................
!     temp3 contains the r.h.s.and temp4 the l.h.s.
!..................................................................

      do 10 i=1,iy
        if ((i.eq.itl .or. i.eq.itu) .and. setup0%cqlpmod.ne."enabled")go to 10
!%OS  do 2 j=1,jx
        do 2 j=2,jx-1
          temp4(i,j)=(f(i,j,k,l_)-f_(i,j,k,l_)) &
            *vptb(i,lr_)*cynt2(i,l_)*cint2(j) &
            *one_
          temp3(i,j)=gfi(i,j,k)*qz(j)
          temp3(i,j)=temp3(i,j)-gfi(i,j-1,k)*qz(j)
          temp3(i,j)=temp3(i,j)+(hfi(i,j,k,l_)-hfi(i-1,j,k,l_))*ry(i,j,l_)
          temp3(i,j)=(temp3(i,j)+vptb(i,lr_) &
            *(cah(i,j)*f(i,j,k,l_)+so(i,j))+spasou(i,j,k,l_)+ &
            cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j)) ) &
            *cynt2(i,l_)*cint2(j)*one_
          temp3(i,j)=temp3(i,j)*dtreff
          sumleft=sumleft+temp4(i,j)
          sumright=sumright+temp3(i,j)
          temp5(i,j)=gfi(i,j,k)
          temp6(i,j)=hfi(i,j,k,l_)
          temp1(i,j)=  qz(j)*(gfi(i,j,k)-gfi(i,j-1,k))
          temp2(i,j)=ry(i,j,l_)*(hfi(i,j,k,l_)-hfi(i-1,j,k,l_))
          zdns(l_)=zdns(l_)+temp3(i,j)/dtreff- &
            vptb(i,lr_)*spasou(i,j,k,l_)*cynt2(i,l_)*cint2(j)
 2      continue
 10   continue

!..................................................................
!     Pass/trapped boundary
!..................................................................

      if (setup0%cqlpmod .eq. "enabled") go to 50

!%OS  should be j=1,jx? may depend on lbdry
      do 3 j=2,jx-1
        temp4(itl,j)=(f(itl,j,k,l_)-f_(itl,j,k,l_)) &
          *vptb(itl,lr_)*cynt2(itl,l_) &
          *cint2(j)*one_
        temp4(itu,j)=temp4(itl,j)
        temp3(itl,j)=gfi(itl,j,k)*qz(j)
        temp3(itl,j)=temp3(itl,j)-gfi(itl,j-1,k)*qz(j)
        temp3(itl,j)= &
          temp3(itl,j)-r2y(j,l_)*(hfi(itl-1,j,k,l_)-2.*hfi(itl,j,k,l_)-hfi(itu,j,k,l_))
        temp3(itl,j)=temp3(itl,j)+(cah(itl,j)*f(itl,j,k,l_) &
          +so(itl,j))*vptb(itl,lr_) + spasou(itl,j,k,l_)
        temp3(itl,j)=temp3(itl,j)*dtreff*cint2(j)*cynt2(itl,l_) &
          *one_
        temp3(itu,j)=temp3(itl,j)
        sumleft=sumleft+2.*temp4(itl,j)
        sumright=sumright+temp3(itl,j)*2.
        temp5(itl,j)=gfi(itl,j,k)
        temp5(itu,j)=gfi(itu,j,k)
        temp6(itl,j)=hfi(itl,j,k,l_)
        temp6(itu,j)=hfi(itu,j,k,l_)
 3    continue

 50   continue

!%OS
!     compute highest errors
      do 99 ii=1,10
        zermx(ii) = 0.0
 99   continue
!
      do 100 i=1,iy
        zsumj3(i) = 0.0
        zsumj4(i) = 0.0
        do 110 j=1,jx
          zsumj3(i) = zsumj3(i) + temp3(i,j)
          zsumj4(i) = zsumj4(i) + temp4(i,j)
          zerabs = abs(temp3(i,j)-temp4(i,j))
          if (zerabs .le. zermx(10)) go to 110
          if (zerabs .gt. zermx(1)) then
            do 111 ii=9,1,-1
              zermx(ii+1) = zermx(ii)
              imx(ii+1) = imx(ii)
              jmx(ii+1) = jmx(ii)
 111        continue
            zermx(1) = zerabs
            imx(1) = i
            jmx(1) = j
            go to 110
          endif
          do 112 int=9,1,-1
            if (zerabs.le.zermx(int) .and. zerabs.gt.zermx(int+1)) &
              then
              do 113 ii=9,int+1,-1
                zermx(ii+1) = zermx(ii)
                imx(ii+1) = imx(ii)
                jmx(ii+1) = jmx(ii)
 113          continue
              zermx(int+1) = zerabs
              imx(int+1) = i
              jmx(int+1) = j
              go to 110
            endif
 112      continue
 110    continue
        zdiffj(i) = zsumj4(i)-zsumj3(i)
 100  continue
      do 120 j=1,jx
        zsumi3 = 0.0
        zsumi4 = 0.0
        do 121 i=1,iy
          zsumi3 = zsumi3 + temp3(i,j)
          zsumi4 = zsumi4 + temp4(i,j)
 121    continue
        zdiffi(j) = zsumi4-zsumi3
 120  continue
!
!%OS  if (l_.eq.1 .and. (n/2)*2..eq.n)
!%OS  + write(6,'(1pe10.2,2i5)') (zermx(ii),imx(ii),jmx(ii),ii=1,10)

!..................................................................
!     Symmetry about pi/2. in trapped region.
!..................................................................

      if (symtrap .eq. "enabled") then
        do 4 i=iyh+1,itu
          ii=iy+1-i
          do 6 j=1,jx
            temp3(i,j)=temp3(ii,j)
 6        continue
 4      continue
      endif

!..................................................................
!     Compute the normalized error.
!..................................................................


      error=(sumleft-sumright)/xlndn(k,lr_)
      if (iactst.eq."abort") then
        if (error.gt.1.e-8) call diagwrng(15)
      endif

!%OS
      if (l_ .eq. lrors) then
        idumy=0
 999    continue
      endif
!%OS
      return
      end subroutine impchk


end module impchk_mod
