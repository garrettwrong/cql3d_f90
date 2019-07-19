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

module diagscal_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use r8subs_mod, only : dcopy
  use r8subs_mod, only : dscal

  !---END USE

!
!

contains

      subroutine diagscal(k)
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dscal, dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine scales the distribution function to maintain
!     FSA constant density, or given time-dependent density.
!     It is used only if lbdry(k)="scale" or "consscal" (which uses
!     conservative differencing at v=0 and scaling of density).
!yup140806: Prior, it was only enabled for lbdry(k)="scale".
!..................................................................

      ratio(k,lr_)=1.
!yup+bh140806: Enable new "consscal" option
      if (lbdry(k).ne."scale" .and. lbdry(k).ne."consscal") return

!..................................................................
!     Compute the new density..
!..................................................................

      call dcopy(iyjx2,f(0:iy+1,0:jx+1,k,l_),1,temp2(0:iy+1,0:jx+1),1)
      call dcopy(iyjx2,temp2(0:iy+1,0:jx+1), 1,temp1(0:iy+1,0:jx+1),1)
      call bcast(tam1,zero,jx)
      call bcast(tam4,zero,jx)

      ! Now integrate over d3v to find n0 and <n>
      do 1 i=1,iy
        do 2 j=1,jx
            if(temp1(i,j).gt.zero) then
               tam1(j)=tam1(j)+temp2(i,j)*cynt2(i,l_) !for n0 (midplane)
               tam4(j)=tam4(j)+vptb(i,lr_)*temp2(i,j)*cynt2(i,l_) !<n> ZOW only
            endif
 2      continue
 1    continue
      gn=0.d0
      sden=0.
      do 3 j=1,jx
          gn=gn+tam1(j)*cint2(j) !  n0 == Midplane density
          sden=sden+tam4(j)*cint2(j)*one_  ! <n> ZOW only
 3    continue

!..................................................................
!     Determine the velocity mesh point j such that for electron
!     runaway calculations the distribution density up to x(j)
!     is held constant. Used only for lbdry(k)="fixed" calculations.
!     This option is currently inoperable.
!..................................................................

      j=jx
!     esfac=1.
!     if (k .ne. kelecg) go to 21
!     evol=esfac*sqrt(2.)*vth(kelec,lr_)/sqrt(abs(eoved))
!     do 20 j=1,jx
!     if (x(j)*vnorm .gt. evol) go to 21
!20   continue
!21   continue

!..................................................................
!     runden is the new density - will be used in denominator
!     of scale factor.
!..................................................................

      if ( j .gt. jx) j=jx
      runden=0.
      do 24 jj=1,j
        runden=runden+tam4(jj)*cint2(jj)*one_ !<n> ZOW only
 24   continue

!      write(*,'(a,2i3,4e12.5)')'k,lr,xlndn00(k,lr),runden,sden,gn=',
!     +                          k,lr_,xlndn00(k,lr_),runden,sden,gn

!..................................................................
!     If time-dependent density, add Maxwellian particles (at present
!     temperature or temp_den.ne.0.) to achieve reden(k,lr_).
!     Else, scale f to original density.
!..................................................................

!YuP,BH180918: enables proper treatment of spline-t profile option,
!YuP,BH180918: adding enein_t
      if (nbctime.gt.0 .and. (redenc(1,k).ne.zero &
            .or.  enein_t(1,k,1).ne.zero) ) then
         if (temp_den.ne.zero) then
            temp_ad=temp_den
         else
            temp_ad=temp(k,lr_)
         endif
         thta=fmass(k)*clite2/(temp_ad*ergtkev)
         do j=1,jx
            swwtemp=exp(-gamm1(j)*thta)
            do i=1,iy
               temp1(i,j)=swwtemp
            enddo
         enddo
         call bcast(tam1,zero,jx)
         do i=1,iy
            do j=1,jx
               tam1(j)=tam1(j)+vptb(i,lr_)*temp1(i,j)*cynt2(i,l_)
            enddo
         enddo

         sden1=0.
         do j=1,jx
            sden1=sden1+tam1(j)*cint2(j)
         enddo

         sden1=sden1/zmaxpsi(lr_)

         ratio1=(reden(k,lr_)-runden/zmaxpsi(lr_))/sden1
         do j=1,jx
            do i=1,iy
               f(i,j,k,l_)=f(i,j,k,l_)+ratio1*temp1(i,j)
            enddo
         enddo ! j
      else  !---> nbctime=0 (not a time-dependent profile)
         ratio(k,lr_)=xlndn00(k,lr_)/runden ! field-line-aver <n>
         call dscal(iyjx2,ratio(k,lr_),f(0:iy+1,0:jx+1,k,l_),1) !rescale f
!         write(*,'(a,3i4,3e13.5)')'diagscal: n,k,lr,xlndn00,runden,gn',&
!                                    n,k,lr_,xlndn00(k,lr_),runden,gn
         write(*,'(a,3i4,e13.5)') &
          'diagscal: n,k,lr=;  f is rescaled by ratio()=', &
                     n,k,lr_,  ratio(k,lr_)
         ! Note: it is not always a good idea to rescale the distr.func.
         ! Example: When there is a large-power NBI source,
         ! comparing to the initial background
         ! set in cqlinput. The value of xlndn00
         ! (in ratio(k,lr_)=xlndn00(k,lr_)/runden) is based
         ! on the initial density, so it does not include particles from NBI.
         ! And the value of runden (or sden) does include all sources.
         ! In such a case, it is better to use lbdry(k)="conserv"
      endif

      return
      end subroutine diagscal

end module diagscal_mod
