module sourcef_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use r8subs_mod, only : dscal
  use soup_mod, only : soup

  !---END USE

!
!

contains

      subroutine sourcef
      use param_mod
      use comm_mod
      use r8subs_mod, only : dscal
      implicit integer (i-n), real*8 (a-h,o-z)
      save


!..................................................................
!     If soucoord.ne."disabled",
!     computes the source m of species k - simple Gaussian profiles.
!     If knockon.ne."disabled",
!     also calc. knock-on source of electrons (Besedin and Pankratov).
!..................................................................

      do 900 k=1,ngen
          if (soucoord .ne. "disabled") then
        do 800 m=1,nso
          if (n .lt. nonso(k,m) .or. n .gt. noffso(k,m)) go to 800
          call bcast(temp2(0:iyjx2-1,0),zero,iyjx2)
          do 700 l=1,lz

!..................................................................
!     Determine the local source at z(l,lr_) in preparation for bounce-av
!..................................................................

            do 650 i=1,imax(l,lr_)
              ii=iy+1-i
              call soup(cosz(i,l,lr_),l,k,m)
              do 630 j=1,jx
                temp1(i,j)=soupp(j,lr_)
 630          continue
              call soup(cosz(ii,l,lr_),l,k,m)
              do 610 j=1,jx
                temp1(ii,j)=soupp(j,lr_)
 610          continue
 650        continue

!..................................................................
!     Do the bounce average of the source..
!..................................................................

            do 500 i=1,imax(l,lr_)
              ii=iy+1-i
              ax=dtau(i,l,lr_)/tau(i,lr_)
              if ((l.ne.lz .and. lmax(i,lr_).eq.l)) then
                ax=ax+dtau(i,l+1,lr_)/tau(i,lr_)
              endif
              do 450 j=1,jx
                temp2(i,j)=temp2(i,j)+ax*temp1(i,j)
                temp2(ii,j)=temp2(ii,j)+ax*temp1(ii,j)
 450          continue
 500        continue
 700      continue  ! On l

!..................................................................
!     Zero out source where less than small fraction of the
!     peak.  This facilitates use of ineg='enabled' option
!     zeroing out f where f<0 occurs, for all j above max of
!     the source.   BH120317.
!..................................................................

          temp2frac=1.d-8
          temp2max=zero
          do j=1,jx
             do i=1,iy
                temp2max=max(temp2max,temp2(i,j))
             enddo
          enddo
          temp2frac=temp2frac*temp2max
          do j=1,jx
             do i=1,iy
                if (temp2(i,j).lt.temp2frac) temp2(i,j)=zero
             enddo
          enddo


!..................................................................
!     Compute the density in preparation for scaling to achieve
!     the desired current.
!..................................................................

          s1=0.
          do 300 i=1,iy
            do 250 j=1,jx
              s1=s1+temp2(i,j)*cynt2(i,l_)*cint2(j)*vptb(i,lr_)
 250        continue
 300      continue

!..................................................................
!     Determine q1 so that the flux surface averaged current
!     will be asor (particles/sec/cc) after scaling.
!..................................................................

          if (s1.ne.zero) q1=asor(k,m,lr_)/s1*zmaxpsi(lr_)/one_

!..................................................................
!     scale the current to be equal to asor
!..................................................................

          call dscal(iyjx2,q1,temp2(0:iyjx2,0),1)
          do 200 j=1,jx
            do 150 i=1,iy
              source(i,j,k,indxlr_)=source(i,j,k,indxlr_)+temp2(i,j)
 150        continue
 200      continue

!..................................................................
!     xlncur is in units of particles/cm**2 (field line density)
!..................................................................

          xlncur(k,lr_)=xlncur(k,lr_)+asor(k,m,lr_)*zmaxpsi(lr_)
          xlncurt(lr_)=xlncurt(lr_)+xlncur(k,lr_)
 800    continue
      endif

 900  continue
      return
      end
end module sourcef_mod
