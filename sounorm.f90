module sounorm_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use soup_mod, only : soup

  !---END USE


!*****************************************************************

!
!

contains

      subroutine sounorm
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine establishes an array of normalization constants
!     which allows the code to force a particular poloidal source
!     profile. It is called only at initialization.
!..................................................................




!..................................................................
!     Set a flag.
!..................................................................

      isounor=1
      do 200 l=1,lz
        do 100 k=1,ngen
          call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)
          do 50 m=1,nso

!..................................................................
!     Call a routine which determines the non-normalized source
!     distribution for a given i and all j.
!..................................................................

            do 40 i=1,iy
              call soup(coss(i,l_),l,k,m)
              do 30 j=1,jx
                temp1(i,j)=soupp(j,lr_)
 30           continue
 40         continue

!..................................................................
!     Set the array sounor such that when the computed source
!     profile is multiplied by sounor the resulting current
!     density is unity.
!..................................................................

            s=0.
            do 10 i=1,iy
              do 20 j=1,jx
                s=s+temp1(i,j)*cynt2(i,l_)*cint2(j)
 20           continue
 10         continue
            if (s.ne.zero) sounor(k,m,l,lr_)=1./(s*one_)
 50       continue
 100    continue
 200  continue

!..................................................................
!     reset the flag (subroutine soup)
!..................................................................

      isounor=0
      return
      end subroutine sounorm
      
      
end module sounorm_mod
