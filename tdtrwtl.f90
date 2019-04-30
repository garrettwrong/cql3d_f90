module tdtrwtl_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine tdtrwtl
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     Computes the Chang-Cooper weights for the transport
!     calculation.
!..............................................................


      do 1 k=1,ngen
        do 15 j=1,jx
          do 10 l=0,lrors-1
            do 20 i=1,iymax
              dl(i,j,k,l)=.5
 20         continue
 10       continue
 15     continue
 1    continue

      return
      end
end module tdtrwtl_mod
