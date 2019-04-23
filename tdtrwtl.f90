module tdtrwtl_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine tdtrwtl
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

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
