module urfwr0c_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine urfwr0c(x,nrayelt,nray,nrayelts)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      complex*16 x,c16i
      dimension x(nrayelts,*),nrayelt(*)
!....................................................................
!     For formatted o/p purposes, values of abs(x) .lt. 1.e-300
!     are set equal to 0. (Actually, .gt.1.e-999 OK for format,
!     but some fortrans restricted to .lt.1.e-37).
!....................................................................

      c16i=(0.d0,1.d0)
      deps=1.d-300
      zero=0.d0
      do 10  iray=1,nray
        do 11  is=1,nrayelt(iray)

           z16r=0.5*(x(is,iray)+conjg(x(is,iray)))
           z16i=-c16i*0.5*(x(is,iray)-conjg(x(is,iray)))
           if (abs(z16r).lt.deps) &
                x(is,iray)=cmplx(zero, z16i)

           if (abs(z16i).lt.deps) &
                x(is,iray)=cmplx(z16r, zero)

           if (abs(z16r).lt.deps .and. abs(z16i).lt.deps) &
                x(is,iray)=cmplx(zero, zero)

 11     continue
 10   continue

      return
      end
end module urfwr0c_mod
