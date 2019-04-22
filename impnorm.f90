module impnorm_mod

!
!

contains

      subroutine impnorm(xnorm,a,rhs,nn)
      implicit integer (i-n), real*8 (a-h,o-z)

!..................................................................
!     This routine normalizes the matrix a so that the maximum
!     coefficient for each equation is of order 1.
!..................................................................
      dimension a(nn)

      xnorm=0.d0
      do 1 i=1,nn
        xnorm=xnorm+dabs(a(i))
 1    continue

      if (xnorm.gt.0.d0) then
      rhs=rhs/xnorm
      do 2 i=1,nn
        a(i)=a(i)/xnorm
 2    continue
      endif

      return
      end
end module impnorm_mod
