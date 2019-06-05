module cfpmodbe_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine cfpmodbe(x,bk1,bk2)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Evaluate Bessel functions for relativistic calculations
!..................................................................

      if(x.ge.2.) goto 1000
      t=(x/3.75)**2
      bi0=1.+t*(3.5156229+t*(3.0899424+t*(1.2067492+ &
        t*(0.2659732+t*(0.0360768+t*0.0045813)))))
      bi1=.5+t*(0.87890594+t*(0.51498869+t*(0.15084934+ &
        t*(0.02658733+t*(0.00301532+t*0.00032411)))))
      bi1=x*bi1
      t=(x*0.5)**2
      bk0=-.57721566+t*(0.42278420+t*(0.23069756+t*(0.03488590+ &
        t*(0.00262698+t*(0.00010750+t*0.00000740)))))
!990131      bk0=-alog(.5*x)*bi0+bk0
      algx=log(.5*x)
      bk0=-algx*bi0+bk0
      bk1=1.+t*(0.15443144+t*(-.67278579+t*(-.18156897+ &
        t*(-.01919402+t*(-.00110404-t*0.00046860)))))
!990131      bk1=alog(.5*x)*bi1+bk1/x
      bk1=algx*bi1+bk1/x
      exp_x=exp(x)
      bk0=exp_x*bk0
      bk1=exp_x*bk1
      goto 2000

!..................................................................
!     These calculate BK0 and BK1 if x>2.
!..................................................................

 1000 continue
      exsx=1./sqrt(x)
      t=2./x
      bk0=1.25331414+t*(-.07832358+t*(0.02189568+t*(-.01062446+ &
        t*(0.00587872+t*(-.0025154+t*0.00053208)))))
      bk0=bk0*exsx
      bk1=1.25331414+t*(0.23498619+t*(-.03655620+t*(0.01504268+ &
        t*(-.00780353+t*(0.00325614-t*0.00068245)))))
      bk1=bk1*exsx
 2000 continue
      bk2=bk0+2.*bk1/x
      return
      end subroutine cfpmodbe
      
end module cfpmodbe_mod
