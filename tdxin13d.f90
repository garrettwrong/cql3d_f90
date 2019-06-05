module tdxin13d_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use profaxis_mod, only : profaxis

  !---END USE


!
!

contains

      subroutine tdxin13d(a,rya,klrz,m,k,expn1,expm1)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!.............................................................
!     this is a utility parabolic "fill in" routine
!.............................................................

      dimension rya(0:klrz), a(m,0:klrz)
      data em90 /1.d-90/
      if (abs(a(k,0)).le. em90) a(k,0)=em90
      dratio=a(k,1)/a(k,0)
      do 1 ll=1,klrz
        call profaxis(rn,expn1,expm1,dratio,rya(ll))
        a(k,ll)=a(k,0)*rn
 1    continue
      return
      end subroutine tdxin13d
      
end module tdxin13d_mod
