module tdxin33d_mod

  !---BEGIN USE

  use profaxis_mod, only : profaxis

  !---END USE

!
!

contains

      subroutine tdxin33d(a,rya,klrz,expn1,expm1)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Fill in input arrays between center and edge parabolically.
!..................................................................



      dimension rya(0:klrz), a(0:klrz)
      data em90 /1.d-90/
      if (abs(a(0)) .le. em90) a(0)=em90
      dratio=a(1)/a(0)
      do 1 ll=1,klrz
        call profaxis(rn,expn1,expm1,dratio,rya(ll))
        a(ll)=a(0)*rn
 1    continue
      return
      end
end module tdxin33d_mod
