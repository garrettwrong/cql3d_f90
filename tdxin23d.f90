module tdxin23d_mod

  !---BEGIN USE

  use profaxis_mod, only : profaxis

  !---END USE

!
!

contains

      subroutine tdxin23d(a,rya,klrz,ngn,nso,k,kk,expn1,expm1)
      use param_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!..................................................................
!     Fill in input arrays between center and edge parabolically.
!..................................................................


      dimension rya(0:klrz), a(ngn,nso,0:klrz)
      data em90 /1.d-90/
      if (abs(a(k,kk,0)) .le. em90) a(k,kk,0)=em90
      dratio=a(k,kk,1)/a(k,kk,0)
      do 1 ll=1,klrz
        call profaxis(rn,expn1,expm1,dratio,rya(ll))
        a(k,kk,ll)=a(k,kk,0)*rn
 1    continue
      return
      end
end module tdxin23d_mod
