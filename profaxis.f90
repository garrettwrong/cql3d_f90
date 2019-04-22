module profaxis_mod

!
!

contains

      subroutine profaxis(rn,expn1,expm1,dratio,rova)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save
!---------------------------------------------------------------------
!     Expands "parabolic" profiles by computing the ratio rn
!     of the local parameter value at normalized radius rova to the central
!     value, given exponents expn1 and expm1 of the parabola, and
!     the ratio dratio of edge to central parameter value.
!---------------------------------------------------------------------

      if(abs(dratio-1.).lt.1.e-12) then
        rn=1.0
      else
        rn=dratio+(1.0-dratio)*(1.-rova**expn1)**expm1
      endif
      return
      end
end module profaxis_mod
