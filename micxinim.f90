module micxinim_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine micxinim
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!......................................................................
!     Creates a single integration coefficient. (if l_.eq.lmdpln_):
!     2*pi*sin(theta0)*|cos(theta0)|*tau*d(theta0), tau with no v.
!......................................................................


      do 10 i=1,iy_(l_)
        vpint(i,lr_)=cynt2(i,l_)*abs(coss(i,l_))*tau(i,lr_)
 10   continue

      return
      end
end module micxinim_mod
