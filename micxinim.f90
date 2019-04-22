module micxinim_mod

!
!

contains

      subroutine micxinim
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

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
