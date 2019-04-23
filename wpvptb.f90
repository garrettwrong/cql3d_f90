module wpvptb_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine wpvptb
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!..................................................................
!     Reset some bounce average quantities for CQLP run
!..................................................................


!.......................................................................
!     1. Set lambda=vptb=cos(theta_0)*tau to 1.0
!.......................................................................

      do 100 i=1,iy
        vptb(i,lr_)=1.0
 100  continue

      return
      end
end module wpvptb_mod
