module wpvptb_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine wpvptb
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

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
      end subroutine wpvptb

end module wpvptb_mod
