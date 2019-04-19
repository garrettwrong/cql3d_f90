c
c
      subroutine wpvptb
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Reset some bounce average quantities for CQLP run
c..................................................................


c.......................................................................
c     1. Set lambda=vptb=cos(theta_0)*tau to 1.0
c.......................................................................

      do 100 i=1,iy
        vptb(i,lr_)=1.0
 100  continue

      return
      end
