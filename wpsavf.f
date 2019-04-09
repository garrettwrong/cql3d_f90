c
c
      subroutine wpsavf
      use param_mod
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     This routine saves various versions of the distribution
c     function.
c..............................................................

      include 'comm.h'
c.......................................................................

      if (relaxtsp.ne."enabled" .or. n.lt.nonavgf .or.
     +  n.gt.nofavgf) go to 200

c.......................................................................
cl    1. relaxtsp=enabled: mix f_n and f_n+1
c.......................................................................

      do 100 k=1,ngen
        do 110 j=0,jx+1
          do 120 l=0,ls+1
            do 130 i=0,iy_(l)+1
              fnp0(i,j,k,l)=0.5*(fnp0(i,j,k,l)+fnp1(i,j,k,l))
              fnp1(i,j,k,l)=fnp0(i,j,k,l)
 130        continue
 120      continue
 110    continue
 100  continue
      call dcopy(iyjx2*ngen*ls,fnp1(0,0,1,1),1,f(0,0,1,1),1)
      call dcopy(iyjx2*ngen*ls,fnp1(0,0,1,1),1,f_(0,0,1,1),1)

      return

c.......................................................................
cl    2.  f=f_n+1
c.......................................................................

 200  continue

      call dcopy(iyjx2*ngen*(ls+2),fnp1(0,0,1,0),1,fnp0(0,0,1,0),1)
      call dcopy(iyjx2*ngen*ls,fnp1(0,0,1,1),1,f(0,0,1,1),1)
      call dcopy(iyjx2*ngen*ls,fnp1(0,0,1,1),1,f_(0,0,1,1),1)

      return
      end
