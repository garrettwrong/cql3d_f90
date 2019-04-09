c
c
      subroutine cfpsymt
      use param_mod
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Symmetrize about pi/2  (f.p. coefficients)
c..................................................................

      save
      include 'comm.h'


      do 1 k=1,ngen
        do 2 i=1,iyh
          ii=iy+1-i
          do 3 j=1,jx
            cal(ii,j,k,l_)=cal(i,j,k,l_)
            cbl(ii,j,k,l_)=cbl(i,j,k,l_)
            cfl(ii,j,k,l_)=cfl(i,j,k,l_)
            ccl(ii,j,k,l_)=-ccl(i,j,k,l_)
            cdl(ii,j,k,l_)=-cdl(i,j,k,l_)
            cel(ii,j,k,l_)=-cel(i,j,k,l_)
 3        continue
 2      continue
 1    continue
      return
      end
