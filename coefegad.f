

c
c
      subroutine coefegad(k)
      use param_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     Adds in coefficient of energy loss term to fp coefficients.
c..................................................................

      include 'comm.h'

      if(k.gt.ngen)  return
      do 10  j=1,jx
        do 11  i=1,iy
          da(i,j)=da(i,j)+egylosa(i,j,k,indxlr_)
 11     continue
 10   continue
      return
      end
