c
c
      subroutine coefsyad(k)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     Adds in contribution of synchrotron radiation to coefficients..
c..................................................................


      if (k .ne. kelecg .or. syncrad .eq. "disabled") return
      do 20 i=1,iy
        do 21 j=1,jx
          da(i,j)=da(i,j)+synca(i,j,indxlr_)
          dd(i,j)=dd(i,j)+syncd(i,j,indxlr_)
 21     continue
 20   continue
      return
      end
