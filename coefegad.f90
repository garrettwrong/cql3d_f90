module coefegad_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE



!
!

contains

      subroutine coefegad(k)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Adds in coefficient of energy loss term to fp coefficients.
!..................................................................


      if(k.gt.ngen)  return
      do 10  j=1,jx
        do 11  i=1,iy
          da(i,j)=da(i,j)+egylosa(i,j,k,indxlr_)
 11     continue
 10   continue
      return
      end subroutine coefegad
      
end module coefegad_mod
