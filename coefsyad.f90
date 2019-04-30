module coefsyad_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine coefsyad(k)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Adds in contribution of synchrotron radiation to coefficients..
!..................................................................


      if (k .ne. kelecg .or. syncrad .eq. "disabled") return
      do 20 i=1,iy
        do 21 j=1,jx
          da(i,j)=da(i,j)+synca(i,j,indxlr_)
          dd(i,j)=dd(i,j)+syncd(i,j,indxlr_)
 21     continue
 20   continue
      return
      end
end module coefsyad_mod
