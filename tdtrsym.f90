module tdtrsym_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine tdtrsym
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      do 5 k=1,ngen
        do 10 l=1,lrors-1
          do 20 i=itl_(l),iyh_(l)
            i_=iy_(l)+1-i
            do 30 j=1,jx
              frn(i,j,k,l)=(frn(i,j,k,l)+frn(i_,j,k,l))*.5
              frn(i_,j,k,l)=frn(i,j,k,l)
 30         continue
 20       continue
 10     continue
 5    continue
      return
      end
end module tdtrsym_mod
