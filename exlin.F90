module exlin_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      real(c_double) function exlin(y1,y2,x1,x2,xx)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!mnt  linear interpolation.
!
      exlin=(y2-y1)/(x2-x1)*(xx-x1)+y1
      return
      end function exlin

end module exlin_mod
