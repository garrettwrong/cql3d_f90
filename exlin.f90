module exlin_mod

!
!

contains

      real*8 function exlin(y1,y2,x1,x2,xx)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
!
!mnt  linear interpolation.
!
      exlin=(y2-y1)/(x2-x1)*(xx-x1)+y1
      return
      end
end module exlin_mod
