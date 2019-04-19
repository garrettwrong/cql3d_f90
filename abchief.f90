c
c
      subroutine abchief
      use param_mod
      use cqcomm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c
c      cputime=second()

c.......................................................................
c     2d code controlled by achief1, called at beginning of tdchief
c.......................................................................
      call tdchief
      return
      end
