c
c
      subroutine cfpcoefc
      use param_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c..................................................................
c     Subroutine to handle coefc subroutines
c..................................................................
      include 'comm.h'
c
c
      if(relativ.eq."fully") then
        call cfpcoefr
      else
        call cfpcoefn
      endif
c
      return
      end
