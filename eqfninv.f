
c
c
      real*8 function eqfninv(psival,scalfct)
      use param_mod
      use cqcomm_mod
      implicit integer (i-n), real*8 (a-h,o-z)



c..................................................................
c     This routine returns E given psival..
c..................................................................

      if (eqmodel.eq."power") then
        eqfninv=(psival/scalfct)**(1./eqpower)
      else
        call eqwrng(6)
      endif
      return
      end
