module eqfninv_mod

  !---BEGIN USE

  use eqwrng_mod, only : eqwrng

  !---END USE


!
!

contains

      real*8 function eqfninv(psival,scalfct)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)



!..................................................................
!     This routine returns E given psival..
!..................................................................

      if (eqmodel.eq."power") then
        eqfninv=(psival/scalfct)**(1./eqpower)
      else
        call eqwrng(6)
      endif
      return
      end
end module eqfninv_mod
