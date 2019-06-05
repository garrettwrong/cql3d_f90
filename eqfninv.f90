module eqfninv_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use eqwrng_mod, only : eqwrng

  !---END USE


!
!

contains

      real(c_double) function eqfninv(psival,scalfct)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)



!..................................................................
!     This routine returns E given psival..
!..................................................................

      if (eqmodel.eq."power") then
        eqfninv=(psival/scalfct)**(1./eqpower)
      else
        call eqwrng(6)
      endif
      return
      end function eqfninv
      
      
end module eqfninv_mod
