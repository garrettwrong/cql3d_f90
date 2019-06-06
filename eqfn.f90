module eqfn_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use eqwrng_mod, only : eqwrng

  !---END USE

!
!

contains

      real(c_double) function eqfn(e,scalfct)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)


!..................................................................
!     This routine returns ad-hoc psi as a function of e.
!..................................................................

      if (eqmodel.eq."power") then
        if (e.eq.zero) then
          eqfn=one
        else
          eqfn=scalfct*e**eqpower
        endif
      else
        call eqwrng(5)
      endif
      return
      end function eqfn

end module eqfn_mod
