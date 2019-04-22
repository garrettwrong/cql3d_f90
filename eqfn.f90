module eqfn_mod

!
!

contains

      real*8 function eqfn(e,scalfct)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)


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
      end
end module eqfn_mod
