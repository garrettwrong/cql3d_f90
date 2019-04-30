module cfpcoefc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use cfpcoefn_mod, only : cfpcoefn
  use cfpcoefr_mod, only : cfpcoefr

  !---END USE

!
!

contains

      subroutine cfpcoefc
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!..................................................................
!     Subroutine to handle coefc subroutines
!..................................................................
!
!
      if(relativ.eq."fully") then
        call cfpcoefr
      else
        call cfpcoefn
      endif
!
      return
      end
end module cfpcoefc_mod
