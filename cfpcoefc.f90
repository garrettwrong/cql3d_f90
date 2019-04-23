module cfpcoefc_mod

  !---BEGIN USE

  use cfpcoefn_mod, only : cfpcoefn
  use cfpcoefr_mod, only : cfpcoefr

  !---END USE

!
!

contains

      subroutine cfpcoefc
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
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
