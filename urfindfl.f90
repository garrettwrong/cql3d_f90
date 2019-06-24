module urfindfl_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfindfl
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine sets defaults for the "urf" module
!..................................................................

!..................................................................
!     Set defaults
!..................................................................

      !hahah go figures, the only unaccounted-for variables

      nrayn=1 !-YuP 101122: added here, instead of nraya in param.h;
              ! will be re-defined in urfsetup.
              ! Will be used to allocate arrays, mostly in urfalloc
      nrayelts=1 !-YuP: added here, instead of nrayelta in param.h;
              ! will be re-defined in urfsetup.
              ! will be used to allocate arrays, mostly in urfalloc.

!..................................................................
!     Setting nmods=nmodsa, for the time being. Generalize later.
!..................................................................

      nmods=nmodsa ! YuP-101220: should be mrfn, but not known yet


      return
      end subroutine urfindfl

end module urfindfl_mod
