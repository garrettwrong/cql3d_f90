module rf_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE


!*****************************************************************

!     FOR UNICOS VERSION RF ROUTINES ARE EXCISED: THEY RESIDE IN
!     CFS UNDER 714/precursr/rf


contains

      subroutine rf(action)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
      character*8 action

      return
      end
end module rf_mod
