module rf_mod

  !---BEGIN USE

  !---END USE


!*****************************************************************

!     FOR UNICOS VERSION RF ROUTINES ARE EXCISED: THEY RESIDE IN
!     CFS UNDER 714/precursr/rf


contains

      subroutine rf(action)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      character*8 action

      return
      end
end module rf_mod
