module ntloop_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

!

contains

      subroutine ntloop
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

      zdttot=1.0
      if (adimeth.eq."enabled" .and. transp.eq."enabled" &
        .and. n.ge.nonadi) zdttot=2.0
      timet=timet+zdttot*dtreff
      time_(l_)=timet

      return
      end subroutine ntloop


end module ntloop_mod
