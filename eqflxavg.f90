module eqflxavg_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use eqorbit_mod, only : eqorbit

  !---END USE

!
!

contains

      subroutine eqflxavg(epsicon_,a,flxavg_,flxavgd_)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      dimension a(*)
!dir$ nobounds
!..................................................................
!     This routine returns the flux surface average of a.
!     It works for both updown and non-up-down symmetric cases.
!..................................................................

      if (eqorb.eq."enabled") then
        call eqorbit(epsicon_)
      endif
      sum1=0.
      sum2=0.

!..................................................................
!     eqdell=dl; eqbpol=B-poloidal; both defined on constant phi
!     flux surface.
!..................................................................
      do 10 l=2,lorbit_
        sum1=sum1+eqdell_(l)/eqbpol_(l)
        sum2=sum2+(a(l)+a(l-1))*.5*eqdell_(l)/eqbpol_(l)
 10   continue
      flxavg_=sum2/sum1
      flxavgd_=sum1
      return
      end subroutine eqflxavg
      
      
end module eqflxavg_mod
