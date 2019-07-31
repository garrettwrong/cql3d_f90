module tdtrcon_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdtrflx_mod, only : tdtrflx

  !---END USE

!
!

contains

      subroutine tdtrcon
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!..............................................................
!     Compute conservation constant.
!..............................................................
      data sgainr /0.0/
! jakub urban 110708: commented out for g95 compiler
! the above blanket save statement should do the job
!      save sgainr

!..............................................................
!     Compute original number of particles in tokamak.
!..............................................................
      if (n.eq.1) then
        total0=0.
        do 10 l=1,setup0%lrz
          ilr=setup0%lrindx(l)
          total0=xlndn0(ilr)/zmaxpsi(ilr)*dvol(ilr)+total0
 10     continue
      endif

!..............................................................
!     Compute total number of particles in device now.
!...............................................................

      total=0.
      do 30 l=1,setup0%lrz
        ilr=setup0%lrindx(l)
        do 40 k=1,ngen
          total=xlndn(k,ilr)/zmaxpsi(ilr)*dvol(ilr)+total
 40     continue
 30   continue

!................................................................
!     Call routine to compute number of particles lost at limiter
!     this time step.
!................................................................

      call tdtrflx
      sgainr=sgainr+flxout
      conserv=(total-total0-sgainr)/(total*.5+total0*.5)
      return
      end subroutine tdtrcon


end module tdtrcon_mod