module urfwrong_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine urfwrong(kerr)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)


!..................................................................
!     Flags errors in urf model - and terminates execution.
!..................................................................

!MPIINSERT_INCLUDE

! print error messages - on mpirank.eq.0 only
!MPIINSERT_IF_RANK_EQ_0
      if (kerr.eq.1) then
        WRITE(*,10010)
      elseif (kerr.eq.2) then
        WRITE(*,10020)
      elseif (kerr.eq.3) then
        WRITE(*,10030)
      elseif (kerr.eq.4) then
        WRITE(*,10040)
      elseif (kerr.eq.5) then
        WRITE(*,10050)
      elseif (kerr.eq.6) then
        WRITE(*,10060)
      elseif(kerr.eq.7) then
        WRITE(*,10070)
      else
        WRITE(*,10990)
      endif
!MPIINSERT_ENDIF_RANK

      stop 'urfwrong:' ! stop at all MPI cores

10010 format("subroutine urfbes: insufficient elements -  Bessel table")
10020 format("harmonic number...nharm(k) must be .le. nharma")
10030 format("subroutine urfb0 - unable to determine poloidal angle.")
10040 format("subroutine urfb0 - cannot locate the resonance region", &
        "for v-perp=0")
10050 format(" ")
10060 format("subroutine urfread: too many ray elements provided by", &
        /, "ray tracing code for some ray. Not enough RAM memory")
10070 format("failure in urfpack, can't find resonance region..")
10990 format("unspecified error originating from urf module.")
      return
      end
end module urfwrong_mod
