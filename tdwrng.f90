module tdwrng_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE



!-----------------------------------------------------------------

contains

      subroutine tdwrng(kerr)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!...
!
!mnt  diagnostic error messages
!
!...

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
      elseif (kerr.eq.7) then
        WRITE(*,10070)
      elseif (kerr.eq.8) then
        WRITE(*,10080)
      elseif (kerr.eq.9) then
        WRITE(*,10090)
      elseif (kerr.eq.10) then
        WRITE(*,10100)
      else
        WRITE(*,10990)
      endif
!MPIINSERT_ENDIF_RANK

      stop ! stop at all MPI cores

!      if (kerr.eq.0) return
10010 format("machine = toroidal only  for 3-d calc. in cqlinput")
10020 format("wdmodel =ech1 only for 3-d calc. in cqlinput")
10030 format("rya(ll).lt.rya(ll-1) for an ll; rzset error (cqlinput)")
10040 format ("subroutine tdsxr - thpol out of bounds")
10050 format("subroutine tdsxr - field line does not intersect viewing", &
        " angle")
10060 format("subroutine tdtodskr - hpdealc returns error message")
10070 format("subroutine ???? - iostatus returns error message")
10080 format("transport model- meshy=fixed_mu.  Insufficient resolution" &
        "in the trapped region for some l. Set tfac larger")
10090 format("transport model - transport mesh falls EXACTLTY on p/t",/, &
        "boundary. Perturb tfac")
10100 format("STOP in tdtry: iactst=abort forces stop. error>1.e-8")
10990 format("unspecified")
      return
      end subroutine tdwrng
      
end module tdwrng_mod
