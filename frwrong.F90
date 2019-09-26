!
!
      subroutine frwrong(kerr)
      use param_mod
      use cqlcomm_mod

      implicit none
      integer kerr ! input

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

! print error messages - on mpirank.eq.0 only
#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      if (kerr.eq.1 .or. kerr.eq.2) then
        WRITE(*,10010)
      elseif (kerr.eq.3) then
        WRITE(*,10030)
      elseif (kerr.eq.4) then
        WRITE(*,10040)
      elseif (kerr.eq.5) then
        WRITE(*,10050)
      else
        WRITE(*,10990)
      endif
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

10010 format("kspeci not properly input for use with fr-module")
10030 format("subroutine freyasou - loop 110")
10040 format("subroutine freyasou - rcon iteration -meth2")
10050 format("subroutine freyasou - luf failure i.gt.iy")
10990 format( "fr module-unspecified")
      stop 'frwrong:' ! stop at all MPI cores
      end subroutine frwrong
