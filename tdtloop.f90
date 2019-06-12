module tdtloop_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aindfpa_mod, only : ainadjnl
  use aindfpa_mod, only : ainadjnl_fsetup_setup0
  use tdeqdsk_mod, only : tdeqdsk
  use tdtscout_mod, only : tdtscout
  use tdwritef_mod, only : tdwritef

  !---END USE

!
!

contains

  subroutine tdtloop
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      use tdeqdsk_mod, only : tdeqdsk
      use aindfpa_mod, only : ainadjnl_fsetup_setup0, ainadjnl
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Determines whether or not it is time to halt execution.
!..................................................................

      include 'name.h'
!MPIINSERT_INCLUDE

      character cptline*80
      integer getpid
      real :: cputime

      if (n .ge. nstop) go to 20 ! Finished; Finalize and Stop program.
      go to 30 ! Not finished yet; Proceed to return/end, next time step

 10   itloop=1
 20   continue
      call tdtscout
      call tdeqdsk

!MPIINSERT_IF_RANK_EQ_0
      ! make plots on mpirank.eq.0 only
      if (setup0%noplots.ne."enabled1") then
         call pgclos  ! PGCLOS
         PRINT *,'PGPLOT CLOSED at time step n=',n
      endif
!MPIINSERT_ENDIF_RANK

!..................................................................
!     Check-point this job, if ichkpnt.ne."disabled" into file ichkpnt.
!     This has been disabled, due to lack of generality of dropfile().
!..................................................................
      if (ichkpnt.ne."disabled") then
!BH080118        i=getpid(0)
        i=9999               !Dummy value
        write(cptline,100) i,ichkpnt
 100    format("chkpnt -p ",i5," -f ",a8)
        print *,cptline
!        i=dropfile(ichkpnt)
      endif

!.......................................................................
!     save distribution fnctn, restore original cqlinput, cputime, STOP
!.......................................................................

      if (setup0%nlwritf.ne."disabled") call tdwritef

!MPIINSERT_IF_RANK_EQ_0
         call ainadjnl(1)
         !restore cqlinput if &FSETUP was renamed to &SETUP0 (earlier)
         call ainadjnl_fsetup_setup0(1)
         !pause
!MPIINSERT_ENDIF_RANK

!MPIINSERT_BARRIER

      call cpu_time(cputime)
!MPIINSERT_IF_RANK_EQ_0
      PRINT *, 'tdtloop: >>> END OF CQL3D PROGRAM <<<'
      PRINT *, 'tdtloop: Execution time (seconds)', cputime
!MPIINSERT_ENDIF_RANK


 30   continue
      return
      end subroutine tdtloop


end module tdtloop_mod
