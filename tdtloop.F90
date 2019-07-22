! Copyright 2019 Garrett Wright, Princeton Plasma Physics Laboratory,
!    contracted by the U.S. Department of Energy (DE-AC02-09CH11466).
!
! This file is part of cql3d_f90. See LICENSE.
!
! cql3d_f90 is free software: you can redistribute it and/or modify it
! under the terms of the GNU Affero General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! cql3d_f90 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with cql3d_f90.  If not, see <https://www.gnu.org/licenses/>.

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

  subroutine tdtloop(nml_file)
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      use tdeqdsk_mod, only : tdeqdsk
      use aindfpa_mod, only : ainadjnl_fsetup_setup0, ainadjnl
      implicit integer (i-n), real(c_double) (a-h,o-z)
      character(len=*), intent(in), optional :: nml_file
      save

!..................................................................
!     Determines whether or not it is time to halt execution.
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

      character cptline*80
      integer getpid
      real :: cputime

      if (n .ge. nstop) go to 20 ! Finished; Finalize and Stop program.
      go to 30 ! Not finished yet; Proceed to return/end, next time step

 10   itloop=1
 20   continue
      call tdtscout
      call tdeqdsk

#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      ! make plots on mpirank.eq.0 only
      if (setup0%noplots.ne."enabled1") then
#ifndef NOPGPLOT
         call pgclos  ! PGCLOS
#endif
         PRINT *,'PGPLOT CLOSED at time step n=',n
      endif
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

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

#ifdef __MPI
      if(mpirank.eq.0) then
#endif
         if(present(nml_file)) call ainadjnl(1, nml_file)
         !restore cqlinput if &FSETUP was renamed to &SETUP0 (earlier)
         if(present(nml_file)) call ainadjnl_fsetup_setup0(1, nml_file)
         !pause
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

#ifdef __MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
#endif

      call cpu_time(cputime)
#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      PRINT *, 'tdtloop: >>> END OF CQL3D PROGRAM <<<'
      PRINT *, 'tdtloop: Execution time (seconds)', cputime
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif


 30   continue
      return
      end subroutine tdtloop


end module tdtloop_mod
