module pltinit_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine pltinit
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!................................................................
!     Initiates PGPLOT
!................................................................

!MPIINSERT_INCLUDE

!      integer pgbeg
      integer pgopen

      INTEGER PG_I, PG_L, PG_C1, PG_C2, PG_NC !YuP
      CHARACTER*16 PG_VAL ! YuP

!MPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return

!      write(*,*) 'PLTINIT-1'
!      Remember, pgbeg should be integer
!      ier=pgbeg(0,'plot.ps/VPS',1,1)
!      Remember, pgopen should be integer
      write(t_,1000) mnemonic(1:length_char(mnemonic))
 1000 format(a,".ps/VCPS") !YuP: was /VPS (vertical black&white)
                           ! Use /VCPS for vertical Color pages
      ier=PGOPEN(t_)
      CALL PGSCI(1)
      CALL PGSLW(lnwidth)
      write(*,*) 'PLTINIT-1 ier=1 is OK: ier=',ier
!      ier=pgbeg(0,'?',1,1)
!      if (ier.ne.1) write(*,*)
!     +              'Problem1 with initiating PGPLOT library'
!      ier=pgbeg(0,'?',1,1)
!      write(*,*) 'ier=',ier

      !YuP[2018-02-07] Added:
      !inquire PGPLOT general information: use PGQINF()
      !---------- First arg is the input:
      !'VERSION'     - version of PGPLOT software in use.
      !'STATE'       - status of PGPLOT ('OPEN' if a graphics device
      !                is open for output, 'CLOSED' otherwise).
      !'USER'        - the username associated with the calling program.
      !'NOW'         - current date and time (e.g., '17-FEB-1986 10:04').
      !'DEVICE'    * - current PGPLOT device or file.
      !'FILE'      * - current PGPLOT device or file.
      !'TYPE'      * - device-type of the current PGPLOT device.
      !'DEV/TYPE'  * - current PGPLOT device and type, in a form which
      !                is acceptable as an argument for PGBEG.
      !'HARDCOPY'  * - is the current device a hardcopy device? ('YES' or
      !                'NO').
      !'TERMINAL'  * - is the current device the user's interactive
      !                terminal? ('YES' or 'NO').
      !'CURSOR'    * - does the current device have a graphics cursor?
      !                ('YES' or 'NO').
      ! Two other arg. are outputs:
      CALL PGQINF('TYPE', PG_VAL, PG_L)
      WRITE (*,*) 'PGPLOT device type: ', PG_VAL(1:PG_L)
      CALL PGQINF('DEVICE', PG_VAL, PG_L)
      WRITE (*,*) 'PGPLOT device: ', PG_VAL(1:PG_L)
      CALL PGQINF('USER', PG_VAL, PG_L)
      WRITE (*,*) 'PGPLOT user: ', PG_VAL(1:PG_L)
      CALL PGQINF('NOW', PG_VAL, PG_L)
      WRITE (*,*) 'PGPLOT time now: ', PG_VAL(1:PG_L)

      !Inquire color index range:
      CALL PGQCIR(PG_C1, PG_C2)
      PG_NC = MAX(0, PG_C2-PG_C1+1)
      WRITE (*,*) 'Number of color indices used for image: ', PG_NC
      ! On Yuri's PC: printed --   PG_NC=240
      IF (PG_NC .LT.8) THEN
         WRITE (*,*) 'Not enough colors available on this device'
         STOP
      ELSE
         WRITE (*,*)
      END IF

      return
      end
end module pltinit_mod
