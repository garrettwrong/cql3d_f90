module ainplt_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine ainplt
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

      REAL RILIN

!....................................................
!     This routine plots out  data read in by code.
!....................................................

!MPIINSERT_INCLUDE

      character*300 line, line_
      character*24 text24
      character*300 text256,text256_
      integer stime
!     text(1:4) is charcter*8, from comm.h


!MPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return

      if (pltinput .eq. "disabled") go to 3
!%OS  call gxglfr(0)
      iunit=2
      ilen=0
      n_count=0
      open(unit=iunit,file='cqlinput',delim='apostrophe',status='old')
!     ctime is a unix library function which
!     returns the time as a character*24 variable.
!     Doesn't work with pathscale compiler.
      WRITE(*,*)
      WRITE(*,*)'ainplt:  If special_calls.eq."enabled" (default nml)'
      WRITE(*,*)'ainplt:    then code will bomb with compilers not'
      WRITE(*,*)'ainplt:    implementing SYSTEM call'
      WRITE(*,*)
!BH111102      if (special_calls.eq.'enabled') then
         call GET_DATE_TIME (text24)
!BH111102      else
!BH111102         text24='special_calls.ne.enabled'
!BH111102      endif

      CALL PGPAGE
      write(t_,2001)
      RILIN=0.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,2002)
      RILIN=RILIN+2.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,2003)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,2004)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,2005)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,2006) text24
      RILIN=RILIN+2.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!     Get, write, and plot machine characteristics
!     by writing a file, reading it, and removing.
!     special_calls.ne.enabled branches around the system calls,
!     which are not enabled for some systems.  Could populate
!     uname_output and pwd_output files extenal to this code
!     using a script which then invokes this code (Ed D'Azevedo).

      if (special_calls.eq.'enabled') then
         call system('uname -a > uname_output')
!_cray      call ishell('uname -a > uname_output')
         open(unit=13,file='uname_output',delim='apostrophe', &
              status='old')
         read(13,100) line
         close(unit=13)
         call system('rm uname_output')
!_cray      call ishell('rm uname_output')

         call system('pwd > pwd_output')
!_cray      call ishell('pwd > pwd_output')
         open(unit=13,file='pwd_output',delim='apostrophe',status='old')
         read(13,101) text256
         close(unit=13)
         call system('rm pwd_output')
!_cray      call ishell('rm pwd_output')

      elseif (special_calls.eq.'external') then
         open(unit=13,file='uname_output',delim='apostrophe', &
              status='old')
         read(13,100) line
         close(unit=13)
         open(unit=13,file='pwd_output',delim='apostrophe',status='old')
         read(13,101) text256
         close(unit=13)

      else
         line='special_calls.ne.enabled'
         text256='special_calls.ne.enabled, possibly use script'

      endif                     ! on special_calls

      WRITE(*,*) ' '
      WRITE(*,2006) text24
      WRITE(*,*) ' '
!     write MACHINE:
      WRITE(*,2007)
      lenmac=len_trim(line)
!      WRITE(*,100) line
      do i=1,(lenmac/60+1)
         write(*,102)line(1+(i-1)*60:i*60)
      enddo
!     write PWD:
      WRITE(*,2008)
      lenpwd=len_trim(text256)
!      write(*,*)'lenpwd=',lenpwd
      do i=1,(lenpwd/60+1)
         write(*,102)text256(1+(i-1)*60:i*60)
      enddo
!      WRITE(*,101) text256

      WRITE(*,2009) version
      WRITE(*,*) ' '
 100  format(a100)
 101  format(a256)
 102  format(a60)



      write(t_,2007)
      RILIN=RILIN+2.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
!      write(line_,1002) line
!      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,line_)
      do i=1,(lenmac/60+1)
         write(line_,1004) line(1+(i-1)*60:i*60)
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,line_)
      enddo
      write(t_,2008)
      RILIN=RILIN+2.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      do i=1,(lenpwd/60+1)
         write(text256_,1004) text256(1+(i-1)*60:i*60)
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,text256_)
      enddo
      write(line,2009) trim(version)
      RILIN=RILIN+2.
      CALL PGMTXT('T',-RILIN,0.,0.,line)

!     PGQINF('VERSION', enquires pgplot version being used.
      CALL PGQINF('VERSION',t_,ilength)
      if (ilength.gt.100) stop 'unlikely ilength problem in ainplt'
      write(line,2010) t_(1:ilength)
      RILIN=RILIN+2.
      CALL PGMTXT('T',-RILIN,0.,0.,line)

      WRITE(*,*) line
!
      CALL PGPAGE


 2001 format("OUTPUT FROM THE FOKKER-PLANCK CODE CQL3D.")
 2002 format("FOR QUESTIONS CONTACT ")
 2003 format("BOB HARVEY - (858)509-2131, bobh@compxco.com")
 2004 format("CQL3D IS A PRODUCT OF NERSC/GA/EPFL/CompX")
 2005 format("COLLABORATION.")
 2006 format("DATE/TIME is ",a)
 2007 format("MACHINE:")
 2008 format("PWD:")
 2009 format("CQL3D VERSION:  ",a)
 2010 format("PGPLOT VERSION: ",a)
!......................................................................
!     transcribe namelist data to plot file if pltinput .ne. "disabled"
!......................................................................

!$$$      call gscpvs(.01,1.)
      RILIN=-1.
      n_count=0
      ilnperpag= 50 ! lines per page           ! 41 originally
 1    continue   ! loop in reading lines in cqlinput file
      !YuP/was: read(iunit,1003,iostat=istat) line ! 80 symbols
      read(iunit,1005,iostat=istat) line !YuP[2018-01-04] 300 symbols
      !YuP: read up to 300 symbols from given line,
      !     then wrap the text in/over 60 symbols.
      !Note: character*300 line, line_
      lenmac=len_trim(line)
      iend=(lenmac-1)/60 +1  ! print up to 60 symbols,
      !then wrap the text over to additional lines.
      !For short lines (1-60 symbols), iend=1.
      if (line.eq."end".or.(istat.ne.0)) go to 3 ! finish reading.
      n_count=n_count+1 !If text is short, it is printed as a single line
      ! Also count extra lines (beyond i=1) to be added
      ! during "wrapping" of a long text:
      n_count=n_count+(iend-1)
      !write(*,*) lenmac,iend,n_count, '  line=',trim(line)
      !pause
      ! Start new page when the number of lines is above ilnperpag
      if (n_count.ge.ilnperpag) then
        ! This textline (possibly wrapped into iend lines)
        ! will be printed at the new page, so it will take iend lines.
        n_count=iend !!1 ! restart counting from new page
        RILIN=-1.
        CALL PGPAGE
      endif
      !The "wrapping" procedure:
      do i=1,iend
         write(line_,1004) line(1+(i-1)*60:i*60)
         !i=1: line(1  : 60) !for short lines (1-60 symbols), iend=1
         !i=2: line(61 : 120)  ! for lines with 61-120 symbols, iend=2
         !i=3: line(121: 180)
         !i=4: line(181: 240)
         !i=5: line(241: 300)  [max possible i, with present setup]
         !Note -- it is declared: character*300 line, line_
         RILIN=RILIN+1. !shift the text position down by one character height
         CALL PGMTXT('T',-RILIN,0.,0.,line_)
      enddo
      !YuP: old version, without wrapping:
      !write(line_,1002) line
      !RILIN=RILIN+1.
      !CALL PGMTXT('T',-RILIN,0.,0.,line_)
      go to 1 ! back to next line
 3    close(unit=iunit)

 1002 format(a80)
 1003 format(a80)
 1004 format(a60)
 1005 format(a300)
      return
      end



      subroutine GET_DATE_TIME (result)

!     get current date and time as a 24-character ASCII string

      implicit none

!BH111102      integer TIME
!BH111102      character*24 CTIME*24, result*24

!BH111102      result=CTIME (TIME () )

!BH111102:  Replacing with f90 intrinsic subroutine

      character*8 date
      character*10 time
      character*5 zone
      character*24 result
      integer values(8)

      call date_and_time(date,time,zone,values)

      result=date(1:4)//'/'//date(5:6)//'/'//date(7:8) &
           //' '//time(1:2)//':'//time(3:4)//' '//time(5:10) &
           //'s'

      return
      end
end module ainplt_mod
