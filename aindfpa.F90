module aindfpa_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE


contains

      subroutine ainadjnl(kopt)
      implicit integer (i-n), real(c_double) (a-h,o-z)
#ifdef __MPI
!MPI >>>
      include 'mpilib.h'
!MPI <<<
#endif
      Save inlmod

!..................................................................
!     kopt=0: Test if two &setup namelist sections in the cqlinput
!             If so, make copy of cqlinput to cqlinput_tmp
!             Modify cqlinput to &setup0 and &setup sections
!     kopt=1: Restore original cqlinput (if modification was made)
!
!BH070414:    This adjustment of namelist structure is to facilitate
!             writes and reads of namelist for the SWIM Integrated
!             Plasma Simulator (IPS), and to maintain backwards
!             compatibility.
!             Previously, there were two namelist sections both
!             designated &setup.  For SWIM, it has been necessary
!             to distinguish these (we use &setup0 and &setup),
!             for purposes of writing out the complete input namelist.
!..................................................................

      character*128  command
      character*8 line
      character*9 line0
      parameter (long_enough=100000)  !To accomodate namelist sections
                                       !written on one line
      character(len=long_enough) :: line1   !Automatic array
      character(len=long_enough+1) :: line2   !Automatic array

!.......................................................................

#ifdef __MPI
!MPI >>>
      if(mpirank.ne.0) return
!MPI <<<
#endif

      if (kopt.eq.0) then

         inlmod=0
         open(unit=4,file="cqlinput",delim='apostrophe',status="old")
 1       read(unit=4,fmt=1003,end=2) line
         if (line.eq." &setup" .or. line.eq." &SETUP" &
             .or. line.eq."&setup" .or. line.eq."&SETUP") then
            inlmod=inlmod+1
            if (inlmod.eq.2) go to 2
         endif
         go to 1
 2       continue
         WRITE(*,*)'inlmod = ',inlmod
         close(4)

 1003    format(a8)

!.......................................................................
!     Make copy if inlmod=2, and adjust nl structure
!.......................................................................

         if (inlmod.le.1) then
            goto 999            !Implies new namelist &setup0/&setup
                                !regimen, i.e., only one &setup.

         else

!     Copy cqlinput to cqlinput_tmp
            open(unit=4,file="cqlinput",delim='apostrophe',status="old")
            open(unit=5,file="cqlinput_tmp",delim='apostrophe', &
                 status="replace")
 3          read(unit=4,fmt='(a)',end=4) line1
            len_line1=len_trim(line1)
            if (len_line1.ge.(long_enough-1)) then
               WRITE(*,*)'len_line1,long_enough',len_line1,long_enough
               STOP 'Adjust long_enough'
            endif
            write(unit=5, fmt='(a)') trim(line1)
            go to 3
 4          continue
            close(4)
            close(5)

!     Reopen cqlinput as new file, and put adjusted namelist into it.
!     There are 2 &setup namelist sections.
            open(unit=4,file="cqlinput",delim='apostrophe', &
                 status="replace")
            open(unit=5,file="cqlinput_tmp",delim='apostrophe', &
                 status="old")
            ifirst=0
 5          read(unit=5,fmt='(a)',end=6) line1
            if ((line1(1:8).eq." &setup" .or. line1(1:8).eq." &SETUP" &
                .or.line1(1:7).eq."&setup" .or. line1(1:7).eq."&SETUP") &
                 .and. (ifirst.eq.0)) then
               ifirst=1
!BH080118               line2(1:9)=" &setup0"
               line2(1:8)=" &setup0"
!BH080118                line2(10:len(line1)+1)=line1(9:len(line1))
               line2(9:len(line1)+1)=line1(9:len(line1))
            else
               line2=line1(1:len(line1))
            endif
            write(unit=4, fmt='(a)') trim(line2)
            go to 5
 6          continue

         endif

         close(4)
         close(5)
         go to 999

      elseif (kopt.eq.1) then

!.......................................................................
!     If cqlinput changed, copy cqlinput_tmp back to cqlinput
!.......................................................................
         if (inlmod.eq.2) then
            open(unit=5,file="cqlinput",delim='apostrophe', &
                 status="replace")
            open(unit=4,file="cqlinput_tmp",delim='apostrophe', &
                 status="old")
 7          read(unit=4,fmt='(a)',end=8) line1
            if (len_trim(line1).ge.(long_enough-1)) &
                 STOP 'Adjust long_enough'
            write(unit=5, fmt='(a)') trim(line1)
            go to 7
 8          continue
            close(4)
            close(5)

         endif

      endif                     !on kopt

 999  return
      end subroutine ainadjnl

!
!
!==================================================================


      subroutine ainadjnl_fsetup_setup0(kopt)
      !YuP[2017]
      implicit integer (i-n), real(c_double) (a-h,o-z)
#ifdef __MPI
!MPI >>>
      include 'mpilib.h'
!MPI <<<
#endif
      Save inlmod

!..................................................................
!     kopt=0: Test if &fsetup namelist sections in the cqlinput
!             is present (instead of setup0, as in svn versions of CQL3D)
!             If so, make copy of cqlinput to cqlinput_tmp
!             Modify/rename $fsetup to &setup0
!     kopt=1: Restore original cqlinput (if modification was made)
!
!     The usage of &fsetup was introduced by MIT and PPPL
!     because of problems with Pathscale compiler
!(cannot read more than five letters in the name of namelist group)
!..................................................................

      character*128  command
      character*8 line
      character*9 line0
      parameter (long_enough=100000)  !To accomodate namelist sections
                                       !written on one line
      character(len=long_enough) :: line1   !Automatic array
      character(len=long_enough+1) :: line2   !Automatic array

!.......................................................................

#ifdef __MPI
!MPI >>>
      if(mpirank.ne.0) return
!MPI <<<
#endif

      if (kopt.eq.0) then

         inlmod=0
         open(unit=4,file="cqlinput",delim='apostrophe',status="old")
 1       read(unit=4,fmt=1003,end=2) line
         if (line.eq." &fsetup" .or. line.eq." &FSETUP" .or. &
             line.eq."&fsetup"  .or. line.eq."&FSETUP"  ) then
             inlmod=inlmod+1
           WRITE(*,*)'FSETUP namelist is present.'
           WRITE(*,*)'cqlinput is adjusted to change FSETUP to SETUP0.'
           WRITE(*,*)'After reading from the adjusted file,'
           WRITE(*,*)'  it will be restored back.'
           if (inlmod.eq.1) go to 2
         endif
         go to 1
 2       continue

         close(4)

 1003    format(a8)

!.......................................................................
!     Make copy if inlmod=1, and adjust nl structure
!.......................................................................

!     Copy cqlinput to cqlinput_tmp

         if (inlmod.le.0) then
!BH180613           WRITE(*,*)'FSETUP namelist is absent (probably using SETUP0)'
           goto 999 !Implies new namelist &setup0 is present (not $fsetup)

         else ! inlmod=1

            !Save the original file as cqlinput_tmp
            open(unit=4,file="cqlinput",delim='apostrophe',status="old")
            open(unit=5,file="cqlinput_tmp",delim='apostrophe', &
                 status="replace")
 3          read(unit=4,fmt='(a)',end=4) line1
            len_line1=len_trim(line1)
            if (len_line1.ge.(long_enough-1)) then
               WRITE(*,*)'len_line1,long_enough',len_line1,long_enough
               STOP 'Adjust long_enough'
            endif
            write(unit=5, fmt='(a)') trim(line1)
            go to 3
 4          continue
            close(4)
            close(5)

!     Reopen cqlinput as new file, and put adjusted namelist into it.
!     There are 2 &setup namelist sections.
            open(unit=4,file="cqlinput",delim='apostrophe', &
                 status="replace")
            open(unit=5,file="cqlinput_tmp",delim='apostrophe', &
                 status="old")
            ifirst=0
 5          read(unit=5,fmt='(a)',end=6) line1
            if ((line1(1:8).eq." &fsetup" .or. &
                line1(1:8).eq." &FSETUP"  .or. &
                line1(1:7).eq."&fsetup"   .or. &
                line1(1:7).eq."&FSETUP") &
                 .and. (ifirst.eq.0)) then
               ifirst=1
               line2(1:8)=" &setup0"
               line2(9:len(line1)+1)=line1(9:len(line1))
            else
               line2=line1(1:len(line1))
            endif
            write(unit=4, fmt='(a)') trim(line2)
            go to 5
 6          continue
            close(4)
            close(5)

         endif ! inlmod

         go to 999

      elseif (kopt.eq.1) then

!.......................................................................
!     If cqlinput changed, copy cqlinput_tmp back to cqlinput
!.......................................................................
         if (inlmod.eq.1) then
            open(unit=5,file="cqlinput",delim='apostrophe', &
                 status="replace")
            open(unit=4,file="cqlinput_tmp",delim='apostrophe', &
                 status="old")
 7          read(unit=4,fmt='(a)',end=8) line1
            if (len_trim(line1).ge.(long_enough-1)) &
                 STOP 'Adjust long_enough'
            write(unit=5, fmt='(a)') trim(line1)
            go to 7
 8          continue
            close(4)
            close(5)
         endif

      endif                     !on kopt

 999  return
      end subroutine ainadjnl_fsetup_setup0

!

!==================================================================


      subroutine ain_transcribe(filename)
      implicit integer (i-n), real(c_double) (a-h,o-z)
#ifdef __MPI
!MPI >>>
      include 'mpilib.h'
!MPI <<<
#endif

!..................................................................
!     Transcibe filename contents to the standard output unit.
!..................................................................

      character(len=*), intent(in) :: filename
      parameter(long_enough=1000000)
      character(len=long_enough) :: line1          !Automatic array
      logical logic1

#ifdef __MPI
!MPI >>>
      if(mpirank.ne.0) return
!MPI <<<
#endif

      WRITE(*,*)
      WRITE(*,*)  'ain_transcribe write of nl to stdout:'
      WRITE(*,*)  '[set nmlstout="disabled" to turn it off]'
      WRITE(*,*)

      max_length=0
      WRITE(*,*)'ain_transcibe: filename =',filename
      inquire(file=filename,iostat=kiostat,opened=logic1,number=inumber)
      WRITE(*,*)'ain_transcribe: inquire on ',filename, &
           ' opened=',logic1,'iostat=',kiostat,'unit=',inumber
      if(.NOT. logic1) then
         iunit = 20
         open(unit=iunit,file=filename,delim='apostrophe',status="old", &
              iostat=kiostat)
       endif
      if (kiostat.ne.0) WRITE(*,*)'ain_transcribe: kiostat=',kiostat
      if (kiostat.ne.0) STOP 'ain_transcribe: prblm with filename'
 3    read(unit=iunit,fmt='(a)',end=4) line1
      len_line1=len_trim(line1)
      if (len_line1.gt.max_length) max_length=len_line1
      WRITE(*,*) trim(line1)
      if (len_line1.ge.(long_enough-1)) then
         WRITE(*,*)'ain_transcribe,long_enough',len_line1,long_enough
         STOP 'Adjust long_enough'
      endif
      go to 3
 4    continue
      WRITE(*,*)'ain_transcribe:  max line length =',max_length
      close(20)
      WRITE(*,*)

      return
      end subroutine ain_transcribe


end module aindfpa_mod
