!
!
      subroutine frset(lrz,noplots,nmlstout)
      use bcast_mod, only : bcast
      use param_mod
      use zfreya_mod, only : random_my

      implicit none
      integer lrz ! input
      
      include 'frcomm.h'
      include 'frname.h' ! contains namelist variables (*,frsetup) 
!MPIINSERT_INCLUDE

      character*1024 t_ ! local
      character*8 noplots,nmlstout ! input
      
      integer j,istat ! local
      ! nion,ibion ! arg of frinitz(), to be found; declared in frcomm.h
      real(c_double) :: dummy,fd ! local

      REAL RILIN

!..................................................................
!     frset initializes non-namelist variables.
! BH081020c     Also, resets namelist variable nprim=1, if not =1.
!..................................................................



!..................................................................
!     Jump out if lrz=1
!..................................................................

      if (lrz.eq.1) then
         ranseed=7**7
         read(2,frsetup)
         if (nmlstout.eq."enabled") write(6,frsetup)
         dummy=RANDOM_my(ranseed)
         goto 300 !-> allocate vars. and return
      endif

      ncorin=75
      do 1710 j=1,kprim
        znipm(j)=0.0d0
 1710 atwpm(j)=0.0d0
      do 1720 j=1,kimp
        iz(j)=0
        atwim(j)=0.0d0
 1720 zniim(j)=0.0d0

      ncont=30

!..................................................................
!     Set defaults for beam smoothing (none)...
!..................................................................

      smooth=.05

!..................................................................


!..................................................................
!     Read in data for NFREYA
!..................................................................

      read(2,frsetup)
      if (nmlstout.eq."enabled") write(6,frsetup)

      if (frmod.eq."disabled") return


!BH081020      if (nprim.ne.1) then
!BH081020         nprim=1
!BH081020         write(*,210)
!BH081020      endif
!BH081020 210  format('WARNING: Resetting nprim=1. Check FREYA if want nprim>1')


!     Check nbeams not too large.
      if (nbeams .gt. kb) then
         WRITE(*,*)
         WRITE(*,*)'Namelist nbeams.gt.kb:  increase parameter kb'
         stop
         write(*,*)
      endif

      if (multiply.ne."disabled".and.multiplyn.eq.0) then
         WRITE(*,214)
 214     format(//,'Must set multiplyn, if multiply.ne.disabled',//)
         stop 'Must set multiplyn'
      endif
!
!     ION PARAMETERS
!--------------------------------------------------------------------
!     nprim       Number of primary ion species
!     nimp        Number of impurity ion species
!     namep(i)    Name of ith primary ion species
!     'h' , protons
!     'd' , deuterons
!     't' , tritons
!     'dt', mixture of d and t
!     'he', thermal alphas
!     namei(i)    Name of ith impurity ion species
!     'he', helium
!     'c' , carbon
!     'o' , oxygen
!     'si', silicon
!     'ar', argon
!     'ti', titanium
!     'cr', chromium
!     'fe', iron
!     'ni', nickel
!     'kr', krypton
!     'mo', molybdenum
!     'w' , tungsten
!     fd          Number fraction of deuterons in d-t mixture
!--------------------------------------------------------------------
      fd=.5d0   !No longer namelist input

!..................................................................
!     Call subroutine linked to CQL3D common blocks for more initializat
!..................................................................

      call frinitz(nprim,nimp,nion,ibion,namep,namei,atw,fd,smooth)


!..................................................................
!     Initiate random number generator for NFREYA
!..................................................................
!     See comments in freyasou.f/subroutine frnnoa.
!990131      call ranset(ranseed)

!BH131015: ranseed is initialized in frinitl to 7**7
!BH131015: This call to RANDOM_my initialized a RN sequence.
      dummy=RANDOM_my(ranseed)

!     Reset npart, for nubeam list case
      if (read_birth_pts.eq."enabled") npart=nbirth_pts


!MPIINSERT_IF_RANK_EQ_0
      ! make plots on mpirank.eq.0 only
      if (noplots .ne. "enabled1") then
#ifndef NOPGPLOT
         write(t_,1000)
 1000 format("FR (freya beam deposition) model parameters:")
      RILIN=11.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,1001)
 1001 format("npart is the number of ions launched")
      RILIN=12.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,1002) npart
 1002 format("====>npart = ",i7)
      RILIN=13.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      endif
!MPIINSERT_ENDIF_RANK


!..................................................................
!     Set constants (for real(c_double))
!..................................................................
 300  continue

      two=2.d0
      half=0.5d0

      allocate(rpts(npart),STAT=istat)
      write(*,*)'frset  rpts: istat=',istat
      if(istat.eq.0) call bcast(rpts,zero,npart)

      allocate(xpts(npart),STAT=istat)
      write(*,*)'frset  xpts: istat=',istat
      if(istat.eq.0) call bcast(xpts,zero,npart)

      allocate(ypts(npart),STAT=istat)
      write(*,*)'frset  ypts: istat=',istat
      if(istat.eq.0) call bcast(ypts,zero,npart)

      allocate(zpts(npart),STAT=istat)
      write(*,*)'frset  zpts: istat=',istat
      if(istat.eq.0) call bcast(zpts,zero,npart)

      allocate(vx(npart),STAT=istat)
      write(*,*)'frset  vx: istat=',istat
      if(istat.eq.0) call bcast(vx,zero,npart)

      allocate(vy(npart),STAT=istat)
      write(*,*)'frset  vy: istat=',istat
      if(istat.eq.0) call bcast(vy,zero,npart)

      allocate(vz(npart),STAT=istat)
      write(*,*)'frset  vz: istat=',istat
      if(istat.eq.0) call bcast(vz,zero,npart)

      return
      end subroutine frset
