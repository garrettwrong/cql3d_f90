module urfinitl_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfinitl
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine does some post-namelist-read initialization
!     and checking of input.
!..................................................................


!MPIINSERT_INCLUDE

!MPIINSERT_IF_RANK_EQ_0
      ! make plots on mpirank.eq.0 only
      if (noplots.ne."enabled1") then
      write(t_,1000)
 1000 format("Urf (lower hybrid, fast wave, ech, ebw...) parameters:")
      CALL PGMTXT('T',-7.,0.,0.,t_)

!      write(t_,1001) nrayn,nrayelts
 1001 format("====>NRAYn =",i5,"     ====>NRAYELTs = ",i5)
!      CALL PGMTXT('T',-8.,0.,0.,t_)

      write(t_,1002) nmodsa
 1002 format("====>NMODSA = ", i3)
      CALL PGMTXT('T',-9.,0.,0.,t_)
      endif
!MPIINSERT_ENDIF_RANK

!..................................................................
!     Return if urfmod.eq."disabled"
!..................................................................

      if(urfmod.eq."disabled") return


!..................................................................
!     If eqsource="tsc", then ensure that xbr [Prior to genray]
!     reads extended eqdsk:
!..................................................................

      if(eqsource.eq."tsc")  ieqbrurf=4

!..................................................................
!     Counter for calculation of urf-diffusion coefficients:
!..................................................................
      nurf=0

!.......................................................................
!     Re-run a few flux surfaces => "do not update delpwr" option
!.......................................................................
      if (urfrstrt .eq. "enabled") then
        nrfpwr=0
        nrfitr1=0
        nrfitr2=0
        nrfitr3=0
      endif

!.......................................................................
!     Check some input values
!.......................................................................

      if (nurftime.gt.nbctimea) stop "nurftime.gt.nbctimea"

!.......................................................................
!     Check that rftype() or lh/fw/ech wave-type designators are used,
!     but not both.
!     irftype=0, older method, lh/fw/ech wage-type designation
!             1, new method, via rftype(1:nmodsa).ne.'notset'
!.......................................................................

      irftype=0
      do i=1,nmodsa
         if (rftype(i).ne."notset" .and. irftype.eq.0) irftype=1
      enddo
      if ( (lh.eq."enabled" .or. fw.eq."enabled" .or. &
           ech.eq."enabled") .and. irftype.eq.1) then
         WRITE(*,*)"STOP: Can't set both lh/fw/ech and rftype()"
         STOP
      endif

!.......................................................................
!     If irftype=1, rfread="netcdf" is only choice
!     (If irftype=0, rfread can be either "text" or "netcdf".)
!.......................................................................

      if (irftype.eq.1 .and. rfread.ne."netcdf") then
         write(*,*)
         WRITE(*,*)"urfinitl STOP: Incompatible rftype and rfread"
         STOP
         write(*,*)
      endif

!.......................................................................
!     Count the number of rf types, mrf
!.......................................................................

      mrf=0
      if ( irftype.eq.0 ) then
         if (lh.eq."enabled") then
            mrf=mrf+1
            rftype(mrf)="lh"
         endif
         if (ech.eq."enabled") then
            mrf=mrf+1
            rftype(mrf)="ech"
         endif
         if (fw.eq."enabled") then
            mrf=mrf+1
            rftype(mrf)="fw"
         endif
      endif

      if ( irftype.eq.1 ) then
         do k=1,nmodsa
            if (rftype(k).ne."notset") mrf=mrf+1
            if (rftype(k).eq."notset") go to 5
         enddo
 5       continue
      endif

!     Make sure not attempting to damp/ql diffuse species which
!     is not present.
      do krf=1,mrf
         if (nrfspecies(krf).gt.ngen) STOP 'nrfspecies().gt.ngen'
      enddo

!...................................................................
!     For each wave type:
!     nharms() is the number of cyclotron harmonics calculated,
!     starting at nharm1().
!     (nharms.gt.1 .and. mrf.gt.1) is now permitted by the
!     storage scheme [BH060314].

!     mrfn is the number of wave "modes", that is, the sum over
!     wave types of the number of harmonics for each wave type,
!     i.e., each harmonic for each wave type is counted as a
!     separate "mode" or interation with the plasma distribution(s).
!
!     NOTE: It is possible that mrfn may be reset to a different
!           value in urfread, based on reading old rf input files
!           which had nharms set = 0. This perhaps justifies some
!           repeat coding in urfsetup [BH080919, but I think this
!           is completely obsolete.]
!...................................................................

      mrfn=0
      do k=1,mrf
         if (nharms(k).eq.0) then
            mrfn=mrfn+1   !  I.E., each rf type implies
                          !  at least 1 harmonic.
         else
            mrfn=mrfn+nharms(k)
         endif
      enddo

      if (mrfn.gt.nmodsa) then
         write(*,*)'urfinitl: mrfn>nmodsa.  mrfn,nmodsa=',mrfn,nmodsa
         STOP 'Increase nmodsa.'
      endif

!.......................................................................
!     Set up URF module file names, for rfread="netcdf".
!     rffile(1:nmodsa) are set by default to "notset",
!     and may be set in the namelist input.

!     If rffile(1)="mnemonic", rffile(1:3) is based
!     on the namelist input character variable mnemonic.
!.......................................................................

      if (rfread.eq."netcdf") then

!     irftype=0 case:
         krf=0
         if (ech.eq."enabled") then
            krf=krf+1
            if (rffile(krf).eq."notset") rffile(krf)="rayech.nc"
         endif
         if (fw.eq."enabled") then
            krf=krf+1
            if (rffile(krf).eq."notset") rffile(krf)="rayfw.nc"
         endif
         if (lh.eq."enabled") then
            krf=krf+1
            if (rffile(krf).eq."notset") rffile(krf)="raylh.nc"
         endif

         if (rffile(1).eq."mnemonic") then
            write(t_,110) mnemonic(1:length_char(mnemonic))
 110        format(a,"_rf.nc")
            rffile(1)=t_

!     Here we use achar(48+1)="1", etc., to name input files
            do i=1,mrf-1
               if (mrf.gt.10) then
                  WRITE(*,*)'urfinitl:  Expand calc of file name'
                  STOP
               endif
               write(t_,111) mnemonic(1:length_char(mnemonic)),i
 111           format(a,"_rf.",i1,".nc")
               rffile(i)=t_

               write(*,*)'urfinitl: i,rffile(i) =',i,rffile(i)
            enddo
         endif
      endif  !On rfread

      return
      end subroutine urfinitl

end module urfinitl_mod
