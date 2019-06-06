module urfread_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use r8subs_mod, only : dcopy
  use urfread__mod, only : urfread_

  !---END USE

!
!

contains

      subroutine urfread
      use param_mod
      use cqlcomm_mod
      use netcdfrf_mod, only : netcdfrf
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This routine reads disk files generated by ray tracing codes,
!     in two possible formats specified by rfread="netcdf" or "text".
!
!     For rf wave types with any of nharms(k),k=1,mrf that are
!     equal to 0, then nharms(k) is set equal to 1, and
!     nharm1(k) is obtained from the ray data file.
!
!     First, ray data is read for each wave type.
!     [Each wave "type" is read from its own ray data file.]
!     Then ray data is copied to additional ray data sets for
!     cases where nharms(k).gt.1.
!.......................................................................

!MPIINSERT_INCLUDE


!MPIINSERT_IF_RANK_EQ_0

      krf=0  !  Counter for wave types.
!     irftype (0 for old input method, 1 for new) is dtrmnd in urfinitl.

      if (irftype.eq.0) then    !Old type input of wave type

      if (ech.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(2,krf)
         else
            open(unit=23,file='rayech',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               call urfread_(krf,23)
               close(unit=23)
            else
               WRITE(*,*)'urfread: check that "rayech" file is present'
               stop
            endif
         endif
         if (nharms(krf).eq.0) then
            nharms(krf)=1
            nharm1(krf)=nharm(krf)   !nharm(krf) read in from ray data.
                                     !This option is for backwards
                                     !compatibility.
         endif
      endif
      if (fw.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(2,krf)
         else
            open(unit=24,file='rayfw',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               call urfread_(krf,24)
               close(unit=24)
            else
               write(*,*)'urfread: check that "rayfw" file is present'
               stop
            endif
         endif
         if (nharms(krf).eq.0) then
            nharms(krf)=1
            nharm1(krf)=nharm(krf)   !nharm(krf) read in from ray data.
         endif
      endif
      if (lh.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(2,krf)
         else
            open(unit=20,file='raylh',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               call urfread_(krf,20)
               close(unit=20)
            else
               write(*,*)'urfread: check that "raylh" file is present'
               stop
            endif
         endif
         if(partner.eq.'lsc') then
            do is=1,nrayelt(iray,krf)
               wnpar(is,iray,krf)=-wnpar(is,iray,krf)
            enddo
         endif
         if (nharms(krf).eq.0) then
            nharms(krf)=1
            nharm1(krf)=nharm(krf)   !nharm(krf) read in from ray data.
         endif
      endif

      endif   ! on irftype.eq.0

      if (irftype.eq.1) then     !New, more general input of wave types.

         do krf=1,mrf
            if (rftype(krf).eq."ec") then
!BH080911         krf=krf+1
!BH080911         Move krf down.
               call netcdfrf(2,krf)
            endif
            if (rftype(krf).eq."fw") then
!BH080911         krf=krf+1
               call netcdfrf(2,krf)
            endif
            if (rftype(krf).eq."lh") then
!BH080911         krf=krf+1
               call netcdfrf(2,krf)
               if(partner.eq.'lsc') then
                  do is=1,nrayelt(iray,krf)
                     wnpar(is,iray,krf)=-wnpar(is,iray,krf)
                  enddo
               endif
            endif
         enddo

      endif  ! on irftype.eq.1

!MPIINSERT_ENDIF_RANK


!MPIINSERT_BCAST_NRAYN
!MPIINSERT_BCAST_NRAYELTS
!MPIINSERT_BCAST_RAYS_DATA
!MPIINSERT_BARRIER


!.......................................................................
!     irftype indicates old (=0) or new (=1) rf data input method.
!     If this is a multispecies run, then irftype must .eq.1
!     for input to more than the first species.
!     For irftype.eq.1 .and .ngen.ge.2, determine if multiple
!     ray types are to be applied separately to each species,
!     or if they are the same rays applied simultaneously ("combined")
!     to multiple species.
!     This issue arises when the same ray data in multiple ray input
!     files is to be applied to different general species (ngen.ge.2),
!     but the damping along the rays is due to the combined effect
!     of the relevant general species.
!     Determine this by whether the input files rffile() are
!     the same file, or not (as determined from the 1st two species).
!     (This logic could be readily generalized to a combination of
!     "separate" and "combined" application of damping and diffusion.)
!.......................................................................

!.......................................................................
      !irffile() is indicator whether input rffile() values are
      !separate (.ne.each other), or duplicates for which power absorption
      !is then combined.
      !irffile()="separate"  !Is default.
                          !Apply damping to ray data files either
                          !as "separate" or combined (damping on
                          !more than one species).
                          !"combine1" is first ray file to be combined,
                          !"combine2" is subsequent identical ray files.
!.......................................................................

      do krf=1,mrf
         irffile(krf)="separate"
      enddo
      if (ngen.ge.2) then  !Assume ray files separate if only 1 species.
      do krf=1,mrf
         if (rffile(krf).ne."notset") then  !should be case for krf=1,mrf

            if (krf.eq.1) then
               irffile(krf)="separate"
               if(mrf.gt.1) then
                  if (rffile(krf).eq.rffile(krf+1)) &
                                 irffile(krf)="combine1"
               endif

            elseif (krf.gt.1 .and. krf.lt. mrf) then
               irffile(krf)="separate"
               if(rffile(krf).eq.rffile(krf-1)) then
                  irffile(krf)="combine2"
               elseif(rffile(krf).eq.rffile(krf+1)) then
                  irffile(krf)="combine1"
               endif

            elseif (krf.eq.mrf) then
               irffile(krf)="separate"
               if(rffile(krf).eq.rffile(krf-1)) &
                                   irffile(krf)="combine2"

            endif  !On krf



         endif                  ! on "notset"
      enddo                     ! on krf
      endif                     ! on ngen

!.......................................................................
!
!     Set up table krfn(1:mrfn), points to wave type index for each mode.
!     Set up table irfn(1:mrf), points to the wave mode index of
!       the lowest harmonic for each wave type.
!     Set up nharm(1:mrfn).
!     REMEMBER:
!     mrf is the number of RF types being calculated.
!     mrfn is the number of wave "modes", that is, the sum over
!     wave types of the number of harmonics for each wave type

!     NOTE:  Cannot set these arrays earlier in code because
!            nharm is possibly read above from ray data files,
!            giving value for nharm1.
!
!.......................................................................

      kk=0
      do krf=1,mrf
         do kkk=1,nharms(krf)
            kk=kk+1
            krfn(kk)=krf
            if (kkk.eq.1) irfn(krf)=kk
            nharm(kk)=nharm1(krf)+(kkk-1)
         enddo
      enddo

!.......................................................................
!     irfm(1:mrf) is total number of modes for each distinct (i.e., not
!     identical) rffile, irffile(krf)="separate".
!     For "combined" cases (i.e., identical rffile files), it gives
!     total number of modes partially summed over the combined case files.
!.......................................................................


      do krf=1,mrf
         if (irffile(krf).eq."separate") then
            irfm(krf)=nharms(krf)
         elseif (irffile(krf).eq."combine1") then
            irfm(krf)=nharms(krf)
            kkrf=krf
         elseif(irffile(krf).eq."combine2") then
            irfm(kkrf)=irfm(kkrf)+nharms(krf)
            irfm(krf)=irfm(kkrf)
         endif
      enddo



      write(*,*)
      write(*,*)'urfread: mrf = ',mrf
      write(*,*)'urfread: mrfn = ',mrfn
      write(*,*)'urfread: irfn = ',(irfn(krf),krf=1,mrf)
      write(*,*)'urfread: irfm = ',(irfm(krf),krf=1,mrf)
      write(*,*)'urfread: irffile = ',(irffile(krf),krf=1,mrf)
      write(*,*)'urfread: krfn = ',(krfn(krf),krf=1,mrfn)
      write(*,*)'urfread: nharm1 = ',(nharm1(krf),krf=1,mrf)
      write(*,*)'urfread: nharms = ',(nharms(krf),krf=1,mrf)
      write(*,*)'urfread: nharm = ',(nharm(krf),krf=1,mrfn)


!     Duplicate ray data in sets nharms(krf) long using the first ray data
!     for each ray type krf, if there is more than one harmonic:
!     Need to do this duplication in reverse order, k=mrfn,1,-1
!     in order not to overwrite data initially stored as krf=1,mrf.

      do k=mrfn,1,-1

         nray(k)=nray(krfn(k))
         freqcy(k)=freqcy(krfn(k))
         omega(k)=omega(krfn(k))

         do 30 iray=1,nray(krfn(k))
            nrayelt(iray,k)=nrayelt(iray,krfn(k))
            istarts(iray,k)=istarts(iray,krfn(k))
            call dcopy(nrayelt(iray,k),ws(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
!BH081106     &         ws(1,iray,k),krfn(k))   Causes big bug when
!BH081106                                krfn(k).ne.1: i.e., neg B0.
                 ws(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),seikon(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 seikon(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),spsi(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 spsi(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),wr(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 wr(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),wphi(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 wphi(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),wz(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 wz(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),wnpar(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 wnpar(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),wnper(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 wnper(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),delpwr(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 delpwr(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),sdpwr(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 sdpwr(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),wdnpar(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 wdnpar(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),fluxn(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 fluxn(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),sbtot(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 sbtot(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),sene(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 sene(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),salphac(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 salphac(1:nrayelt(iray,k),iray,k),1)
            call dcopy(nrayelt(iray,k),salphal(1:nrayelt(iray,k), &
                 iray,krfn(k)),1, &
                 salphal(1:nrayelt(iray,k),iray,k),1)

            do is=1,nrayelt(iray,k)
               cwexde(is,iray,k)=cwexde(is,iray,krfn(k))
               cweyde(is,iray,k)=cweyde(is,iray,krfn(k))
               cwezde(is,iray,k)=cwezde(is,iray,krfn(k))
            enddo

 30      continue
      enddo


!     Temporary printout of some input ray data
!     Dimensioning is nrayelt(nrayn,mrfn),
!                     ws(nrayelts,nrayn,mrfn),etc.
!$$$      write(*,*)'****************************************************'
!$$$      do krf=1,mrfn
!$$$      write(*,*)'urfread:  mode krf = ',krf
!$$$      write(*,*)'nray(krf),nharms(krfn(krf)),freqcy(krf),nharm(krf)',
!$$$     1           nray(krf),nharms(krfn(krf)),freqcy(krf),nharm(krf)
!$$$      write(*,*)'urfread:nrayelt(*,krf)',(nrayelt(i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:ws(1,*,krf)',(ws(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:spsi(1,*,krf)',(spsi(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wr(1,*,krf)',(wr(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wz(1,*,krf)',(wz(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wphi(1,*,krf)',(wphi(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wnpar(1,*,krf)',(wnpar(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wnper(1,*,krf)',(wnper(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:delpwr(1,*,krf)',
!$$$     1     (delpwr(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:sdpwr(1,*,krf)',(sdpwr(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wdnpar(1,*,krf)',
!$$$     1     (wdnpar(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:cwexde(1,*,krf)',
!$$$     1     (cwexde(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:cweyde(1,*,krf)',
!$$$     1     (cweyde(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:cwezde(1,*,krf)',
!$$$     1     (cwezde(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:fluxn(1,*,krf)',(fluxn(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:sbtot(1,*,krf)',(sbtot(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:sene(1,*,krf)',(sene(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:salphac(1,*,krf)',
!$$$     1     (salphac(1,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:salphal(1,*,krf)',
!$$$     1     (salphal(1,i,krf),i=1,nray(krf))
!$$$
!$$$      write(*,*)'urfread:ws(2,*,krf)',(ws(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:spsi(2,*,krf)',(spsi(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wr(2,*,krf)',(wr(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wz(2,*,krf)',(wz(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wphi(2,*,krf)',(wphi(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wnpar(2,*,krf)',(wnpar(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wnper(2,*,krf)',(wnper(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:delpwr(2,*,krf)',
!$$$     1     (delpwr(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:sdpwr(2,*,krf)',(sdpwr(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:wdnpar(2,*,krf)',
!$$$     1     (wdnpar(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:cwexde(2,*,krf)',
!$$$     1     (cwexde(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:cweyde(2,*,krf)',
!$$$     1     (cweyde(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:cwezde(2,*,krf)',
!$$$     1     (cwezde(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:fluxn(2,*,krf)',(fluxn(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:sbtot(2,*,krf)',(sbtot(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:sene(2,*,krf)',(sene(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:salphac(2,*,krf)',
!$$$     1     (salphac(2,i,krf),i=1,nray(krf))
!$$$      write(*,*)'urfread:salphal(2,*,krf)',
!$$$     1     (salphal(2,i,krf),i=1,nray(krf))
!$$$      enddo
!$$$      write(*,*)'****************************************************'

      return
      end subroutine urfread


end module urfread_mod
