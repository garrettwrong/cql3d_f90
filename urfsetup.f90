module urfsetup_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast
  use urfalloc_mod, only : urfalloc
  use urfread__mod, only : urfread_i

  !---END USE

!
!

contains

      subroutine urfsetup
      use param_mod
      use cqlcomm_mod
      use netcdfrf_mod, only : netcdfrf
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine sets up arrays unique to urf module
!..................................................................

!MPIINSERT_INCLUDE

!..................................................................
!     Return if urfmod.ne."enabled"
!..................................................................

      if (urfmod.ne."enabled") return

!...................................................................
!     mrf is the number of RF types  being calculated.
!     nharms is the number of cyclotron harmonics calculated,
!     starting at nharm1.
!     (nharms.gt.1 .and. mrf.gt.1) is now permitted by  the
!     storage scheme [BH060314].
!...................................................................

      mrf=0
      mrfn=0

!     Check rftype method of input, irftype [Set in urfinitl]
!     irftype=0 is older wave type specification
!     irftype=1 gives direct specification through rftype(1:nmodsa).
!     Count number of wave types, mrf.
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


!     mrfn is the number of wave "modes", that is, the sum over
!     wave types of the number of harmonics for each wave type
      mrfn=0
      do k=1,mrf
         if (nharms(k).eq.0) then
            mrfn=mrfn+1
         else
            mrfn=mrfn+nharms(k)
         endif
      enddo

      if (mrfn.gt.nmodsa) print 100, mrfn,nmodsa
 100  format('urfsetup:   mrfn.gt.nmodsa',i5,'>',i3)


      if (lr_.eq.lrindx(lrz) .and. l_.eq.lrors) then

!-YuP 101122: This part is added to read data on number of rays---------
! It repeats part of urfread.
! If netcdf file is read, use netcdf(kopt,krf) with kopt=3 and then =4.
! kopt=3 is to only read the values of nray(krf),nharm(krf),freqcy(krf).
! kopt=4 is to read nrayelt(iray,krf).
! Other data cannot be read because arrays are not allocated yet.
 1    format(2i5,1pe16.9)

!MPIINSERT_IF_RANK_EQ_0
      krf=0  !  Counter for wave types.
!     irftype is determined in urfinitl.
      if (irftype.eq.0) then    !Old type input of wave type
      if (ech.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(3,krf)
         else
            open(unit=23,file='rayech',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               !-YuP: Part of call urfread_(krf,23)
               read(23,1) nray(krf),nharm(krf),freqcy(krf)
               close(unit=23)
            else
               WRITE(*,*)'urfread: check that "rayech" file is present'
               stop
            endif
         endif
      endif ! ech
      if (fw.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(3,krf)
         else
            open(unit=24,file='rayfw',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               !-YuP: Part of call urfread_(krf,24)
               read(24,1) nray(krf),nharm(krf),freqcy(krf)
               close(unit=24)
            else
               WRITE(*,*)'urfread: check that "rayfw" file is present'
               stop
            endif
         endif
      endif ! fw
      if (lh.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(3,krf)
         else
            open(unit=20,file='raylh',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               !-YuP: Part of call urfread_(krf,20)
               read(20,1) nray(krf),nharm(krf),freqcy(krf)
               close(unit=20)
            else
               WRITE(*,*)'urfread: check that "raylh" file is present'
               stop
            endif
         endif
      endif ! lh
      endif   ! on irftype.eq.0

      if (irftype.eq.1) then     !New, more general input of wave types.
         do krf=1,mrf
            call netcdfrf(3,krf) !get the value of nray(krf)
         enddo
      endif  ! on irftype.eq.1

! Now the values of nray(krf),nharm(krf),freqcy(krf) are known.
! Check the largest nray(krf) and set nrayn=max(nray(1:mrf))
! nrayn will be used to allocate arrays, mostly in urfalloc.
      nrayn=1
      do krf=1,mrf
         nrayn= max(nrayn,nray(krf))
        WRITE(*,'(a,2i9,e12.4)')'URFSETUP: nray(krf),nharm(krf),freqcy', &
          nray(krf),nharm(krf),freqcy(krf)
      enddo
      WRITE(*,*)'URFSETUP: nrayn===',nrayn
!MPIINSERT_ENDIF_RANK

!MPIINSERT_BCAST_NRAYN
!MPIINSERT_BARRIER

!---> Allocation of arrays with (nrayn,mrfn) size
!---> is moved here from urfalloc.
!---> They will be needed to read the value of nrayelt(nrayn,mrfn)
      istat_tot=0
      allocate(nrayelt(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(nrayelt,0,SIZE(nrayelt))
      allocate(jslofas(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(jslofas,0,SIZE(jslofas))
      allocate(nurefls(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(nurefls,0,SIZE(nurefls))
      allocate(keiks(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(keiks,0,SIZE(keiks))
      allocate(jpes(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(jpes,0,SIZE(jpes))
      allocate(jpis(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(jpis,0,SIZE(jpis))
      allocate(istarts(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(istarts,0,SIZE(istarts))
      allocate(iprmt5(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(iprmt5,0,SIZE(iprmt5))
      allocate(jhlfs(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(jhlfs,0,SIZE(jhlfs))
      allocate(sxxrt(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call bcast(sxxrt,zero,SIZE(sxxrt))
      allocate(skpsi(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call bcast(skpsi,zero,SIZE(skpsi))
      allocate(skth(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call bcast(skth,zero,SIZE(skth))
      allocate(skphi(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call bcast(skphi,zero,SIZE(skphi))
      allocate(lrayelt(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(lrayelt,0,SIZE(lrayelt))
      allocate(delpwr0(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call bcast(delpwr0,zero,SIZE(delpwr0))
      allocate(nrayelt0(nrayn,mrfn),STAT=istat)
      istat_tot=istat_tot+istat
      if(istat.eq.0) call ibcast(nrayelt0,0,SIZE(nrayelt0))

!     Check that allocations were OK
!    (should not happen because (nrayn,mrfn)-size is quite small)
      if (istat_tot.ne.0) then
         write(*,*)'urfsetup.f:  Problem with allocation'
         write(*,*)'urfsetup.f:  Stopping. istat_tot=', istat_tot
         STOP
      endif


!MPIINSERT_IF_RANK_EQ_0
!---> Read rays data file again to get nrayelt(1:nrayn,1:mrfn)
      krf=0  !  Counter for wave types.
      if (irftype.eq.0) then    !Old type input of wave type
      if (ech.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(4,krf)
         else
            open(unit=23,file='rayech',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               call urfread_i(krf,23)
               close(unit=23)
            else
               WRITE(*,*)'urfread: check that "rayech" file is present'
               stop
            endif
         endif
      endif ! ech
      if (fw.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(4,krf)
         else
            open(unit=24,file='rayfw',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               call urfread_i(krf,24)
               close(unit=24)
            else
               WRITE(*,*)'urfread: check that "rayfw" file is present'
               stop
            endif
         endif
      endif ! fw
      if (lh.eq."enabled") then
         krf=krf+1
         if (rfread.eq."netcdf") then
            call netcdfrf(4,krf)
         else
            open(unit=20,file='raylh',status='old',IOSTAT=io1)
            if (io1.eq.0) then
               call urfread_i(krf,20)
               close(unit=20)
            else
               WRITE(*,*)'urfread: check that "raylh" file is present'
               stop
            endif
         endif
      endif ! lh
      endif   ! on irftype.eq.0

      if (irftype.eq.1) then     !New, more general input of wave types.
         do krf=1,mrf
            call netcdfrf(4,krf) !get the values of nrayelt(:,krf)
         enddo
      endif  ! on irftype.eq.1

! Now the values of nrayelt(iray,krf) are also known.
! Check the largest nrayelt(:,:) and set nrayelts=max(nrayelt(:,:))
! nrayelts will be used to allocate arrays, mostly in urfalloc.
      nrayelts=1
      do krf=1,mrf
      do iray=1,nray(krf)
         nrayelts= max(nrayelts,nrayelt(iray,krf))
         WRITE(*,'(a,3i9)')'URFSETUP: iray,krf,nrayelt(iray,krf)=', &
                                      iray,krf,nrayelt(iray,krf)
      enddo
      enddo
      WRITE(*,*)'URFSETUP: nrayelts===',nrayelts
!MPIINSERT_ENDIF_RANK

!MPIINSERT_BCAST_NRAYELTS
!MPIINSERT_BARRIER

!-YuP 101122-end--------------------------------------------------------




!..................................................................
!     Allocate arrays
!..................................................................
!       ibytes=8/machinea, a parameter.   machinea=2 works with
!       present 32- and 64-bit machines, with 32-bit intergers.
        jjx=((jx-1)/ibytes)*ibytes+ibytes

        call urfalloc
!MPIINSERT_BARRIER

        do 30 j=1,jx
          sx(j)=x(j)/gamma(j)
          xmdx(j)=xmidpt(j)**2
 30     continue

!-YuP: commented-out because it's done in urfalloc
!-YuP        do 40 k=1,mrfn
!-YuP          call ibcast(lrayelt(1,k),0,nrayn)
!-YuP 40     continue

      endif !   if (lr_.eq.lrindx(lrz) .and. l_.eq.lrors) then


!..................................................................
!     "l" refers to poloidal position z(l,lr_)
!     lr_ is the flux surface label.
!..................................................................

      do 20 l=1,lz
        cosmz(iy,l,lr_)=-cosz(iy,l,lr_)
        do 10 i=1,iy-1
          cosmz(i,l,lr_)=-.5*(cosz(i+1,l,lr_)+cosz(i,l,lr_))
 10     continue
 20   continue
      return
      end subroutine urfsetup
      
      
end module urfsetup_mod
