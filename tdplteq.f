c
c
      subroutine tdplteq(krf)
      use param_mod
      use cqcomm_mod
      use frplteq_mod, only : micfrplt, textt
      implicit integer (i-n), real*8 (a-h,o-z)
CMPIINSERT_INCLUDE

      character*8 pltrays
      character*8 plteq

      REAL rbot,zbot,rtop,ztop
      REAL RTAB1(LFIELDA),RTAB2(LFIELDA)
      REAL PGER1,PGERNNR,PGZERO,PGZMIN,PGZMAX,PGEZNNZ
      REAL wk_ray_r(nrayelts), wk_ray_z(nrayelts)

c..................................................................
c     This routine plots out the contours (flux surfaces) on
c     which the calculations are done.
c..................................................................

CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only
 
      pltrays="enabled" ! for plotting rays over cross-sectional view.
                        ! Can be moved to namelist later.

      if (noplots.eq."enabled1") return
      plteq="enabled"
      if (plteq.eq."disabled") return

c     Set up textt to contain number 1 to 200:
      call micfrplt
c     Would be nice to label the flux surfaces, but
c     it isn't done at the moment....

      if (eqsym.ne."none") then ! Up-dn symm
         ztop=.75*ez(nnz)/(er(nnr)-er(1))+.05
      else ! Not up-dn symm
         ztop=0.95
      endif
      
      ! YuP[03-2016][07-2017] Added plotting rays in cross-sectional view
      if (urfmod.ne."disabled" .and. pltrays.eq.'enabled') then 
         ! over-write setting for page size - plot whole cross-section.
         ! For up-dn symmetrical case only half of surfaces are plotted
         ! but rays could be in the other hemisphere.
         ztop=.95
      endif

      
      ! YuP: min and max value of Z-coord over all flux surfaces:
      solz_min=MINVAL(solz) 
      solz_max=MAXVAL(solz)
      PGZMIN=real(solz_min) 
      PGZMAX=real(solz_max) 
      
CPGPLT      CALL PGPAGE
CPGPLT      CALL PGSVP(.15,.85,.15,ztop)

      PGER1=er(1)
      PGERNNR=er(nnr)
      PGZERO=0.
      PGEZNNZ=ez(nnz)
      if( (eqsym.ne."none") .and.
     +    (urfmod.eq."disabled" .or. pltrays.eq.'disabled') ) then 
         !plot half only
c-YuP         CALL PGSWIN(PGER1,PGERNNR,PGZERO,PGEZNNZ)
c-YuP         CALL PGWNAD(PGER1,PGERNNR,PGZERO,PGEZNNZ)
CPGPLT         CALL PGSWIN(PGER1,PGERNNR,PGZMIN,PGZMAX)
CPGPLT         CALL PGWNAD(PGER1,PGERNNR,PGZMIN,PGZMAX)
      else !eqsym=none; and/or urfmod='enabled',
         ! plot upper and lower halves:
CPGPLT         CALL PGSWIN(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
CPGPLT         CALL PGWNAD(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
      endif
     
CPGPLT      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
CPGPLT      if ( (urfmod.ne."disabled") .and. (pltrays.eq.'enabled') 
CPGPLT     +                            .and. (krf.gt.0)            ) then
CPGPLT         CALL PGLAB('Major radius (cms)','Vert height (cms)',
CPGPLT     +        'Fokker-Planck Flux Surfaces + Rays')
CPGPLT      else ! urfmod='disabled', or krf=0
CPGPLT         CALL PGLAB('Major radius (cms)','Vert height (cms)',
CPGPLT     +        'Fokker-Planck Flux Surfaces')
CPGPLT      endif

      IF (LRZMAX.GT.200) STOP 'TDPLTEQ: CHECK DIM OF TEXTT'

      do 10 l=1,lrzmax
         IF (LORBIT(L).GT.LFIELDA) STOP'TDPLTEQ: CHECK DIM OF RTAB1/2'
        do 20 j=1,lorbit(l)
           RTAB1(j)=solr(lorbit(l)+1-j,l)
           RTAB2(j)=solz(lorbit(l)+1-j,l)
 20     continue
        text(1)=textt(l)
CPGPLT        CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
        ! YuP[03-2016] Added plotting rays in cross-sectional view
        if ( (urfmod.ne."disabled") .and. (pltrays.eq.'enabled')
     +       .and. (eqsym.ne.'none')) then 
         ! Add surfaces in whole cross-section.
         ! For up-dn symmetrical case only half of surfaces are plotted
         ! but rays could be in the other hemisphere, 
         ! so plot the other half:
CPGPLT         CALL PGLINE(LORBIT(L),RTAB1,-RTAB2)
        endif
 10   continue

      if((eqsym.eq.'none').or. 
     + ((urfmod.ne."disabled") .and. (pltrays.eq.'enabled')) ) then
     
        if(ncontr.gt.1) then
          ! YuP[2015/05/03] Add LCFS, if available
          ncontr_= min(ncontr,LFIELDA)
          do ilim=1,ncontr_
             RTAB1(ilim)=rcontr(ilim)
             RTAB2(ilim)=zcontr(ilim)
          enddo
CPGPLT          CALL PGLINE(ncontr_,RTAB1,RTAB2)
        endif
      
        if(nlimiter.gt.1) then
          ! YuP[2019-02-22] Add "last surface" (plasma border), if available
          nline= min(nlimiter,LFIELDA)
          do ilim=1,nline
             RTAB1(ilim)=rlimiter(ilim)
             RTAB2(ilim)=zlimiter(ilim)
          enddo
CPGPLT          CALL PGSLW(lnwidth*2) ! bold
CPGPLT          CALL PGLINE(nline,RTAB1,RTAB2)
          if(machine.eq."mirror") then
CPGPLT          CALL PGLINE(nline,-RTAB1,RTAB2) !mirror area to the left of Z-axis
          endif
CPGPLT          CALL PGSLW(lnwidth) ! restore
        endif
        
      endif !eqsym.eq.'none' or (urfmod.ne."disabled")&(pltrays.eq.'enabled')
      
      
      
c..................................................................
c YuP[03-2016] Added plotting rays in cross-sectional view
      if ((urfmod.ne."disabled") .and. (pltrays.eq.'enabled') .and. 
     +    (krf.gt.0)  ) then
        !do krf=1,mrfn ! YuP[2019-02-22] separate page for each krf mode now
        do iray=1,nray(krf)  !Loop over rays
           nrayelt00=nrayelt(iray,krf)
c           write(*,*)'tdplteq: krf, iray, nrayelt00=',
c     +                         krf, iray, nrayelt00
c           write(*,'(a,3i6)')
c     +      'tdplteq: iray,lloc(nrayelt00),llray(nrayelt00)=',
c     +      iray,lloc(nrayelt00,iray,krf),llray(nrayelt00,iray,krf) !local to ray element point.
           wk_ray_r=0.0 ! reset
           wk_ray_z=0.0 ! reset
           do is=1,nrayelt00
             wk_ray_r(is)=wr(is,iray,krf)
             wk_ray_z(is)=wz(is,iray,krf)
           enddo
c           write(*,*)'tdplteq: krf, minval(wk_ray_r),maxval(wk_ray_r)=',
c     +                         krf, minval(wk_ray_r),maxval(wk_ray_r)
CPGPLT           CALL PGLINE(nrayelt00,wk_ray_r,wk_ray_z)
        enddo  ! iray 
        !Add some info on rays:
        write(t_,'(a,i3, a,2i3, a,e12.5)') 'krf=',krf, 
     +                        '  nharm,nharms=', nharm(krf),nharms(krf),
     +                        '    f[Hz]=', freqcy(krf)
CPGPLT        CALL PGMTXT('T',0.9,0.0,0.,t_) ! 0.9=just outside of viewport
        !enddo
      endif  
c..................................................................

      return
      end
