
c
      subroutine frplteq(xpts,ypts,zpts,rpts,vx,vy,vz,ipts,curdep,
     1     nfrplt,frplt)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE


      character*8 textt,frplt
      common /plttextt/ textt(200)

      dimension xpts(*),ypts(*),zpts(*),rpts(*),vx(*),vy(*),vz(*)

      REAL*4 RBOT,RTOP,ZBOT,ZTOP
      REAL*4 RTAB1(LFIELDA),RTAB2(LFIELDA)
      REAL*4 RPG1,RPG2, xyplotmax

      data nconskp /2/

c..................................................................
c     frplt="enabled":
c     This routine plots out the contours (flux surfaces) on
c     which the calculations are done and it plots out some of the
c     birth points of the neutral beam source.
c     frplt="plotwrit":
c     Also outputs points to ascii file freya_points.txt.
c     frplt="write_op":
c     No plotting, just output to ascii file freya_points.txt.
c..................................................................


CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return

      if (frplt.eq."disabled") return
      
      if (frplt.eq."write_op") goto 200 !Skip plotting but save data into freya_points.txt
      
      call micfrplt
      
c---- PLOTS in (R,Z) ---------------------------------------------------   
      rmincon1=rmincon
      if(machine.eq."mirror") rmincon1=-rmaxcon

      if (zmaxcon.gt..5*(rmaxcon-rmincon1)) then
        delr=(rmaxcon-rmincon1)*.9/(2.*zmaxcon)
        zbot=.01 !.05
        ztop=.99 !.95
        rbot=.01+(.9-delr)/2.
        rtop=.99-(.9-delr)/2.
      else
        delz=(2.*zmaxcon)*.9/(rmaxcon-rmincon1)
        rbot=.01
        rtop=.99
        zbot=.01+(.9-delz)/2.
        ztop=.99-(.9-delz)/2.
      endif

      CALL PGPAGE
      CALL PGSVP(rbot,rtop,zbot,ztop)
      RBOT=rmincon
      RTOP=rmaxcon*1.2 ! give 20% more outside of last surface
      if(machine.eq."mirror") then
         RBOT=-RTOP ! to plot left and right sides of flux surfaces
      endif
      ZBOT=zmincon
      ZTOP=zmaxcon
      CALL PGSWIN(rbot,rtop,zbot,ztop)
      CALL PGWNAD(rbot,rtop,zbot,ztop)  ! limits
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
      if(machine.eq."mirror") then      
      CALL PGLAB('X (cms)','Z (cms)', 'NBI Deposition')
      else
      CALL PGLAB('Major radius (cms)','Vert height (cms)',
     +           'NBI Deposition')
      endif
      xyplotmax=0. ! to set limits in (X,Y) plots
      do 10 l=1,lrz,nconskp
        l1=l
        if (l1.gt.lrz-nconskp) l1=lrz
        call tdnflxs(l1)
        call dcopy(lfield,solr(1,lr_),1,tlorb1,1)
        call dcopy(lfield,solz(1,lr_),1,tlorb2,1)

        do 20 j=1,lorbit(lr_)
          solr(j,lr_)=tlorb1(lorbit(lr_)+1-j)
          solz(j,lr_)=tlorb2(lorbit(lr_)+1-j)
 20     continue
        do j=1,lorbit(lr_)
           RTAB1(j)=solr(j,lr_)
           RTAB2(j)=solz(j,lr_)
        enddo
        CALL PGLINE(LORBIT(LR_),RTAB1,RTAB2)
        if(machine.eq."mirror") then
        CALL PGLINE(LORBIT(LR_),-RTAB1,RTAB2) !mirror area to the left of Z-axis
        endif

c       if eqsym.ne."none", still need to plot lower flux surface
        if (eqsym.ne."none") then
           do 30 j=1,lorbit(lr_)
              solz(j,lr_)=-solz(j,lr_)
 30        continue
           DO J=1,LORBIT(LR_)
              RTAB2(J)=SOLZ(J,LR_)
           ENDDO
           CALL PGLINE(LORBIT(LR_),RTAB1,RTAB2)
           if(machine.eq."mirror") then
           CALL PGLINE(LORBIT(LR_),-RTAB1,RTAB2) !mirror area to the left of Z-axis
           endif
        endif

        call dcopy(lfield,tlorb1,1,solr(1,lr_),1)
        call dcopy(lfield,tlorb2,1,solz(1,lr_),1)
 10   continue ! l=1,lrz,nconskp
 
      if(ncontr.gt.1) then
        ! YuP[2015/05/03] Add LCFS, if available
        ncontr_= min(ncontr,LFIELDA)
        r_surf= MAXVAL(rcontr)
        xyplotmax= max(xyplotmax,r_surf)
        do ilim=1,ncontr_
           RTAB1(ilim)=rcontr(ilim)
           RTAB2(ilim)=zcontr(ilim)
        enddo
        CALL PGSLS(2) ! 2-> dashed
        CALL PGLINE(ncontr_,RTAB1,RTAB2)
        if(machine.eq."mirror") then
        CALL PGLINE(ncontr_,-RTAB1,RTAB2) !mirror area to the left of Z-axis
        endif
        CALL PGSLS(1) ! 1-> restore solid line
      endif
      if(nlimiter.gt.1) then
        ! YuP[2016] Add "last surface" (plasma border), if available
        nline= min(nlimiter,LFIELDA)
        r_surf= MAXVAL(rlimiter)
        xyplotmax= max(xyplotmax,r_surf)
        do ilim=1,nline
           RTAB1(ilim)=rlimiter(ilim)
           RTAB2(ilim)=zlimiter(ilim)
        enddo
        CALL PGSLW(lnwidth*2) ! bold
        CALL PGLINE(nline,RTAB1,RTAB2)
        if(machine.eq."mirror") then
        CALL PGLINE(nline,-RTAB1,RTAB2) !mirror area to the left of Z-axis
        endif
        CALL PGSLW(lnwidth) ! restore
      endif

c..................................................................
c     Plot nfrplt birth points. These will all be projected onto
c     one poloidal cross-section.
c..................................................................

      iskip=1+ipts/nfrplt
      write(*,*)'frplteq: ipts,nfrplt,iskip',ipts,nfrplt,iskip
      if (iskip .eq. 0) then
         write(*,*) 'frplteq: iskip=0, beam missing plasma? '
      else
         do i=1,ipts,iskip
            RPG1=RPTS(I)
            if(machine.eq."mirror") then
              RPG1=XPTS(I)
              !area to the left of Z-axis is included, 
              ! so the horizontal axis is X
            endif
            RPG2=ZPTS(I)
            CALL PGPT1(RPG1,RPG2,17)
            xyplotmax= max(xyplotmax,RPG1) ! limits for plots in (X,Y)
         enddo
      endif



c---- PLOTS in (X,Y) (top view) -------------------------------------------   
      rbot=.1
      rtop=.9
      CALL PGPAGE
      CALL PGSVP(rbot,rtop,rbot,rtop)
      RBOT=-rmaxcon
      RTOP= rmaxcon
      CALL PGSWIN(rbot,rtop,rbot,rtop)
      CALL PGWNAD(-xyplotmax,xyplotmax,-xyplotmax,xyplotmax) ! limits 
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
      CALL PGLAB('X (cms)','Y (cms)',
     +           'NBI Deposition')
      ! Plot circles for the largest and smallest FP surfaces.
      nline=LFIELDA ! could be other number, but RTAB1 has LFIELDA size
      r_surf=rpcon(lrz)  ! R radius of largest FP surf, outboard
      do iline=1,nline
         tora= (iline-1)*twopi/(nline-1)
         RTAB1(iline)= r_surf*cos(tora)
         RTAB2(iline)= r_surf*sin(tora)
      enddo
      CALL PGLINE(nline,RTAB1,RTAB2)
      r_surf=rmcon(lrz)  ! R radius of largest FP surf, inboard
      do iline=1,nline
         tora= (iline-1)*twopi/(nline-1)
         RTAB1(iline)= r_surf*cos(tora)
         RTAB2(iline)= r_surf*sin(tora)
      enddo
      CALL PGLINE(nline,RTAB1,RTAB2)
      r_surf=rpcon(1)  ! R radius of smallest FP surf, outboard 
      do iline=1,nline
         tora= (iline-1)*twopi/(nline-1)
         RTAB1(iline)= r_surf*cos(tora)
         RTAB2(iline)= r_surf*sin(tora)
      enddo
      CALL PGLINE(nline,RTAB1,RTAB2)
      r_surf=rmcon(1)  ! R radius of smallest FP surf, inboard 
      do iline=1,nline
         tora= (iline-1)*twopi/(nline-1)
         RTAB1(iline)= r_surf*cos(tora)
         RTAB2(iline)= r_surf*sin(tora)
      enddo
      CALL PGLINE(nline,RTAB1,RTAB2)
 
      if(ncontr.gt.1) then
        ! YuP[2016] Add "last surface" (plasma border), if available
        nline=LFIELDA ! could be other number, but RTAB1 has LFIELDA size
        r_surf= MAXVAL(rcontr)
        do iline=1,nline
          tora= (iline-1)*twopi/(nline-1)
          RTAB1(iline)= r_surf*cos(tora)
          RTAB2(iline)= r_surf*sin(tora)
        enddo
        CALL PGSLS(2) ! 2-> dashed
        CALL PGLINE(nline,RTAB1,RTAB2)
        CALL PGSLS(1) ! 1 - solid
      endif
      if(nlimiter.gt.1) then
        ! YuP[2016] Add "last surface" (plasma border), if available
        nline=LFIELDA ! could be other number, but RTAB1 has LFIELDA size
        r_surf= MAXVAL(rlimiter)
        do iline=1,nline
          tora= (iline-1)*twopi/(nline-1)
          RTAB1(iline)= r_surf*cos(tora)
          RTAB2(iline)= r_surf*sin(tora)
        enddo
        CALL PGSLW(lnwidth*2) ! bold
        CALL PGLINE(nline,RTAB1,RTAB2)
        CALL PGSLW(lnwidth) ! restore
      endif

c..................................................................
c     Plot nfrplt birth points. These will all be projected onto
c     one poloidal cross-section.
c..................................................................

      iskip=1+ipts/nfrplt
      if (iskip .eq. 0) then
         write(*,*) 'frplteq: iskip=0, beam missing plasma? '
      else
         do i=1,ipts,iskip
            RPG1=XPTS(I)
            RPG2=YPTS(I)
            CALL PGPT1(RPG1,RPG2,17)
         enddo
      endif



c--------------------------------------------------------------------      

 200  continue ! to skip plots
 
      if (frplt.eq."plotwrit" .or. frplt.eq."write_op") then

         open(unit=19,file="freya_points.txt",status="replace")
         write(19,1000) 
     1        'Freya birth points: pnt number,x,y,Z,R,vx,vy,vz (cgs)'
         do i=1,ipts
            write(19,1001) i,xpts(i),ypts(i),zpts(i),rpts(i),
     +       vx(i),vy(i),vz(i)
         enddo
 1000    format(a)
 1001    format(i7,1x,7ES12.4E2)
         close(19)

      endif  ! On frplt
      
      
c      STOP  ! TEMP
      return
      end
