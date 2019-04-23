module pltstrml_mod

  !---BEGIN USE

  use aminmx_mod, only : aminmx
  use bcast_mod, only : bcast
  use coefmidt_mod, only : coefmidt
  use coefstup_mod, only : coefstup
  use pack21_mod, only : pack21
  use pltmain_mod, only : gxglfr
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine pltstrml
      use param_mod
      use comm_mod
      use advnce_mod
      use pltdf_mod, only : cont, tempcntr, nconta
      use pltdf_mod, only : wx, wy, IIY, JXQ
      use pltmain_mod, only : gxglfr
      use r8subs_mod, only : luf, dcopy
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real*8 (a-h,o-z)
!

!     This routine plots streamlines of the steady state phase
!     flow. Since the steady state solution does not represent
!     the divergence of the gradient of a function throughout the
!     domain but is sliced into thirds by the pass/trapped
!     boundary, the problem is done three times, once in the
!     passing region, then the co-passing region and then the
!     trapped region. The contours of the resulting functions
!     will represent the stream lines of the flow at steady
!     state only in each of the regions. The contours lines
!     should represent the tangent to the vector field plotted
!     in routine pltvec when all physical processes are
!     included and when the problem is at steady state.
!
!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!     YuP[2018-01-04] Adjusted, to make plots of streamlines.
!
!MPIINSERT_INCLUDE

      integer pltcase
      character*64 tt_
      character*64 tx_,ty_

      REAL RILIN !-> For PGPLOT (text output positioning)

!     PASSING ARRAYS TO PGFUNC1, FOR PGPLOT PGCONX:
      REAL xpt,ypt
      REAL RCONT,RXMAXQ,RTEMP1,RXPTS,RYPTS
      DIMENSION RCONT(NCONTA),RTEMP1(iy,jx),RXPTS(2),RYPTS(2)
!     wx IS V-NORM ARRAY, wy IS THETA ARRAY.  TYPE REAL.
      real*4 RTAB1(iy),RTAB2(iy) ! local
      EXTERNAL PGFUNC1

!MPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return

      !mcont=ncont ! ncont is set in cqlinput (default is 25)
      mcont=60 !YuP: looks like, from setup below,
               !     it should be at least 32

      if(ASSOCIATED(wx)) then
        ! wx and wy are already allocated => do nothing
      else ! Not allocated yet
        allocate(wx(jx))
        allocate(wy(iy))
      endif

      if (mcont.gt.nconta) stop 'in pltcont'

!     Streamlines are plotted in x (=u/vnorm)-space.
!     pltlim and pltlimm are used to limit region of plot to x.lt.1.
!      xmaxq=pltlimm !default is pltlim=disabled; then xmaxq is set here.
!      if (pltlim.ne."disabled") then
!         if (pltlim.eq.'x') then
!            pltlimmm=pltlimm
!         elseif (pltlim.eq.'u/c') then
!            pltlimmm=pltlimm*cnorm
!         elseif (pltlim.eq.'energy') then
!            pltlimmm=sqrt((1.+pltlimm/restmkev)**2-1.)*cnorm
!         endif
!	 xmaxq=pltlimmm
!      endif

      write(t_,5000)
 5000 format("Stream Function of Steady State Phase Flow")
      tt_=trim(t_) ! for the title, above plot



      do 500 k=1,ngen !----------------------------------------------

!-----YuP[2018-01-08] revised to match cqlinput_help:
!**    pltlim= "disabled",  plots versus x (u/vnorm) from
!**                    x=0. to 1. (default:pltlim="disabled",pltlimm=1.)
!**            "x",    plot 1d and 2d plots versus x
!**                    from 0. to pltlimm.
!**            "u/c",  plot 1d and 2d plots versus u/c
!**                    from 0. to pltlimm.
!**            "energy", plot 1d plots verus energy (kev)
!**                    from 0. to pltlimm (kev).
!yup                   BUT, for 2d plots, use u/c units, not keV
         if (pltlim.eq."disabled") then ! whole range in x(j)
            jxq=jx
            xmaxq=x(jxq)
            iyjxq=iy*jxq
            tx_='v_parallel (u/vnorm)'
            ty_='v_perp (u/vnorm)'
         elseif (pltlim.eq.'u/c' .or. pltlim.eq.'energy') then
            if (pltlim.eq.'u/c') then
               pltlimmm=pltlimm
            else ! pltlim.eq.'energy'
               pltlimmm=sqrt((1.+pltlimm/restmkev)**2-1.)
            endif
            jxq=min(luf(pltlimmm,uoc,jx),jx)
            xmaxq=uoc(jxq)
            iyjxq=iy*jxq
            tx_='v_parallel (u/c)'
            ty_='v_perp (u/c)'
         elseif (tandem.eq."enabled" .and. fmass(k).gt.1.e-27) then
            jxq=jlwr
            xmaxq=xlwr
            iyjxq=iy*jlwr
            tx_='v_parallel (u/vnorm)'
            ty_='v_perp (u/vnorm)'
         else ! 'x'
            pltlimmm=pltlimm
            jxq=min(luf(pltlimmm,x,jx),jx)
            xmaxq=x(jxq)
            iyjxq=iy*jxq
            tx_='v_parallel (u/vnorm)'
            ty_='v_perp (u/vnorm)'
         endif

         if (pltlim.eq.'u/c' .or. pltlim.eq.'energy') then
            do j=1,jxq
            tam1(j)=x(j)/cnorm
            enddo
         else
            do j=1,jxq
            tam1(j)=x(j)
            enddo
         endif



         call coefstup(k)
!
!     The coefficients of the equation are currently defined on the
!     same mesh as the distribution function f. The fluxes are best
!     defined (from the point of view of differencing and particle
!     conservation) on mid meshpoints. We undertake here to
!     interpolate the coefficients as needed to either (i,j+1/2)
!     (velocity flux) or to (i+1/2,j) (theta flux).
!     Finally to enforce boundary conditions (zero flux in general
!     except at the pass/trapped boundary) certain coefficients
!     are zeroed out or suitably averaged at specific mesh points.
!     The numbers 1,2,3 appearing in the calls below signify
!     which coefficient is being treated.
!
!     the theta flux..
!
        call coefmidt(dd,1)
        call coefmidt(de,2)
        call coefmidt(df,3)
        call bcast(temp5(0,0),zero,iyjx2)
!
!     In the case implct .eq. "disabled" copy the former values
!     of the distribution function into temporary arrays.
!
        if (implct .eq. "disabled") then
          call dcopy(iyjx2,fxsp(0:iyjx2-1,0,k,l_),1, &
                temp1(0:iyjx2-1,0),1)
          call dcopy(iyjx2,f(0:iyjx2-1,0,k,l_),1,temp2(0:iyjx2-1,0),1)
        endif
!
!     Now proceed with the integration over H
!
        if (implct .eq. "enabled") then
!
!     initialization..
!
          do 15 i=1,iy
 15       temp4(i,1)=(hfi(i,1)+hfi(i-1,1))*.5*dxp5(1)
!
!     Now complete the integration for all x(j)
!
          do 20 j=2,jx-1
            do 31 i=1,iy
              temp4(i,j)=(hfi(i,j)+hfi(i-1,j))*.5*dxp5(j) &
                +temp4(i,j-1)
 31         continue
 20       continue
        else
!
!     initialization..
!
          do 19 i=1,iy
 19       temp4(i,1)=(hfu(i,1)+hfu(i-1,1))*.5*dxp5(1)
!
!     Now complete the integration for all x(j)
!
          do 21 j=2,jx-1
            do 32 i=1,iy
              temp4(i,j)=(hfu(i,j)+hfu(i-1,j))*.5*dxp5(j) &
                +temp4(i,j-1)
 32         continue
 21       continue
        endif
!
!     Patch in values at the pass/trapped boundary to keep
!     the contour plotter happy..
!
        do 40 j=1,jx
          temp4(itl,j)=temp4(itl-1,j)
 40     temp4(itu,j)=temp4(itu+1,j)
!
!     j=jx
!
        do 250 i=1,iy
          temp4(i,jx)=temp4(i,jx-1)
 250    continue
!
!     patch in a value at x=0
!
        do 50 i=1,iy
          temp4(i,1)=temp4(1,1)
 50     continue
        do 210 i=iyh+1,itu-1
          do 211 j=1,jx
            temp4(i,j)=temp4(iy+1-i,j)
 211      continue
 210    continue
!
!     This completes the definition of the function.
!     Determine maximum an minimum values in each of the
!     three regions.
!
        cn1=temp4(1,1)
        cx1=temp4(1,1)
        cx2=cx1
        cx3=cx1
        cn2=cn1
        cn3=cn1

        do 70 j=2,jx
          call aminmx(temp4(1:itl,j),1,itl,1,swwmin,swwmax,kmin,kmax)
!990131          cn1=amin1(cn1,swwmin)
!990131          cx1=amax1(cx1,swwmax)
          cn1=min(cn1,swwmin)
          cx1=max(cx1,swwmax)
!-sww do 60 i=1,itl
!-sww if (cn1 .gt. temp4(i,j)) cn1=temp4(i,j)
!-sww if (cx1 .lt. temp4(i,j)) cx1=temp4(i,j)
!-sww60continue
          do 61 i=itl+1,iyh
            if (cn2 .gt. temp4(i,j)) cn2=temp4(i,j)
            if (cx2 .lt. temp4(i,j)) cx2=temp4(i,j)
 61       continue
          call aminmx(temp4(itu+1:iy-itu,j),1,iy-itu, &
               1,swwmin,swwmax,kmin,kmax)
!990131          cn3=amin1(cn3,swwmin)
!990131          cx3=amax1(cx3,swwmax)
          cn3=min(cn3,swwmin)
          cx3=max(cx3,swwmax)
!-sww do 62 i=itu+1,iy
!-sww if (cn3 .gt. temp4(i,j)) cn3=temp4(i,j)
!-sww if (cx3 .lt. temp4(i,j)) cx3=temp4(i,j)
!-sww62continue
 70     continue
!
!     Scale the distribution fn. in each region
!     so that it's maximum is 1.
!
        do 80 j=1,jx
          do 81 i=1,itl
 81       temp4(i,j)=1.-(temp4(i,j)-cn1)/(cx1-cn1)
 80     continue

        do 82 j=1,jx
          do 83 i=itl+1,itu-1
 83       temp4(i,j)=1.-(temp4(i,j)-cn2)/(cx2-cn2)
          do 84 i=itu,iy
 84       temp4(i,j)=1.-(temp4(i,j)-cn3)/(cx3-cn3)
 82     continue
!
!     determine x-parallel and x-perp mesh
!
!        do 110 j=1,jx
!          do 111 i=1,iy
!            cf(i,j)=x(j)*sinn(i,l_) !-> wx now
!            cd(i,j)=x(j)*coss(i,l_) !-> wy now
!            ! NOT NEEDED?
! 111      continue
! 110    continue

        !call pack21(cf,1,iy,1,jx,tem3,iy,jx) ! NOT NEEDED?
        !call pack21(cd,1,iy,1,jx,tem6,iy,jx) ! NOT NEEDED?
!BobH990608:Is there a problem here with getting desired data into tem5?
!           I.E., is temp4(0,*) temp4(*,0) wanted?
        call pack21(temp4,0,iyp1,0,jxp1,tem5,iy,jx)
!
!
!     determine the contour step size

!990131        dmaxr=alog(.15)
!990131        dminr=alog(constr)
        p15=.15
        dmaxr=log(p15)
        dminr=log(constr)
        dcontr=(dmaxr-dminr)/30.
        cont(1)=dmaxr
        do 800 m=2,30
 800    cont(m)=cont(m-1)-dcontr
        do 801 m=1,30
 801    cont(m)=1.-exp(cont(m))
!990131        smaxr=alog(.845)
!990131        sminr=alog(constr)
        p845=.845
        smaxr=log(p845)
        sminr=log(constr)
        dcontr=(smaxr-sminr)/float(mcont-30)
        cont(31)=smaxr
        do 200 m=32,mcont
          cont(m)=cont(m-1)-dcontr
 200    continue
        do 201 ku=31,mcont
          cont(ku)=exp(cont(ku))
 201    continue


        if (pi .gt. 3.) go to 603
        mau=mcont/3
        if (mau .gt. jx) mau=jx
        inc=jx/mau
        cont(1)=temp5(itl+2,3)
        do 600 nc=2,mau
          is=3+(nc-1)*inc
          if (is .gt. jx) go to 600
          cont(nc)=temp5(itl+2,is)
 600    continue
        cont(1)=temp5(iyh-1,3)
        do 601 nc=2,mau
          is=3+(nc-1)*inc
          if (is .gt. jx) go to 601
          cont(nc)=temp5(iyh-1,is)
 601    continue
 603    continue

        dmin=0.d0
        dmax=0.d0
        DO J=1,JXQ
         DO I=1,iy
            DMIN=MIN(temp1(I,J),DMIN)
            DMAX=MAX(temp1(I,J),DMAX)
            RTEMP1(I,J)=temp4(I,J) ! the cont() levels of this function
            ! will be plotted.
         ENDDO
        ENDDO
        admin=abs(dmin)

        DO J=1,JXQ
         wx(J)=TAM1(J)
        ENDDO
        DO I=1,iy
         wy(I)=y(I,L_)
        ENDDO
        DO JS=1,mcont
         RCONT(JS)=CONT(JS)
         !write(*,*)'pltcont: lr_,JS,CONT(JS)=',lr_,JS,CONT(JS)
        ENDDO
        IIY=iy
        RXMAXQ=XMAXQ


        call GXGLFR(0) ! new page for each k
        CALL PGSVP(.2,.8,.65,.9)
        IF ( RXMAXQ.eq.0. ) THEN
           RXMAXQ=1.
        ENDIF
        CALL PGSWIN(-RXMAXQ,RXMAXQ,0.,RXMAXQ)
        CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
        CALL PGLAB(tx_,ty_,tt_)

        t0t=sin(thb(l_))/cos(thb(l_))  ! PLOT t-p boundary (ZOW cone)
        if (t0t .lt. 1.) then
          RXPTS(1)=0.
          RYPTS(1)=0.
          RXPTS(2)=XMAXQ
          RYPTS(2)=XMAXQ*T0T
          CALL PGLINE(2,RXPTS,RYPTS)
          RXPTS(2)=-XMAXQ
          CALL PGLINE(2,RXPTS,RYPTS)
        else
          RXPTS(1)=0.
          RYPTS(1)=0.
          RXPTS(2)=XMAXQ/T0T
          RYPTS(2)=XMAXQ
          CALL PGLINE(2,RXPTS,RYPTS)
          RXPTS(2)=-XMAXQ/T0T
          CALL PGLINE(2,RXPTS,RYPTS)
        endif

        !plot v=vnorm line
        if (pltlim.eq.'u/c' .or. pltlim.eq.'energy') then
          do i=1,iy
          RTAB1(i)= coss(i,lr_)/cnorm ! v_par/c  (cnorm is c/vnorm)
          RTAB2(i)= sinn(i,lr_)/cnorm ! v_perp/c
          enddo
        else ! v/vnorm units
          do i=1,iy
          RTAB1(i)= coss(i,lr_) ! v_par/vnorm
          RTAB2(i)= sinn(i,lr_) ! v_perp/vnorm
          enddo
        endif
        CALL PGLINE(iy,RTAB1,RTAB2) ! v=vnorm line (or v=vnorm/c)

        !plot v=vth line, for the k-th gen. species
!..................................................................
!     Note: the do loop below uses vth(),
!     vth is the thermal velocity =sqrt(T/m) (at t=0 defined in ainpla).
!     But, T==temp(k,lr) can be changed in profiles.f,
!     in case of iprote (or iproti) equal to "prbola-t" or "spline-t"
!..................................................................
        if (pltlim.eq.'u/c' .or. pltlim.eq.'energy') then
          do i=1,iy
          RTAB1(i)= (vth(k,lr_)/clight)*coss(i,lr_) ! vth_par/c
          RTAB2(i)= (vth(k,lr_)/clight)*sinn(i,lr_) ! vth_perp/c
          enddo
        else ! v/vnorm units
          do i=1,iy
          RTAB1(i)= (vth(k,lr_)/vnorm)*coss(i,lr_) ! vth_par/vnorm
          RTAB2(i)= (vth(k,lr_)/vnorm)*sinn(i,lr_) ! vth_perp/vnorm
          enddo
        endif
        ! Five different line styles are available:
        ! 1 (full line), 2 (dashed), 3 (dot-dash-dot-dash), 4 (dotted),
        CALL PGSLS(4)
        CALL PGLINE(iy,RTAB1,RTAB2) ! v=vth line
        CALL PGSLS(1) ! 1-> restore solid line
        CALL PGSLW(lnwidth) !lnwidth=3 line width in units of 0.005
        !--------------------------------------------------------
        CALL PGCONX(RTEMP1,iy,jx,1,iy,1,JXQ,RCONT,mcont,PGFUNC1)
        !subr.PGFUNC1(VISBLE,yplt,xplt,zplt) uses /PGLOCAL1/wx,wy,IIY,JXQ
        !--------------------------------------------------------
        !Add some text on the plot:
        CALL PGSCH(1.0) ! set character size; default is 1.
        write(t_,5001) k
 5001   format("Species number k=",i3)
        CALL PGMTXT('B',6.,0.,0.,t_)
        write(t_,150) n,timet
 150    format("time step n=",i5,5x," time=",1pe10.2," secs")
        CALL PGMTXT('B',7.,0.,0.,t_)
        rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
        write(t_,151) rovera(lr_),rr
 151    format( "r/a=",1pe10.3,5x," radial position (R)=",1pe12.4," cm")
        CALL PGMTXT('B',8.,0.,0.,t_)
        write(t_,153) rya(lr_), rpcon(lr_), lr_
 153    format( "rya=",1pe10.3,5x," R=rpcon=",1pe10.3," cm,  Surf#",i4)
        CALL PGMTXT('B',9.,0.,0.,t_)
        ! print contour values under the plot:
        do js=1,mcont
           tempcntr(js)=cont(js)
        enddo
        write(t_,560)
 560    format("Contour values:")
        RILIN=11.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)
        do jcs=1,mcont,4
          write(t_,570) (tempcntr(jc),jc=jcs,min(jcs+3,mcont))
          if ((mcont/4)*4.ne.mcont .and. mcont-jcs.le.2) then
            icend=4 * 16 + 1
            t_(icend:icend)="$"
          endif
          RILIN=RILIN+1.
          CALL PGMTXT('B',RILIN,-.2,0.,t_)
        enddo

        CALL PGSLS(1) ! restore: solid line
        CALL PGSLW(lnwidth) ! restore linewidth
        CALL PGSCH(1.0) ! recover default 1.0 fontsize

 500  continue ! k species ------------------------------------------

 570  format(4(1pe16.6))

      return
      end
end module pltstrml_mod
