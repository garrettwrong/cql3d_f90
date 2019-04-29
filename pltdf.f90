!
!
module pltdf_mod

  !---BEGIN USE

  use lookup_mod, only : lookup_tdf
  use r8subs_mod, only : dcopy

  !---END USE

  use iso_c_binding, only : c_double
  integer, parameter, public :: nconta=100
  ! these are shared with other routines (outside the contains)
  real(c_double), public :: cont(nconta),tempcntr(nconta)
  ! PGLOCAL1 from pltcont.f
  real(c_double), pointer, public :: wx(:), wy(:)
  integer, public :: IIX, IIY, JXQ

  save

contains

  subroutine pltdf
    use param_mod
    use comm_mod
    use r8subs_mod, only : dcopy
    implicit integer (i-n), real*8 (a-h,o-z)

    !..................................................................
    !     if (pltd.eq."enabled" or pltd.eq."color") then
    !     subroutine pltdf contour plots the distribution function, f
    !     if (pltd.eq."df"  or "df_color") then
    !     subroutine pltdf also plots the difference between f
    !     at the current time and f at the previous time
    !     step.
    !..................................................................


    REAL RILIN

    if (noplots.eq."enabled1") return

    if (pltd.eq."disabled") return


    do 10 k=1,ngen

       ! This part is plotted for any pltd.ne.'disabled'

       call dcopy(iyjx2,f(0:iyjx2-1,0,k,l_),1,temp1(0:iyjx2-1,0),1)
       write(t_,550) k
550    format(1x,"Species ",i2," Distribution Function Contour Plot")
               CALL PGPAGE
       itype=1 ! means: plots are made for distr.func f
       call pltcont(k,2,t_,itype) ! for f()
       write(t_,560)
560    format("Contour values:")
       RILIN=10.
               CALL PGMTXT('B',RILIN,-.2,0.,t_)


       do 11 jcs=1,ncont,4
          write(t_,570) (tempcntr(jc),jc=jcs,min(jcs+3,ncont))
          if ((ncont/4)*4.ne.ncont .and. ncont-jcs.le.2) then
             icend=4 * 16 + 1
             t_(icend:icend)="$"
          endif
          RILIN=RILIN+1.
                    CALL PGMTXT('B',RILIN,-.2,0.,t_)
11     end do


       if (n.eq.0) goto 10

       ! Additionally, plot f(n+1)-f(n)
       if (pltd.eq."df" .or. pltd.eq."df_color")then

          do 20 i=1,iy
             do 21 j=1,jx
                temp1(i,j)=f(i,j,k,l_)-f_(i,j,k,l_)
21           end do
20        end do
          write(t_,530) k,n
530       format(1x, "Contours of df/dt for species",i2,1x,"during timestep",i5)
                  CALL PGPAGE
          itype=2 ! means: plots are made for df
                  call pltcont(k,1,t_,itype) ! for df
          RILIN=10.
          write(t_,560)
                  CALL PGMTXT('B',RILIN,-.2,0.,t_)

          do 12 jcs=1,ncont,4
             write(t_,570) (tempcntr(jc),jc=jcs,min(jcs+3,ncont))
             if ((ncont/4)*4.ne.ncont .and. ncont-jcs.le.2) then
                icend=4 * 16 + 1
                t_(icend:icend)="$"
             endif
             RILIN=RILIN+1.
                       CALL PGMTXT('B',RILIN,-.2,0.,t_)
12        end do

       endif !  pltd.eq."df" .or. pltd.eq."df_color"


10  end do

570 format(4(1pe16.6))

    return

!bug end do
  end subroutine pltdf

      subroutine pltcont(k,pltcase,tt_,itype)
      use param_mod
      use comm_mod
      use r8subs_mod, only : luf
      !YuP[2018-02-07] added input itype, to identify what is plotted.
      ! itype=1 for plots of f(),
      !      =2 for df
      !      =3 for source
      !      =4 for urfb
      !      =5 for vlfb
      !      =6 for vlhb
      !      =7 for rdcb
      !      =8 for loss region (from pltlosc)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
!MPIINSERT_INCLUDE

      integer pltcase
      character*(*) tt_
      character*64 tx_,ty_
      REAL wk_tam(jx)
      real*8 wkd(jx)
!...
!mnt  This routine performs contour plots of distributions
!mnt  given in temp1(0:iyp1,0,jxp1) as a function of v,theta,
!mnt  as specified by x,y.
!mnt    pltcase=1: geometric progression of plot contours,
!mnt    pltcase.ne.1: contours equispaced for max. at temp(k,lr_).
!mnt    k gives species index.
!mnt    tt_ is contour plot heading.
!mnt    Additional annotation of the plot can be added from
!mnt    the calling routine.
!...
!
!     PASSING ARRAYS TO PGFUNC1, FOR PGPLOT PGCONX:
      REAL xpt,ypt
      REAL RCONT,RXMAXQ,RTEMP1,RXPTS,RYPTS
      DIMENSION RCONT(NCONTA),RTEMP1(iy,jx),RXPTS(2),RYPTS(2)
!     wx IS V-NORM ARRAY, wy IS THETA ARRAY.  TYPE REAL.
      real*4 RTAB1(iy),RTAB2(iy) ! local

      !YuP[2018-01-27] Local, for PGPLOT:
      parameter(npar=200, nprp=100)
      real*8 vpar(npar),vprp(nprp) ! rectangular grid for plots
      real*8 fparprp(npar,nprp) !f(i,j) will be interpolated to this grid
      real*4 BRIGHT,CONTRA,FMIN,FMAX,RVMIN,RVMAX, &
             TRPG(6),RTEMP2(npar,nprp)
      REAL*4 RCONTLOG(NCONTA),FLOGMIN,FLOGMAX



!MPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return

      if(ASSOCIATED(wx)) then
        ! wx and wy are already allocated => do nothing
      else ! Not allocated yet
        allocate(wx(jx))
        allocate(wy(iy))
      endif

      if (ncont.gt.nconta) stop 'in pltcont'

      dmin=0.
      dmax=0.

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
         else
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
            tam1(j)=uoc(j) !==x(j)/cnorm
         enddo
      else ! 'disabled', 'x', or tandem
         do j=1,jxq
            tam1(j)=x(j)
         enddo
      endif

!

      DO J=1,JXQ
         DO I=1,iy
            DMIN=MIN(temp1(I,J),DMIN)
            DMAX=MAX(temp1(I,J),DMAX)
            RTEMP1(I,J)=temp1(I,J)
         ENDDO
      ENDDO


      admin=abs(dmin)

      if( (itype.eq.1 .and. pltd.eq.'color')   .or. &
          (itype.eq.1 .and. pltd.eq.'df_color')      .or. &
          (itype.eq.2 .and. pltd.eq.'df_color'   )   .or. &
          (itype.eq.3 .and. pltso.eq.'color')  .or. &
          (itype.eq.3 .and. pltso.eq.'first_cl')  .or. &
          (itype.eq.4 .and. plturfb.eq.'color')   .or. &
          (itype.eq.7 .and. pltrdc.eq.'onecolor')   .or. &
          (itype.eq.7 .and. pltrdc.eq.'allcolor')     ) then
          !Other itype can be added later

      !YuP[2018-01-27] Added: interpolate f(i,j) to fparprp(npar,nprp)
      ! over the rectangular (vpar,vprp) grid.
      ! Set rect. grid, then interpolate f(i,j) to this grid
      f_zero=1.d-100
      vmax= xmaxq  ! either u/c or u/vnorm units
      vmin=-xmaxq  ! either u/c or u/vnorm units
      dvpar=(vmax-vmin)/(npar-1)
      do ipar=1,npar
        vpar(ipar)= vmin+dvpar*(ipar-1) ! either u/c or u/vnorm units
      enddo
      dvprp=(vmax-0.d0)/(nprp-1)
      do iprp=1,nprp
        vprp(iprp)= 0.d0+dvprp*(iprp-1) ! either u/c or u/vnorm units
      enddo
      ! Interpolate f(i,j) to fparprp(npar,nprp)
      ! ADJUST f(i,j) - to eliminate zero values.
      ! This is important if we want to consider LOG10(f).
      ! For now, simply set it to f_zero
      fparprp= f_zero ! ~zero level ! initialize
      dxlast=tam1(jxq)-tam1(jxq-1) ! will be used for (vpar,vprp) points
                           ! outside of tam1(jxq) semi-circle
      do iprp=1,nprp
      do ipar=1,npar
         ! for a given (vpar,vprp) find the four nearest points
         ! in (i,j) grid
         !(y(i)= pitch angle in radians)
         vloc= sqrt(vpar(ipar)**2 + vprp(iprp)**2) !either u/c or u/vnorm units
         ploc= atan2(vprp(iprp),vpar(ipar)) ! pitch angle [rad]
         !If vprp=0, the result is 0 (if vpar >0) or pi (if vpar <0).
         !-1-> For the given vloc, find the nearest j index in tam1(j) grid:
         call lookup_tdf(vloc,tam1,jxq,rweightu,rweightl,jloc)
         ! vloc is between tam1(jloc-1) and tam1(jloc) grid points
         ! rweightu,rweightl are the weight factors for lin. interpolation.
         !-2-> For the given ploc, find the nearest i index
         !     in y0pi(i) pitch grid.
         call lookup_tdf(ploc,y,iy,pweightu,pweightl,iloc)
         !write(*,*)ploc,iloc
         ! ploc is between y(iloc-1) and y(iloc).
         !-3-> Interpolate values of f(i,j) from four points
         !     to (vpar,vprp) point, with ~zero outside of tam1(jxq) border
         if(vloc.gt.tam1(jxq)+dxlast)then
           ! Far outside of tam1(jxq) semi-circle
           fparprp(ipar,iprp)=f_zero ! basically zero (not defined)
         elseif(vloc.gt.tam1(jxq))then
           ! Outside of tam1(jxq) semi-circle,
           ! but close (vloc is less than tam1(jxq)+dxlast).
           ! Make a gradual drop to ~zero level
           fmm= temp1(iloc-1,jxq)
           f0m= temp1(iloc,  jxq)
           fm0= f_zero ! ~zero level
           f00= f_zero
           floc= rweightl*( pweightl*fmm + pweightu*f0m ) &
                +rweightu*( pweightl*fm0 + pweightu*f00 )
           fparprp(ipar,iprp)=floc !! set
         else ! vloc.le.tam1(jxq)i.e. interior to tam1(jxq) semi-circle
           fmm= temp1(iloc-1,jloc-1)
           f0m= temp1(iloc,  jloc-1)
           fm0= temp1(iloc-1,jloc)
           f00= temp1(iloc,  jloc)
           floc= rweightl*( pweightl*fmm + pweightu*f0m ) &
                +rweightu*( pweightl*fm0 + pweightu*f00 )
           fparprp(ipar,iprp)=floc !! set
         endif
         ! adjust - to eliminate neg. values:
         fparprp(ipar,iprp)=max(fparprp(ipar,iprp),f_zero)
         ! Because for PGIMAG, we need LOG10() scale.
      enddo
      enddo

      endif ! on color option

!**bh
!
      if(dmax.le.0)  go to 999
!**bh
!
!     In the case abs(dmin).le.contrmin*dmax, then plot contours
!     occur at levels:
!     cont(j)=dmax * contrmin**(1-(j-.5)/ncont),  j=1,ncont
!     This can be described as a geometric progression of values
!     from (near) dmax down to contrmin*dmax.
!
      if (dmin.ge.0. .or. admin.lt.contrmin*dmax) then
         k2=1
         if(pltcase.eq.1) then
            smin=log(contrmin*dmax)
            if (admin/dmax .gt. contrmin) smin=log(admin)
            smax=log(dmax)
            dcont=(smax-smin)/dfloat(ncont)
            cont(1)=smin+.5*dcont
            do 20 kc=2,ncont
               cont(kc)=cont(kc-1)+dcont
 20         continue
            do 30 kc=1,ncont
               cont(kc)=exp(cont(kc))
 30         continue
         else

!     Modifying above case, which is usual situation for a function
!     for which dmin .ge.0. .or. not very  negative, to contours which
!     will be equispaced for a Maxwellian distribution with temperature
!     equal to that defined for the distribution:
            emax=-temp(k,lr_)*log(contrmin)
            if (emax.gt.enorm) emax=enorm
            gammax=1.+emax/restmkev
            uocmax=sqrt(gammax**2-1.)
            do j=1,ncont
               cont(j)=j/dfloat(ncont)*uocmax
            enddo
            do j=1,ncont
               cont(j)=dmax* &
                    exp(-restmkev*(sqrt(1.+cont(j)**2)-1.)/temp(k,lr_))
            enddo
         endif

      else
        if (dmax .gt. 0.) then
          k2=ncont*.1+1
          ncontp=ncont-k2+1
          ncontm=k2-1
          smaxp=log(dmax)
          sminp=log(contrmin*dmax)
        else
          ncontm=ncont
          ncontp=1
          k2=1
        endif
        sminm=log(-contrmin*dmin)
        if (dmax/dmin.gt.contrmin) sminm=log(-dmax)
        smaxm=log(-dmin)
        dcontp=(smaxp-sminp)/dfloat(ncontp)
        dcontm=(smaxm-sminm)/dfloat(ncontm)
        cont(1)=smaxm-.5*dcontm
        do 40 kc=2,ncontm
          cont(kc)=cont(kc-1)-dcontm
 40     continue
        do 50 kc=1,ncontm
          cont(kc)=-exp(cont(kc))
 50     continue
        if (dmax .gt. 0.) then
          cont(k2)=sminp-.5*dcontp
          do 60 kc=k2+1,ncont
            cont(kc)=cont(kc-1)+dcontp
 60       continue
          do 70 kc=k2,ncont
            cont(kc)=exp(cont(kc))
 70       continue
        endif
      endif
      if (k2.eq.0) k2=1 ! YuP: bug? was if(k.eq.0)
      do 71 js=1,ncont
        tempcntr(js)=cont(js)
 71   continue

      DO J=1,JXQ
         wx(J)=TAM1(J)
      ENDDO
      DO I=1,iy
         wy(I)=y(I,L_)
      ENDDO
      DO JS=1,NCONT
         RCONT(JS)=CONT(JS)
         !write(*,*)'pltcont: t_,lr_,JS,CONT(JS)=',lr_,JS,CONT(JS)
         if(CONT(JS).gt.0.d0)then
           RCONTLOG(JS)=LOG10(CONT(JS)) !Can be used for (as an option):
           !CALL PGCONT(RTEMP2,npar,nprp,1,npar,1,nprp,RCONTLOG,-NCONT,TR)
           ! where RTEMP2(ipar,iprp)=log10(fparprp(ipar,iprp))
           ! The LOG10 scale is needed for PGIMAG()
         else
           RCONTLOG(JS)=-100.
         endif
      ENDDO
      IIY=iy
      RXMAXQ=XMAXQ

      CALL PGSVP(.2,.8,.65,.9)
        IF ( RXMAXQ.eq.0. ) THEN
           RXMAXQ=1.
        ENDIF
      CALL PGSWIN(-RXMAXQ,RXMAXQ,0.,RXMAXQ)
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
      CALL PGLAB(tx_,ty_,tt_)

      if( (itype.eq.1 .and. pltd.eq.'color')   .or. &
          (itype.eq.1 .and. pltd.eq.'df_color')      .or. &
          (itype.eq.2 .and. pltd.eq.'df_color'   )   .or. &
          (itype.eq.3 .and. pltso.eq.'color')  .or. &
          (itype.eq.3 .and. pltso.eq.'first_cl')  .or. &
          (itype.eq.4 .and. plturfb.eq.'color')   .or. &
          (itype.eq.7 .and. pltrdc.eq.'onecolor')   .or. &
          (itype.eq.7 .and. pltrdc.eq.'allcolor')     ) then
          !Other itype can be added later

      !BH,YuP[2018-01-26] Version 2 for color map of distr.func.
      !Set up the color map.
      BRIGHT=0.5 ! 0.8 gives yellow_low -- red_upper (no blue)
      CONTRA=0.8 !0.5 gives light-blue background
                 !1.0 gives dark-blue(almost black)
      CALL PALETT(2, CONTRA, BRIGHT)
      !First arg:
      ! 1- gray scale
      ! 2- rainbow (ok)
      ! 3- heat    (bad: gives black background)
      ! 4- weird IRAF ( really weird: black background and random colors)
      ! 5- AIPS ( not so bad, but 2 probably is the best)
      FMIN=max(f_zero,dmin) ! or could be a fraction of it.
      FMAX=dmax ! or could be a fraction of it.
      FLOGMAX=LOG10(dmax)
      FLOGMIN=LOG10(dmax*contrmin) ! => FLOGMIN=FLOGMAX-LOG10R
      !write(*,*)'FLOGMIN,FLOGMAX=',FLOGMIN,FLOGMAX
      !dmin, dmax are found above, for TEMP1(I,J)=f(I,J) (lin. scale).
      !      but RCONT() levels are log()

      !Set the coordinate transformation matrix TR(1:6):
      !world coordinate = pixel number.
      ! From PGPLOT manual:
      !The transformation matrix TR is used to calculate the world
      !coordinates of the center of the "cell" that represents each
      !array element. The world coordinates of the center of the cell
      !corresponding to array element A(I,J) are given by:
      !          X = TR(1) + TR(2)*I + TR(3)*J
      !          Y = TR(4) + TR(5)*I + TR(6)*J
      !
      ! Based on our X==Vpar= RVMIN  + dvpar*(ipar-1)
      !              Y==Vprp=  0     + dvprp*(iprp-1)
      ! we set:
      RVMIN=vmin ! vmin=-xmaxq  ! either u/c or u/vnorm units
      TRPG(1) = RVMIN-dvpar
      TRPG(2) = dvpar
      TRPG(3) = 0.0
      TRPG(4) = 0.0-dvprp
      TRPG(5) = 0.0
      TRPG(6) = dvprp
      !FLOGMIN=FLOGMAX
      do ipar=1,npar
      do iprp=1,nprp
         RTEMP2(ipar,iprp)=log10(fparprp(ipar,iprp)) ! to REAL*4
         ! LOG scale - better for color plots?
         !FLOGMIN=MIN(FLOGMIN,RTEMP2(ipar,iprp))
      enddo
      enddo

      ! Draw the map with PGIMAG.
      ! Valid only for RTEMP2(ipar,iprp) over rectangular grid!
      CALL PGIMAG(RTEMP2,npar,nprp,1,npar,1,nprp,FLOGMIN,FLOGMAX,TRPG)
      ! Colorbar:
      CALL pgwedg('RI', 2.0, 5.0, FLOGMIN,FLOGMAX, 'log10(..)')

      endif ! on color option

      ! Overlay contours, with black color:
      CALL PGSCI(1) ! 1==black

      ! Plot original f(i,j) in (v,pitch) coord
      ! (lin.scale of f, but log scale of RCONT)
      CALL PGCONX(RTEMP1,iy,jx,1,iy,1,JXQ,RCONT,NCONT,PGFUNC1)

      ! Or plot the interpolated fparprp(ipar,iprp)
      ! over (vpar,vprp) rectangular grid:
      !CALL PGSLW(1) !lnwidth=1 thin line
      !CALL PGCONT(RTEMP2,npar,nprp,1,npar,1,nprp,RCONTLOG,-NCONT,TR)
      CALL PGSLW(3) !restore lnwidth=3 normal line width



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


      !YuP[03-2016] Added: plot v=vnorm line
      !plot v=vth line, for the k-th gen. species
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
      CALL PGLINE(iy,RTAB1,RTAB2)

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
      CALL PGLINE(iy,RTAB1,RTAB2)
      CALL PGSLS(1) ! 1-> restore solid line
      CALL PGSLW(lnwidth) !lnwidth=3 line width in units of 0.005


      if (k.eq.0) return
      rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
      write(t_,150) n,timet
        CALL PGMTXT('B',6.,0.,0.,t_)
      write(t_,151) rovera(lr_),rr
        CALL PGMTXT('B',7.,0.,0.,t_)
      write(t_,153) rya(lr_), rpcon(lr_), lr_
        CALL PGMTXT('B',8.,0.,0.,t_)
 150  format("time step n=",i5,5x,"time=",1pe10.2," secs")
 151  format( "r/a=",1pe10.3,5x,"radial position (R)=",1pe12.4," cm")
 153  format( "rya=",1pe10.3,5x,"R=rpcon=",1pe12.4," cm,  Surf#",i4)
 999  return
      end
!
!
      subroutine PGFUNC1(VISBLE,yplt,xplt,zplt)
      use param_mod
      use comm_mod
      INTEGER VISBLE
      REAL xplt,yplt,zplt

!
!

!     wx IS V-NORM ARRAY, wy IS THETA ARRAY.  TYPE REAL.
!     xplt (yplt) IS FRACTIONAL INDEX IN V-NORM (THETA) ARRAY.
!     THIS SUBROUTINE MOVES PEN TO NORMALIZED V_PARALLEL,V_PERP COORDS.

      IIX=INT(xplt)
      IF (IIX.GE.1 .AND. IIX.LT.jx) THEN
         XX=wx(IIX) + (xplt-IIX)*(wx(IIX+1)-wx(IIX))
      ELSEIF (IIX.LE.1) THEN
         XX=wx(1)
      ELSE
         XX=wx(jx)
      ENDIF

      IIY=INT(yplt)
      IF (IIY.GE.1 .AND. IIY.LT.iy) THEN
         YY=wy(IIY) + (yplt-IIY)*(wy(IIY+1)-wy(IIY))
      ELSEIF (IIY.LE.1) THEN
         YY=wy(1)
      ELSE
         YY=wy(iy)
      ENDIF

      XWORLD=XX*COS(YY)
      YWORLD=XX*SIN(YY)

!      write(*,*) 'visble,x,y,z,xworld,yworld',visble,x,y,z,xworld,yworld

      IF (VISBLE.EQ.0) THEN
         CALL PGMOVE(XWORLD,YWORLD)
      ELSE
         CALL PGDRAW(XWORLD,YWORLD)
      ENDIF

      RETURN
      END




!====================================================================
!====================================================================
      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
      use iso_c_binding, only : c_float
!-----------------------------------------------------------------------
! Set a "palette" of colors in the range of color indices used by
! PGIMAG.
! From pgdemo4.f
!-----------------------------------------------------------------------
      INTEGER TYPE
      REAL(c_float) CONTRA, BRIGHT
!
      REAL(c_float) GL(2), GR(2), GG(2), GB(2)
      REAL(c_float) RL(9), RR(9), RG(9), RB(9)
      REAL(c_float) HL(5), HR(5), HG(5), HB(5)
      REAL(c_float) WL(10), WR(10), WG(10), WB(10)
      REAL(c_float) AL(20), AR(20), AG(20), AB(20)
!
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
!
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
               0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
               0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
               0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      IF (TYPE.EQ.1) THEN
!        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
!        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
!        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
!        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
!        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      END


end module pltdf_mod
