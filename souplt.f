c
c
      subroutine souplt
      use param_mod
      use cqcomm_mod
      use pltdf_mod, only : cont, tempcntr, nconta, JXQ
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     contour plots the ion source.
c..................................................................


      REAL RILIN
      REAL RTAM1(jx),RTAM2(jx)
      REAL REMAX,REMIN

      if (pltso.eq."disabled") return
      
      
      do 10 k=1,ngen
      
         temp1=0.d0 ! initialize for each k species

         if (xlncur(k,lr_).lt.1.e-10) goto 10
cBH171231         if(frmodp.eq.'enabled')then ! NBI source
           call dcopy(iyjx2,source(0:iyjx2-1,0,k,indxlr_),1,
     +        temp1(0:iyjx2-1,0),1)
cBH171231         endif
         write(t_,550) k
 550     format(1x,"Species ",i2,
     +        " Source Function (units: dist. f/sec)")
         CALL PGPAGE
         itype=3 ! means: plots are made for source
         call pltcont(k,1,t_,itype) ! for source
cBH171231         crnt_nbi=xlncur(k,lr_)*zmaxpsii(lr_) ! [ptcl/sec/cm^3]
         crnt=xlncur(k,lr_)*zmaxpsii(lr_) ! [ptcl/sec/cm^3]

cBH171231         write(t_,540) crnt_nbi
cBH171231 540     format("NBI source rate=",1pe11.4," ptcls/cc/sec")
         write(t_,540) crnt
 540     format("Particle source rate=",1pe11.4," ptcls/cc/sec")
         RILIN=10.
         CALL PGMTXT('B',RILIN,-.2,0.,t_)
         write(t_,542) entr(k,5,l_)
 542     format("Total source power [entr(..5..)]=",1pe11.4," W/cc")
         RILIN=RILIN+2.
         CALL PGMTXT('B',RILIN,-.2,0.,t_)

         write(t_,560)
 560     format("Contour values:")
         RILIN=RILIN+2.
         CALL PGMTXT('B',RILIN,-.2,0.,t_)

         do  jcs=1,ncont,4
            write(t_,570) (tempcntr(jc),jc=jcs,min(jcs+3,ncont))
            if ((ncont/4)*4.ne.ncont .and. ncont-jcs.le.2) then
               icend=4 * 16 + 1
               t_(icend:icend)="$"
            endif
            RILIN=RILIN+1.
            CALL PGMTXT('B',RILIN,-.2,0.,t_)

         enddo
         
 10   continue ! k species

 570  format(4(1pe16.4))
         
c     Plot pitch angle integrated source:
      call pltsofvv

c     Plot the speed-integrated source,
c     as a function of pitch angle theta0 at the midplane    
cyup      call pltso_theta  !YuP[06-2016]
   
      return
      end
c
c=======================================================================
      subroutine pltsofvv
      use param_mod
      use cqcomm_mod
      use r8subs_mod, only : luf
      use r8subs_mod, only : dcopy
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real*8 (a-h,o-z)
c
c     Calculates the pitch angle integrated source
c**********************
c

      REAL RILIN, RXMAXQ
      REAL RTAM1(jx),RTAM2(jx), wk_tam(jx)
      REAL REMAX,REMIN
      real*8 wkd(jx) 

      character*8 target
      character*8 tx_

      
      do 20 k=1,ngen
      
        temp3=0.d0 ! initialize for each k species
        if (xlncur(k,lr_).lt.1.e-10) goto 20
cBH171231        if(frmodp.eq.'enabled')then ! NBI source
          call dcopy(iyjx2,source(0:iyjx2-1,0,k,indxlr_),1,
     +       temp3(0:iyjx2-1,0),1) ! temp3
cBH171231        endif

c-----YuP[2018-01-08] revised to match cqlinput_help:
c**    pltlim= "disabled",  plots versus x (u/vnorm) from
c**                    x=0. to 1. (default:pltlim="disabled",pltlimm=1.)
c**            "x",    plot 1d and 2d plots versus x 
c**                    from 0. to pltlimm.
c**            "u/c",  plot 1d and 2d plots versus u/c
c**                    from 0. to pltlimm.
c**            "energy", plot 1d and 2d plots verus energy (kev)
c**                    from 0. to pltlimm (kev).
         if (pltlim.eq."disabled") then
            target="mainmesh"
            jxq=jx
            xmaxq=x(jxq)
            tx_='u/vnorm' ! or 'u/u\dnorm\u'
            do j=1,jxq
               !tam1(j)=x(j)
               wk_tam(j)=x(j)
            enddo
         endif
         
         if (tandem.eq."enabled" .and. fmass(k).gt.1.e-27) then
            target="ionmesh"
            jxq=jlwr 
            xmaxq=xlwr ! xlwr is set in cqlinput
            !iyjxq=iy*jlwr
            tx_='u/vnorm'
            do j=1,jxq
               !tam1(j)=x(j)
               wk_tam(j)=x(j)
            enddo
            ! If pltlim .ne. "disabled", plot versus
            ! 'x', 'u/c', or 'energy', up to maximum pltlimm.
         elseif (pltlim.eq.'x') then
            target="mainmesh"
            pltlimmm=pltlimm
            jxq=min(luf(pltlimmm,x,jx),jx)
            xmaxq=x(jxq) ! (upper limit)
            !iyjxq=iy*jxq
            tx_='u/vnorm'
            do j=1,jxq
               !tam1(j)=x(j)
               wk_tam(j)=x(j)
            enddo
         elseif (pltlim.eq.'u/c') then
            target="mainmesh"
            pltlimmm=pltlimm
            jxq=min(luf(pltlimmm,uoc,jx),jx)
            xmaxq=uoc(jxq) ! (upper limit); uoc(j)=x(j)/cnorm
            !iyjxq=iy*jxq
            tx_='u/c' ! or TX_='u/c\dlight\u'
            do j=1,jxq
               !tam1(j)=uoc(j)
               wk_tam(j)=uoc(j)
            enddo
         elseif (pltlim.eq.'energy') then
            target="mainmesh"
            pltlimmm=pltlimm
            wkd(1:jx)=enerkev(1:jx,k)
            jxq=min(luf(pltlimmm,wkd,jx),jx)
            xmaxq=enerkev(jxq,k) !YuP[2018-01-08] added 2nd index (k)
            !iyjxq=iy*jxq
            tx_='Energy (keV)'        
            do j=1,jxq
               !tam1(j)=enerkev(j,k) !YuP[2018-01-08] added 2nd index (k)
               wk_tam(j)=enerkev(j,k)
            enddo
         endif
         RXMAXQ=xmaxq
c-----YuP[2018-01-02] done



c       Obtain integrated distribution in tam1
        call fofv(target,"nonorm") ! uses temp3; out: tam1(j)
        !In case of isotropized distribution: 
        !  Int(tam1* 2*pi * vnorm**3 * x**2 *dx) is density
        !YuP/note: It seems the value of "target" has no effect?
        !The integration[summation] is done over all j=1:jx and i=1:iy

        ! Vertical axis:
        call aminmx(tam1,1,jxq,1,fmin,fmax,kmin,kmax)
        !write(*,*)'pltsofvv: FMIN,FMAX',FMIN,FMAX ! can be 0
        
        if(fmax.lt.em90)then
          fmax=em90
        endif
        
        fmin=1.e-08*fmax
        do jj=1,jxq
          if (tam1(jj) .lt. fmin ) tam1(jj)=fmin
        enddo



        DO J=1,JXQ
           RTAM1(J)=wk_tam(J) ! either u/c or u/vnorm
           TAM2(J)=ABS(TAM1(J))
           TAM2(J)=MAX(em300,TAM2(j))
           RTAM2(J)=LOG10(TAM2(J))
        ENDDO

        REMIN=LOG10(fmin)
        REMAX=LOG10(fmax)

        CALL PGPAGE
        CALL PGSVP(.2,.8,.45,.9)
        IF ( Remax-Remin .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           Remax= Remin+1.e-16
        ENDIF
        CALL PGSWIN(Rtam1(1),RXMAXQ,Remin,Remax)
        CALL PGBOX('BCNST',0.,0,'BCNSTL',0.,0)
        CALL PGSAVE
        CALL PGSCH(1.44)
        CALL PGLAB(tx_, 'Source', 'Pitch Angle Avg Source vs. u')
        CALL PGUNSA
        CALL PGLINE(JXQ,RTAM1,RTAM2)



        write(t_,10040) k
10040   format("Particle source integrated over theta0 for species",i3)

        RILIN=8.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)

        write(t_,10041) 
10041   format("(normed so int(0,1)*2pi*x**2*dx=mid-plane source)")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)

        write(t_,10042) vnorm
10042   format("vnorm=",1x,1pe12.4," cm/s")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)


        rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
        write(t_,10030) n,timet
10030   format("time step (n) is",i5,5x,"time=",1pe12.4," secs")
        RILIN=RILIN+2.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)

        write(t_,10031) rovera(lr_),rr
10031   format("r/a=",1pe12.4,5x,"radial position (R) =",1pe12.4," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)


 20   continue ! k species


      return
      end


c
c=======================================================================
      subroutine pltso_theta
      use param_mod
      use cqcomm_mod
      use r8subs_mod, only : dcopy
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real*8 (a-h,o-z)
c     YuP[06-2016]
c     Calculates and plots the speed-integrated source,
c     as a function of pitch angle theta0 at the midplane 
c     Called for every flux surface lr_ (lr_ is stored in comm.h)
c**********************
c
CMPIINSERT_INCLUDE

      REAL RILIN
      REAL RTAM1(iy),RTAM2(iy)
      REAL REMAX,REMIN
      real*8 wk_so(iy)
      character*8 vert_scale ! 'log10' or 'lin'
      
      vert_scale='log10' !'lin'
      do i=1,iy
        RTAM1(i)=y(i,lr_)*180.0/pi  ! horizontal.axis: theta0 (degree)
      enddo
      
      do 20 k=1,ngen ! sources for each general sp. are plotted
      
        temp3=0.d0 ! initialize for each k species (i,j)
        if (xlncur(k,lr_).lt.1.e-10) goto 20
        
        if(frmodp.eq.'enabled')then ! NBI source
          call dcopy(iyjx2,source(0:iyjx2-1,0,k,indxlr_),1,
     +          temp3(0:iyjx2-1,0),1) ! temp3
        endif
         

c       Obtain integrated distribution into wk_so(1:iy)
        do i=1,iy
           wk_so(i)=0. ! initialize
        do j=1,jx
           wk_so(i)= wk_so(i)+ temp3(i,j)*cint2(j)
           !cint2= x**2 *dx  (and remember that temp3 includes vnorm^3)
           !cynt2= 2pi*sin(theta0)*dtheta0
        enddo
        WRITE(*,'(a,2i5,2e13.5)')
     +   'pltso_theta: lr_,i,y(i,lr_)*180.0/pi,wk_so(i)=',
     +        lr_,i, y(i,lr_)*180.0/pi, wk_so(i)
        enddo 
        !In case of isotropized distribution: 
        !  Int(wk_so(i)*cynt2(i)) is density (per sec.)

        call aminmx(wk_so,1,iy,1,fmin,fmax,kmin,kmax)
        if(fmax.lt.em90) goto 20 !-> Almost no source , next k species
        
        if(vert_scale.eq.'log10')then
          fmin=1.e-03*fmax ! limit the lower range (for log scale)
          do i=1,iy
            if (wk_so(i) .lt. fmin ) wk_so(i)=fmin
            RTAM2(i)=LOG10(wk_so(i))
          enddo
          REMIN=LOG10(fmin)
          REMAX=LOG10(fmax)
        else ! lin scale
          IF ( fmax-fmin .le. 1.e-16 ) THEN 
           ! fmax~fmin => extend plot limits a bit
           fmax= fmax+1.e-16
           fmin= fmin-1.e-16
          ENDIF
          RTAM2(1:iy)=wk_so(1:iy) ! Units: vnorm^3 *(reactions/sec)/cm^3 
          ! Should be divided by vnorm^3 to obtain physical units 
          ! (reactions/sec)/cm^3 /(cm/sec)^3
          REMIN=fmin*0.99
          REMAX=fmax*1.01
        endif

        CALL PGPAGE
        CALL PGSVP(.2,.8,.45,.9)
        CALL PGSWIN(Rtam1(1),Rtam1(iy),Remin,Remax)
        if(vert_scale.eq.'log10')then
        CALL PGBOX('BCNST',0.,0,'BCNSTL',0.,0)
        else
        CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
        endif
        CALL PGSAVE
        CALL PGSCH(1.44)
        CALL PGLAB('theta0 (degree)','S0(theta0)','v-integrated Source')
        CALL PGUNSA
        CALL PGLINE(iy,RTAM1,RTAM2)

        write(t_,10040) k
10040   format("Particle source integrated over v for species",i2)

        RILIN=8.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)

        write(t_,10041) 
10041   format("(int(0,pi)*S0*2pi*sin(theta0)*dtheta0= ptcls/sec)")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)

        rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
        write(t_,10030) n,timet
10030   format("time step (n) is",i5,5x,"time=",1pe12.4," secs")
        RILIN=RILIN+2.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)

        write(t_,10031) rovera(lr_),rr
10031   format("r/a=",1pe12.4,5x,"radial position (R)=",1pe12.4," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,-.2,0.,t_)


 20   continue ! k species


      return
      end
c

