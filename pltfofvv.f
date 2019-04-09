c
c
      subroutine pltfofvv
      use param_mod
      use pltmain_mod, only : gpcv2d, gsvp2d, gswd2d, gxglfr
      implicit integer (i-n), real*8 (a-h,o-z)
c
**SCC10/7/94
c   prppru calculates the pitch angle integrated distribution function
c**********************
c
c
c     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
c     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
c
      include 'comm.h'

      REAL RILIN !-> For PGPLOT (text output positioning)

      character*8 target
      do 20 k=1,ngen
        call GXGLFR(0) ! new page(s)
        call dcopy(iyjx2,f(0,0,k,l_),1,temp3(0,0),1)
        if (tandem.eq."enabled" .and. k.eq.kionn) then
          target="ionmesh"
          jxq=jlwr
        else
          target="mainmesh"
          jxq=jx
        endif

c       Obtain integrated distribution in tam1
        call fofv(target,"nonorm")
        !YuP/note: It seems the value of "target" has no effect?
        !The integration[summation] is done over all j=1:jx and i=1:iy

        call aminmx(tam1,1,jxq,1,fmin,fmax,kmin,kmax)
        fmin=1.d-08*fmax
        call GSVP2D(.2,.8,.25,.95)
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlog$",x(1),x(jxq),fmin,fmax)
        do jj=1,jxq
          if (tam1(jj) .lt. fmin ) tam1(jj)=fmin
        enddo
        
        call GPCV2D(x,tam1,jxq)
        
        write(t_,10040) k
10040 format("distribution integrated over theta0 for species",2x,i5)
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10041) 
10041 format("(normed so that int(0,1)*2pi*x**2*dx=mid-plane ne)")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10042) vnorm
10042 format("vnorm=",1x,e14.6,"cm/s")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)
        
        write(t_,10020)
10020 format( "(log plot)")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)
        
        rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
        write(t_,10030) n,timet
10030 format("time step (n) is",i5,5x,"time=",e14.6," secs")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10031) rovera(lr_)
10031 format("r/a=",e14.6)
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10032) rr
10032 format("radial position (R) =",e14.6," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        CALL PGSCH(1.0) ! recover default 1.0 fontsize
        
 20   continue


      return
      end


c
c
      subroutine fofv(target,action)
      use param_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      character*(*) target,action
      save

c...............................................................
cmnt  This routine takes data stored in temp3 in (y,x) coordinates
c     and integrates over theta0 to get tam1 which is, in this case,
c     the isotropized distribution; Int(tam1*2*pi*vnorm**3*x**2*dx) is
c     the density.
c     SCChiu  10/7/94
c      xul,xuu: lower and upper values of x (normalized u)
c      ytl,ytu:  "     "   "      "    "  y (theta0)
c...............................................................

      include 'comm.h'

      logical trnsfm
      trnsfm=(target.eq."velocity".and. relativ .ne. "disabled")

        do jp=1,jx
          tam2(jp)=0.
          do ip=1,iy-1
            ip1=ip+1
            tam2(jp)=tam2(jp)+cynt2(ip,l_)*(temp3(ip,jp)
     1               +temp3(ip1,jp))
          enddo
          tam1(jp)=tam2(jp)/twopi
        enddo

      return
      end
