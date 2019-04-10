c
c
      module pltmain_mod
      use iso_c_binding, only : c_float
      character*7, private :: scale
      real(c_float), private ::  X, Y, ANGLE, FJUST
            
      save

      contains

      subroutine pltmain
      use param_mod
      use cqcomm_mod
      use pltdf_mod, only: pltdf
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
c     This routine controls driver plots
c
c
c     Modified some Graflib to pgplot calls by Yuri Petrov, 090727,
c     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
c
      REAL RILIN
      
CMPIINSERT_INCLUDE

CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only
      
      if (noplots.eq."enabled1") return
      if (n.eq.0 .and. lrzmax.gt.1) return
      if (mplot(l_).eq."disabled") return
      
      rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
      if (pltend.eq."notplts") goto 10
      if (pltend.eq."last" .and. n.lt.nstop) goto 10

      CALL PGPAGE
      CALL PGSCH(1.0) ! restore to default font size
      !(sometimes font is too big from previous plot)

      RILIN=0.
      if (cqlpmod .ne. "enabled") then
         CALL PGMTXT('T',-RILIN,0.,0.,"LOCAL RADIAL QUANTITIES")
      else 
         CALL PGMTXT('T',-RILIN,0.,0.,"LOCAL PARALLEL QUANTITIES")
      endif
      RILIN=RILIN+1.

      write(t_,150) n,timet
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,1501) lr_,lrz
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,151) rovera(lr_),rr
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,153) rya(lr_),rpcon(lr_)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
 150  format("time step n=",i5,","5x,"time=",1pe12.4," secs")
 1501 format("flux surf=",i3,5x,"total flux surfs=",i3)
 151  format("r/a=",1pe10.3,5x,"radial position (R)=",1pe12.4," cms")
 153  format("rya=",1pe10.3,5x,"R=rpcon=",1pe10.3," cm")
 
      if (cqlpmod .eq. "enabled") then
        write(t_,152) l_,sz(l_)
        RILIN=RILIN+1.
        CALL PGMTXT('T',-RILIN,0.,0.,t_)
 152    format("orbit at s(",i5,") = ",1pe10.2)
      endif

      vnormdc=vnorm/clight
      vtedc=vthe(lr_)/clight
      vtdvnorm=vthe(lr_)/vnorm
      if (tandem.eq."enabled") then ! YuP[08-2017] added, 
        ! for the case of ngen=2 tandem i+e runs:
        !In this case, enorm=enorme
        !and  xlwr=sqrt(enormi*fmass(kelecg)/(enorme*fmass(kionn)))
        write(t_,'(a,2f11.3)') ' enormi, enorme(=enorm) (kev) =',
     +                           enormi, enorme
      else
        write(t_,'(a,f11.3)') ' enorm (kev) =' ,enorm
      endif
        RILIN=RILIN+1.
        CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,161)  vnormdc
        RILIN=RILIN+1.
        CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,162)  vtedc
        RILIN=RILIN+1.
        CALL PGMTXT('T',-RILIN,0.,0.,t_)
      write(t_,163)  vtdvnorm
        RILIN=RILIN+1.
        CALL PGMTXT('T',-RILIN,0.,0.,t_)
      do k=1,ntotal
         vtdvnorm= vth(k,lr_)/vnorm
         !YuP/note: For time-dependent profiles, 
         ! temp() can evolve in time, and so vth(k,*) can evolve, too.
         ! See profiles.f.
         write(t_,'(a,i2,a,f15.7)') "k=",k, "  vth(k)/vnorm =", vtdvnorm
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)
      enddo

 160  format("enorm (kev) = ",f11.3)
 161  format("vnorm/c = ",f15.7)
 162  format("vthe (sqrt(te/me))/c = ",f15.7)
 163  format("vthe/vnorm = ",f15.7)
 
      if (cqlpmod .eq. "enabled") then
        zvthes=vth(kelec,l_)/clight
        zvtheon=vth(kelec,l_)/vnorm
        write(t_,164) zvthes
        write(t_,165) zvtheon
        RILIN=RILIN+1.
        CALL PGMTXT('T',-RILIN,0.,0.,t_)
 164    format(";","vthe(s) (sqrt(te/me))/c = ",f15.7)
 165    format("vthe(s)/vnorm = ",f15.7)
      endif

      ncplt=ncplt+1
 10   continue

c     Plot energy, density, toroidal current, conservation diagnostic vs
c     time
c
      if (pltend.ne."disabled" .and. nch(l_).ge.2) then
        if (cqlpmod .ne. "enabled") call pltendn
        if (cqlpmod .eq. "enabled") call pltends
      endif
c
c     Plot ion source if marker is engaged.
c
      isouplt=0
      if (pltso.eq."enabled" .or. pltso.eq."color") then
        call souplt
      elseif ( (pltso.eq."first" .or. pltso.eq."first_cl") .and.
     +         isouplt.eq.0 .and.
     +         n.ge.nonso(1,1) .and. n .le. noffso(1,1)) then
        call souplt
        isouplt=1
      endif
c
c     Plot electron resistivity and related quantities.
c
      if (abs(elecfld(lr_)).gt.1.e-9 .and. n.ge.nonel) then
        if (pltrst.eq."enabled" .and. nch(l_).ge.2 .and. kelecg.gt.0)
     +    then
          call pltrstv
        endif
      endif
c
c     Plot normalized v-flux..
c
      if (abs(elecfld(lr_)).gt.1.e-9 .and. n.ge.nonel) then
cBH        if (pltvflu.eq."enabled" .and. n.gt.1 .and. kelecg.gt.0) then
        if (pltvflu.eq."enabled" .and. n.ge.1 .and. kelecg.gt.0) then
          call pltvflux
        endif
      endif
c
c     Plot power deposited in plasma by various mechanisms vs.
c     time.
c
      if (pltpowe.eq."enabled" .and. nch(l_).ge.2) then
        call pltpower
      elseif (pltpowe.eq."last" .and. n.eq.nstop) then
        call pltpower
      endif
c...  
cmnt  Distribution f slices at constant pitch and/or pitch angle 
cmnt    averaged f,  vs velocity or energy coordinate.
c...  
      if (pltfvs.eq."enabled".or.pltfofv.eq."enabled") call pltfvsv
c...  
c...  
cmnt  Distribution fluxes vs velocity for some values of theta
c...  
      if (pltflux.eq."enabled") call pltfluxs


c...  
cmnt  Contour the loss region.
c...  
      if (pltlos.ne."disabled") call pltlosc
c...

  
cmnt  Plot the reduced parallel distribution (f integrated on vpp)
c...  
      if (pltprpp.eq."enabled") call pltprppr
c
c     Plot the density as a function of poloidal angle for a
c     set of energy ranges..
c
      if (n.ne.0 .and. pltdn.ne."disabled" .and. cqlpmod.ne."enabled") 
     +  call pltdnz
c
c     plot contours of df/dt next..
c
      if (pltd.ne."disabled") call pltdf
c
c     vector flux plots follow..
c
      if (n.eq.0 .or. pltvecal.eq."disabled") goto 670
      call pltvec(4)
      if (pltvecrf.ne."disabled") call pltvec(3)
      if (pltvece.ne."disabled") call pltvec(2)
      if (pltvecc.ne."disabled") call pltvec(1)
 670  continue
c
c     Plot the stream lines of the steady state flux
c
      if (pltstrm.ne."disabled" .and. n.ge.1) call pltstrml

      return
      end





C==== CONVERT SOME GRAFLIB ROUTINES to PGPLOT ========================
C  Yuri Petrov, 090727
c---------------------------------------------------------------------
      subroutine gxglfr(n) 
      integer n       
        CALL PGPAGE
      return 
      end
c---------------------------------------------------------------------
      subroutine gsvp2d(xmin,xmax,ymin,ymax) ! called explicitly with real*4 args
      ! xmin,...defines where on the frame the data is plotted. 
      REAL*8 xmin,xmax,ymin,ymax
      REAL*4 PGxmin,PGxmax,PGymin,PGymax ! PGPLOT uses REAL*4
        PGxmin = xmin
        PGxmax = xmax
        PGymin = ymin
        PGymax = ymax
        CALL PGSVP(PGxmin,PGxmax,PGymin,PGymax)
c PGSVP (XLEFT, XRIGHT, YBOT, YTOP)
c XLEFT  (input)  : x-coordinate of left hand edge of viewport, in NDC.
c XRIGHT (input)  : x-coordinate of right hand edge of viewport,in NDC.
c YBOT   (input)  : y-coordinate of bottom edge of viewport, in NDC.
c YTOP   (input)  : y-coordinate of top  edge of viewport, in NDC.        
      return 
      end
c---------------------------------------------------------------------
      subroutine gswd2d(scales,xmin,xmax,ymin,ymax) ! NOT USED ANYMORE?
      use r8subs_mod, only : rbound
      ! xmin,... defines the user coordinate system.            
      implicit integer (i-n), real*8 (a-h,o-z)
      REAL*4 PGxmin,PGxmax,PGymin,PGymax, RPG1,RPG2 ! PGPLOT uses REAL*4
      !REAL*4 RBOUND
      character*7 scales ! = "linlin$" or "linlog$","loglin$","loglog$"
        scale  = scales ! To gpcv2d -> PGLINE
        PGxmin = RBOUND(xmin)
        PGxmax = RBOUND(xmax)
        PGymin = RBOUND(ymin)
        PGymax = RBOUND(ymax)
        IF ( PGymax-PGymin .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           PGymax= PGymin+1.e-16
        ENDIF
        CALL PGSCH(1.) ! set character size; default is 1.              
        if(scale.eq."linlin$") then
          CALL PGSWIN(PGxmin,PGxmax,PGymin,PGymax)
          CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
        endif
        if(scale.eq."loglin$") then 
          CALL PGSWIN(log10(PGxmin),log10(PGxmax),PGymin,PGymax)
          CALL PGBOX('BCNSTL',0.,0,'BCNST',0.,0)
        endif
        !----------------------------
        PGymin= max(PGymin,1.e-32) ! cannot be negative
        PGymax= max(PGymax,1.e-32) ! cannot be negative
        RPG1= log10(PGymin)
        RPG2= log10(PGymax)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        if(scale.eq."linlog$") then 
          CALL PGSWIN(PGxmin,PGxmax,RPG1,RPG2)
          CALL PGBOX('BCNST',0.,0,'BCNSTL',0.,0)
        endif
        if(scale.eq."loglog$") then
          CALL PGSWIN(log10(PGxmin),log10(PGxmax),RPG1,RPG2)
          CALL PGBOX('BCNSTL',0.,0,'BCNSTL',0.,0)
        endif
      return 
      end
c---------------------------------------------------------------------
      subroutine gpcv2d(xarray,yarray,length)   
      ! Plot a line (xarray,yarray) of length n.    
      implicit integer (i-n), real*8 (a-h,o-z)
      real*8 xarray(length),yarray(length)
      integer length
      REAL*4 PGxarray(length),PGyarray(length)
      small_p = EPSILON(1.0) !a positive number that is almost negligible
        if(scale.eq."linlin$") then
          do n=1,length
             PGxarray(n)= xarray(n) ! Convert to REAL*4 for PGPLOT
             PGyarray(n)= yarray(n) ! Convert to REAL*4 for PGPLOT
          enddo
        endif        
        if(scale.eq."linlog$") then
          do n=1,length
             PGxarray(n)= xarray(n) ! Convert to REAL*4 for PGPLOT
             PGyarray(n)= log10( max(small_p,abs(yarray(n))) ) 
          enddo
        endif
        if(scale.eq."loglin$") then
          do n=1,length
             PGxarray(n)= log10( max(small_p,abs(xarray(n))) )
             PGyarray(n)= yarray(n) ! Convert to REAL*4 for PGPLOT
          enddo
        end if
        if(scale.eq."loglog$") then
          do n=1,length
             PGxarray(n)= log10( max(small_p,abs(xarray(n))) ) 
             PGyarray(n)= log10( max(small_p,abs(yarray(n))) )
          enddo
        endif
        CALL PGLINE(length,PGxarray,PGyarray) 
c Primitive routine to draw a Polyline. A polyline is one or more
c connected straight-line segments.  The polyline is drawn using
c the current setting of attributes color-index, line-style, and
c line-width. The polyline is clipped at the edge of the window.
c  N=length (input): number of points defining the line; the line
c                    consists of (N-1) straight-line segments.
c                    N should be greater than 1 (if it is 1 or less,
c                    nothing will be drawn).
c  X=xarray (input): world x-coordinates of the points.
c  Y=yarray (input): world y-coordinates of the points.
c The dimension of arrays X and Y must be greater than or equal to N.
c The "pen position" is changed to (X(N),Y(N)) in world coordinates
c (if N > 1).
      return 
      end
c---------------------------------------------------------------------
      subroutine gpln2d(x1, x2, y1, y2) 
      ! draw line between two points.
      ! x1,x2,y1,y2 defines the two points to draw a line between. 
      implicit integer (i-n), real*8 (a-h,o-z)
      REAL*4 PGx1, PGy1, PGx2, PGy2
        PGx1 = x1 ! Convert to REAL*4
        PGx2 = x2 ! Convert to REAL*4
        PGy1 = y1 ! Convert to REAL*4
        PGy2 = y2 ! Convert to REAL*4
        CALL PGMOVE (PGx1, PGy1)
c Move the "pen" to the point with world
c coordinates (X,Y). No line is drawn.
        CALL PGDRAW (PGx2, PGy2)
c Draw a line from the current pen position to the point
c with world-coordinates (X,Y). The line is clipped at the edge of the
c current window. The new pen position is (X,Y) in world coordinates.
      return 
      end
c---------------------------------------------------------------------
      subroutine gslnsz(size) ! called explicitly with real*4 args
      ! set line size (width), where 0. is the default.
      real*8 size
      INTEGER  LW
      LW = int(size*10. + 1.) !-YuP: Not sure if this conversion is ok
        CALL PGSLW(LW)
      ! Set the line-width attribute. This attribute affects lines, graph
      ! markers, and text. The line width is specified in units of 1/200 
      ! (0.005) inch (about 0.13 mm) and must be an integer in the range
      ! 1-201. On some devices, thick lines are generated by tracing each
      ! line with multiple strokes offset in the direction perpendicular 
      ! to the line.
      return 
      end
c---------------------------------------------------------------------
      subroutine gslnst(LS) 
      ! sets line style: 1-solid, 2-dashed, 3-dotted, 4-dash-dotted, etc.
      implicit integer (i-n), real*8 (a-h,o-z)
      INTEGER  LS
        CALL PGSLS(LS)
c Set the line style attribute for subsequent plotting. This
c attribute affects line primitives only; it does not affect graph
c markers, text, or area fill.
c Five different line styles are available, with the following codes:
c 1 (full line), 2 (dashed), 3 (dot-dash-dot-dash), 4 (dotted),
c 5 (dash-dot-dot-dot). The default is 1 (normal full line).
      return 
      end
c---------------------------------------------------------------------
      subroutine  gscpvs(gl_x,gl_y) ! called explicitly with real*4 args
      ! set current position for text
      real*4 gl_x,gl_y
      REAL*4 X,Y,ANGLE,FJUST
      common/GRAFLIB_PGPLOT_text/ X,Y,ANGLE,FJUST
        X = gl_x  ! To PGPTXT(X,Y,ANGLE,FJUST,TEXT)
        Y = gl_y
      return 
      end
c---------------------------------------------------------------------
      subroutine  gstxan(gl_angle) ! called explicitly with real*4 args
      ! set angle to plot text
      real*4 gl_angle
      REAL*4 X,Y,ANGLE,FJUST
      common/GRAFLIB_PGPLOT_text/ X,Y,ANGLE,FJUST
        ANGLE = gl_angle ! To PGPTXT(X,Y,ANGLE,FJUST,TEXT) 
      return 
      end
c---------------------------------------------------------------------
      subroutine  gstxjf(just1,just2)
      ! set justification of string
      character*(*) just1 !can be 'left', 'right', or 'center'
      character*(*) just2 !can be 'top', 'bottom', or 'center'
      REAL*4 X,Y,ANGLE,FJUST
      common/GRAFLIB_PGPLOT_text/ X,Y,ANGLE,FJUST
         FJUST = 0.0
         if(just1.eq."left")   FJUST = 0.0
         if(just1.eq."center") FJUST = 0.5
         if(just1.eq."right")  FJUST = 1.0
         ! To PGPTXT(X,Y,ANGLE,FJUST,TEXT)
      return 
      end
c---------------------------------------------------------------------
      subroutine  gstxno(size) ! called explicitly with real*4 args
      ! set number of characters per line, i.e., character width.
      ! Default: size=90.  (not sure)
      real*4 size
      REAL*4 PGsize 
         PGsize = size ! Convert to REAL*4
         CALL PGSCH(90./PGsize) ! set character size; default is 1.
      return 
      end
c---------------------------------------------------------------------
      subroutine  gptx2d(text) 
      ! Plot Text
      character*(*) text
      REAL*4 X,Y,ANGLE,FJUST
      common/GRAFLIB_PGPLOT_text/ X,Y,ANGLE,FJUST
        call PGPTXT(X, Y, ANGLE, FJUST, text)
c Primitive routine for drawing text. The text may be drawn at any
c angle with the horizontal, and may be centered or left- or right-
c justified at a specified position.  Routine PGTEXT provides a
c simple interface to PGPTXT for horizontal strings. Text is drawn
c using the current values of attributes color-index, line-width,
c character-height, and character-font.  Text is NOT subject to
c clipping at the edge of the window.
c X      (input)  : world x-coordinate.
c Y      (input)  : world y-coordinate. The string is drawn with the
c                   baseline of all the characters passing through
c                   point (X,Y); the positioning of the string along
c                   this line is controlled by argument FJUST.
c ANGLE  (input)  : angle, in degrees, that the baseline is to make
c                   with the horizontal, increasing counter-clockwise
c                   (0.0 is horizontal).
c FJUST  (input)  : controls horizontal justification of the string.
c                   If FJUST = 0.0, the string will be left-justified
c                   at the point (X,Y); if FJUST = 0.5, it will be
c                   centered, and if FJUST = 1.0, it will be right
c                   justified. [Other values of FJUST give other
c                   justifications.]
c TEXT   (input)  : the character string to be plotted.
      return 
      end
c---------------------------------------------------------------------
      

      end module pltmain_mod
