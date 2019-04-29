module pltmain_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use bcast_mod, only : bcast
  use coefefad_mod, only : coefefad
  use coefegad_mod, only : coefegad
  use coeffpad_mod, only : coeffpad
  use coefload_mod, only : coefload
  use coefmidt_mod, only : coefmidt
  use coefmidv_mod, only : coefmidv
  use coefrfad_mod, only : coefrfad
  use coefstup_mod, only : coefstup
  use coefsyad_mod, only : coefsyad
  use diagentr_mod, only : gfi
  use fle_mod, only : fle_fsa
  use fle_mod, only : fle_pol
  !XXXX use pack21_mod, only : pack21
  use pltdf_mod, only : pltcont
  use pltdf_mod, only : pltdf
  use pltdf_mod, only : pgfunc1
  use pltdnz_mod, only : pltdnz
  use pltendn_mod, only : pltendn
  use pltfvsv_mod, only : pltfvsv
  use pltlosc_mod, only : pltlosc
  use pltpower_mod, only : pltpower
  use pltvec_mod, only : pltvec
  use pltvflux_mod, only : pltvflux
  use r8subs_mod, only : dcopy

  !---END USE

  use iso_c_binding, only : c_float

  !XXXX
  external pack21

  character*7, private :: scale
  real(c_float), private ::  X, Y, ANGLE, FJUST
  save

contains

  subroutine pltmain
    use param_mod
    use comm_mod
    use pltdf_mod, only: pltdf
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save
    !
    !     This routine controls driver plots
    !
    !
    !     Modified some Graflib to pgplot calls by Yuri Petrov, 090727,
    !     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
    !
    REAL RILIN

    !MPIINSERT_INCLUDE

    !MPIINSERT_IF_RANK_NE_0_RETURN
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
150 format("time step n=",i5,","5x,"time=",1pe12.4," secs")
1501 format("flux surf=",i3,5x,"total flux surfs=",i3)
151 format("r/a=",1pe10.3,5x,"radial position (R)=",1pe12.4," cms")
153 format("rya=",1pe10.3,5x,"R=rpcon=",1pe10.3," cm")

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
       write(t_,'(a,2f11.3)') ' enormi, enorme(=enorm) (kev) =',enormi, enorme
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

160 format("enorm (kev) = ",f11.3)
161 format("vnorm/c = ",f15.7)
162 format("vthe (sqrt(te/me))/c = ",f15.7)
163 format("vthe/vnorm = ",f15.7)

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
10  continue

    !     Plot energy, density, toroidal current, conservation diagnostic vs
    !     time
    !
    if (pltend.ne."disabled" .and. nch(l_).ge.2) then
       if (cqlpmod .ne. "enabled") call pltendn
       if (cqlpmod .eq. "enabled") call pltends
    endif
    !
    !     Plot ion source if marker is engaged.
    !
    isouplt=0
    if (pltso.eq."enabled" .or. pltso.eq."color") then
       call souplt
    elseif ( (pltso.eq."first" .or. pltso.eq."first_cl") .and. &
         isouplt.eq.0 .and. &
         n.ge.nonso(1,1) .and. n .le. noffso(1,1)) then
       call souplt
       isouplt=1
    endif
    !
    !     Plot electron resistivity and related quantities.
    !
    if (abs(elecfld(lr_)).gt.1.e-9 .and. n.ge.nonel) then
       if (pltrst.eq."enabled" .and. nch(l_).ge.2 .and. kelecg.gt.0) then
          call pltrstv
       endif
    endif
    !
    !     Plot normalized v-flux..
    !
    if (abs(elecfld(lr_)).gt.1.e-9 .and. n.ge.nonel) then
       !BH        if (pltvflu.eq."enabled" .and. n.gt.1 .and. kelecg.gt.0) then
       if (pltvflu.eq."enabled" .and. n.ge.1 .and. kelecg.gt.0) then
          call pltvflux
       endif
    endif
    !
    !     Plot power deposited in plasma by various mechanisms vs.
    !     time.
    !
    if (pltpowe.eq."enabled" .and. nch(l_).ge.2) then
       call pltpower
    elseif (pltpowe.eq."last" .and. n.eq.nstop) then
       call pltpower
    endif
    !...
    !mnt  Distribution f slices at constant pitch and/or pitch angle
    !mnt    averaged f,  vs velocity or energy coordinate.
    !...
    if (pltfvs.eq."enabled".or.pltfofv.eq."enabled") call pltfvsv
    !...
    !...
    !mnt  Distribution fluxes vs velocity for some values of theta
    !...
    if (pltflux.eq."enabled") call pltfluxs


    !...
    !mnt  Contour the loss region.
    !...
    if (pltlos.ne."disabled") call pltlosc
    !...


    !mnt  Plot the reduced parallel distribution (f integrated on vpp)
    !...
    if (pltprpp.eq."enabled") call pltprppr
    !
    !     Plot the density as a function of poloidal angle for a
    !     set of energy ranges..
    !
    if (n.ne.0 .and. pltdn.ne."disabled" .and. cqlpmod.ne."enabled") call pltdnz
    !
    !     plot contours of df/dt next..
    !
    if (pltd.ne."disabled") call pltdf
    !
    !     vector flux plots follow..
    !
    if (n.eq.0 .or. pltvecal.eq."disabled") goto 670
    call pltvec(4)
    if (pltvecrf.ne."disabled") call pltvec(3)
    if (pltvece.ne."disabled") call pltvec(2)
    if (pltvecc.ne."disabled") call pltvec(1)
670 continue
    !
    !     Plot the stream lines of the steady state flux
    !
    if (pltstrm.ne."disabled" .and. n.ge.1) call pltstrml

    return
  end subroutine pltmain





  !==== CONVERT SOME GRAFLIB ROUTINES to PGPLOT ========================
  !  Yuri Petrov, 090727
  !---------------------------------------------------------------------
  subroutine gxglfr(n)
    integer n
            CALL PGPAGE
    return
  end subroutine gxglfr
  !---------------------------------------------------------------------
  subroutine gsvp2d(xmin,xmax,ymin,ymax)
    ! xmin,...defines where on the frame the data is plotted.
    REAL xmin,xmax,ymin,ymax ! let these take the compiler default...
    real(c_float) PGxmin,PGxmax,PGymin,PGymax ! PGPLOT uses real(c_float)
    PGxmin = xmin
    PGxmax = xmax
    PGymin = ymin
    PGymax = ymax
            CALL PGSVP(PGxmin,PGxmax,PGymin,PGymax)
    ! PGSVP (XLEFT, XRIGHT, YBOT, YTOP)
    ! XLEFT  (input)  : x-coordinate of left hand edge of viewport, in NDC.
    ! XRIGHT (input)  : x-coordinate of right hand edge of viewport,in NDC.
    ! YBOT   (input)  : y-coordinate of bottom edge of viewport, in NDC.
    ! YTOP   (input)  : y-coordinate of top  edge of viewport, in NDC.
    return
  end subroutine gsvp2d

  !---------------------------------------------------------------------
  subroutine gswd2d(scales,xmin,xmax,ymin,ymax) ! NOT USED ANYMORE?
    use r8subs_mod, only : rbound
    ! xmin,... defines the user coordinate system.
    implicit integer (i-n), real(c_double) (a-h,o-z)
    real(c_float) PGxmin,PGxmax,PGymin,PGymax, RPG1,RPG2 ! PGPLOT uses real(c_float)
    !real(c_float) RBOUND
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
  end subroutine gswd2d
  !---------------------------------------------------------------------
  subroutine gpcv2d(xarray,yarray,length)
    ! Plot a line (xarray,yarray) of length n.
    implicit integer (i-n), real(c_double) (a-h,o-z)
    real(c_double) xarray(length),yarray(length)
    integer length
    real(c_float) PGxarray(length),PGyarray(length)
    small_p = EPSILON(1.0) !a positive number that is almost negligible
    if(scale.eq."linlin$") then
       do n=1,length
          PGxarray(n)= xarray(n) ! Convert to real(c_float) for PGPLOT
          PGyarray(n)= yarray(n) ! Convert to real(c_float) for PGPLOT
       enddo
    endif
    if(scale.eq."linlog$") then
       do n=1,length
          PGxarray(n)= xarray(n) ! Convert to real(c_float) for PGPLOT
          PGyarray(n)= log10( max(small_p,abs(yarray(n))) )
       enddo
    endif
    if(scale.eq."loglin$") then
       do n=1,length
          PGxarray(n)= log10( max(small_p,abs(xarray(n))) )
          PGyarray(n)= yarray(n) ! Convert to real(c_float) for PGPLOT
       enddo
    end if
    if(scale.eq."loglog$") then
       do n=1,length
          PGxarray(n)= log10( max(small_p,abs(xarray(n))) )
          PGyarray(n)= log10( max(small_p,abs(yarray(n))) )
       enddo
    endif
            CALL PGLINE(length,PGxarray,PGyarray)
    ! Primitive routine to draw a Polyline. A polyline is one or more
    ! connected straight-line segments.  The polyline is drawn using
    ! the current setting of attributes color-index, line-style, and
    ! line-width. The polyline is clipped at the edge of the window.
    !  N=length (input): number of points defining the line; the line
    !                    consists of (N-1) straight-line segments.
    !                    N should be greater than 1 (if it is 1 or less,
    !                    nothing will be drawn).
    !  X=xarray (input): world x-coordinates of the points.
    !  Y=yarray (input): world y-coordinates of the points.
    ! The dimension of arrays X and Y must be greater than or equal to N.
    ! The "pen position" is changed to (X(N),Y(N)) in world coordinates
    ! (if N > 1).
    return
  end subroutine gpcv2d
  !---------------------------------------------------------------------
  subroutine gpln2d(x1, x2, y1, y2)
    ! draw line between two points.
    ! x1,x2,y1,y2 defines the two points to draw a line between.
    implicit integer (i-n), real(c_double) (a-h,o-z)
    real(c_float) PGx1, PGy1, PGx2, PGy2
    PGx1 = x1 ! Convert to real(c_float)
    PGx2 = x2 ! Convert to real(c_float)
    PGy1 = y1 ! Convert to real(c_float)
    PGy2 = y2 ! Convert to real(c_float)
            CALL PGMOVE (PGx1, PGy1)
    ! Move the "pen" to the point with world
    ! coordinates (X,Y). No line is drawn.
            CALL PGDRAW (PGx2, PGy2)
    ! Draw a line from the current pen position to the point
    ! with world-coordinates (X,Y). The line is clipped at the edge of the
    ! current window. The new pen position is (X,Y) in world coordinates.
    return
  end subroutine gpln2d
  !---------------------------------------------------------------------
  subroutine gslnsz(size) ! called explicitly with real(c_float) args
    ! set line size (width), where 0. is the default.
    real size ! take the compiler default
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
  end subroutine gslnsz
  !---------------------------------------------------------------------
  subroutine gslnst(LS)
    ! sets line style: 1-solid, 2-dashed, 3-dotted, 4-dash-dotted, etc.
    implicit integer (i-n), real(c_double) (a-h,o-z)
    INTEGER  LS
            CALL PGSLS(LS)
    ! Set the line style attribute for subsequent plotting. This
    ! attribute affects line primitives only; it does not affect graph
    ! markers, text, or area fill.
    ! Five different line styles are available, with the following codes:
    ! 1 (full line), 2 (dashed), 3 (dot-dash-dot-dash), 4 (dotted),
    ! 5 (dash-dot-dot-dot). The default is 1 (normal full line).
    return
  end subroutine gslnst
  !---------------------------------------------------------------------
  subroutine  gscpvs(gl_x,gl_y) ! called explicitly with real(c_float) args
    ! set current position for text
    real(c_float) gl_x,gl_y
    X = gl_x  ! To PGPTXT(X,Y,ANGLE,FJUST,TEXT)
    Y = gl_y
    return
  end subroutine gscpvs
  !---------------------------------------------------------------------
  subroutine  gstxan(gl_angle) ! called explicitly with real(c_float) args
    ! set angle to plot text
    real(c_float) gl_angle
    ANGLE = gl_angle ! To PGPTXT(X,Y,ANGLE,FJUST,TEXT)
    return
  end subroutine gstxan
  !---------------------------------------------------------------------
  subroutine  gstxjf(just1,just2)
    ! set justification of string
    character*(*) just1 !can be 'left', 'right', or 'center'
    character*(*) just2 !can be 'top', 'bottom', or 'center'

    FJUST = 0.0
    if(just1.eq."left")   FJUST = 0.0
    if(just1.eq."center") FJUST = 0.5
    if(just1.eq."right")  FJUST = 1.0
    ! To PGPTXT(X,Y,ANGLE,FJUST,TEXT)
    return
  end subroutine gstxjf
  !---------------------------------------------------------------------
  subroutine  gstxno(size) ! called explicitly with real(c_float) args
    ! set number of characters per line, i.e., character width.
    ! Default: size=90.  (not sure)
    real(c_float) size
    real(c_float) PGsize
    PGsize = size ! Convert to real(c_float)
             CALL PGSCH(90./PGsize) ! set character size; default is 1.
    return
  end subroutine gstxno
  !---------------------------------------------------------------------
  subroutine  gptx2d(text)
    ! Plot Text
    character*(*) text
            call PGPTXT(X, Y, ANGLE, FJUST, text)
    ! Primitive routine for drawing text. The text may be drawn at any
    ! angle with the horizontal, and may be centered or left- or right-
    ! justified at a specified position.  Routine PGTEXT provides a
    ! simple interface to PGPTXT for horizontal strings. Text is drawn
    ! using the current values of attributes color-index, line-width,
    ! character-height, and character-font.  Text is NOT subject to
    ! clipping at the edge of the window.
    ! X      (input)  : world x-coordinate.
    ! Y      (input)  : world y-coordinate. The string is drawn with the
    !                   baseline of all the characters passing through
    !                   point (X,Y); the positioning of the string along
    !                   this line is controlled by argument FJUST.
    ! ANGLE  (input)  : angle, in degrees, that the baseline is to make
    !                   with the horizontal, increasing counter-clockwise
    !                   (0.0 is horizontal).
    ! FJUST  (input)  : controls horizontal justification of the string.
    !                   If FJUST = 0.0, the string will be left-justified
    !                   at the point (X,Y); if FJUST = 0.5, it will be
    !                   centered, and if FJUST = 1.0, it will be right
    !                   justified. [Other values of FJUST give other
    !                   justifications.]
    ! TEXT   (input)  : the character string to be plotted.
    return
  end subroutine gptx2d
  !---------------------------------------------------------------------

      subroutine pltends
      use param_mod
      use comm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!     Plot energy, density, parallel current and density conservation
!     constant vs. time., at a given s, distance along the magnetic field.
!     A cqlp.eq.'enabled' routine.
!
!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!

      REAL RILIN !-> For PGPLOT (text output positioning)
      dimension wk_nch(nonch)
!
      if (noplots.eq."enabled1") return
      if (pltend.eq."disabled") return
      dgts=1.e-8
      rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
      do 220 k=1,ngen
        if (pltend.eq."notplts") then
          goto 100
        elseif (pltend.eq."last") then
          if (n.lt.nstop) goto 100
        endif
!...
!mnt  Generate plot "endn"
!...
        call GXGLFR(0)
        call aminmx(pdens(1:nch(l_),k,l_),1,nch(l_), &
             1,emin,emax,kmin,kmax)
        if (abs(emin-emax).lt.emax*dgts) emax=emin+.001*abs(emin)
        if(emax.gt.0.) emax=emax*1.05 ! extend the upper range
        !write(*,*) 'pltends-1:', ptime(1,l_),ptime(nch(l_),l_),emin,emax
        call GSVP2D(.2,.8,.6,.9) !---------------> 1st subplot
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),emin, &
             emax*1.05)
        call GPCV2D(ptime(1:nch(l_),l_),pdens(1:nch(l_),k,l_),nch(l_)) ! density(time)
        CALL PGLAB('time (sec)',' ',' ')
        !call GSCPVS(.08,.75) ! set current position for text
        write(t_,10120) k
10120   format("density(s) of species",i2)
        RILIN=0.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('LV',RILIN,0.,0.,t_) ! Left-Vertical

        call aminmx(pengy(1:nch(l_),k,l_),1,nch(l_), &
                1,emin,emax,kmin,kmax)
        if (abs(emin-emax).lt.emax*dgts) emax=emin+.001*abs(emin)
        if(emax.gt.0.) emax=emax*1.05 ! extend the upper range
        !write(*,*) 'pltends-2:', ptime(1,l_),ptime(nch(l_),l_),emin,emax
        call GSVP2D(.2,.8,.2,.5) !---------------> 2nd subplot
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),emin, &
             emax*1.05)
        call GPCV2D(ptime(1:nch(l_),l_), pengy(1:nch(l),k,l_), nch(l_)) ! energy(time)
        CALL PGLAB('time (sec)',' ',' ')
        !call GSCPVS(.08,.35) ! set current position for text
        write(t_,10110) k
10110   format("energy(s) of species",i2)
        RILIN=0.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('LV',RILIN,0.,0.,t_) ! Left-Vertical

        write(t_,10140) denpar(k,ls_),enrgypa(k,ls_)
10140   format("local density(s) (/cm**2) = ",e15.6, &
          ";  energy(s) (kev) =",e15.6)
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

        write(t_,10150) n,timet,k
10150   format("time step (n) is",i5,5x,"time=",1pe14.6," secs", &
               "   Species k=",i2)
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

        write(t_,10152) rovera(lr_),rr
10152   format("r/a=",e14.6,5x,"radial position (R)=",e14.6," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

        write(t_,10153) sz(l_)
10153   format("parallel position (s) =",e14.6," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

!...
!mnt  Generate plot "curr"
!...
        call GXGLFR(0) ! new page
        call aminmx(pcurr(1:nch(l_),k,l_), &
             1,nch(l_),1,emin,emax,kmin,kmax)
        if (abs(emin-emax).lt.emax*dgts) emax=emin+.001*abs(emin)
        if(emax.gt.0.) emax=emax*1.05 ! extend the upper range
        !write(*,*) 'pltends-3:', ptime(1,l_),ptime(nch(l_),l_),emin,emax
        call GSVP2D(.2,.8,.6,.9)
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),emin,emax)
        call GPCV2D(ptime(1:nch(l_),l_),pcurr(1:nch(l_),k,l_),nch(l_)) ! current_dens(time)
        CALL PGLAB('time (sec)',' ',' ')
        curramp=currm(k,l_)/3.e9
        write(t_,10160) k,curramp
10160   format("Local in s current density",";", &
          "of species ",i2," =",e14.5,";", &
          "units are Amps/cm**2")
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom


 100    continue

        if (pltend.eq."tplts") goto 220
!...
!mnt  Generate plot "currv"
!...
        call GXGLFR(0) ! new page
        call aminmx(currv(1:jx,k,l_),1,jx, &
            1,fnmin,fnmax,kmin,kmax)
        if (abs(fnmin-fnmax).lt.fnmax*dgts) fnmax=fnmin+.001*abs(fnmin)
        !write(*,*) 'pltends-4:', x(1),x(jx),fnmin,fnmax
        call GSVP2D(.2,.8,.6,.9) !---------------> 1st subplot
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",x(1),x(jx),fnmin,fnmax)
        call GPCV2D(x,currv(1:jx,k,l_),jx)
        CALL PGLAB('v/vnorm','current: J(v)',' ') !

        call aminmx(currvs(1:jx,k),1,jx, &
           1,fnmin,fnmax,kmin,kmax)
        if (abs(fnmin-fnmax).lt.fnmax*dgts) fnmax=fnmin+.001*abs(fnmin)
        !write(*,*) 'pltends-5:', x(1),x(jx),fnmin,fnmax
        call GSVP2D(.2,.8,.2,.5) !---------------> 2nd subplot
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",x(1),x(jx),fnmin,fnmax)
        call GPCV2D(x,currvs(1:jx,k),jx)
        CALL PGLAB('v/vnorm','Integral from 0 to v of J(v)dv',' ') !

        write(t_,10183) currvs(jx,k)
10183   format("current =",e10.4,"Amps/Watt")
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

 220  continue ! k=1,ngen
!
!     Plot the density conservation diagnostic vs time
!
      if (pltend.eq."notplts") then
        if (n.lt.nstop) return
      elseif (pltend.eq."last") then
        if (n.lt.nstop) return
      endif
!...
!mnt  Generate plot "consn(l_)"
!...
      call GXGLFR(0) ! new page
      call aminmx(consnp(1:nch(l_),l_),1,nch(l_), &
        1,emin,emax,kmin,kmax)
      !write(*,*) 'pltends-6:', ptime(1,l_),ptime(nch(l_),l_),emin,emax
      call GSVP2D(.2,.8,.6,.9)
      CALL PGSCH(1.) ! set character size; default is 1.
      call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),emin,emax)
      call GPCV2D(ptime(1:nch(l_),l_),consnp(1:nch(l_),l_),nch(l_))
      CALL PGLAB('time (sec)','consn(l_) conservation diag',' ')

      write(t_,10250) consn(l_)
10250 format("consn(l_)=",1pe12.4)
      RILIN=5.
      CALL PGSCH(0.8) ! set character size; default is 1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)

      write(t_,10251)
10251 format("Perfect conservation should yield  machine accuracy,")
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)

      write(t_,10252)
10252 format("or about 1.e-14:")
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)

      write(t_,10150) n,timet
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)

      write(t_,10152) rovera(lr_),rr
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)

      write(t_,10153) sz(l_)
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,0.,0.,t_)


      do 280 k=1,ngen
        call GXGLFR(0) ! new page(s)
!$$$    Possibly write t_ greater than present dimension character*512:
!$$$        write(t_,10260) k,(sgaint(i,k,l_),i=1,8)
!$$$        call gptx2d(t_)
!$$$10260   format("total gain (+) or loss (-) to date for species",i3,";",
!$$$     &    "(in particles/cm**2 - units of line density)",";",
!$$$     &    "forcing nonnegative f x-sweep (implct=disabled)",e12.5,";",
!$$$     &    "forcing nonnegative f y-sweep (implct=disabled)",e12.5,";",
!$$$     &    "due to particle source=",e12.5,";","due to runaway=",e12.5,
!$$$     +    ";",
!$$$     &    "due to lossmode(k)=",e12.5,";","due to torloss(k)=",e12.5,
!$$$     +    ";",
!$$$     &    "due to fusion losses=",e12.5,";",
!$$$     &    "forcing nonnegative distribution (implct=enabled)=",e12.4,
!$$$     +    ";","$")
 280  continue
        write(t_,10150) n,timet
        RILIN=0.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('T',RILIN,0.,0.,t_) ! Top

        write(t_,10152) rovera(lr_),rr
        RILIN=RILIN-1.
        CALL PGMTXT('T',RILIN,0.,0.,t_) ! Top

        write(t_,10153) sz(l_)
        RILIN=RILIN-1.
        CALL PGMTXT('T',RILIN,0.,0.,t_) ! Top

      return
      end

  !From pltfluxs
      subroutine pltfluxs
      use param_mod
      use comm_mod
      use advnce_mod
      use pltdf_mod, only : JXQ
      use r8subs_mod, only : luf
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)

!...................................................................
!     Plots velocity (momentum-per-mass)  fluxes versus u,
!     for selected values of theta.
!     Do combined fluxes, and individual fluxes, per pltflux1 vector.
!...................................................................
!
!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!

      save


      REAL RILIN !-> For PGPLOT (text output positioning)
      real(c_double) wkd(jx)
      CHARACTER*64 TX_

!     pltflux1(1)=1. ==> sum of fluxes
!              2         collisions
!              3         parallel electric field
!              4         rf
!              5         synchrotron
!              6         Bremssstrahlung (+ phenomenological energy loss)
!              7         Krook operator slowing down

!...................................................................
!     Outer loop is over the various plots, as enabled by pltflux1
!...................................................................

      do 500 kk=1,7

      if (pltflux1(kk).ne.1.) go to 500

      do 400 k=1,ngen

!     Do velocity and theta fluxes separately:
      do 410 kkk=1,2
!
        call bcast(da,zero,iyjxp1)
        call bcast(db,zero,iyjxp1)
        call bcast(dc,zero,iyjxp1)
        call bcast(dd,zero,iyp1jx)
        call bcast(de,zero,iyp1jx)
        call bcast(df,zero,iyp1jx)

        go to (10,20,30,40,50,60,70),kk
 10     call coefstup(k)
        write(t_,910) k
        go to 100
 20     call coeffpad(k)
        write(t_,920) k
        go to 100
 30     if (elecfld(lr_)  .lt. 1.e-09) go to 500
        call coefefad(k)
        write(t_,930) k
        go to 100
 40     continue
        if (urfmod.eq."disabled" .and. vlfmod.eq."disabled" &
             .and. vlhmod.eq."disabled") go to 500
        xrf=0.
        if (n .lt. nonrf(k) .or. n .ge. noffrf(k)) go to 41

        call coefrfad(k,xrf)
          write(*,'(a,2i6,e12.2)') &
          'pltfluxs->coefrfad: n,lr_,sum(dbb)=', &
          n,lr_,sum(dbb)

 41     continue
        if (xrf.gt.0.) then
          write(t_,940) k
        endif
        go to 100

 50     continue
        if (syncrad.eq."disabled") go to 500
        call coefsyad(k)
        write(t_,950) k
        go to 100

 60     continue
        if (bremsrad.eq."disabled" &
             .and. torloss(k).eq."disabled") go to 500
        call coefegad(k)
        write(t_,960) k
        go to 100

 70     call coefload(k)
        write(t_,970) k



 100    continue

!...................................................................
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
!     first the velocity flux..
!...................................................................

        call coefmidv(da,1)
        call coefmidv(db,2)
        call coefmidv(dc,3)

!...................................................................
!     the theta flux..
!...................................................................

        call coefmidt(dd,1)
        call coefmidt(de,2)
        call coefmidt(df,3)


!...................................................................
!     Fluxes
!...................................................................


        if (kkk.eq.1) then
        do 200 j=1,jx
           do 199 i=1,iy
              temp1(i,j)=gfi(i,j,k)
 199       continue
 200    continue
        else
         do 202 j=1,jx
           do 201 i=1,iy
              temp1(i,j)=hfi(i,j)
 201       continue
 202    continue
        endif



!...
!mnt  Plot fluxes as a function of velocity for various
!mnt  angles.
!...
      iii=(itl+iyh)/2
      trpmd=y(iii,l_)
      do 210 i=1,iyh
        if (y(i,l_) .ge. trpmd) goto 220
 210   continue
 220   continue
      midtrp=i
      imsh(1)=1
      imsh(2)=itl
      imsh(3)=midtrp
      imsh(4)=iyh
      imsh(5)=iy

!        if (tandem.eq."enabled" .and. fmass(k).gt.1.e-27) then
!          jxq=jlwr
!          xmaxq=xlwr
!          iyjxq=iy*jlwr
!        else
!          jxq=jx
!          xmaxq=xmax
!          iyjxq=iyjx
!        endif

        if (pltlim.eq."disabled") then
           jxq=jx
           xmaxq=x(jxq)
           iyjxq=iyjx
           tx_='u/vnorm'
           do j=1,jxq
              tam1(j)=x(j)
           enddo
        endif

        if (tandem.eq."enabled" .and. fmass(k).gt.1.e-27) then
           jxq=jlwr
           xmaxq=xlwr
           iyjxq=iy*jlwr
           pltlim='x' ! YuP: reset?
           pltlimm=xmaxq
           tx_='u/vnorm'
           do j=1,jxq
              tam1(j)=x(j)
           enddo
!       If pltlim .ne. "disabled", plot versus
!       'x', 'u/c', or 'energy', up to maximum pltlimm.
        elseif (pltlim.eq.'x') then
           jxq=min(luf(pltlimm,x,jx),jx)
           xmaxq=x(jxq)
           iyjxq=iy*jxq
           tx_='u/vnorm'
           do j=1,jxq
              tam1(j)=x(j)
           enddo
        elseif (pltlim.eq.'u/c') then
           jxq=min(luf(pltlimm,uoc,jx),jx)
           xmaxq=uoc(jxq)
           iyjxq=iy*jxq
           tx_='u/c'
           do j=1,jxq
              tam1(j)=uoc(j)
           enddo
        elseif (pltlim.eq.'energy') then
           pltlimmm=pltlimm
           wkd(1:jx)=enerkev(1:jx,k)
           jxq=min(luf(pltlimmm,wkd,jx),jx)
           xmaxq=enerkev(jxq,k) !YuP[2018-01-08] added 2nd index (k)
           iyjxq=iy*jxq
           tx_='Energy (keV)'
           do j=1,jxq
              tam1(j)=enerkev(j,k) !YuP[2018-01-08] added 2nd index (k)
           enddo
        endif


        call bcast(tam4,zero,jxq)
        do 240 i=1,iy
          do 230 j=1,jxq
            tam4(j)=tam4(j)+temp1(i,j)*cynt2(i,l_)/twoint(l_)
 230       continue
 240     continue
        call aminmx(tam4,1,jxq,1,emin,emax,kmin,kmax)
        do 260 iu=1,5
          i=imsh(iu)
          do 250 j=1,jxq
            tam2(j)=temp1(i,j)
 250       continue
          call aminmx(tam2,1,jxq,1,fmin,fmax,kmin,kmax)
          if (fmin .lt. emin) emin=fmin
          if (fmax .gt. emax) emax=fmax
 260     continue
        emax=emax*1.03
        emin=emax/1.e+12

        call GXGLFR(0) ! new page
        call GSVP2D(.2,.8,.3,.9)
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlog$",tam1(1),tam1(jxq),emin,emax)
        !-YuP:   call GSCVLB(1)
        !-YuP:   call GSCVFT(0.)
        do 310 iu=1,6
          if (iu .eq. 6) then
            do 270 j=1,jxq
              tam2(j)=tam4(j)
 270         continue
            text(1)="avg"
            goto 290
          endif
          do 280 j=1,jxq
            tam2(j)=f(imsh(iu),j,k,l_)
 280       continue
          if (iu .eq. 1) text(1)="pll"
          if (iu .eq. 2) text(1)="trp/ps"
          if (iu .eq. 3) text(1)="midtrp"
          if (iu .eq. 4) text(1)="perp"
          if (iu .eq. 5) text(1)="pll-pi"
 290       continue
          !-YuP:   call GSCVTX(loc(text))
          do 300 j=1,jxq
            if(tam2(j).le.emin) tam2(j)=emin
 300      continue
          call GPCV2D(tam1,tam2,jxq)
          xu=float(iu)
          !-YuP:   call GSCVFT(xu/6.)
 310    continue
        !-YuP:   call GSCVLB(0)

       CALL PGLAB(tx_,' ',' ') ! YuP/added: horizontal axis label

!     Write previously set title
        RILIN=1.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('T',RILIN,0.,0.,t_) ! Top

        rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
        write(t_,10020) k
10020 format("Flux vs. velocity for some angles; species number = ",i3)
        RILIN=3.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

        write(t_,10010) n,timet
10010 format("time step (n) is",i5,5x,"time=",e14.6," secs")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10011) rovera(lr_),rr
10011 format("r/a=",e14.6,5x,"radial position (R) =",e14.6," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10030)
10030 format("pll    ---- theta = 0 radians")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10031)
10031 format("pll-pi ---- theta = pi radians")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10032) y(itl,l_)
10032 format("trp/ps ---- theta = ",e13.5," radians")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10033) trpmd
10033 format("midtrp ---- theta = ",e13.5," radians")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10034) y(iyh,l_)
10034 format("perp   ---- theta = ",e13.5," radians")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10035)
10035 format("avg    ---- theta averaged over pi radians")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        CALL PGSCH(1.0) ! recover default 1.0 fontsize

 410    continue
 400    continue ! k species
 500  continue ! skipping handle

      CALL PGSCH(1.0) ! recover default 1.0 fontsize

 910  format("species no.",i2,5x,"combined velocity space fluxes")
 920  format("species no.",i2,5x,"Fokker-Planck velocity space flux")
 930  format("species no.",i2,5x,"electric field velocity space flux")
 940  format("species no.",i2,5x,"RF velocity space flux")
 950  format("species no.",i2,5x,"synchrotron velocity space flux")
 960  format("species no.",i2,5x,"Brems+phenom velocity space flux")
 970  format("species no.",i2,5x,"Krook velocity space flux")

      return
      end

  !from pltprppr

  subroutine pltprppr
      use param_mod
      use comm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!...
!mnt  This routine plots the parallel distribution function.
!...
!
!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!
!     Need to re-work this a bit, for plot of pitch angle averaged
!     reduced distributions (knockon=.ne."disabled" case).
!

      REAL RILIN !-> For PGPLOT (text output positioning)

      character*8 target
      if (noplots.eq."enabled1") return

!     Return, if using a theta average distribution
!     (which may conflict with call fle_fsa or fle_pol, below).
      if (flemodel.ne."pol"  .and. flemodel.ne."fsa") return

      do 20 k=1,ngen
        if (tandem.eq."enabled" .and. k.eq.kionn) then
          xll=-xlwr
          xlu=xlwr
          xpl=0.
          xpu=xlwr
          target="ionmesh"
        else
          xll=-xmax
          xlu=xmax
          xpl=0.
          xpu=xmax
          if (xprpmax.ne.1.) xpu=xprpmax
          target="mainmesh"
        endif
!       Obtain equatorial plane pitch-avg'd distribution for ko model
!       case.  (Else setup interfere's with subroutine sourceko).
        if (knockon.ne."disabled") then
!           call fle("setup",0)
!           call fle("calc",1)
        elseif (lrz.eq.1) then
           call fle_pol("setup",0)
           call fle_pol("calc",1)
        else
           call fle_fsa("setup")
           call fle_fsa("calc")
        endif

!     plot (fll(+xpar)-fll(-xpar)), and xpar*(fll(+xpar)-fll(-xpar))
        jpxyh=(jfl+1)/2
        jpxyhm=jpxyh-1
        tem1(1)=0.0
        tem2(1)=0.0
        do 30  jp=1,jpxyhm
          tem1(jp)=fl(jpxyh+jp-1)-fl(jpxyh-jp)
          tem2(jp)=xlm(jpxyh+jp-1)*tem1(jp)
 30     continue
        call aminmx(tem1,1,jpxyhm,1,fmin1,fmax1,kmin,kmax)
        call aminmx(tem2,1,jpxyhm,1,fmin2,fmax2,kmin,kmax)
!990131        fmin=amin1(fmin1,fmin2)
!990131        fmax=amax1(fmax1,fmax2)
        fmin=min(fmin1,fmin2)
        fmax=max(fmax1,fmax2)

        call GXGLFR(0) ! new page(s)
        call GSVP2D(.2,.8,.25,.95) ! (XLEFT, XRIGHT, YBOT, YTOP)
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",0.d0,xlu,fmin,fmax)
        call GPCV2D(xlm(jpxyh:jpxyh+jpxyh), tem1, jpxyhm)
        call GPCV2D(xlm(jpxyh:jpxyh+jpxyh), tem2, jpxyhm)

        write(t_,10011) k
10011 format("asymmetric cmpt of f_par, and xpar*cmpt, species:",1x,i5)
        RILIN=5.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10012)
10012 format("(f_par normed so int(-1,+1)=equatorial ne)")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        call GXGLFR(0) ! new page(s)
        call aminmx(fl(1:2*jpxyhm),1,2*jpxyhm,1,fmin,fmax,kmin,kmax)
        fmin=1.d-08*fmax
        call GSVP2D(.2,.8,.25,.95)
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlog$",xll,xlu,fmin,fmax)
        do 10 jj=1,2*jpxyhm
          if (fl(jj) .lt. fmin ) fl(jj)=fmin
 10     continue
      !-YuP:   call GSVTCL("on$")
        call GPCV2D(xlm,fl,2*jpxyhm)
      !-YuP:   call GSVTCL("off$")

        write(t_,10013) k
10013 format("parallel distribution function for species:",1x,i5)
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10014)
10014 format("(normed so int(-1,+1)=equatorial ne)")
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

        write(t_,10031) rovera(lr_),rr
10031 format("r/a=",e14.6,5x,"radial position (R) =",e14.6," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        CALL PGSCH(1.0) ! recover default 1.0 fontsize

 20   continue

      return
      end

   !from pltrstv

   subroutine pltrstv
      use param_mod
      use comm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!     Plot electron resistivity and related quantities.
!

!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!

      REAL RILIN !-> For PGPLOT (text output positioning)

!
      if (noplots.eq."enabled1") return

      if (kelecg .eq. 0 .or. abs(elecfld(lr_)) .lt. 1.e-10) go to 190

      call GXGLFR(0)
      call aminmx(sptzrp(2:nch(l_),lmdpln_),1,nch(l_)-1, &
           1,emin,emax,kmin,kmax)
      call aminmx(restp(2:nch(l_),lr_),1,nch(l_)-1, &
           1,fmin,fmax,kmin,kmax)
      if (fmin .lt. emin) emin=fmin
      if (fmax .gt. emax) emax=fmax
      call GSVP2D(.2,.8,.6,.9) !---------------> 1st subplot
      CALL PGSCH(1.) ! set character size; default is 1.
      call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_), &
       .95d0*emin,1.05d0*emax)
      text(1)="spitzer"
      !-YuP:   call GSCVLB(1)
      !-YuP:   call GSCVTX(loc(text))

      call GPCV2D(ptime(2:nch(l_)-1,l_), sptzrp(2:nch(l_)-1,lmdpln_), &
       nch(l_)-1)
      text(1)="rstvty"
      call GPCV2D(ptime(2:nch(l_)-1,l_), restp(2:nch(l_)-1,lr_), &
       nch(l_)-1)
      !-YuP:   call GSCVLB(0)
      write(t_,170)
 170  format("upper graph - flux avg and spitzer resistivities")
      RILIN=2.
      CALL PGSCH(0.9) ! set character size; default is 1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      write(t_,171)
 171  format("lower graph - ratio of resist to spitzer or neo resist")
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      illeff=lr_
      if (cqlpmod .eq. "enabled") illeff=ls_
      call aminmx(rovsp(2:nch(l_),illeff),1,nch(l_)-1, &
            1,emin,emax,kmin,kmax)

      call GSVP2D(.2,.8,.2,.5) !---------------> 2nd subplot
      call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_), &
       .95d0*emin,1.05d0*emax)
      CALL PGLAB('time (sec)',' ',' ')
      call GPCV2D(ptime(2:nch(l_)-1,l_), rovsp(2:nch(l_)-1,illeff), &
       nch(l_)-1)

      if (efswtchn.eq."neo_hh") then
         write(t_,10164)
10164 format("(efswtchn=neo_hh)")
         RILIN=RILIN-1.
         CALL PGMTXT('T',RILIN,0.,0.,t_)
      endif



      call GXGLFR(0) ! new page
      call GSVP2D(0.2,0.8,0.2,0.7) ! (XLEFT, XRIGHT, YBOT, YTOP)
      write(t_,70)
 70   format("--calculated resistivity and other related quantities--")
      RILIN=2.
      CALL PGSCH(0.9) ! set character size; default is 1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,80) sptzr(l_)
 80   format("spitzer restvty= ",1pe14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,81) resist
 81   format("toroidal restvty=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,82) resistn
 82   format("neoclass restvty=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,83) rovs(lr_)
 83   format("ratio of resistivities=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,84) rovsf
 84   format("O(epsilon**.5) expansion for resistivity ratio=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,85) elecr(lr_)
 85   format("E-Dreicer=",e14.5,"vlts/cm")
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,86) rovsloc(l_)
 86   format("local resistivity over spitzer=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,90) rovscf
 90   format("Small eps(lr_) fla for resist ratio (connor)=",e16.6)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,91) rovsc(lr_)
 91   format("gen. epsilon fla for resist ratio (connor)=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,92) xconn
 92   format("^\i(connor)=",e16.6)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,93) elecfld(lr_)
 93   format("electric field=",e16.6,"vlts/cm")
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,94) eovedd
 94   format("E-parallel/E-Dreicer=",e16.6)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,95) tauee(lr_)
 95   format("tauee(lr_)=",e16.6,"secs")
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
!
!     add some other relevant quantities
!
      write(t_,96)(btor0(lr_)/bmod0(lr_))**2
 96   format("b_phi/b at outer midpplane=",1pe14.4)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,97) onovrp(2,lr_)*rpcon(lr_)**2
 97   format("R(z=0)**2 * <1/R**2>      =",  e14.4)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,98) psiavg(2,lr_)
 98   format("<(B(z)/B(0))**2>          =",  e14.4)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)


 190  continue
      return
      end

  ! pltstrml
      subroutine pltstrml
      use param_mod
      use comm_mod
      use advnce_mod
      use pltdf_mod, only : cont, tempcntr, nconta
      use pltdf_mod, only : wx, wy, IIY, JXQ
      use r8subs_mod, only : luf, dcopy
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
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
      real(c_float) RTAB1(iy),RTAB2(iy) ! local

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
!        xmaxq=pltlimmm
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
        call bcast(temp5(0:iyjx2,0),zero,iyjx2)
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

  ! from pltfofvv
  subroutine pltfofvv
      use param_mod
      use comm_mod
      use r8subs_mod, only : dcopy
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!**SCC10/7/94
!   prppru calculates the pitch angle integrated distribution function
!**********************
!
!
!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!

      REAL RILIN !-> For PGPLOT (text output positioning)

      character*8 target
      do 20 k=1,ngen
        call GXGLFR(0) ! new page(s)
        call dcopy(iyjx2,f(0:iyjx2-1,0,k,l_),1,temp3(0:iyjx2-1,0),1)
        if (tandem.eq."enabled" .and. k.eq.kionn) then
          target="ionmesh"
          jxq=jlwr
        else
          target="mainmesh"
          jxq=jx
        endif

!       Obtain integrated distribution in tam1
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


!
!
      subroutine fofv(target,action)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      character*(*) target,action
      save

!...............................................................
!mnt  This routine takes data stored in temp3 in (y,x) coordinates
!     and integrates over theta0 to get tam1 which is, in this case,
!     the isotropized distribution; Int(tam1*2*pi*vnorm**3*x**2*dx) is
!     the density.
!     SCChiu  10/7/94
!      xul,xuu: lower and upper values of x (normalized u)
!      ytl,ytu:  "     "   "      "    "  y (theta0)
!...............................................................


      logical trnsfm
      trnsfm=(target.eq."velocity".and. relativ .ne. "disabled")

        do jp=1,jx
          tam2(jp)=0.
          do ip=1,iy-1
            ip1=ip+1
            tam2(jp)=tam2(jp)+cynt2(ip,l_)*(temp3(ip,jp) &
                     +temp3(ip1,jp))
          enddo
          tam1(jp)=tam2(jp)/twopi
        enddo

      return
      end

  ! from souplt
      subroutine souplt
      use param_mod
      use comm_mod
      use pltdf_mod, only : cont, tempcntr, nconta, JXQ
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     contour plots the ion source.
!..................................................................


      REAL RILIN
      REAL RTAM1(jx),RTAM2(jx)
      REAL REMAX,REMIN

      if (pltso.eq."disabled") return


      do 10 k=1,ngen

         temp1=0.d0 ! initialize for each k species

         if (xlncur(k,lr_).lt.1.e-10) goto 10
!BH171231         if(frmodp.eq.'enabled')then ! NBI source
           call dcopy(iyjx2,source(0:iyjx2-1,0,k,indxlr_),1, &
              temp1(0:iyjx2-1,0),1)
!BH171231         endif
         write(t_,550) k
 550     format(1x,"Species ",i2, &
              " Source Function (units: dist. f/sec)")
         CALL PGPAGE
         itype=3 ! means: plots are made for source
         call pltcont(k,1,t_,itype) ! for source
!BH171231         crnt_nbi=xlncur(k,lr_)*zmaxpsii(lr_) ! [ptcl/sec/cm^3]
         crnt=xlncur(k,lr_)*zmaxpsii(lr_) ! [ptcl/sec/cm^3]

!BH171231         write(t_,540) crnt_nbi
!BH171231 540     format("NBI source rate=",1pe11.4," ptcls/cc/sec")
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

!     Plot pitch angle integrated source:
      call pltsofvv

!     Plot the speed-integrated source,
!     as a function of pitch angle theta0 at the midplane
!yup      call pltso_theta  !YuP[06-2016]

      return
      end
!
!=======================================================================
      subroutine pltsofvv
      use param_mod
      use comm_mod
      use r8subs_mod, only : luf
      use r8subs_mod, only : dcopy
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!     Calculates the pitch angle integrated source
!**********************
!

      REAL RILIN, RXMAXQ
      REAL RTAM1(jx),RTAM2(jx), wk_tam(jx)
      REAL REMAX,REMIN
      real(c_double) wkd(jx)

      character*8 target
      character*8 tx_


      do 20 k=1,ngen

        temp3=0.d0 ! initialize for each k species
        if (xlncur(k,lr_).lt.1.e-10) goto 20
!BH171231        if(frmodp.eq.'enabled')then ! NBI source
          call dcopy(iyjx2,source(0:iyjx2-1,0,k,indxlr_),1, &
             temp3(0:iyjx2-1,0),1) ! temp3
!BH171231        endif

!-----YuP[2018-01-08] revised to match cqlinput_help:
!**    pltlim= "disabled",  plots versus x (u/vnorm) from
!**                    x=0. to 1. (default:pltlim="disabled",pltlimm=1.)
!**            "x",    plot 1d and 2d plots versus x
!**                    from 0. to pltlimm.
!**            "u/c",  plot 1d and 2d plots versus u/c
!**                    from 0. to pltlimm.
!**            "energy", plot 1d and 2d plots verus energy (kev)
!**                    from 0. to pltlimm (kev).
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
!-----YuP[2018-01-02] done



!       Obtain integrated distribution in tam1
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


!
!=======================================================================
      subroutine pltso_theta
      use param_mod
      use comm_mod
      use r8subs_mod, only : dcopy
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!     YuP[06-2016]
!     Calculates and plots the speed-integrated source,
!     as a function of pitch angle theta0 at the midplane
!     Called for every flux surface lr_ (lr_ is stored in comm.h)
!**********************
!
!MPIINSERT_INCLUDE

      REAL RILIN
      REAL RTAM1(iy),RTAM2(iy)
      REAL REMAX,REMIN
      real(c_double) wk_so(iy)
      character*8 vert_scale ! 'log10' or 'lin'

      vert_scale='log10' !'lin'
      do i=1,iy
        RTAM1(i)=y(i,lr_)*180.0/pi  ! horizontal.axis: theta0 (degree)
      enddo

      do 20 k=1,ngen ! sources for each general sp. are plotted

        temp3=0.d0 ! initialize for each k species (i,j)
        if (xlncur(k,lr_).lt.1.e-10) goto 20

        if(frmodp.eq.'enabled')then ! NBI source
          call dcopy(iyjx2,source(0:iyjx2-1,0,k,indxlr_),1, &
                temp3(0:iyjx2-1,0),1) ! temp3
        endif


!       Obtain integrated distribution into wk_so(1:iy)
        do i=1,iy
           wk_so(i)=0. ! initialize
        do j=1,jx
           wk_so(i)= wk_so(i)+ temp3(i,j)*cint2(j)
           !cint2= x**2 *dx  (and remember that temp3 includes vnorm^3)
           !cynt2= 2pi*sin(theta0)*dtheta0
        enddo
        WRITE(*,'(a,2i5,2e13.5)') &
         'pltso_theta: lr_,i,y(i,lr_)*180.0/pi,wk_so(i)=', &
              lr_,i, y(i,lr_)*180.0/pi, wk_so(i)
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
!


end module pltmain_mod
