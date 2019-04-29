module impavnc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use coefmidt_mod, only : coefmidt
  use coefmidv_mod, only : coefmidv
  use coefstup_mod, only : coefstup
  use coefwti_mod, only : coefwti
  use coefwtj_mod, only : coefwtj
  use diagimpd_mod, only : diagimpd
  use esefld_mod, only : fluxpar
  use impchk_mod, only : impchk
  use impnorm_mod, only : impnorm
  use r8subs_mod, only : dcopy
  use r8subs_mod, only : dgbtrf
  use r8subs_mod, only : dgbtrs
  use tdtrvsou_mod, only : tdtrvsou

  !---END USE

!
!

contains

      subroutine impavnc
      use param_mod
      use comm_mod
      use advnce_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.................................................................
!     This routine advances the Fokker-Planck equation one time step
!     for the implicit case (implct=enabled)
!     (IMPAVNC calls IMPCHK which differentiates the solution and
!     compares it with the r.h.s of the advanced F.P. equation to
!     check for accuracy)
!.................................................................
!BH170802:  It uses a different treatment of the f(v=0) BC, and is no
!BH170802:  longer called (for many years).


      parameter(md1abd=3*(iya+3)+1)

      character*8 alloc, factor, dalloc
      dimension ipivot(iyjxa), zmatcont(12), janelt(12)
      character*1 transpose

!.................................................................
!     matrix declared as pointer
!     need inewjmax in common to use with pointer
!.................................................................
      dimension abd(:,:)
      pointer abd
      common/gauss_/ inewjmax,abd,ialloc

      save inewmax, mlmax, mumax

!MPIINSERT_INCLUDE
!.......................................................................
!     Indicate progress of code
!.......................................................................

      write(*,*)'impavnc: n,lr_',n,lr_

!.......................................................................
!     nadjoint is set = 0, to remove effect of Olivier Sauter's adjoint
!        operator changes for this subroutine.
!.......................................................................

       nadjoint=0

!.......................................................................
!     Need to install Sauter's 030509 bug fix.  See impanvce0.f
!     BH050804:   Bug fix only applies to lbdry0="enabled" case.
!.......................................................................

!       write(*,*) 'impavnc  NEED TO INSTALL SAUTERS 030509 BUGFIX'
!       STOP

!.......................................................................
!     ifirst is a counter for eseswtch.eq."enabled" case
!.......................................................................

       ifirst=1

!.......................................................................
!     define effective number of unknowns in the theta direction and
!     resulting band widths
!     see below for comments on inew, etc.
!.......................................................................

      inew=iy/2 + itl - 1
      if (symtrap .ne. "enabled") inew=iy
      if (nadjoint.eq.1 .and. symtrap.eq."enabled") inew=2*itl-1

!MPIREPLACE_IFFIRST
      if (l_.eq.1 .and. n.eq.1) then
        inewmax=inew
        if (meshy.eq."fixed_mu".or. manymat.eq."enabled") then
          do 15 ll=2,lrors
            inew1=iy_(ll)/2+itl_(ll)-1
            if (symtrap .ne. "enabled") inew1=iy_(ll)
            if (nadjoint.eq.1 .and. symtrap.eq."enabled") &
              inew1=2*itl_(ll)-1
            if (inew1.gt.inewmax) inewmax=inew1
 15       continue
        endif
        mlmax=inewmax+2
        mumax=inewmax+2
        inewjmax=inewmax*jx
      endif

      ml=inew+2
      mu=inew+2
      if (symtrap.ne."enabled" .and. cqlpmod.ne."enabled") then
!     itu+1 contributes to the equation for itl
        ml=inew+3
        mu=inew+3
        mlmax=inewmax+3
        mumax=inewmax+3
      endif

      ibandto=ml+mu+1
      inewjx=inew*jx
      navnc=navnc+1
      iyjx9=iyjx*9

!.................................................................
!     Copy current distribution into f_..
!.................................................................

!      call scopy(iyjx2a*ngen,f(0,0,1,l_),1,f_(0,0,1,l_),1)
      call dcopy(iyjx2a*ngen,f(0,0,1,l_),1,f_(0,0,1,l_),1)

!..................................................................
!     Logic for allocation..
!..................................................................

      alloc="disabled"
!     allocate once for all time-step
!MPIREPLACE_ALLOC
      if (l_.eq.1 .and. (nflux.eq.1.or.irstart.eq.1) .and. n.eq.1) &
        then
         alloc="enabled"
         ialloc=1
      endif
!     Special case for multiple calls at n=nstop
      if (n.eq.nstop.and.ialloc.eq.0) alloc="enabled"

!..................................................................
!     Logic for deallocation..
!..................................................................

      dalloc="disabled"
!     deallocate once at last time-step
!MPIREPLACE_DEALLOC
      if (l_.eq.lrors .and. n.ge.nstop .and. nflux.ge.nrstrt) then
         dalloc="enabled"
         ialloc=0
      endif

!..................................................................
!     Allocate space.
!..................................................................

      lenabd=md1abd*inewjmax
      if (alloc.eq."enabled") then
!     will preset the values to zero in loop k=1,ngen

!_cray   call hpalloc(abdptr,lenabd,ierr,1)
!_pc     Using unix library libU77.a
!_pc     abdptr=malloc(8*lenabd)
!_pc     if (abdptr.eq.0) stop 'failed malloc in impavnc0'
         allocate(abd(md1abd,inewjmax),STAT=istat)
         call bcast(abd,zero,size(abd))


      endif

!.................................................................
!     loop over all time advanced (general) species..
!.................................................................

      do 600 k=1,ngen

!..................................................................
!     if (colmodl.eq.4) f(i,j,ngen,l_) contains a drifting Maxwellian
!     whose momentum corresponds to that of the general electron
!     species contained in f(i,j,kelecg,l_);
!     the purpose of this species is to cool the general electrons
!     without absorbing net momentum;
!     the species is treated as a background or fixed field and is
!     updated only with respect to its momentum if the general
!     electron species momentum changes during the time step.
!..................................................................

        if (colmodl.eq.4 .and. k.eq.ngen) goto 600

!.................................................................
!     normalized time step..
!.................................................................

        rbgn=1./dtreff

!.................................................................
!     The routine coefstup(k) does the following:
!     Preset all coefficients to zero. da, db, and dc are the three
!     coefficients associated with the velocity (v) flux while
!     dd, de, df are associated with the theta flux. The contribution
!     from the Krook operator (cah(i,j)*f(i,j,l_)) is stored in cah.
!     Note cah does not contribute to da through df. These six
!     coefficients incorporate only wave-particle, collisional and
!     d.c. electric field effects. The ion source is represented
!     through so(i,j).
!     The routine coefstup(k) adds into the above arrays all the
!     relevant physics.
!.................................................................

        call coefstup(k)

!.................................................................
!     Determine whether matrix factorization and related space
!     allocation and de-allocation are to be done.
!     To save memory, if ngen (the number of general species) is
!     greater than 1, factorization is done every time step.
!.................................................................

        impchg=impcoef+imprf+impelec+impcah+irstart

!..................................................................
!     Logic for factorization..
!..................................................................

        factor="disabled"
        if (ngen.gt.1 .or. nflux.eq.1) then
          factor="enabled"
        endif
        if (impchg.gt.0) factor="enabled"

!.................................................................
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

!     first the velocity flux..
!.................................................................

        call coefmidv(da,temp1,1,vptb(1,lr_))
        call coefmidv(db,temp1,2,vptb(1,lr_))
        call coefmidv(dc,temp1,3,vptb(1,lr_))

!.................................................................
!     the theta flux..
!.................................................................

        call coefmidt(dd,temp1,1)
        call coefmidt(de,temp1,2)
        call coefmidt(df,temp1,3)

!.................................................................
!     The differencing is done ala Chang and Cooper as generalized
!     by Karney to 2-D. This requires defining a couple of arrays
!     di(i,j,l_) and dj(i,j,l_) which are used to define f at mid mesh
!     points: to wit f(i,j+1/2,l_)=f(i,j+1,l_)*(1.-dj(i,j,l_))+ f(i,j,l_
!     dj(i,j,l_). Similarly for theta.
!     The routines coefwtj and coefwti are used to
!     determine the dj and di arrays.
!.................................................................
        call coefwtj(k)
        call coefwti(k)

!.................................................................
!     Loop-back point for eseswtch.eq.enabled calculation
!.................................................................

 10     continue

!.................................................................
!     Proceed to call the Gaussian elimination routine. i indexes
!     theta, j indexes speed and ieq keeps a running total. It
!     is necessary to separate out the various boundary equations
!     from the the inner 9-point equation. Boundary equations can
!     involve from 1 to 12 points (or unknowns) depending on location
!     in velocity space.

!     ieq will keep track of which equation is being considered
!     below. This defines the row of the matrix being inverted. The
!     theta mesh points (index i) are incremented in the inner
!     loop.
!.................................................................

        ieq=0

!.................................................................
!     abd is the band matrix in a format expected by the sgbtrf routine.
!     That is, if a(i,j) is the matrix to invert, we have:
!     abd(iband+i-j,j) = a(i,j), with iband = ml+mu+1=total bandwidth
!     As if each column is moved up such that each diagonal element is
!     at abd(iband,j).  Note:  i,j in this comment designate the row
!     and column, respectively, of the general matrix a(,) to be
!     inverted, not the velocity space locations.
!
!     if factor .eq. "enabled" then we have to recompute the entire
!     matrix in preparation for factorization. Otherwise all
!     that has to be done is to recompute the r.h.s.  - The
!     code has to be called with factor .eq. "enabled" before it
!     can be called with factor="disabled".
!
!     See below for comments on unknowns numbering
!.................................................................

        if (factor.eq."disabled") then
          ipassbnd=iyh+1-itl
          if (symtrap.ne."enabled" .or. nadjoint.eq.1) ipassbnd=0
          do 2002 j=1,jx
            ihalf=0
            icntsww=1
            if (symtrap .ne. "enabled") icntsww=itl-iyh
            if (nadjoint.eq.1 .and. symtrap.eq."enabled") then
              icntsww=0
              ihalf=1
            endif
            do 2003 indexsww=1,inew
              ieq=ieq+1
              if (j .eq. 1 .and. lbdry(k).ne."conserv") then
                rhs(ieq)=f_(iyh,1,k,l_)
              else
                if(indexsww.le.ipassbnd) then
                  i=iyh+1-indexsww
                elseif(ihalf.eq.0) then
                  i = itl - icntsww
                  ihalf = 1
                else
                  i = itu + icntsww
                  ihalf = 0
                  icntsww = icntsww + 1
                endif
                neq=ieq+(k-1)*iyjx
!.......................................................................
!               Modification for eseswtch.eq."enabled"
!               This code section only to be used for factor.eq.disabled
!.......................................................................

                rhs(ieq)=z00(i,j,k)
                if (cqlpmod.eq."enabled".and.eseswtch.eq."enabled") then
                   cnst=-bnumb(k)/fmass(k)*charge/vnorm
                   if (j.ne.1.or.j.ne.jx) then
                      dfdvz=half*(f(i,j+1,k,l_)-f(i,j-1,k,l_)) &
                               *coss(i,l_)/dx(j)
                      if (i.ne.1 .or. i.ne.iy) then
                         dfdvz=dfdvz &
                               -half*(f(i+1,j,k,l_)-f(i-1,j,k,l_)) &
                               *sinn(i,l_)/x(j)/dy(i,l_)
                      endif
                   elseif (j.eq.1) then
                      dfdvz=zero
                   elseif (j.eq.jx) then
                      dfdvz=(f(i,jx,k,l_)-f(i,jx-1,k,l_)) &
                            *coss(i,l_)/dx(jx)
                      dfdvz=dfdvz &
                            -half*(f(i+1,j,k,l_)-f(i-1,j,k,l_)) &
                            *sinn(i,l_)/x(j)/dy(i,l_)
                   endif
                   rhs(ieq)=cnst*dfdvz
                endif
                rhs(ieq)=rhs(ieq)*scal(neq,l_)
                if (nadjoint.eq.1 .and. symtrap.eq."enabled" &
                  .and. i.eq.itu) rhs(ieq)=0.0
              endif
 2003       continue
 2002     continue
        else  ! factor.ne."disabled", i.e. FACTORIZE

!.......................................................................
!
!     reinitialize to zero matrix and rhs
!
          call bcast(abd,zero,lenabd)
          do i=1,inewjx
            rhs(i)=0.0
            ipivot(i)=0.0
          end do

!..................................................................
!     The next set of variables are defined for use in code for
!     performing the second compression.  The value I is defined as
!     in the past and IEQ is incremented by unity each iteration.
!     However, due to the remapping of equation (node) number to
!     (i,j) pair of coordinates, the value of I "jumps" around.
!
!     As the index IEQ increases (for a given value of J), we start
!     at I=IYH (iy/2) and move clockwise to ITL+1. Then we alternate
!     between the left and right sides, moving clockwise on the right
!     side and CCW on the left side.  As an example, for the case
!     of iy=16, itl=5, a sequence of inew=12 (12=iyh + itl-1) terms is
!     repeated for each value of J (n=12*(j-1))
!
!     ieq=n+1, n+2,  n+3,  n+4, n+5, n+6 , n+7 , n+8 , n+9 , n+10, n+11, n+12
!     i =iyh,iyh-1,itl+1,itl,itl-1,itu+1,itl-2,itu+2,itl-3,itu+3,itl-4,itu+4
!     (in this case, iyh-2=6=itl+1, itu=iy+1-itl=12)
!
!     idistl - the distance, from point (i,j) to point (i+1,j)
!
!     idistr - the distance, from point (i,j) to point (i-1,j)
!
!     IHALF is used during the mapping of ieq (node number) to
!     values of I (mesh point I,J). IHALF is 0 if I<=iyh, otherwise
!     IHALF=1.
!     ICNTSWW - keeps track of the amount to add to ITL (IHALF=0)
!     or ITU (IHALF=1) to obtain current I value.
!
!     if (symtrap.ne."enabled") => no symmetry in trapping region =>
!     compute all the iy mesh points:
!     iyh,iyh+1,iyh-1,iyh+2,...,1,iy (=> assumes iyh=iy/2)
!
!     nadjoint=1: Solve the adjoint equation. As f_adj=0 in trapped region,
!     solve in passing region only (i<itl, i>itu) and  at i=itu, setting
!     f_adj(itu)=0, and f(itl)=f(itu) as boundary cond.
!     The first point is i=itu, then itl-1,itu+1,... =>inew=2*itl-1
!..................................................................

          ipassbnd=iyh+1-itl
          if (symtrap.ne."enabled" .or. nadjoint.eq.1) ipassbnd=0
          ieq=0
          iflag=0

!%OS
          zzzto=0.0
          zzz12=second()
!%OS

!..................................................................
!     do loop 12 is the major loop over j (momentum).
!..................................................................

          do 12 j=1,jx

            icntsww=1
            if (symtrap .ne. "enabled") icntsww=itl-iyh
            ihalf=0
            if (nadjoint.eq.1 .and. symtrap.eq."enabled") then
              icntsww = 0
              ihalf = 1
            endif
            idistl = -1
            idistr = 1

!..................................................................
!     indexsww is the major loop over I (theta).
!..................................................................

            do 11 indexsww=1,inew
              if(indexsww.le.ipassbnd) then
                i=iyh+1-indexsww
              elseif(ihalf.eq.0) then
                i = itl - icntsww
                ihalf = 1
                idistr = 2
                idistl = -2
                if (symtrap .ne. "enabled") then
                  if (i .eq. iyh) idistl=+1
                else
                  if (i .eq. itl-1) idistl=-1
                endif
              else
                i = itu + icntsww
                ihalf = 0
                icntsww = icntsww + 1
                idistr = -2
                idistl = 2
                if (symtrap.ne."enabled" .and. i.eq.iyh+1) idistr=-1
              endif

!.................................................................
!     Separate out various boundary equations:
!     itl (itu) is the lower (upper) p/t boundary.
!     iyh is the mesh point index just below pi/2.
!     iy is the mesh point index of pi and
!     jx is the maximum velocity mesh point index.
!     For the current equation being considered (equation number
!     ieq) define the columns of the sparse matrix of order
!     iyjx where non-zero elements are to be found. This number
!     will vary from 1 to 12 elements in a given row depending
!     on whether the equation is an interior equation (9) or
!     a pass/trapped equation (12 in general) or another
!     boundary equation (1 or 2 or 4 or 6 or 8).
!..................................................................

              ieq=ieq+1
              ibandpieq=ibandto + ieq
              if (i.eq.itu .and. cqlpmod.ne."enabled") go to 31
              if (i.eq.itu .and. nadjoint.eq.1.and.symtrap.eq."enabled") &
                go to 31
              if (j .eq. jx) go to 1
              if (j .eq. 1) go to 2
              if (i .eq. 1) go to 3
              if (i .eq. iy) go to 4
              if (i.eq.itl .and. cqlpmod.ne."enabled") go to 5
              if (i.eq.iyh .and. symtrap.eq."enabled") goto 998

!.................................................................
!     Define the r.h.s.
!.................................................................

              rhs(ieq)=z00(i,j,k)

!.................................................................
!     Non-boundary region, nine point scheme
!..................................................................

              janelt(1)=ieq-inew+idistr
              janelt(2)=ieq-inew
              janelt(3)=ieq-inew+idistl
              janelt(4)=ieq+idistr
              janelt(5)=ieq
              janelt(6)=ieq+idistl
              janelt(7)=ieq+inew+idistr
              janelt(8)=ieq+inew
              janelt(9)=ieq+inew+idistl
              nvar=9

!.................................................................
!     Define the 9 coefficients of the left hand side
!     for example xmm(i,j)*f(i-1,j-1,l_)+ x0m(i,j-1)*f(i,j-1,l_) +..
!     + seven other terms on left side = z00(i,j,k)
!     xmm, xpm, t0m, z00 etc are statement functions defined in
!
!     abd stores the transposed of the band matrix in the LAPACK sgbtrf format
!.................................................................

              zmatcont(1)=xmm(i,j)
              zmatcont(2)=x0m(i,j)
              zmatcont(3)=xpm(i,j)
              zmatcont(4)=xm0(i,j)
              zmatcont(5)=x00(i,j)
              zmatcont(6)=xp0(i,j)
              zmatcont(7)=xmp(i,j)
              zmatcont(8)=x0p(i,j)
              zmatcont(9)=xpp(i,j)

!     Normalize equation
              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

 998          continue

!.................................................................
!     Added for the case i=iyh - using symmetry the "+1" points
!     are the same as the ith point - therefore, the coeffients are
!     added and one element removed.
!     (if symtrap=enabled)
!.................................................................

              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq-inew+idistr
              janelt(2)=ieq-inew
              janelt(3)=ieq+idistr
              janelt(4)=ieq
              janelt(5)=ieq+inew+idistr
              janelt(6)=ieq+inew
              nvar=6
!
              zmatcont(1)=xmm(i,j)
              zmatcont(2)=x0m(i,j)+xpm(i,j)
              zmatcont(3)=xm0(i,j)
              zmatcont(4)=x00(i,j)+xp0(i,j)
              zmatcont(5)=xmp(i,j)
              zmatcont(6)=x0p(i,j)+xpp(i,j)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

!.................................................................
!     j=1 (means velocity=0 boundary condition)
!     First consider case where f is fixed at v=0 (lbdry.ne."conserv")
!     Do not fill in matrix, will be done after "5000 continue"
!..................................................................

 2            continue
              if (lbdry(k).eq."scale" .or. lbdry(k).eq."fixed") then
                janelt(1)=ieq
                nvar=1
                zmatcont(1)=1.
                rhs(ieq)=f_(iyh,1,k,l_)
                xnorm=1.
                go to 5000
              endif

!..................................................................
!     Now the conservative boundary condition. i=1 first...
!.................................................................
              if (i .eq. 1) then
                rhs(ieq)=z00(i,j,k)

                janelt(1)=ieq
                janelt(2)=ieq+idistl
                janelt(3)=ieq+inew
                janelt(4)=ieq+inew+idistl
                nvar=4
!
                zmatcont(1)=x00(i,j)
                zmatcont(2)=xp0(i,j)
                zmatcont(3)=x0p(i,j)
                zmatcont(4)=xpp(i,j)

                call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

                go to 5000

!.................................................................
!     Now for the case j=1 and i=iy (theta=pi).
!.................................................................

              elseif (i .eq. iy) then
                rhs(ieq)=z00(i,j,k)

                janelt(1)=ieq+idistr
                janelt(2)=ieq
                janelt(3)=ieq+inew+idistr
                janelt(4)=ieq+inew
                nvar=4
!
                zmatcont(1)=xm0(i,j)
                zmatcont(2)=x00(i,j)
                zmatcont(3)=xmp(i,j)
                zmatcont(4)=x0p(i,j)

                call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

                go to 5000

!.................................................................
!     j=1 ,i=itl
!.................................................................

              elseif (i.eq.itl .and. cqlpmod.ne."enabled") then
                rhs(ieq)=z00(i,j,k)

                janelt(1)=ieq+idistr
                janelt(2)=ieq+idistr+inew
                janelt(3)=ieq
                janelt(4)=ieq+inew
                janelt(5)=ieq+idistl
                janelt(6)=ieq+idistl+inew
!     itu+1
                if (symtrap .eq. "enabled") then
                  janelt(7)=ieq+2
                  janelt(8)=ieq+inew+2
                else
                  janelt(7)=ieq+3
                  janelt(8)=ieq+inew+3
                endif

                nvar=8
!
                zmatcont(1)=tm0(j)
                zmatcont(2)=tmp(j)
                zmatcont(3)=t00(j)
                zmatcont(4)=t0p(j)
                zmatcont(5)=tp0(j)
                zmatcont(6)=tpp(j)
                zmatcont(7)=tu0(j)
                zmatcont(8)=tup(j)

                call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

                go to 5000

!.................................................................
!     j=1, i=iyh
!.................................................................

              elseif(i.eq.iyh .and. symtrap.eq."enabled") then
                rhs(ieq)=z00(i,j,k)

                janelt(1)=ieq+idistr
                janelt(2)=ieq
                janelt(3)=ieq+inew+idistr
                janelt(4)=ieq+inew
                nvar=4
!
                zmatcont(1)=xm0(i,j)
                zmatcont(2)=x00(i,j)+xp0(i,j)
                zmatcont(3)=xmp(i,j)
                zmatcont(4)=x0p(i,j)+xpp(i,j)

                call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

                go to 5000

!.................................................................
!     case j=1, i not equal 1 or iy or ...
!.................................................................

              else
                rhs(ieq)=z00(i,j,k)

                janelt(1)=ieq+idistr
                janelt(2)=ieq
                janelt(3)=ieq+idistl
                janelt(4)=ieq+inew+idistr
                janelt(5)=ieq+inew
                janelt(6)=ieq+inew+idistl
                nvar=6
!
                zmatcont(1)=xm0(i,j)
                zmatcont(2)=x00(i,j)
                zmatcont(3)=xp0(i,j)
                zmatcont(4)=xmp(i,j)
                zmatcont(5)=x0p(i,j)
                zmatcont(6)=xpp(i,j)
                call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

                go to 5000
              endif

!.................................................................
!     i=1 case (theta=0) : j=2,3,4.....,jx-1
!.................................................................

 3            continue
              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq-inew
              janelt(2)=ieq-inew+idistl
              janelt(3)=ieq
              janelt(4)=ieq+idistl
              janelt(5)=ieq+inew
              janelt(6)=ieq+inew+idistl
              nvar=6
!
              zmatcont(1)=x0m(i,j)
              zmatcont(2)=xpm(i,j)
              zmatcont(3)=x00(i,j)
              zmatcont(4)=xp0(i,j)
              zmatcont(5)=x0p(i,j)
              zmatcont(6)=xpp(i,j)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000
 4            continue

!.................................................................
!     i=iy case (theta=pi) : j=2,3,4,......,jx-1
!.................................................................

              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq-inew+idistr
              janelt(2)=ieq-inew
              janelt(3)=ieq+idistr
              janelt(4)=ieq
              janelt(5)=ieq+inew+idistr
              janelt(6)=ieq+inew
              nvar=6
!
              zmatcont(1)=xmm(i,j)
              zmatcont(2)=x0m(i,j)
              zmatcont(3)=xm0(i,j)
              zmatcont(4)=x00(i,j)
              zmatcont(5)=xmp(i,j)
              zmatcont(6)=x0p(i,j)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

 5            continue

!.................................................................
!     Pass-trapped boundary (lower): all velocities except v=0.
!     (if symtrap=enabled)
!.................................................................

              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq-inew+idistr
              janelt(2)=ieq+idistr
              janelt(3)=ieq+inew+idistr
              janelt(4)=ieq-inew
              janelt(5)=ieq
              janelt(6)=ieq+inew
              janelt(7)=ieq-inew+idistl
              janelt(8)=ieq+idistl
              janelt(9)=ieq+inew+idistl
!     itu+1
              if (symtrap .eq. "enabled") then
                janelt(10)=ieq+2-inew
                janelt(11)=ieq+2
                janelt(12)=ieq+inew+2
              else
                janelt(10)=ieq+3-inew
                janelt(11)=ieq+3
                janelt(12)=ieq+inew+3
              endif
              nvar=12
!
              zmatcont(1)=tmm(j)
              zmatcont(2)=tm0(j)
              zmatcont(3)=tmp(j)
              zmatcont(4)=t0m(j)
              zmatcont(5)=t00(j)
              zmatcont(6)=t0p(j)
              zmatcont(7)=tpm(j)
              zmatcont(8)=tp0(j)
              zmatcont(9)=tpp(j)
              zmatcont(10)=tum(j)
              zmatcont(11)=tu0(j)
              zmatcont(12)=tup(j)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

 1            continue

!.................................................................
!     Upper boundary j=jx...
!     if pass/trapped boundary (and j=jx) then go to 7
!.................................................................

              if (i.eq.itl .and. cqlpmod.ne."enabled") go to 7

!.................................................................
!     if theta=pi go to 14; if theta=0 go to 13
!.................................................................

              if (i .eq. iy) go to 14
              if (i .eq. 1) go to 13
              if(i.eq.iyh .and. symtrap.eq."enabled") goto 991

!.................................................................
!     inner theta points at v=vmax (hyperbolic)
!.................................................................

              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq-inew+idistr
              janelt(2)=ieq+idistr
              janelt(3)=ieq-inew
              janelt(4)=ieq
              janelt(5)=ieq+idistl-inew
              janelt(6)=ieq+idistl
              nvar=6
!
              zmatcont(1)=xmm(i,j)
              zmatcont(2)=xm0(i,j)
              zmatcont(3)=x0m(i,j)
              zmatcont(4)=x00(i,j)
              zmatcont(5)=xpm(i,j)
              zmatcont(6)=xp0(i,j)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

 991          continue

!.......................................................................
!     i=iyh and symtrap=enabled
!.......................................................................
              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq-inew+idistr
              janelt(2)=ieq+idistr
              janelt(3)=ieq-inew
              janelt(4)=ieq
              nvar=4
!
              zmatcont(1)=xmm(i,j)
              zmatcont(2)=xm0(i,j)
              zmatcont(3)=x0m(i,j)+xpm(i,j)
              zmatcont(4)=x00(i,j)+xp0(i,j)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

 7            continue

!.................................................................
!     v=vmax ; theta=pass/trapped bndry
!     (if symtrap=enabled)
!.................................................................

              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq-inew+idistr
              janelt(2)=ieq+idistr
              janelt(3)=ieq-inew
              janelt(4)=ieq
              janelt(5)=ieq+idistl-inew
              janelt(6)=ieq+idistl
!     itu+1
              if (symtrap .eq. "enabled") then
                janelt(7)=ieq+2-inew
                janelt(8)=ieq+2
              else
                janelt(7)=ieq+3-inew
                janelt(8)=ieq+3
              endif

              nvar=8
!
              zmatcont(1)=tmm(jx)
              zmatcont(2)=tm0(jx)
              zmatcont(3)=t0m(jx)
              zmatcont(4)=t00(jx)
              zmatcont(5)=tpm(jx)
              zmatcont(6)=tp0(jx)
              zmatcont(7)=tum(jx)
              zmatcont(8)=tu0(jx)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

 14           continue

!..................................................................
!     v=vmax , theta=pi
!.................................................................

              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq
              janelt(2)=ieq+idistr
              janelt(3)=ieq-inew
              janelt(4)=ieq-inew+idistr
              nvar=4
!
              zmatcont(1)=x00(i,j)
              zmatcont(2)=xm0(i,j)
              zmatcont(3)=x0m(i,j)
              zmatcont(4)=xmm(i,j)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000
 13           continue

!.................................................................
!     v=vmax, theta=0
!.................................................................

              rhs(ieq)=z00(i,j,k)

              janelt(1)=ieq
              janelt(2)=ieq+idistl
              janelt(3)=ieq+idistl-inew
              janelt(4)=ieq-inew

              nvar=4
!
              zmatcont(1)=x00(i,j)
              zmatcont(2)=xp0(i,j)
              zmatcont(3)=xpm(i,j)
              zmatcont(4)=x0m(i,j)

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

 31           continue

!.................................................................
!     i=itu
!     if symtrap.ne.enabled .or. nadjoint.eq.1
!.................................................................

              rhs(ieq)=z00(i,j,k)

              if (nadjoint.eq.1 .and. symtrap.eq."enabled") then
!     set f(itu)=0
                janelt(1)=ieq
                nvar=1
                zmatcont(1)=1.0
                rhs(ieq)=0.0

              else if (j .eq. 1) then

                janelt(1)=ieq+idistr
                janelt(2)=ieq+inew+idistr
                janelt(3)=ieq
                janelt(4)=ieq+inew
                janelt(5)=ieq+idistl
                janelt(6)=ieq+inew+idistl
                janelt(7)=ieq+1
                janelt(8)=ieq+inew+1
                nvar=8
!
                zmatcont(1)=tp0(j)
                zmatcont(2)=tpp(j)
                zmatcont(3)=t00(j)
                zmatcont(4)=t0p(j)
                zmatcont(5)=tu0(j)
                zmatcont(6)=tup(j)
                zmatcont(7)=tm0(j)
                zmatcont(8)=tmp(j)
              else if (j .eq. jx) then

                janelt(1)=ieq-inew+idistr
                janelt(2)=ieq+idistr
                janelt(3)=ieq-inew
                janelt(4)=ieq
                janelt(5)=ieq-inew+idistl
                janelt(6)=ieq+idistl
                janelt(7)=ieq+1-inew
                janelt(8)=ieq+1
                nvar=8
!
                zmatcont(1)=tpm(jx)
                zmatcont(2)=tp0(jx)
                zmatcont(3)=t0m(jx)
                zmatcont(4)=t00(jx)
                zmatcont(5)=tum(jx)
                zmatcont(6)=tu0(jx)
                zmatcont(7)=tmm(jx)
                zmatcont(8)=tm0(jx)

              else

                janelt(1)=ieq-inew+idistr
                janelt(2)=ieq+idistr
                janelt(3)=ieq+inew+idistr
                janelt(4)=ieq-inew
                janelt(5)=ieq
                janelt(6)=ieq+inew
                janelt(7)=ieq-inew+idistl
                janelt(8)=ieq+idistl
                janelt(9)=ieq+inew+idistl
                janelt(10)=ieq+1-inew
                janelt(11)=ieq+1
                janelt(12)=ieq+inew+1
                nvar=12
!
                zmatcont(1)=tpm(j)
                zmatcont(2)=tp0(j)
                zmatcont(3)=tpp(j)
                zmatcont(4)=t0m(j)
                zmatcont(5)=t00(j)
                zmatcont(6)=t0p(j)
                zmatcont(7)=tum(j)
                zmatcont(8)=tu0(j)
                zmatcont(9)=tup(j)
                zmatcont(10)=tmm(j)
                zmatcont(11)=tm0(j)
                zmatcont(12)=tmp(j)

              endif

              call impnorm(xnorm,zmatcont,rhs(ieq),nvar)

              go to 5000

!..................................................................
!     After an equation has been set, code jumps to 5000
!..................................................................


 5000         continue

!.......................................................................
!     put normalized contributions into matrix
!.......................................................................

              do jcont=1,nvar
                abd(ibandpieq-janelt(jcont),janelt(jcont))= &
                  zmatcont(jcont)
              end do

!.......................................................................

              kku=ieq+(k-1)*iyjx
              scal(kku,l_)=1./xnorm

!..................................................................
!     End of major j and I loops..
!..................................................................

 11         continue
 12       continue

!%OS
!%OS  print *,' time used for constructing matrix= ',second()-zzz12
!%OS
!..................................................................
!     next endif is the end for "if(factor.eq."disabled") then
!     thus matrix and rhs are defined
!..................................................................

        endif

!..................................................................
!     Call direct solve Gaussian elimination algorithm.
!..................................................................

!
!%OS        zzzt1=second()

!     factorize matrix
!_cray  64-bit-compiler uses sgbtrf (from lapack library).
!_pc    32-bit-compiler uses dgbtrf (from lapack library).
        if (factor .ne. "disabled") then
!          call sgbtrf(inewjx,inewjx,ml,mu,abd,md1abd,ipivot,info)
          call dgbtrf(inewjx,inewjx,ml,mu,abd,md1abd,ipivot,info)
          if (info .ne. 0) then
            print *,' warning after sgbtrf in impavnc: info = ',info
            stop 'impavnc 1'
          endif
        endif

!     solve system
!_cray  64-bit-compiler uses sgbtrs (from lapack library).
!_pc    32-bit-compiler uses dgbtrs (from lapack library).
        inbrhs = 1
        transpose = 'n'
!        call sgbtrs(transpose,inewjx,ml,mu,inbrhs,abd,md1abd,ipivot
        call dgbtrs(transpose,inewjx,ml,mu,inbrhs,abd,md1abd,ipivot &
          ,rhs,inewjx,info)
        if (info .ne. 0) then
          print *,' warning after sgbtrs in impavnc: info = ',info
          stop 'impavnc 2'
        endif

!%OS
!%OS        print *,' sol in sgbtr= ',second()-zzzt1

!.................................................................
!     The next small block expands the compressed f matrix - and replicates
!     the mirror-imaged values in the trapping region. (i.e. copy solution rhs to f)
!..................................................................

        ix=iy-itu
        do 966 j=1,jx
          i1=(j-1)*inew

!..................................................................
!     For each j (velocity) value, take the first ipassbnd # of elements
!     from rhs and place them in the right half of the pass-trapped
!     region starting at iyh.
!..................................................................

          ihalf=0

          if (symtrap.eq."enabled".and. nadjoint.ne.1) then
            do 958 ii=1,ipassbnd
              f(iyh+1-ii,j,k,l_)=rhs(i1+ii)
 958        continue
            icntsww=1
          else if (nadjoint.eq.1 .and. symtrap.eq."enabled") then
            ipassbnd=0
            icntsww=0
            ihalf=1
          else
            ipassbnd=0
            icntsww=itl-iyh
          endif

!..................................................................
!     Alternating between the left and right sides, place the remaining
!     values (for this value of J) below ITL and above ITU.
!..................................................................

          do 959 ii=ipassbnd+1,inew
            if(ihalf.eq.0) then
              f(itl-icntsww,j,k,l_)=rhs(i1+ii)
              ihalf=1
            else
              f(itu+icntsww,j,k,l_)=rhs(i1+ii)
              ihalf=0
              icntsww=icntsww+1
            endif
 959      continue

!..................................................................
!     Mirror the right half of the pass-trapped region into the left hal
!..................................................................

          if (symtrap .eq. "enabled") then
            if (nadjoint .eq. 1) f(itl,j,k,l_)=f(itu,j,k,l_)
            do 957 i=iyh+1,itu
              ii=iy+1-i
              f(i,j,k,l_)=f(ii,j,k,l_)
 957        continue
          endif

 966    continue


!.................................................................
!     Differentiate the solution and compare with r.h.s.
!.................................................................

        if (iactst.ne."disabled") call impchk(k)

!.................................................................
!     Redefine f at v=0 so it is unique. (should be before check?)
!.................................................................

!%OS
        if (noffelpr .eq. 9) then

!          call scopy(iyjx2a,f(0,0,k,l_),1,fxsp(0,0,k,l_),1)
          call dcopy(iyjx2a,f(0,0,k,l_),1,fxsp(0,0,k,l_),1)
          s=0.
          t=0.
          do 21 i=1,iy
            s=s+vptb(i,lr_)*cynt2(i,l_)
            t=t+vptb(i,lr_)*cynt2(i,l_)*f(i,1,k,l_)
 21       continue
          do 22 i=1,iy
            f(i,1,k,l_)=t/s
 22       continue

!%OS
        endif

!.......................................................................
!     compute velocity source term needed in transport equ. for ADI method
!.......................................................................

        if (adimeth.eq."enabled" .and. transp.eq."enabled" &
          .and. n.ge.nonadi) call tdtrvsou(k)

!.................................................................
!     Call a routine which computes density gains and losses
!.................................................................

        call diagimpd(k)

!%OS
        if (noffelpr .ne. 9) then

!          call scopy(iyjx2a,f(0,0,k,l_),1,fxsp(0,0,k,l_),1)
          call dcopy(iyjx2a,f(0,0,k,l_),1,fxsp(0,0,k,l_),1)
          s=0.
          t=0.
          do 2100 i=1,iy
            s=s+vptb(i,lr_)*cynt2(i,l_)
            t=t+vptb(i,lr_)*cynt2(i,l_)*f(i,1,k,l_)
 2100     continue
          do 2200 i=1,iy
            f(i,1,k,l_)=t/s
 2200     continue

        endif
!%OS

!.......................................................................
!     eseswtch="enabled" section, to control and calculate
!     fluxes for quasi-neutral electrostatic field.
!.......................................................................

        if (cqlpmod.eq."enabled" .and. eseswtch.eq."enabled") then

           if (ifirst.eq.1) then

!             Store new f into fh
!              call scopy(iyjx2a*ngen,f(0,0,1,l_),1,fh(0,0,1,l_),1)
              call dcopy(iyjx2a*ngen,f(0,0,1,l_),1,fh(0,0,1,l_),1)
              factor="disabled"
              ifirst=2
              go to 10

           elseif (ifirst.eq.2) then

!             g function is in f.
!             Calculate fluxes, restore f with f_, and return.

              flux1(l_)=fluxpar( &
                   1,x,coss(1,l_),cynt2(1,l_),cint2,temp1,iy,jx,iya)
              flux2(l_)=fluxpar( &
                   1,x,coss(1,l_),cynt2(1,l_),cint2,temp2,iy,jx,iya)
!              call scopy(iyjx2a*ngen,f_(0,0,1,l_),1,f(0,0,1,l_),1)
              call dcopy(iyjx2a*ngen,f_(0,0,1,l_),1,f(0,0,1,l_),1)

              go to 999
           endif

        endif

!..................................................................
!     End of "k" (species) loop.
!..................................................................

 600  continue

      if(dalloc.eq."enabled") then

!_cray   call hpdeallc(abdptr,lenabd)
!_pc     call free(abdptr)
         deallocate(abd,STAT=istat)

      endif

      irstart=0
 999  continue
      return

      end subroutine impavnc

!      CONTAINS !removed this construct for franklin.nerc.gov: pg compiler

         real(c_double) function z00(i,j,k)
         use param_mod
      use comm_mod
         implicit integer (i-n), real(c_double) (a-h,o-z)
         save


!.......................................................................
!     z00 is the right hand side of the equation, and holds the
!     explicit-in-time rhs of the FP difference equations.
!
!     The terms involving the factors bsl, bsu , x**_ and t0**_
!     are related to calculation of the bootstrap effect.
!     We assume virtually that the distribution is skewed asymetrically
!     in the trapped region...that is we assume (virtually) that
!     f(itl).ne.f(itu) and that the difference is driven by
!     a df/dr term through bsl and bsu. Since this term involves f at
!     different radial positions, it cannot figure into the solution
!     implicitly, that is, it is differenced explicitly. The resulting
!     contributions appear below. There will be contributions from
!     i=itl-1, itu+1, itl and itu only.
!     All contributions are zero elsewhere, and are zero everywhere
!     if bootcalc= "disabled".   (Refer to Harvey et al, 1993 Sherwood
!     Theory Mtg; E. Westerhof and A.G. Peters, Computer Physics Comm.,
!     Vol. 95, p. 131-138 (1996).)
!.......................................................................


!     statement functions for itl or itu [depracated in f95]

      t0ml_(j)=qz(j)*( &
       cl(itl-1,j-1)*dj(itl,j-1,l_)*eym5(itl,l_)) &
       +r2y(j)*(-de(itl-1,j)*(1.-di(itl-1,j-1,l_))) &
       /(2.*dx(j))

      t00l_(j)= &
       +qz(j)*( &
       -cl(itl-1,j)*dj(itl,j,l_)*eym5(itl,l_) &
       +cl(itl-1,j-1)*(1.-dj(itl,j-1,l_))*eym5(itl,l_)) &
       +r2y(j)*(dd(itl-1,j)*(1.-di(itl-1,j,l_)) &
       +df(itl-1,j)*eym5(itl,l_))

      t0pl_(j)=qz(j)*( &
       -cl(itl-1,j)*eym5(itl,l_)*(1.-dj(itl,j,l_))) &
       +r2y(j)*de(itl-1,j)/(2.*dx(j))*(1.-di(itl-1,j+1,l_))


      t0mu_(j)=qz(j)*( &
       -cl(itu+1,j-1)*dj(itu,j-1,l_)*eyp5(itu,l_)) &
       +r2y(j)*( &
       +de(itu,j)*di(itu,j-1,l_))/(2.*dx(j))

      t00u_(j)= &
       +qz(j)*( &
       +cl(itu+1,j)*dj(itu,j,l_)*eyp5(itu,l_) &
       -cl(itu+1,j-1)*(1.-dj(itu,j-1,l_))*eyp5(itu,l_)) &
       +r2y(j)*( &
       -dd(itu,j) &
       *di(itu,j,l_) &
       +df(itu,j)*eyp5(itu,l_))

      t0pu_(j)=qz(j)*( &
       +cl(itu+1,j)*(1.-dj(itu,j,l_))*eyp5(itu,l_)) &
       +r2y(j)*( &
       -de(itu,j)*di(itu,j+1,l_)/ &
       (2.*dx(j)))



!  Test if bootstrap calc is needed (giving iboot=1)
         iboot=0
         if (bootcalc.ne."disabled" &
              .and. (i.eq.(itl-1).or.i.eq.itl.or. &
                     i.eq.itu.or.i.eq.(itu+1))) iboot=1

         z00f=vptb(i,lr_)*(f_(i,j,k,l_)/dtreff+so(i,j)) + &
              spasou(i,j,k,l_)

         if (iboot.eq.1) then
            if (i.eq.(itl-1)) then
!              itl-1 case:
               z00itl1=z00f &
                    -xpm(i,j)*bsl(j-1,k,l_)-xp0(i,j)*bsl(j,k,l_) &
                    -xpp(i,j)*bsl(j+1,k,l_)
               z00t=z00itl1
            else
!              itu+1 case:
               if (i.eq.(itu+1))then
                  z00itu1=z00f &
                       -xmm(i,j)*bsu(j-1,k,l_)-xm0(i,j)*bsu(j,k,l_) &
                       -xmp(i,j)*bsu(j+1,k,l_)
                  z00ff=z00itu1
!              itl case:
               else
                  z00itl=z00f &
                       -(t0ml_(j)*bsl(j-1,k,l_)+t00l_(j)*bsl(j,k,l_) &
                       +t0pl_(j)*bsl(j+1,k,l_)+t0mu_(j)*bsu(j-1,k,l_) &
                       +t00u_(j)*bsu(j,k,l_)+t0pu_(j)*bsu(j+1,k,l_))
                  z00ff=z00itl
               endif
               z00t=z00ff
            endif
            z00=z00t

         else
            z00=z00f
         endif

         end function z00



!      end subroutine impavnc
end module impavnc_mod
