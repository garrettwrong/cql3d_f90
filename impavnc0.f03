!
!
module impavnc0_mod
  use iso_c_binding, only : c_double
  ! XXX this is probably a bug in the original code.
  use advnce_mod, only :  k => advnce_k ! use advnce_mod's advnce_k
  use bsu_mod, only : bsu
  use bsl_mod, only : bsl
  ! gauss_, these are used outside module in it3dalloc.f
  integer, public ::  inewjmax, ialloc      
  ! XXX this was an integer somewhere?
  real(c_double), pointer, public ::  abd(:,:)
  integer, pointer, public :: ipivot(:)
  ! ampf, these are used outside module in it3dalloc.f
  real(c_double), pointer, public :: ampfda(:,:), ampfdd(:,:)      
  save

contains

  subroutine impavnc0(kopt)
    use param_mod
    use cqcomm_mod
    use advnce_mod
    use r8subs_mod, only : dgbtrf, dgbtrs, dcopy
    implicit integer (i-n), real*8 (a-h,o-z)

    !.................................................................
    !     Impose unicity point at u=0 as boundary condition
    !
    !     kopt=0:
    !     This routine advances the Fokker-Planck equation one time step
    !     for the implicit case (implct=enabled)
    !     (IMPAVNC0 calls IMPCHK which differentiates the solution and
    !     compares it with the r.h.s of the advanced F.P. equation to 
    !     check for accuracy)
    !     kopt=3:
    !     Calc h and g functions for iterative Amp-Far eqn solve
    !.................................................................


    !MPIINSERT_INCLUDE     

    !BH080303    real*4 etime,tm1,tm(2) !Dclrtns for lf95 ETIME(). See man etime.
    !BH080303    Using f95 intrinsic subroutine, cpu_time
    real*4 tm1,tm(2)    !Dclrtns for lf95 ETIME(). See man etime.

    character*8 alloc, factor, dalloc
    dimension zavarj1(iy+4),ijaj1(iy+4)
    dimension  zmatcont(12), janelt(12)
    character*1 transpose

    integer icount_imp   !Counter for calls to impavnc0
    data icount_imp/0/
    save inewmax, mlmax, mumax, icount_imp, icoeff, ieq, &
         ieq_tot, icoeff_est, icoeff_est_tot

    twopi2=twopi*twopi
    !.......................................................................
    !     Indicate progress of code
    !.......................................................................

    !     icount_imp counts total calls to impavnc0 (from achiefn).
    !     At each time step: first call to impavnc0 is with lr_=1.
    !     ifirst_lr=1 indicates lr_=1 call (otherwise ifirst_lr=0).
    !     ilast_lr=1 indicates lr_=lrz call (otherwise ilast_lr=0).

    icount_imp=icount_imp+1
    ifirst_lr=0
    ilast_lr=0
    if (l_.eq.1 .or. lrz.eq.1) ifirst_lr=1
    if (l_.eq.lrz) ilast_lr=1
    !      write(*,'(a,4i5)')'impavnc0: n,icount_imp,ifirst_lr,ilast_lr',
    !     +                     n,icount_imp,ifirst_lr,ilast_lr

    !.......................................................................
    !     nadjoint is set = 0, to remove effect of Olivier Sauter's adjoint
    !        operator changes for this subroutine.   This is not checked
    !        out for iterative solve.  Sauter has used it with the
    !        soln_methed="direct" solve, in bootstrap current calculations
    !        (O.Sauter, C. Angioni, and Y.R. Lin-Liu, Pop 66, 2834 (1999).
    !.......................................................................

    nadjoint=0

    !.......................................................................
    !     icount_ese/icount_ampf are counters for eseswtch or ampfmod cases
    !.......................................................................

    icount_ese=1
    icount_ampf=1

    !.......................................................................
    !     define effective number of unknowns in the theta direction and
    !     resulting band widths
    !     see below for comments on inew, etc.
    !.......................................................................

    itl=itl_(l_)
    itu=itu_(l_)
    !BH071029:  During work on 3d fully-implicit solve, found a problem
    !BH071029:  with iyh=itl cases.  Fixing this case, and moving ipassbnd
    !BH071029:  setting to here [See further notes below for ipassbnd]:
    !BH071029:  [But subsequently found had problems with iyh=itl elsewhere
    !BH071029:  in the code, so suggest bigger rya(1) or iy for time being.]
    if (itl.lt.iyh) then
       inew=iy/2 + itl - 1
       ipassbnd=iyh-itl+1
    else  !i.e., itl=iyh
       inew=iy
       ipassbnd=0  !i.e., no trapped region for iyh=itl, 
       !as in following case.
    endif
    if (symtrap .ne. "enabled") then
       inew=iy
       ipassbnd=0
    endif
    if (nadjoint.eq.1 .and. symtrap.eq."enabled") inew=2*itl-1
    if (symtrap.ne."enabled" .or. nadjoint.eq.1) ipassbnd=0

    !     On first call to impanvc0, set up some max dimensions
    !cc      if(icount_imp.eq.1) then 
    inewmax=inew
    if (meshy.eq."fixed_mu".or. manymat.eq."enabled") then
       do 15 ll=2,lrors
          inew1=iy_(ll)/2+itl_(ll)-1
          if (symtrap .ne. "enabled") inew1=iy_(ll)
          if (nadjoint.eq.1 .and. symtrap.eq."enabled") &
               inew1=2*itl_(ll)-1
          if (inew1.gt.inewmax) inewmax=inew1
15     end do
    endif
    mlmax=inewmax+2 ! YuP: not used?
    mumax=inewmax+2 ! YuP: not used?
    inewjmax=inewmax*jx
    !cc      endif

    ml=inew+2
    mu=inew+2
    if (symtrap.ne."enabled" .and. cqlpmod.ne."enabled") then
       !     itu+1 contributes to the equation for itl
       ml=inew+3
       mu=inew+3
       mlmax=inewmax+3 ! YuP: not used?
       mumax=inewmax+3 ! YuP: not used?
    endif

    ibandto=ml+mu+1
    inewjx=inew*jx
    navnc=navnc+1
    md1abd= 3*(iy+3)+1 ! Should be at least 2*ml+mu+1 (see below)


    !YuP[07-2017] Moved this part from micxinit and tdinitl, 
    ! as it may depend on k in a multi-species run
    ieq_tot=0
    do ll=1,lrors
       if(symtrap.eq."enabled") then
          !Number of indep theta pts, for symmetric trapped region
          inew_(ll)=iy_(ll)/2+itl_(ll)-1   
       else
          inew_(ll)=iy_(ll) !general case (not symmetric trapped reg.)
       endif
       inewjx_(ll)=inew_(ll)*jx   !Number eqns for ll
       ieq_tot=ieq_tot+inew_(ll)*jx 
       if (ll.eq.1) then
          ieq_(ll)=1
       else  ! Eqn no. at beginning of each flux surface:
          ieq_(ll)=ieq_(ll-1)+inew_(ll-1)*jx  
       endif
    enddo
    ieq_(lrors+1)=ieq_tot
    ieq_tot=ieq_(lrors+1)

    !BH070419:   ITSOL Add, down to ....
    !.................................................................
    !     Set size of LAPACK and SPARSKIT2 matrices, and allocate space
    !     Probably should move these calcs to it3dalloc.f
    !.................................................................

    !      write(*,*)'impavnc0: iy,iyh,itl,itu,ml,mu,inew,jx,inewjx',
    !     +     iy,iyh,itl,itu,ml,mu,inew,jx,inewjx

    iwk_ilu=25000000  !sub ilut will give error if not suff lrg
    !                          !NEED bigger than single flux surface solve
    !                          !values, for given iy,jx.
    krylov=50  !Krylov subspace size, see below

    if ( soln_method.eq.'itsol' ) then
       icoeff_est=(iy+3)*3+(jx-1)*9*inew  !Number of coeffs + abit
       icsrij=icoeff_est
       icsrip=inewjx+1
       icsri2=2*inewjx
       icsrikry=(krylov+1)*inewjx
    endif  !  on soln_method.eq.'itsol'

    if ( soln_method.eq.'itsol1' ) then
       lapacki=ml+mu+1
       lapackj=inewjx
       icsrij=lapacki*lapackj
       icsrip=lapackj+1
       icsri2=2*lapackj
       icsrikry=(krylov+1)*lapackj
    endif  !  on soln_method.eq.'itsol1'

    if ( soln_method.eq.'it3dv' .and. ifirst_lr.eq.1 ) then
       ! no transport here
       icoeff_est_tot=0
       do ll=1,lrors
          icoeff_est=(iy_(ll)+2)*3+(jx-1)*9*inew_(ll) !Number coeffs + abit
          icoeff_est_tot=icoeff_est_tot+icoeff_est
       enddo
       icsrij=icoeff_est_tot
       icsrip=ieq_tot+1
       icsri2=2*ieq_tot
       icsrikry=(krylov+1)*ieq_tot
    endif

    if ( soln_method.eq.'it3drv' .and. ifirst_lr.eq.1 ) then
       icoeff_est_tot=0
       do ll=1,lrors
          icoeff_est=(iy_(ll)+2)*3+(jx-1)*9*inew_(ll) & !Number coeffs + abit &
          + 2*inew_(ll)*jx  !  for additional Drr coeffs
          icoeff_est_tot=icoeff_est_tot+icoeff_est
       enddo
       icsrij=icoeff_est_tot
       icsrip=ieq_tot+1
       icsri2=2*ieq_tot ! for jw_ilu(icsri2) can be just ieq_tot ?
       icsrikry=(krylov+1)*ieq_tot
       ! Above gives storage for the velocity space coeffs.
       ! Additional storage is needed for the radial coeff matrix
       ! This data will be stored in ar_csr, jar_csr, iar_csr
       icsrijr=ieq_tot*3   ! roughly .le.3 coeffs per equation
       ! Also require larger matrix to contain the CSR addition of
       ! the velocity and radial space coeff matrices.
       icsrijc=icsrij+2*ieq_tot !sum of a_csr and off-diagon. ar_csr
    endif

    iwk_ilu=max(icsrijc*5,50000000)
    !iwk_ilu=25000000 !sub ilut will give error (-2,-3) if not suff lrg
    !NEED bigger than single flux surface solve
    !values, for given iy,jx.

    !     Allocate it3d related storage [No-op for soln_method='direct']

    !-YuP      if (ifirst_lr.eq.1 .and. n.eq.1 .and. nefiter.eq.1) call it3dalloc
    if(soln_method.ne.'direct') then 
       if(ASSOCIATED(rhs0)) then  !-YuP  New logic for allocation
          ! rhs0 and other arrays are already allocated => do nothing
       else ! Not allocated yet
          call it3dalloc
       endif
    endif

    !     Just to check:
    !      write(*,*)'impavnc0:size(abd_lapack,1),size(abd_lapack,2)=',
    !     +     size(abd_lapack,1),size(abd_lapack,2)
    !      write(*,*)'impavnc0:size of a_csr,ja_csr,ia_csr,alu,jlu,ju',
    !     +     size(a_csr),size(ja_csr),size(ia_csr),size(alu),
    !     +     size(jlu),size(ju)
    !      write(*,*)'impavnc0:size of jw_ilu,w_ilu_rhs0,sol,vv',
    !     +     size(jw_ilu),size(w_ilu),size(rhs0),size(sol),size(vv)



    !     So, no confusion:
    if (soln_method.eq.'itsol1')  &
         call bcast(abd_lapack(1,1),zero,icsrij)

    !      write(*,*)'impavnc0: size(ja_csr),icsrip=',size(ja_csr),icsrij
    !      write(*,*)'impavnc0: size(ia_csr),icsrip=',size(ia_csr),icsrip

    if (soln_method.eq.'itsol' .or. soln_method.eq.'itsol1') then
       call bcast(a_csr(1),zero,icsrij)
       call ibcast(ja_csr(1),0,icsrij)
       call ibcast(ia_csr(1),0,icsrip)
    elseif (soln_method.eq.'it3dv') then
       if (ifirst_lr.eq.1) then
          call bcast(a_csr(1),zero,icsrij)
          call ibcast(ja_csr(1),0,icsrij)
          call ibcast(ia_csr(1),0,icsrip)
       endif
    elseif (soln_method.eq.'it3drv') then
       if (ifirst_lr.eq.1) then
          !BH080429   Problem:  ipofi only calculated once in code, in tdtranspn.
          !BH080429            call ibcast(ipofi(1,1),0,iy*lrz)
          call bcast(a_csr(1),zero,icsrij)
          call ibcast(ja_csr(1),0,icsrij)
          call ibcast(ia_csr(1),0,icsrip)
          call bcast(ar_csr(1),zero,icsrijr)
          call ibcast(jar_csr(1),0,icsrijr)
          call ibcast(iar_csr(1),0,icsrip)
          call bcast(ac_csr(1),zero,icsrijc)
          call ibcast(jac_csr(1),0,icsrijc)
          call ibcast(iac_csr(1),0,icsrip)
       endif
    endif

    !     Copy current distribution into f_..   YuP: moved to tdchief:
    !cc      call dcopy(iyjx2*ngen,f(0,0,1,l_),1,f_(0,0,1,l_),1) 

    !      write(*,*) 'impavnc0, l_,elecfld:', l_ ,elecfld(l_)  
    !..................................................................
    !     Logic for allocation/deallocation
    !..................................................................

    !-YuP      alloc="disabled"
    !-YuP      dalloc="disabled"  !-YuP: not used
    !     allocate once for all time-step
    !-YuP      if ( soln_method.eq.'direct' .or. soln_method.eq.'itsol1') then

    !-YuP  C--MPIREPLACE_ALLOC  
    !-YuP: Now allocation is based on first call to impavnc0
    !-YuP      if (l_.eq.1  .and. n.eq.1 .and. nefiter.eq.1)  then 
    !-YuP         alloc="enabled"
    !-YuP         ialloc=1
    !-YuP      endif
    !     Special case for multiple calls at n=nstop
    !-YuP      if (n.eq.nstop.and.ialloc.eq.0) then
    !-YuP          alloc="enabled"
    !-YuP  Problem: called during first iteration of eflditer;
    ! Should be after last nefiter iteration (which is unknown).
    ! Not really necessary to allocate at last time step?         
    !-YuP      endif

    !     deallocate once at last time-step
    !-YuP C--MPIREPLACE_DEALLOC !-YuP: no deallocation anymore
    !-YuP      if (l_.eq.lrors .and. n.ge.nstop) then
    !-YuP         dalloc="enabled" !-YuP: not used
    !-YuP         ialloc=0 !-YuP: not used
    !-YuP      endif

    !-YuP      endif  ! On soln_method

    !..................................................................
    !     Allocate space for abd and ipivot. 
    ! YuP-101231  New logic for allocation: 
    ! simply check whether the array is already allocated or not.
    !..................................................................
    if (soln_method.eq.'direct'.or.soln_method.eq.'itsol1') then      
       lenabd=md1abd*inewjmax
       !-YuP      if (alloc.eq."enabled") then
       if(ASSOCIATED(ipivot)) then  
          ! abd and ipivot are already allocated => do nothing 
       else ! Not allocated yet
          write(*,*)'impavnc0 Before allocate; Size req. for abd:',lenabd
          allocate(abd(md1abd,inewjmax),STAT=istat) !Big size ~3*jx*iy^2
          write(*,*)'impavnc0 abd: Allocate/istat=',istat
          !XXXcall bcast(abd,0,lenabd)
          abd = 0.
          allocate(ipivot(iyjx),STAT=istat)
          call ibcast(ipivot,0,iyjx)
          ! will preset the values to zero in loop k=1,ngen
       endif
    endif ! on soln_method= direct or itsol1

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
       !     Logic for LU factorization..
       !..................................................................

       factor="disabled"
       if (ngen.gt.1 ) then
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

       call coefmidv(da,1)
       call coefmidv(db,2)
       call coefmidv(dc,3)

       !.................................................................
       !     the theta flux..
       !.................................................................

       call coefmidt(dd,1)
       call coefmidt(de,2)
       call coefmidt(df,3)

       !.................................................................
       !     The differencing is done ala Chang and Cooper as generalized
       !     by Karney to 2-D. This requires defining a couple of arrays
       !     di(i,j,k,l_) and dj(i,j,k,l_) which are used to define f 
       !     at mid mesh points: 
       !     to wit f(i,j+1/2,k,l_)=f(i,j+1,k,l_)*(1.-dj(i,j,k,l_))+ 
       !     f(i,j,k,l_)*dj(i,j,k,l_). Similarly for theta.
       !     The routines coefwtj and coefwti are used to
       !     determine the dj and di arrays.
       !.................................................................

       call coefwtj(k)
       call coefwti(k)

       !.................................................................
       !     Loop-back point for eseswtch or ampfmod calculations
       !.................................................................

10     continue

       !.................................................................
       !     Proceed to setup coeffs and call the Gaussian elimination routine.
       !     i indexes theta, j indexes speed and ieq keeps a running total. It
       !     is necessary to separate out the various boundary equations
       !     from the the inner 9-point equation. Boundary equations can
       !     involve from 1 to 12 points (or unknowns) depending on location
       !     in velocity space.  In addition, lbdry(k)='conserv'.or.'consscal'
       !     is a special case at j=1 involving inew unknowns.

       !     ieq will keep track of which equation is being considered
       !     below. This defines the row of the matrix being inverted. The
       !     theta mesh points (index i) are incremented in the inner
       !     loop.
       !.................................................................

       if (soln_method.eq."itsol" .or. soln_method.eq."itsol1" &
            .or. soln_method.eq."direct") then
          ieq=0
          icoeff=0         !  CSR coefficient storage counter (ITSOL)
       elseif(soln_method.eq.'it3dv' .or. soln_method.eq.'it3drv') then
          ! reset counters only at first flux surface on new time-step
          if ( ifirst_lr.eq.1 ) then
             ieq=0
             icoeff=0      !  CSR coefficient storage counter 
          endif
       endif
       if (itl.lt.iyh) then
          inew=iy/2 + itl - 1
          ipassbnd=iyh-itl+1
       else                    !i.e., itl=iyh
          inew=iy
          ipassbnd=0           !i.e., no trapped region for iyh=itl, 
          !      as in following case.
       endif


       !cc        if(ieq. eq. ieq_(l_+1)-1)  ieq=ieq_(l_)-1     
       !.................................................................
       !     abd is the band matrix in a format expected by the sgbtrf routine.
       !     That is, if a(i,j) is the matrix to invert, we have:
       !     abd(iband+i-j,j) = a(i,j), with iband = ml+mu+1=total bandwidth
       !     As if each column is moved up/down such that each diagonal element
       !     is at abd(iband,j).  NOTE:  i,j in this comment designate the row
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

       if (factor.eq."disabled") then   !Using prior factorization...
          !endif on factor at line 1834
          if (.not.(ampfmod.eq."enabled" .and. kopt.eq.3 &
               .and.cqlpmod.ne."enabled" .and. k.eq.kelec)) then  

             write(*,*)'impavnc0:  factor.eq."disabled", OK?'
             ipassbnd=iyh+1-itl
             if (symtrap.ne."enabled" .or. nadjoint.eq.1) ipassbnd=0
             do 2002 j=1,jx
                ihalf=0
                icntsww=1
                if (iyh.eq.itl) icntsww=0
                if (symtrap .ne. "enabled") icntsww=itl-iyh
                if (nadjoint.eq.1 .and. symtrap.eq."enabled") then
                   icntsww=0
                   ihalf=1
                endif
                do 2003 indexsww=1,inew
                   ieq=ieq+1
                   if (j .eq. 1 .and. (lbdry(k).ne."conserv" .and.  &
                        lbdry(k).ne."consscal")) then
                      rhs(ieq)=f_(iyh,1,k,l_) ! YuP: Why iyh?  BH: all same
                   else ! includes j=1 when lbdry(k).eq."conserv"
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
                      !For scal(neq,l_) :
                      neq=ieq+(k-1)*iyjx   ! for multiple species (k.gt.1)

                      alambda= vptb(i,lr_) ! ZOW or Hybrid-FOW

                      rhs(ieq)=z00(i,j,k)
                      !.......................................................................
                      !               Modification for eseswtch.eq."enabled" to obtain Kupfer
                      !               g function.
                      !               This code section only to be used for factor.eq.disabled
                      !.......................................................................

                      rhs(ieq)=z00(i,j,k)
                      if (cqlpmod.eq."enabled".and.eseswtch.eq."enabled") then
                         cnst=-bnumb(k)/fmass(k)*charge/vnorm
                         if (j.ne.1.or.j.ne.jx) then
                            dfdvz=half*(f(i,j+1,k,l_)-f(i,j-1,k,l_)) &
                                 *coss(i,l_)*dxi(j)
                            if (i.ne.1 .or. i.ne.iy) then
                               dfdvz=dfdvz  &
                                    -half*(f(i+1,j,k,l_)-f(i-1,j,k,l_)) &
                                    *sinn(i,l_)*xi(j)*dyi(i,l_)
                            endif
                         elseif (j.eq.1) then
                            dfdvz=zero
                         elseif (j.eq.jx) then
                            dfdvz=(f(i,jx,k,l_)-f(i,jx-1,k,l_)) &
                                 *coss(i,l_)*dxi(jx)
                            dfdvz=dfdvz  &
                                 -half*(f(i+1,j,k,l_)-f(i-1,j,k,l_)) &
                                 *sinn(i,l_)*xi(j)*dyi(i,l_)
                         endif
                         rhs(ieq)=cnst*dfdvz
                      endif  !  On cqlpmod/eseswtch

                      rhs(ieq)=rhs(ieq)*scal(neq,l_)

                      if (nadjoint.eq.1 .and. symtrap.eq."enabled" &
                           .and. i.eq.itu) rhs(ieq)=0.0

                   endif  ! On if j.eq.1.and.lbdtry(k).ne.conserv;  else

2003            end do ! On indexsww
2002         end do  ! On j
          endif  ! On .not.(ampfmod....


          !.......................................................................
          !               Modification for ampfmod.eq."enabled" to obtain Kupfer
          !               g funtion.
          !               This code section only to be used for factor.eq.disabled
          !.......................................................................

          if (ampfmod.eq."enabled" .and. kopt.eq.3 .and.icount_ampf.eq.2 &
               .and.cqlpmod.ne."enabled" .and. k.eq.kelec) then  
             !$$$                   cnst=-bnumb(k)/fmass(k)*charge/vnorm
             !$$$                   if (j.ne.1.or.j.ne.jx) then
             !$$$                      dfdvz=half*(f(i,j+1,k,l_)-f(i,j-1,k,l_))
             !$$$     +                         *coss(i,l_)*dxi(j)
             !$$$                      if (i.ne.1 .or. i.ne.iy) then
             !$$$                         dfdvz=dfdvz 
             !$$$     +                         -half*(f(i+1,j,k,l_)-f(i-1,j,k,l_))
             !$$$     +                         *sinn(i,l_)*xi(j)*dyi(i,l_)
             !$$$                      endif
             !$$$                   elseif (j.eq.1) then
             !$$$                      dfdvz=zero
             !$$$                   elseif (j.eq.jx) then
             !$$$                      dfdvz=(f(i,jx,k,l_)-f(i,jx-1,k,l_))
             !$$$     +                      *coss(i,l_)*dxi(jx)
             !$$$                      dfdvz=dfdvz 
             !$$$     +                      -half*(f(i+1,j,k,l_)-f(i-1,j,k,l_))
             !$$$     +                      *sinn(i,l_)*xi(j)*dyi(i,l_)
             !$$$                   endif
             !$$$                   rhs(ieq)=cnst*dfdvz
             !
             !BH140101: RHS from BA coeffs for tor E field
             !
             if(ASSOCIATED(ampfda)) then
                continue
             else               ! Not allocated yet
                allocate(ampfda(iy,0:jx),STAT=istat)
                if (istat.ne.0) write(*,*)'impanvc0: ampfda prblm'
                allocate(ampfdd(0:iy,jx),STAT=istat)
                if (istat.ne.0) write(*,*)'impanvc0: ampfdd prblm'
             endif
             call bcast(ampfda,zero,SIZE(ampfda))
             call bcast(ampfdd,zero,SIZE(ampfdd))

             !   Set up the ampfda/ampfdd coeffs, and shift to bin bndries
             !   (Following coefedad, and above.  See, also, 131231 notes.)  
             !   (Could dimension by 1:lrz also, and set up once for all l_.)
             do j=1,jx
                do i=1,iy
                   !YuP170227                   ampfda(i,j)=bnumb(k)/fmass(k)*cex(i,j,l_)
                   !YuP170227                   ampfdd(i,j)=bnumb(k)/fmass(k)*cet(i,j,l_)
                   !BH170228                   ampfda(i,j)=bnumb(k)/fmass(k)*cex(i,j,l_)*300.d0
                   !BH170228                   ampfdd(i,j)=bnumb(k)/fmass(k)*cet(i,j,l_)*300.d0
                   !YuP170228                   ampfda(i,j)=bnumb(k)/fmass(k)*cex(i,j,l_)
                   !YuP170228                   ampfdd(i,j)=bnumb(k)/fmass(k)*cet(i,j,l_)
                   ampfda(i,j)=bnumb(k)/fmass(k)*cex(i,j,l_)*300.d0
                   ampfdd(i,j)=bnumb(k)/fmass(k)*cet(i,j,l_)*300.d0
                   ! Note: cex and cet coeffs are 0 in the trapped cone.
                enddo
             enddo
             call coefmidv(ampfda,1)
             call coefmidt(ampfdd,1)


             ipassbnd=iyh+1-itl
             if (symtrap.ne."enabled") ipassbnd=0
             do j=1,jx
                ihalf=0
                icntsww=1
                if (iyh.eq.itl) icntsww=0
                if (symtrap .ne. "enabled") icntsww=itl-iyh
                do indexsww=1,inew
                   ieq=ieq+1
                   !yup                   if (j .eq. 1 .and. (lbdry(k).ne."conserv" 
                   !yup     +                  .and. lbdry(k).ne."consscal")) then
                   if(j.eq.1)then !YuP: for 'g' function, should be true for any lbdry
                      !BH170228                      rhs(ieq)=f_(iyh,1,k,l_) ! YuP: Why iyh?  BH: all same
                      !rhs_old=rhs(ieq)
                      rhs(ieq)=0.d0 !for 'g' function, it should be true for any lbdry
                      !write(*,*)'l_,j,i,rhs=',l_,j,indexsww,rhs_old
                   else ! YuP: j>1, any lbdry
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

                      !For scal(neq,l_) :
                      neq=ieq+(k-1)*iyjx   ! for multiple species (k.gt.1)

                      alambda= vptb(i,lr_) ! ZOW or Hybrid-FOW

                      !rhs_old=rhs(ieq)

                      fh0p=fh(i,j+1,1,l_)+fh(i,j,1,l_)
                      fh0m=fh(i,j,1,l_)+fh(i,j-1,1,l_)
                      fhp0=fh(i+1,j,1,l_)+fh(i,j,1,l_)
                      fhm0=fh(i,j,1,l_)+fh(i-1,j,1,l_)
                      !BH170407  Following had very small effect on a test (see a_change.h),
                      !BH170407  and have concern it could affect stability of the solution.
                      !BH170407                      !BH170406
                      !BH170407                      !Adding nonlinear (NL) term using previous iteration fg:
                      !BH170407                      if (it_ampf.gt.1) then
                      !BH170407                         fh0p=fh0p+delecfld0n(l_,n,it_ampf-1)*
                      !BH170407    1                        (fg(i,j+1,1,l_)+fg(i,j,1,l_))
                      !BH170407                         fh0m=fh0m+delecfld0n(l_,n,it_ampf-1)*
                      !BH170407     1                        (fg(i,j,1,l_)+fg(i,j-1,1,l_))
                      !BH170407                         fhp0=fhp0+delecfld0n(l_,n,it_ampf-1)*
                      !BH170407     1                        (fg(i+1,j,1,l_)+fg(i,j,1,l_))
                      !BH170407                         fhm0=fh0m+delecfld0n(l_,n,it_ampf-1)*
                      !BH170407     1                        (fg(i,j,1,l_)+fg(i-1,j,1,l_))
                      !BH170407                      endif

                      rhs(ieq)=half*( (ampfda(i,j)*fh0p - ampfda(i,j-1)*fh0m) /cint2(j) &
                      + (ampfdd(i,j)*fhp0 - ampfdd(i-1,j)*fhm0) &
                                !YuP170227     +                     /(xsq(j)*sinn(i,l_)*dy(i,l_))) &
                           /(xsq(j)*cynt2(i,l_)/twopi))

                      rhs(ieq)=rhs(ieq)*scal(neq,l_)
                      !                      write(*,'(a,3i4,4e12.3)')'l_,j,i,rhs=',
                      !     +                l_,j,i,rhs_old,rhs(ieq),ampfda(i,j),ampfdd(i,j)

                   endif        ! On if j.eq.1
                enddo           ! On indexsww
             enddo              ! On j
          endif                 ! On ampfmod/kopt

       else  ! factor.ne."disabled", i.e. FACTORIZE

          !.......................................................................
          !
          !     reinitialize to zero matrix and rhs
          !
          !          write(*,*)'impavnc0:size(rhs) ',size(rhs)
          if (soln_method.eq.'direct' .or. soln_method.eq.'itsol1') then
             !XXX call bcast(abd,0,lenabd)
             abd = 0.
             do i=1,inewjx_(l_)
                rhs(i)=0.0
                ipivot(i)=0
             enddo
          elseif (soln_method.eq.'itsol') then
             do i=1,inewjx_(l_)
                rhs(i)=0.0
             enddo
          elseif (ifirst_lr.eq.1) then  ! for it3dv and it3drv
             do i=1,ieq_tot
                rhs(i)=0.0
             enddo
          endif

          !..................................................................
          !     The next set of variables are defined for use in code for
          !     performing the second compression.  The value I is defined as
          !     in the past and IEQ is incremented by unity each iteration.
          !     However, due to the remapping of equation (node) number to
          !     (i,j) pair of coordinates, the value of I "jumps" around.
          !
          !     As the index IEQ increases (for a given value of J), we start
          !     at I=IYH (iy/2) and move clockwise to I=itl-1. Then we alternate
          !     between the left and right sides, moving clockwise on the right
          !     side and CCW on the left side.  As an example, for the case
          !     of iy=16, itl=5, a sequence of inew=12 (12=iyh + itl-1) terms is
          !     repeated for each value of J (n=12*(j-1))
          !
          !     ieq=n+1, n+2,  n+3,  n+4, n+5, n+6 , n+7 , n+8 , n+9 , n+10, n+11, n+12
          !     i =iyh,iyh-1,itl+1,itl,itl-1,itu+1,itl-2,itu+2,itl-3,itu+3,itl-4,itu+4
          !     (in this case, iyh-2=6=itl+1, itu=iy+1-itl=12)
          !
          !     ipassbnd is the number of pitch-angle-points from iyh to itl (i.e.,
          !     number of points to trapped-passing boundary at itl, including iyh).
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

          !BH071029           ipassbnd=iyh+1-itl
          !BH071029           if (symtrap.ne."enabled" .or. nadjoint.eq.1) ipassbnd=0
          !BH070527          ieq=0
          !          write(*,*)'impavnc0, before j and i loops: ieq,icoeff=',
          !     +                                               ieq,icoeff

          !%OS  
          !          zzzto=0.0
          !          zzz12=second()
          !%OS  

          !..................................................................
          !     do loop 12 is the major loop over j (momentum).
          !..................................................................

          zrhsj1=0.0
          !BH040719          call bcast(zavarj1,zero,mu+1)
          !BH040719          call ibcast(ijaj1,0,mu+1)
          call bcast(zavarj1,zero,iy+4)
          call ibcast(ijaj1,0,iy+4)

          do 12 j=1,jx

             icntsww=1
             if (itl.eq.iyh) icntsww=0
             if (symtrap .ne. "enabled") icntsww=itl-iyh
             ihalf=0
             if (nadjoint.eq.1 .and. symtrap.eq."enabled") then
                icntsww = 0
                ihalf = 1
             endif
             idistl = -1
             idistr = 1
             if (symtrap .eq. "enabled") then
                idist_itl= 2
                idist_itu= 0 ! not used in this case         
             elseif (symtrap .eq. "disabled") then
                idist_itl= 3
                idist_itu= 1
             else ! "fow" (not happening here)
                idist_itl= itu-itl+1
                idist_itu=-idist_itl
             endif

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
                !     For lbdry(k).eq. ("conserv" or "consscal"),
                !     the j=1 conservation equations give larger iy bandwidth,
                !     linking j=1,i=inew to all j=2 distribution values (i=1,inew).
                !..................................................................
                alambda= vptb(i,lr_) ! ZOW

                ieq=ieq+1
                ibandpieq=ibandto + ieq

                if(lbdry(k).eq."fixed" .or.  &
                     lbdry(k).eq."scale" .or. lbdry(k).eq."scale_n0")then
                   ! YuP[07-24-2014] Put j=1 case BEFORE itl/itu (goto31)
                   ! There is no trap-pass bndry for j=1
                   if (j .eq. 1) go to 2
                endif

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
                !     Non-boundary region, nine point scheme.
                !     Most of (i,j) points (internal)
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

998             continue

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
                !     Impose all points equal to f(iy,1,k,l_) and add equations
                !     at all i onto eq. for i=iy.
                !     First consider case where f is fixed at v=0 for
                !     (lbdry.ne.("conserv" or "consscal")):
                !     Do not fill in matrix, will be done after "5000 continue"
                !..................................................................

2               continue ! j=1 (v=0) point (totally iy, or "inew" points)
                if (lbdry(k).eq."fixed") then
                   janelt(1)=ieq
                   nvar=1
                   zmatcont(1)=1.
                   rhs(ieq)=f_(iyh,1,k,l_) ! f_ at previous t-step
                   xnorm=1.
                   go to 5000
                endif

                if(lbdry(k).eq."scale")then
                   ! Set f(j=1) to mean of f_(j=2,i) (mean over i).
                   ! Notice: variation of f() along pitch-angle
                   ! can be factor of 2.
                   ! This is consistent with assumption that the new f() 
                   ! should be nearly Maxwellian at v~0, so that df/dv=0 at v=0.
                   ! The overall magnitude of distribution will be set
                   ! by rescaling to a target density (see diagscal)
                   !-> First, find the mean value of f_(i) at j=2:
                   !f_max=0.d0
                   !f_min=1.d100
                   !f_mean=0.d0
                   !i_mean=0
                   !do im=1,iy_(l_)
                   !   f_max=  max(f_max,f_(im,2,k,l_))
                   !   f_min=  min(f_min,f_(im,2,k,l_))
                   !   f_mean= f_mean + f_(im,2,k,l_)
                   !   i_mean= i_mean+1
                   !enddo
                   !f_mean=f_mean/float(i_mean)
                   ! Add to the matrix: f(j=1)=rhs=f_mean(j=2) :
                   janelt(1)=ieq
                   nvar=1
                   zmatcont(1)=1.
                   xnorm=1.
                   rhs(ieq)=f_(iyh,1,k,l_) !original (same as mean(f_(:,j=1,k,l_))
                   !rhs(ieq)=f_(iyh,2,k,l_) ! j=2 at prev.step ! gives neg.dens.
                   !rhs(ieq)= f_max !unstable
                   !rhs(ieq)= f_min ! results in neg.dens.
                   !rhs(ieq)= f_mean !unstable: jiggles in n(rho), enrgy(rho)

                   ! Another version: impose df/dv =0 (for each i).
                   ! Differencing corresponds to f(j=1,i)-f(j=2,i) =0.
                   !                rhs(ieq)=0.d0
                   !                janelt(1)=ieq      !corresponding to f(j=1) (for given i)
                   !                janelt(2)=ieq+inew !corresponding to f(j=2) (and same i)
                   !                nvar=2
                   !                zmatcont(1)=+0.5
                   !                zmatcont(2)=-0.5
                   !                call impnorm(xnorm,zmatcont,rhs(ieq),nvar)
                   ! Tests: not stable in runs with lossmode+bootstrap.
                   ! In runs with only lossmode (no bootstrap) - ok.

                   go to 5000
                endif ! "scale"

                !..................................................................
                ! "conserv"   Now the conservative boundary condition. i=1 first...
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
                   !                if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !                write(*,*)'HERE0.1: i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                   !                if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !                write(*,*)'HERE0.2:i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                   !                if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !                write(*,*)'HERE0.3: i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                   !                if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !                write(*,*)'HERE0.4: i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                   !                if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !                write(*,*)'HERE0.5: i,j,icntsww,xnorm',i,j,icntsww,xnorm
                   go to 5000
                endif

                !.................................................................
                !     i=1 case (theta=0) : j=2,3,4.....,jx-1
                !.................................................................

3               continue
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
4               continue

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

5               continue

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

1               continue

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

991             continue

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

7               continue

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

14              continue

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
13              continue

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

31              continue

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

5000            continue

                !.......................................................................
                !     SOLN_METHOD='direct'/'itsol1': put normalized contribs into matrix
                !.......................................................................

                if (soln_method.eq.'direct'.or.soln_method.eq.'itsol1') then

                   !BH050804 bugfix, adding following if clause:
                   if(lbdry(k).ne."conserv" .and. lbdry(k).ne."consscal")then

                      do jcont=1,nvar
                         abd(ibandpieq-janelt(jcont),janelt(jcont))= &
                              zmatcont(jcont)
                      end do

                   else  !i.e.,lbdry(k).eq. ("conserv" or "consscal")


                      if (j .ne. 1) then

                         do jcont=1,nvar
                            abd(ibandpieq-janelt(jcont),janelt(jcont))= &
                                 zmatcont(jcont)
                         end do

                         !BH050804 bugfix:              else
                         !BH070419:             elseif (lbdry(k).eq."conserv") then
                         !BH070419: Wonder what I was thinking?

                      else  !i.e., j=1

                         !.......................................................................
                         !     for j=1 impose f(i)=cst condition
                         !.......................................................................

                         do 5100 icon=1,nvar
                            if (janelt(icon) .le. inew) then
                               !     (i,j=1) contributions -> (iy,1) point
                               zavarj1(1)=zavarj1(1)+zmatcont(icon)*xnorm
                            else
                               !     (i,j=2) contributions -> (i,2) point
                               zavarj1(janelt(icon)-inew+1)= &
                                    zavarj1(janelt(icon)-inew+1)+zmatcont(icon)*xnorm
                               ijaj1(janelt(icon)-inew+1)=janelt(icon)
                            endif
5100                     end do
                         zrhsj1=zrhsj1+rhs(ieq)*xnorm

                         !     modify matrix and right-hand-side
                         !     assumes point iy is the last of the j=1 sequence
                         !
                         if (ieq .lt. inew) then
                            !     i<iy,j=1: f(i)-f(iy)=0.
                            nvar=2
                            zmatcont(1)=0.5
                            zmatcont(2)=-0.5
                            janelt(1)=ieq
                            janelt(2)=inew
                            rhs(ieq)=0.0
                            xnorm=2.0
                            do jcont=1,nvar
                               abd(ibandpieq-janelt(jcont),janelt(jcont))= &
                                    zmatcont(jcont)
                            end do
                         else
                            !     i=iy
                            !   Replace equation by sum of equations for i=1,iy, 
                            !   using f(i)=f(iy), i=1,iy-1.
                            !   Thus, number of contributions should be at most mu to the  
                            !   right of (i=iy,j=1) point => inew+1.
                            !   Bug fix, Olivier Sauter, 030509: right of (i=iy,j=1) point => mu+1.
                            ijaj1(1)=ieq
                            !   Bug fix, Olivier Sauter, 030509:     nvar=mu + 1
                            nvar=inew + 1
                            rhs(ieq)=zrhsj1
                            call impnorm(xnorm,zavarj1,rhs(ieq),nvar)
                            do jcont=1,nvar
                               abd(ibandpieq-ijaj1(jcont),ijaj1(jcont))= &
                                    zavarj1(jcont)
                            end do
                         endif  ! on ieq.lt.new/else

                      endif  ! on j.ne.1/else

                   endif  ! on lbdry(k)

                endif  ! on soln_method.eq.'direct'
                !   .or.soln_method.eq.'itsol1'

                !.......................................................................

                !.......................................................................
                !      SOLN_METHOD='itsol',  put normalized contribs into matrix
                !.......................................................................

                if (soln_method.eq.'itsol') then

                   ia_csr(ieq)=icoeff+1   ! Coeff counter at beginning of row

                   !BH050804 bugfix, adding following if clause:
                   if(lbdry(k).ne."conserv" .and. lbdry(k).ne."consscal")then

                      do jcont=1,nvar
                         icoeff=icoeff+1
                         a_csr(icoeff)=zmatcont(jcont)
                         ja_csr(icoeff)=janelt(jcont)
                      end do

                   else   !i.e., lbdry(k).eq.('conserv' or 'consscal')


                      if (j .ne. 1) then

                         do jcont=1,nvar
                            icoeff=icoeff+1
                            a_csr(icoeff)=zmatcont(jcont)
                            ja_csr(icoeff)=janelt(jcont)
                         end do

                         !BH050804 bugfix:              else
                         !BH070419:             elseif (lbdry(k).eq."conserv") then
                         !BH070419: Wonder what I was thinking?

                      else   !  i.e., j=1
                         !.......................................................................
                         !     for j=1 impose f(i)=cst condition
                         !.......................................................................

                         do 5200 icon=1,nvar
                            if (janelt(icon) .le. inew) then
                               !     (i,j=1) contributions -> (iy,1) point
                               zavarj1(1)=zavarj1(1)+zmatcont(icon)*xnorm
                            else
                               !     (i,j=2) contributions -> (i,2) point
                               zavarj1(janelt(icon)-inew+1)= &
                                    zavarj1(janelt(icon)-inew+1)+zmatcont(icon)*xnorm
                               ijaj1(janelt(icon)-inew+1)=janelt(icon)
                            endif
5200                     end do
                         zrhsj1=zrhsj1+rhs(ieq)*xnorm

                         !     modify matrix and right-hand-side
                         !     assumes point iy is the last of the j=1 sequence
                         !
                         if (ieq .lt. inew) then
                            !     i<iy,j=1: f(i)-f(iy)=0.
                            nvar=2
                            zmatcont(1)=0.5
                            zmatcont(2)=-0.5
                            janelt(1)=ieq
                            janelt(2)=inew
                            rhs(ieq)=0.0
                            xnorm=2.0
                            do jcont=1,nvar
                               icoeff=icoeff+1
                               a_csr(icoeff)=zmatcont(jcont)
                               ja_csr(icoeff)=janelt(jcont)
                            end do
                         else
                            !     i=iy
                            !   Replace equation by sum of equations for i=1,iy, 
                            !   using f(i)=f(iy), i=1,iy-1.
                            !   Thus, number of contributions should be at most mu to the  
                            !   right of (i=iy,j=1) point => inew+1.
                            !   Bug fix, Olivier Sauter, 030509: right of (i=iy,j=1) point => mu+1.
                            ijaj1(1)=ieq
                            !   Bug fix, Olivier Sauter, 030509:     nvar=mu + 1
                            nvar=inew + 1
                            rhs(ieq)=zrhsj1
                            call impnorm(xnorm,zavarj1,rhs(ieq),nvar)
                            do jcont=1,nvar
                               icoeff=icoeff+1
                               a_csr(icoeff)=zavarj1(jcont)
                               ja_csr(icoeff)=ijaj1(jcont)
                            end do
                         endif  ! on ieq

                      endif  ! on j

                   endif  ! on lbdry(k)

                endif  ! on soln_method.eq.'itsol'

                !.......................................................................
                !      SOLN_METHOD='it3dv' or 'it3drv', put normed contribs into matrix
                !.......................................................................

                if ( soln_method.eq.'it3dv' .or.  &
                     soln_method.eq.'it3drv') then

                   ia_csr(ieq)=icoeff+1   ! Coeff counter at beginning of row

                   !BH050804 bugfix, adding following if clause:
                   if(lbdry(k).ne."conserv" .and. lbdry(k).ne."consscal")then

                      do jcont=1,nvar
                         icoeff=icoeff+1
                         a_csr(icoeff)=zmatcont(jcont)
                         ja_csr(icoeff)=janelt(jcont)
                      end do

                   else   !i.e., lbdry(k).eq.("conserv".or."consscal")


                      if (j .ne. 1) then

                         do jcont=1,nvar
                            icoeff=icoeff+1
                            a_csr(icoeff)=zmatcont(jcont)
                            ja_csr(icoeff)=janelt(jcont)
                         end do

                         !BH050804 bugfix:              else
                         !BH070419:             elseif (lbdry(k).eq."conserv") then
                         !BH070419: Wonder what I was thinking?

                      else   !  i.e., j=1
                         !.......................................................................
                         !     for j=1 impose f(i)=cst condition
                         !.......................................................................

                         do 5300 icon=1,nvar
                            if ((janelt(icon)-(ieq_(l_)-1)) .le. inew) then
                               !     (i,j=1) contributions -> (iy,1) point
                               zavarj1(1)=zavarj1(1)+zmatcont(icon)*xnorm
                            else
                               !     (i,j=2) contributions -> (i,2) point
                               zavarj1((janelt(icon)-(ieq_(l_)-1))-inew+1)= &
                                    zavarj1((janelt(icon)-(ieq_(l_)-1))-inew+1)+ &
                                    zmatcont(icon)*xnorm
                               ijaj1((janelt(icon)-(ieq_(l_)-1))-inew+1)= &
                                    janelt(icon)-(ieq_(l_)-1)
                            endif
5300                     end do
                         zrhsj1=zrhsj1+rhs(ieq)*xnorm

                         !     modify matrix and right-hand-side
                         !     assumes point iy is the last of the j=1 sequence
                         !
                         if ((ieq-(ieq_(l_)-1)).lt. inew) then
                            !     i<iy,j=1: f(i)-f(iy)=0.
                            nvar=2
                            zmatcont(1)=0.5
                            zmatcont(2)=-0.5
                            janelt(1)=ieq
                            janelt(2)=inew+(ieq_(l_)-1)
                            rhs(ieq)=0.0
                            xnorm=2.0
                            do jcont=1,nvar
                               icoeff=icoeff+1
                               a_csr(icoeff)=zmatcont(jcont)
                               ja_csr(icoeff)=janelt(jcont)
                            end do
                         else
                            !     i=iy
                            !   Replace equation by sum of equations for i=1,iy, 
                            !   using f(i)=f(iy), i=1,iy-1.
                            !   Thus, number of contributions should be at most mu to the  
                            !   right of (i=iy,j=1) point => inew+1.
                            !   Bug fix, Olivier Sauter, 030509: right of (i=iy,j=1) point => mu+1.
                            ijaj1(1)=ieq - (ieq_(l_)-1)
                            !   Bug fix, Olivier Sauter, 030509:     nvar=mu + 1
                            nvar=inew + 1
                            rhs(ieq)=zrhsj1
                            call impnorm(xnorm,zavarj1,rhs(ieq),nvar)
                            do jcont=1,nvar
                               icoeff=icoeff+1
                               a_csr(icoeff)=zavarj1(jcont)
                               ja_csr(icoeff)=ijaj1(jcont) + (ieq_(l_)-1)
                            end do
                         endif  ! on ieq

                      endif  ! on j

                   endif  ! on lbdry(k)

                endif  ! on soln_method.eq.'it3dv' or 'it3drv'

                !.......................................................................
                !     Save scaling for reuse by no-factorize solution or radial transp
                !.......................................................................
                if (soln_method.eq."itsol" .or. soln_method.eq."itsol1" &
                     .or. soln_method.eq."direct") then

                   kku=ieq+(k-1)*iyjx  !-YuP-> version from impavnc
                else ! for soln_method it3dv or it3drv
                   kku=ieq-(ieq_(l_)-1)+(k-1)*iyjx 
                endif
                scal(kku,l_)=1./xnorm
                !              if(kku.le.0 .or. kku.gt.iyjx*ngen) then
                !                write(*,*)'impavnc0: k,l_,ieq,ieq_(l_),kku,ieq_tot=',
                !     ~                          k,l_,ieq,ieq_(l_),kku,ieq_tot
                !                pause
                !              endif
                
                !..................................................................
                !     End of major j and I loops..
                !..................................................................
                !              write(*,*)'ij-loop-end: ieq,rhs(ieq) ',ieq,rhs(ieq)
11              continue
12              continue
                
                !..................................................................
                !     Final element for CSR storage
                !..................................................................
                
                if ( soln_method.ne.'direct' ) then
                   ia_csr(ieq+1)=icoeff+1  !CSR storage includes coeff count+1
                   !recorded at last eqn+1.
                   
                   !         Estimating number of coeffs per flux surface
                   !         to compare with actual number:
                   !          write(*,*)'impanvc0: Number of CSR coeffs: icoeff=',icoeff
                   !         icoeff_est=(iy+2)*3+(jx-1)*9*((iyh-itl)+2*itl)  Calc'd above
                   !          write(*,*)'impavnc0: Estimate of of icoeff=', icoeff_est
                   !         Check coeff storage:
                   if (icoeff.gt.size(a_csr)) then
                      WRITE(*,*)'impavnc0:icoeff.gt.size(a_csr)'
                      STOP
                   endif
                endif  ! on soln_method.ne.'direct'
                
                !     Check number of equations:
                !        if (soln_method.eq.'itsol' .or. soln_method.eq.'itsol1') then
                !           write(*,*)'impavnc0: ieq should equal inewjx: ',ieq,inewjx
                !        elseif(soln_method.eq.'it3dv' .or. soln_method.eq.'it3drv') then
                !           write(*,*)'impavnc0: ieq should equal ieq_tot: ',ieq,ieq_tot
                !           pause
                !        endif
                
                
                if (soln_method.eq.'itsol') then
                   iunit=39
                   i1=1
                   i2=icoeff
                   !yup      if (n.eq.1 .and. ifirst_lr.eq.1) open(iunit)
                   !yup      write(iunit,*)'a_csr(i1:i2), l_,i1,i2 ',l_,i1,i2
                   !yup      write(iunit,110) a_csr(i1:i2)
                   !yup      write(iunit,*)'ja_csr(i1:i2), l_,i1,i2 ',l_,i1,i2
                   !yup      write(iunit,111) ja_csr(i1:i2)
                   !yup      write(iunit,*)'ia_csr(1:ieq), l_,ieq ',l_,ieq
                   !yup      write(iunit,111) ia_csr(1:ieq)
                   !yup      write(iunit,*)'rhs(1:ieq), l_,ieq ',l_,ieq 
                   !yup      write(iunit,110) rhs(1:ieq)
                   !yup      if (n.eq.1. and. ilast_lr.eq.1) close(iunit)
                endif  ! on solution_meth=itsol
                
                !..................................................................
                !     Next endif is the end for "if(factor.eq."disabled") then".
                !     Thus matrix and rhs are defined
                !..................................................................
                
             endif     ! on factor
             
             !        write(*,*)'impavnc0:  factor =',factor
             
             
             
             if (soln_method.eq.'itsol' .or. soln_method.eq.'itsol1') then
                n_rows_A=inewjx
             elseif (soln_method.eq.'it3dv'.or. soln_method.eq.'it3drv') then
                n_rows_A=ieq_tot
             endif
             
             
             !..................................................................
             !     SOLN_METHOD = 'it3drv':
             !     Calculate the CSR matrix of coeffs
             !     Add the coefficients into the above velocity space coeffs
             !..................................................................
             
             if (soln_method.eq.'it3drv' .and. ilast_lr.eq.1) then
                call tdtranspn
                
                ! Use SPARSKIT (in BLASSM/blassm.f) CSR addition routine aplb,
                !      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
                !     *     c,jc,ic,nzmax,iw,ierr)
                ! Computes matrix sum:    C = A+B, in CSR format
                ! nzmax	= integer. The  length of the arrays c and jc.
                !         aplb will stop if the result matrix C  has a number 
                !         of elements that exceeds nzmax. See ierr.
                ! ierr	= integer. serving as error message. 
                !         ierr = 0 means normal return,
                !         ierr .gt. 0 means that aplb stopped while computing the
                !         i-th row  of C with i=ierr, because the number 
                !         of elements in C exceeds nzmax.
                ! work arrays:
                ! iw	= integer work array of length equal to the number of
                !         columns in A.
                
                ieqp=54
                i1=1
                i2=iar_csr(ieqp)
                
                call aplb(ieq_tot,ieq_tot,1,a_csr,ja_csr,ia_csr, &
                     ar_csr,jar_csr,iar_csr,ac_csr,jac_csr,iac_csr, &
                     icsrijc,jw_ilu,ierr)
                
                !MPIINSERT_IF_RANK_EQ_0
                WRITE(*,*)'impavnc0 aft.aplb: ieq_tot, iyjx*lrz, ierr', &
                     ieq_tot,iyjx*lrz,ierr
                !MPIINSERT_ENDIF_RANK
                
                ieqp=54
                i1=1
                i2=iac_csr(ieqp)
                
                !BH070419:   Add for iterative sparse matrix solve, down to
                !BH070419:   Call direct solve Gaussian..... 
                !..................................................................
                !     This is for soln_method='itsol1'.
                !     [This method uses abd_lapack which is ml+mu+1 wide, and
                !     uses much more storage than necessary.   The preferred
                !     method is soln_method='itsol' which directly populates 
                !     the CSR coeff matrix.]
                !
                !     Put matrix into SPARSKIT format, i.e., compressed sparse row
                !     storage, CSR.  This is done with SPARSKIT2 routine bndscr,
                !     which converts from LAPACK storage format.
                !     In the present case, their are 1+ml leading rows in abd(,)
                !     which have no input values.   LAPACK storage format
                !     does not include these dummy rows.
                !     Create Harwell/Boeing formatted matrix for plotting, and
                !     call SPARSKIT2 routine which plots outline of matrix entries.
                !..................................................................
                
                !NOTE:  There are a lot of 0.0d0 entries in abd_lapack
                !NOTE:  which we can eliminate.   Best would be to set up
                !NOTE:  a_csr,ja_csr,ia_csr directly from the above
                !NOTE:  coefficient calculation.   Then would have a much
                !NOTE:  smaller sparse coeff set, and consequently much
                !NOTE:  faster execution.
                !       [Now done with soln_method='itsol']
                
                !       i_orig,j_orig refer to original A(,) coeff array
                !       before adjustment to lapack order.
                
             endif  ! on soln_method.eq.'it3drv' and ilast_lr
             
             
             if (soln_method.eq.'itsol1') then
                
                icount=0
                do j=1,inewjx
                   do i=1,ml+mu+1
                      abd_lapack(i,j)=abd(ml+i,j)
                   enddo
                enddo
                
                nabd=ml+mu+1
                lowd=ml+mu+1
                ml_a=ml
                mu_a=mu
                call bndcsr(n_rows_A,abd_lapack(1,1),nabd,lowd,ml_a,mu_a, &
                     a_csr(1),ja_csr(1),ia_csr(1),icsrij,ierr_csr)
                if (ierr_csr.ne.0) then
                   WRITE(*,*)'impavnc0/bndcsr: STOP ierr_csr=',ierr_csr
                   stop
                endif
                
             endif  !  on soln_method.eq.'itsol1'
             
             
             if (soln_method.eq.'itsol' .or. soln_method.eq.'itsol1') then 
110             format(10(1pe10.2))
111             format(10i10)
             endif  ! on soln_method itsol or itsol1
             
             
             if ( (soln_method.eq.'it3dv')  .and.  &
                  ilast_lr.eq.1 .and. n.eq.1  ) then
                !        write out coeffs in a format which can be compared directly 
                !        with soln_method=itsol output.
                iunit=39
                i1=0  !beginning of coeff print for given surface
                i2=0  !end of coeff print for given surface
                l1=1  !eqn number at beginning of given flux surface
                l2=1  !eqn number at beginning of next flux surface
                do ll=1,lrz
                   i1=i2+1
                   l1=l2
                   if (ll.lt.lrz) then
                      l2=ieq_(ll+1)
                      i2=ia_csr(l2)-1
                   else
                      l2=ieq_tot+1
                      i2=ia_csr(l2)-1
                   endif
                enddo
             endif  ! on soln_method .eq. it3dv
             
             if (soln_method.eq.'itsol' .or. soln_method.eq.'itsol1' .or. &
                  soln_method.eq.'it3dv'.or.soln_method.eq.'it3drv') then 
                
                !     Call ilut preconditioner, following SPARSKIT:rilut.f test
                !      subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
                
                ierr=0
                lfil0=min(lfil,n_rows_A)
                
                !-YuP        write (*,*) ' ++ ILUT Preconditioner ++++ '
                !-YuP        write (*,*) ' ++ droptol =, lfil.le.n_rows_A  ++ ',droptol,lfil0
                
                call cpu_time(tm1)
                !BH070523:  Change array input to be first element,
                !BH070523:  per compiler results for lf95 by Ershov:       
                !        call ilut (n_rows_A,a_csr,ja_csr,ia_csr,lfil0,droptol,
                !     1       alu,jlu,ju,iwk_ilu,w_ilu,jw_ilu,ierr)
                ! On return:
                ! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
                !           the L and U factors together.
                ! ju      = integer array of length n containing the pointers to
                !           the beginning of each row of U in the matrix alu,jlu.
                ! ierr    = integer. Error messages.  0 is succesful return.
                !
                !
                !  For soln_method .eq. 'it3dv' or 'it3drv', only call at end
                !  of coeff setup over ll=1,lrz flux surfaces.
                ! Additional arguments used in ilutp:
                permtol= 0.5d0 ! 0.d0 -> means never permute
                mbloc=   n_rows_A
                !iperm_ilu(2*n) is output.
                ! If permtol=0. and mbloc=n_rows_A, 
                ! the result from ilutp() is identical to result from ilut()
                ! 
                
                if ( soln_method.eq."itsol" .or. soln_method.eq."itsol1" ) then
                   call ilut (n_rows_A,a_csr(1),ja_csr(1),ia_csr(1),lfil0, &
                        droptol,alu(1),jlu(1),ju(1),iwk_ilu,w_ilu(1), &
                        jw_ilu(1),ierr)
                   !Setup up above, since vv size dependency
                   do i=1,n_rows_A
                      rhs0(i)=rhs(i)    !Copy rhs, for input to pgmres since
                      !pgmres modifies this input.
                   enddo
                elseif(soln_method.eq.'it3dv') then
                   ! perform soln only at last flux surface (l_=lrz)
                   if ( ilast_lr.eq.1 ) then
                      call ilut (n_rows_A,a_csr(1),ja_csr(1),ia_csr(1),lfil0, &
                           droptol,alu(1),jlu(1),ju(1),iwk_ilu,w_ilu(1), &
                           jw_ilu(1),ierr)
                      !              call ilutp(n_rows_A,a_csr(1),ja_csr(1),ia_csr(1),
                      !     +             lfil0,droptol,permtol,mbloc,
                      !     +             alu(1),jlu(1),ju(1),iwk_ilu,w_ilu(1),jw_ilu(1),
                      !     +             iperm_ilu(1),ierr)
                      !MPIINSERT_IF_RANK_EQ_0
                      !       WRITE(*,*)'impavnc0 aft.ilut:  l_,ierr',l_,ierr,soln_method
                      !MPIINSERT_ENDIF_RANK
                      do i=1,n_rows_A
                         rhs0(i)=rhs(i) !Copy rhs, for input to pgmres since
                         !this input is modified.
                      enddo
                   endif
                elseif(soln_method.eq.'it3drv') then
                   ! perform soln only at last flux surface (l_=lrz)
                   if ( ilast_lr.eq.1 ) then
                      call ilut (n_rows_A,ac_csr(1),jac_csr(1),iac_csr(1),lfil0, &
                           droptol,alu(1),jlu(1),ju(1),iwk_ilu,w_ilu(1), &
                           jw_ilu(1),ierr)
                      !              call ilutp(n_rows_A,ac_csr(1),jac_csr(1),iac_csr(1),
                      !     +             lfil0,droptol,permtol,mbloc,
                      !     +             alu(1),jlu(1),ju(1),iwk_ilu,w_ilu(1),jw_ilu(1),
                      !     +             iperm_ilu(1),ierr)
                      !MPIINSERT_IF_RANK_EQ_0
                      !       WRITE(*,*)'impavnc0 aft.ilut:  l_,ierr',l_,ierr,soln_method
                      !MPIINSERT_ENDIF_RANK
                      do i=1,n_rows_A
                         rhs0(i)=rhs(i) !Copy rhs, for input to pgmres since
                         !this input is modified.
                      enddo
                   endif
                endif

                call cpu_time(tm(1))
                tm1= tm(1) - tm1
                !        write(*,*)'impavnc0: icount_imp, time for ilut tm1=',
                !     +                        icount_imp,tm1
                if (ierr.ne.0) then
                   WRITE(*,*)'impavnc0 after ilut: ierr=',ierr
                   STOP 'if ierr=-2 or -3, reduce lfil or increase iwk_ilu'
                   ! Try to reduce lfil (say, to 10) or increase iwk_ilu
                endif


                !      Iterative solve using:
                !       subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout,
                !     *                    aa, ja, ia, alu, jlu, ju, ierr)
                
                krylov1=krylov        !Size of kylov subspace, whatever that is?
                
                !.......................................................................
                !     Initial soln sol for the iterative method is 
                !     soln at the start of the time step.  At the  first time
                !     step it will be either a Maxwellian or a restart distribution.
                !     
                !     For each j (velocity) value, take the first ipassbnd # of elements
                !     from f and place them in the right half of trapped
                !     region starting at iyh.
                !.......................................................................
                
                if (soln_method.eq."itsol" .or. soln_method.eq."itsol1") then
                   
                   l1=l_
                   l2=l_
                   
                elseif(soln_method.eq.'it3dv' .or. soln_method.eq.'it3drv') then
                   ! perform soln
                   ! only at last flux surface (l_=lrz).
                   if ( ilast_lr.eq.1 ) then
                      l1=1
                      l2=lrz
                   else
                      l1=0
                      l2=-1              !i.e., following ll do-loop is no-op
                   endif
                endif
                
                iii=0                   !index to point to end of eqn/soln 
                !number for ll-flux surface.
                do ll=l1,l2
                   do j=1,jx
                      i1=(j-1)*inew_(ll) + iii  !eqn/soln index at end of
                      !previous j,ll
                      ihalf=0
                      
                      if (symtrap.eq."enabled".and. nadjoint.ne.1) then
                         if (itl_(ll).lt.iyh_(ll)) then
                            inew=iyh_(ll) + itl_(ll) - 1
                            ipassbnd=iyh_(ll)-itl_(ll)+1
                         else            !ie., itl_(ll)=iyh_(ll)
                            inew=iy_(ll)
                            ipassbnd=0   !i.e., no trapped particles
                            !as in following case.
                         endif
                         if (symtrap .ne. "enabled") then
                            inew=iy_(ll)
                            ipassbnd=0
                         endif
                         
                         do  ii=1,ipassbnd
                            sol(i1+ii)=f(iyh_(ll)+1-ii,j,k,ll)
                         enddo
                         icntsww=1
                         if (itl_(ll).eq.iyh_(ll)) icntsww=0
                      else if (nadjoint.eq.1 .and. symtrap.eq."enabled") then
                         ipassbnd=0
                         icntsww=0
                         ihalf=1
                      else
                         ipassbnd=0
                         icntsww=itl_(ll)-iyh_(ll)
                      endif
                      
                      !.......................................................................
                      !     Alternating between the left and right sides, place the remaining
                      !     values (for this value of J) below ITL and above ITU.
                      !.......................................................................
                      
                      do ii=ipassbnd+1,inew_(ll)
                         if(ihalf.eq.0) then
                            sol(i1+ii)=f(itl_(ll)-icntsww,j,k,ll)
                            ihalf=1
                         else
                            sol(i1+ii)=f(itu_(ll)+icntsww,j,k,ll)
                            ihalf=0
                            icntsww=icntsww+1
                         endif
                      enddo              !  On ii
                      
                   enddo                 !  On j
                   
                   iii=iii+inewjx_(ll)  !update index to point to end of eqn/soln 
                   !number for ll-flux surface.
                   
                enddo                    !  On ll       
                
                !     Set prmres error tolerance and max iterations internally:
                epsilon=1.d-5
                maxits=400 !200
                iout=0 !38 ! 0 means no writing to fort.iout file.
                ! Writing into this file causes problems for MPI (I/O).
                ! Keep it 0, to avoid slow runs ! YuP[04-2017]
                ! Alternatively, set in your batch script:
                ! export FORT_BUFFERED=1 
                !(or the csh equivalent - setenv FORT_BUFFERED 1 ) 
                ! At NERSC, it was the default in before 2017
                ! but it caused problems with some
                ! newer compilers. 
                ! It seems its removal causes the I/O to take longer.
                
                
                call cpu_time(tm1)
                
                !--------------------------------------------------------------
                !     call PGMRES
                !--------------------------------------------------------------
                !     Put sol into rhs vector (i.e., as in soln from direct solve)
                if ( soln_method.eq."itsol" .or. soln_method.eq."itsol1" ) then
                   call pgmres(n_rows_A,krylov1,rhs0(1),sol(1),vv(1),epsilon, &
                        maxits,iout,a_csr(1),ja_csr(1),ia_csr(1),alu(1), &
                        jlu(1),ju(1),ierr)
                   rhs(1:n_rows_A)=sol(1:n_rows_A)
                elseif(soln_method.eq.'it3dv') then
                   ! perform soln 
                   ! only at last flux surface (l_=lrz).
                   if ( ilast_lr.eq.1 ) then
                      !MPIINSERT_IF_RANK_EQ_0
                      !       WRITE(*,*)'impavnc0 before pgmres:  ierr',ierr,soln_method
                      !MPIINSERT_ENDIF_RANK
                      call pgmres(n_rows_A,krylov1,rhs0(1),sol(1),vv(1),epsilon, &
                           maxits,iout,a_csr(1),ja_csr(1),ia_csr(1),alu(1), &
                           jlu(1),ju(1),ierr)
                      !MPIINSERT_IF_RANK_EQ_0
                      !       WRITE(*,*)'impavnc0  after pgmres:  ierr',ierr
                      !MPIINSERT_ENDIF_RANK
                      do inr = 1,n_rows_A
                         rhs(inr)=sol(inr) 
                         !-> rhs(1:n_rows_A)=sol(1:n_rows_A) results in stack overflow
                      enddo
                   endif
                elseif(soln_method.eq.'it3drv') then
                   ! perform soln 
                   ! only at last flux surface (l_=lrz).
                   if ( ilast_lr.eq.1 ) then
                      
                      !MPIINSERT_IF_RANK_EQ_0
                      !       WRITE(*,*)'impavnc0 before pgmres:  ierr',ierr,soln_method
                      !MPIINSERT_ENDIF_RANK
                      call pgmres(n_rows_A,krylov1,rhs0(1),sol(1),vv(1),epsilon, &
                           maxits,iout,ac_csr(1),jac_csr(1),iac_csr(1),alu(1), &
                           jlu(1),ju(1),ierr)
                      !MPIINSERT_IF_RANK_EQ_0
                      !       WRITE(*,*)'impavnc0  after pgmres:  ierr',ierr
                      !MPIINSERT_ENDIF_RANK
                      
                      do inr = 1,n_rows_A
                         rhs(inr)=sol(inr) 
                         !-> rhs(1:n_rows_A)=sol(1:n_rows_A) results in stack overflow
                      enddo
                   endif
                endif
                
                call cpu_time(tm(1))
                tm1= tm(1) - tm1
                !       write(*,*)'impavnc0: time for pgmres tm1=',tm1
                if (ierr.ne.0) then
                   !MPIINSERT_IF_RANK_EQ_0
                   WRITE(*,*)'impavnc0 after pgmres, ierr=',ierr
                   WRITE(*,*)'ierr=1 converg. not achieved in itmax iterations'
                   WRITE(*,*)' -1 Initial guess seems to be the exact solution'
                   !MPIINSERT_ENDIF_RANK
                   STOP 'ierr.ne.0, from pgmres'
                endif
                
             endif  ! on   soln_method.eq.'itsol' .or. soln_method.eq.'itsol1'
             ! .or. soln_method.eq.'it3dv' .or. soln_method.eq.'it3drv'
             
             !..................................................................
             !     Call direct solve Gaussian elimination algorithm.
             !..................................................................
             
             if (soln_method.eq.'direct') then
                !
                !%OS        zzzt1=second()
                !BH080303        tm1 = etime(tm)       
                call cpu_time(tm1)
                
                !     factorize matrix
                !_cray  64-bit-compiler uses sgbtrf (from lapack library).
                !_pc    32-bit-compiler uses dgbtrf (from lapack library).
                if (factor .ne. "disabled") then
                   !          call sgbtrf(inewjx,inewjx,ml,mu,abd,md1abd,ipivot,info)
                   !          write(*,*)'impavnc0 before dgbtrf, l_,inewjx,ml=',l_,inewjx,ml
                   ! XXX
                   call dgbtrf(inewjx,inewjx,ml,mu,abd,md1abd,ipivot,info)
                   if (info .ne. 0) then
                      print *,' warning after sgbtrf in impavnc0: info = ',info
                      stop 'impavnc0 1'
                   endif
                endif
                !       tm1 = etime(tm) - tm1
                !       write(*,*)'impavnc0: time for dgbtrf tm1=',tm1
                
                !     solve system
                !_cray  64-bit-compiler uses sgbtrs (from lapack library).
                !_pc    32-bit-compiler uses dgbtrs (from lapack library).
                
                !BH080303        tm1 = etime(tm)
                call cpu_time(tm1)
                
                inbrhs = 1
                transpose = 'n'
                !        call sgbtrs(transpose,inewjx,ml,mu,inbrhs,abd,md1abd,ipivot
                ! XXXX
                call dgbtrs(transpose,inewjx,ml,mu,inbrhs,REAL(abd),md1abd, &
                     ipivot,rhs,inewjx,info)
                
                if (info .ne. 0) then
                   print *,' warning after sgbtrs in impavnc0: info = ',info
                   stop 'impavnc0 2'
                endif
                !       tm1 = etime(tm) - tm1
                !       write(*,*)'impavnc0: time for dgbtrs tm1=',tm1
                
             endif  !on soln_method.eq.'direct'
             
             
             !     Deallocate it3d related storage
             !-YuP       if (n.eq.nstop .and. ilast_lr.eq.1) then
             !-YuP           call it3ddalloc 
             !-YuP  Problem: called during first iteration of eflditer;
             ! Should be after last nefiter iteration (which is unknown).
             ! Not really necessary to deallocate?
             !-YuP       endif
             
             !.................................................................
             !     The next small block expands the compressed solution in rhs()
             !     and replicates the mirror-imaged values in the trapping region 
             !     (i.e. copy solution rhs to f)
             !..................................................................
             
             if (soln_method.eq."itsol" .or. soln_method.eq."itsol1" &
                  .or. soln_method.eq."direct") then
                l1=l_
                l2=l_
                
             elseif(soln_method.eq.'it3dv' .or. soln_method.eq.'it3drv') then
                ! soln performed only at last flux surface
                !          write(*,*)'impavnc0: n,ilast_lr or rhs=>f',n,ilast_lr
                if ( ilast_lr.eq.1 ) then
                   l1=1
                   l2=lrz
                else
                   l1=0
                   l2=-1           !i.e., following ll do-loop is no-op
                endif
             endif  ! On soln_method
             
             iii=0                   !index to point to end of eqn/soln 
             !number for ll-flux surface.
             
             
             ! Now update f() (solution is obtained => ilast_lr=1)
             
             !       jstart = 1 ! For lbdry="conserv" or lbdry="consscal"
             jstart = 1 ! For lbdry .ne. "fixed"
             if(lbdry(k).eq."fixed") jstart=2 
             !For "fixed", f(j=1) is same as at the previous time step.
             
             do ll=l1,l2  
                do 966 j=jstart,jx  ! YuP: jstart=2 for lbdry="fixed"
                   i1=(j-1)*inew_(ll) + iii  !eqn/soln index at end of
                   !previous j,ll
                   
                   !..................................................................
                   !     For each j (velocity) value, take the first ipassbnd # of elements
                   !     from rhs and place them in the right half of the pass-trapped
                   !     region starting at iyh.
                   !..................................................................
                   
                   ihalf=0
                   
                   if (symtrap.eq."enabled".and. nadjoint.ne.1) then
                      if (itl_(ll).lt.iyh_(ll)) then
                         inew=iyh_(ll) + itl_(ll) - 1
                         ipassbnd=iyh_(ll)-itl_(ll)+1
                      else               !ie., itl_(ll)=iyh_(ll)
                         inew=iy_(ll)
                         ipassbnd=0      !i.e., no trapped particles
                         !as in following case.
                      endif
                      if (symtrap .ne. "enabled") then
                         inew=iy_(ll)
                         ipassbnd=0
                      endif
                      
                      do ii=1,ipassbnd
                         f(iyh_(ll)+1-ii,j,k,ll)=rhs(i1+ii)
                      enddo
                      
                      icntsww=1
                      if (itl_(ll).eq.iyh_(ll)) icntsww=0
                   else if (nadjoint.eq.1 .and. symtrap.eq."enabled") then
                      ipassbnd=0
                      icntsww=0
                      ihalf=1
                   else
                      ipassbnd=0
                      icntsww=itl_(ll)-iyh_(ll)
                   endif
                   
                   !..................................................................
                   !     Alternating between the left and right sides, place the remaining
                   !     values (for this value of J) below ITL and above ITU.
                   !..................................................................
                   
                   do 959 ii=ipassbnd+1,inew_(ll)
                      if(ihalf.eq.0) then
                         f(itl_(ll)-icntsww,j,k,ll)=rhs(i1+ii)
                         ihalf=1
                      else
                         f(itu_(ll)+icntsww,j,k,ll)=rhs(i1+ii)
                         ihalf=0
                         icntsww=icntsww+1
                      endif
959                end do
                   
                   !..................................................................
                   !     Mirror the right half of the pass-trapped region into the left hal
                   !..................................................................
                   
                   if (symtrap .eq. "enabled") then
                      if (nadjoint .eq. 1) f(itl_(ll),j,k,ll)=f(itu_(ll),j,k,ll)
                      do 957 i=iyh_(ll)+1,itu_(ll)
                         ii=iy_(ll)+1-i
                         f(i,j,k,ll)=f(ii,j,k,ll)
957                   end do
                   endif
                   
966             end do    !   On j
                
                !..................................................................
                !     Following if-endif fixes a nipple which grows on f at v=0
                !     in lbdry()="scale" conditions.  If namelist lbdry(k)="scale"
                !     we reset it to "consscal'.  "scale" is no longer operable.
                !..................................................................
                
                !BH-YuP: 170623
                !BH-YuP: Following Syun'ichi Shiraiwa, 170623, this adds to f(v=0),
                !BH-YuP: and gives growing distribution function and current.
                !BH-YuP: The problem evidently arose when jstart
                !BH-YuP: Syun'inchi: Probably avg is better.
                !BH-YuP: We will use pitch angle phase space avg
                !BH-YuP:         if(lbdry(k).eq."scale")then
                !BH-YuP:            ! copy MAX[f(j=2)] to f(j=1), YuP150410
                !BH-YuP:            f_j2_max= MAXVAL(f(1:iy,2,k,ll)) 
                !BH-YuP:            f(1:iy,1,k,ll)= f_j2_max
                !BH-YuP:         endif ! "scale"
                
                if (  ampfmod.eq."enabled" .and. kopt.eq.3 &
                     .and. cqlpmod.ne."enabled") then
                   if(icount_ampf.eq.2)then ! this is fg=='g' function
                      f(1:iy,1,k,ll)=0.d0 !this is 'g' function: should be 0 at j=1
                      f(1:iy,0,k,ll)=0.d0 !this is 'g' function: should be 0 at j=0
                   endif
                endif
                
                
                iii=iii+inewjx_(ll)  !update index to point to end of eqn/soln 
                !number for ll-flux surface.	
                
             enddo  !  On ll ; finished with updating f
             
             WRITE(*,'(a,2i6)') &
                  'impavnc[ZOW]: sol.found, f is updated. k,l_=', &
                  k,l_
             ! Note, with MPI, and soln_method.eq.'direct',
             ! parallelization is done over ll (l_) index.
             ! mpirank=0 is not doing calculations, only receiving data.
             
             !.......................................................................
             !     eseswtch="enabled" section, to control and calculate
             !     fluxes for quasi-neutral electrostatic field.
             !.......................................................................
             
             if (cqlpmod.eq."enabled" .and. eseswtch.eq."enabled") then
                
                if (icount_ese.eq.1) then
                   
                   !             Store new f into fh        
                   call dcopy(iyjx2*ngen,f(0:iyjx2*ngen-1,0,1,l_),1, &
                        fh(0:iyjx2*ngen-1,0,1,l_),1)
                   factor="disabled"
                   icount_ese=2
                   go to 10
                   
                elseif (icount_ese.eq.2) then
                   
                   !             g function is in f.
                   !             Calculate fluxes, restore f with f_, and return.
                   
                   flux1(l_)=fluxpar( &
                        1,x,coss(1,l_),cynt2(1,l_),cint2,temp1,iy,jx)
                   flux2(l_)=fluxpar( &
                        1,x,coss(1,l_),cynt2(1,l_),cint2,temp2,iy,jx)
                   call dcopy(iyjx2*ngen,f_(0:iyjx2*ngen-1,0,1,l_),1, &
                        f(0:iyjx2*ngen-1,0,1,l_),1)
                   
                   go to 999
                endif
                
             endif
             
             
             !.......................................................................
             !     ampfmod="enabled"/kopt=3 section, to calculate time-advance
             !     toroidal electric field according to Ampere_Faraday equations.
             !.......................................................................
             
             if (  ampfmod.eq."enabled" .and. kopt.eq.3 &
                  .and. cqlpmod.ne."enabled") then
                
                if (icount_ampf.eq.1) then
                   
                   !             Store new f into fh        
                   call dcopy(iyjx2,f(0:iyjx2-1,0,kelec,l_),1, &
                        fh(0:iyjx2-1,0,1,l_),1)
                   !YuP[05-2017] corrected dcopy in the above line:
                   !was iyjx2*ngen, but we only copy one species (kelec)
                   factor="disabled"
                   icount_ampf=2 !next, re-define rhs() for the 'g' function equations
                   go to 10
                   
                elseif (icount_ampf.eq.2) then
                   
                   !             g function is in f, store into fg
                   call dcopy(iyjx2,f(0:iyjx2-1,0,kelec,l_),1, &
                        fg(0:iyjx2-1,0,1,l_),1)
                   !YuP[05-2017] corrected dcopy in the above line:
                   !was iyjx2*ngen, but we only copy one species (kelec)
                   !BH170312:  Checking effect of zeroing correction part of distn:
                   !BH170312:              call bcast(fg(0,0,1,l_),zero,iyjx2*ngen)
                   !             Restore f with f_, and return.
                   call dcopy(iyjx2,f_(0:iyjx2-1,0,kelec,l_),1, &
                        f(0:iyjx2-1,0,kelec,l_),1)
                   !YuP[05-2017] corrected dcopy in the above line:
                   !was iyjx2*ngen, but we only copy one species (kelec)
                   go to 999
                endif
                
             endif
             
             !.................................................................
             !     Differentiate the solution and compare with r.h.s.
             !.................................................................
             if (soln_method.ne.'it3dv' .or. soln_method.ne.'it3drv') then
                if (iactst.ne."disabled") call impchk(k)
             endif
             
             !.......................................................................
             !     compute vel source term needed in transport eq. for ADI method
             !.......................................................................
             if (adimeth.eq."enabled" .and. transp.eq."enabled" &
                  .and. n.ge.nonadi) call tdtrvsou(k)
             
             ! YuP:        call diagimpd(k) !---> moved to tdchief.f
             
             !..................................................................
             !     End of "k" (species) loop.
             !..................................................................
600       end do ! k
          
          
          
          !-YuP      if(dalloc.eq."enabled") then
          !-YuP  deallocate(abd,STAT=istat)
          !-YuP  deallocate(ipivot,STAT=istat)
          !-YuP  Problem: called during first iteration of eflditer;
          ! Should be after last nefiter iteration (which is unknown).
          ! Not really necessary to deallocate?
          !-YuP      endif
          
          irstart=0
999       continue
          
          
          
          return
  end subroutine impavnc0
        
  
  
  
  !=======================================================================
  !=======================================================================      
  real*8 function z00(i,j,k)
    use param_mod
    use cqcomm_mod
    use advnce_mod
    implicit integer (i-n), real*8 (a-h,o-z)
    !cc         save  ! YuP: Not really needed

    !      CONTAINS     PG90 at GA couldn't accept this construct,
    !                   complaining that functions in advnce.h were 
    !                   being redefined.

    !     removed this construct for franklin.nerc.gov: pg compiler

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
         cl(itl-1,j-1)*dj(itl,j-1,k,l_)*eym5(itl,l_)) &
         +r2y(j)*(-de(itl-1,j)*(1.-di(itl-1,j-1,k,l_))) &
         /(2.*dx(j))

    t00l_(j)= &
         +qz(j)*( &
         -cl(itl-1,j)*dj(itl,j,k,l_)*eym5(itl,l_) &
         +cl(itl-1,j-1)*(1.-dj(itl,j-1,k,l_))*eym5(itl,l_)) &
         +r2y(j)*(dd(itl-1,j)*(1.-di(itl-1,j,k,l_)) &
         +df(itl-1,j)*eym5(itl,l_))

    t0pl_(j)=qz(j)*( &
         -cl(itl-1,j)*eym5(itl,l_)*(1.-dj(itl,j,k,l_))) &
         +r2y(j)*de(itl-1,j)/(2.*dx(j))*(1.-di(itl-1,j+1,k,l_))


    t0mu_(j)=qz(j)*( &
         -cl(itu+1,j-1)*dj(itu,j-1,k,l_)*eyp5(itu,l_)) &
         +r2y(j)*( &
         +de(itu,j)*di(itu,j-1,k,l_))/(2.*dx(j))

    t00u_(j)= &
         +qz(j)*( &
         +cl(itu+1,j)*dj(itu,j,k,l_)*eyp5(itu,l_) &
         -cl(itu+1,j-1)*(1.-dj(itu,j-1,k,l_))*eyp5(itu,l_)) &
         +r2y(j)*( &
         -dd(itu,j) &
         *di(itu,j,k,l_) &
         +df(itu,j)*eyp5(itu,l_))

    t0pu_(j)=qz(j)*( &
         +cl(itu+1,j)*(1.-dj(itu,j,k,l_))*eyp5(itu,l_)) &
         +r2y(j)*( &
         -de(itu,j)*di(itu,j+1,k,l_)/ &
         (2.*dx(j)))



    !  Test if bootstrap calc is needed (giving iboot=1)
    iboot=0
    if (bootcalc.ne."disabled" .and. &
         (i.eq.(itl-1).or.i.eq.itl.or. &
         i.eq.itu.or.i.eq.(itu+1))) iboot=1

    z00f=vptb(i,lr_)*(f_(i,j,k,l_)/dtreff+so(i,j)) + &
         spasou(i,j,k,l_)

    if (iboot.eq.1) then
       if (i.eq.(itl-1)) then
          !              itl-1 case:
          z00itl1=z00f -xpm(i,j)*bsl(j-1,k,l_)-xp0(i,j)*bsl(j,k,l_) &
               -xpp(i,j)*bsl(j+1,k,l_)
          z00t=z00itl1
       else
          !              itu+1 case:
          if (i.eq.(itu+1))then
             z00itu1=z00f -xmm(i,j)*bsu(j-1,k,l_)-xm0(i,j)*bsu(j,k,l_) &
                  -xmp(i,j)*bsu(j+1,k,l_)
             z00t=z00itu1
             !              itl (or itu)case:
          else
             z00itl=z00f - (t0ml_(j)*bsl(j-1,k,l_)+t00l_(j)*bsl(j,k,l_) &
                  +t0pl_(j)*bsl(j+1,k,l_)+t0mu_(j)*bsu(j-1,k,l_) &
                  +t00u_(j)*bsu(j,k,l_)+t0pu_(j)*bsu(j+1,k,l_))
             z00t=z00itl
          endif
       endif
       z00=z00t

    else ! iboot=0
       z00=z00f
    endif

  end function z00

  !      end of file impavnc0.f

end module impavnc0_mod
