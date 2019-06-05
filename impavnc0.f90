!
!
module impavnc0_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast
  use bsu_mod, only : bsu
  use bsl_mod, only : bsl
  use coefmidt_mod, only : coefmidt
  use coefmidv_mod, only : coefmidv
  use coefstup_mod, only : coefstup
  use coefwti_mod, only : coefwti
  use coefwtj_mod, only : coefwtj
  use cqlcomm_mod
  use esefld_mod, only : fluxpar
  use ilut_mod, only : aplb
  use ilut_mod, only : bndcsr
  use ilut_mod, only : ilut
  use ilut_mod, only : pgmres
  use impchk_mod, only : impchk
  use impnorm_mod, only : impnorm
  use r8subs_mod, only : dcopy
  use r8subs_mod, only : dgbtrf
  use r8subs_mod, only : dgbtrs
  use tdtranspn_mod, only : tdtranspn
  use tdtrvsou_mod, only : tdtrvsou

  !---END USE
  ! XXX this is probably a bug in the original code.  YuP: added 'k' as argument in "statement functions", see advnce.f90
  !use advnce_mod, only :  k => advnce_k ! use advnce_mod's advnce_k
  !use advnce_mod  !here: in impavnc0_mod. To get k?
  !integer,pointer  :: k => advnce_k

  ! gauss_, these are used in it3dalloc
  integer, public ::  inewjmax, ialloc
  ! XXX this was an integer somewhere?  ! YuP: ialloc is not used anymore.
  ! XXX vnc vs vnc0 abd etc?   YuP: impavnc should be removed from compilation; only impavnc0 is used
  real(c_double), pointer, private ::  abd(:,:)
  integer, pointer, private :: ipivot(:)
  ! ampf, these are used in it3dalloc
  real(c_double), pointer, private :: ampfda(:,:)
  real(c_double), pointer, private :: ampfdd(:,:)
  save

contains

  subroutine impavnc0(kopt)
    use param_mod
    use cqlcomm_mod
    use advnce_mod !here: in impavnc0(). To get many functions
    use r8subs_mod, only : dgbtrf, dgbtrs, dcopy
    implicit none

    real(c_double) :: alambda
    real(c_double) :: cnst
    real(c_double) :: dfdvz
    real(c_double) :: epsilon
    real(c_double) :: fh0m
    real(c_double) :: fh0p
    real(c_double) :: fhm0
    real(c_double) :: fhp0
    integer :: i
    integer :: i1
    integer :: i2
    integer :: ibandpieq
    integer :: ibandto
    integer :: icntsww
    integer :: icon
    integer :: icount
    integer :: icount_ampf
    integer :: icount_ese
    integer :: idist_itl
    integer :: idist_itu
    integer :: idistl
    integer :: idistr
    integer :: ieqp
    integer :: ierr
    integer :: ierr_csr
    integer :: ifirst_lr
    integer :: ihalf
    integer :: ii
    integer :: iii
    integer :: ilast_lr
    integer :: impchg
    integer :: inbrhs
    integer :: indexsww
    integer :: inew
    integer :: inew1
    integer :: inewjx
    integer :: info
    integer :: iout
    integer :: istat
    integer :: inr
    integer :: ipassbnd
    integer :: iunit
    integer :: j
    integer :: jcont
    integer :: jstart
    integer :: k ! species index (internal loop in k, see below)
    integer :: kku
    integer :: kopt
    integer :: krylov1
    integer :: l1
    integer :: l2
    integer :: lenabd
    integer :: lfil0
    integer :: ll
    integer :: lowd
    integer :: maxit
    integer :: maxits
    integer :: mbloc
    integer :: md1abd
    integer :: ml
    integer :: ml_a
    integer :: mu
    integer :: mu_a
    integer :: n_rows_A
    integer :: nabd
    integer :: nadjoint
    integer :: neq
    integer :: nvar
    real(c_double) :: permtol
    real(c_double) :: xnorm != 1.  !XXXX should this be initialized? YuP: xnorm is output of impnorm() 
    real(c_double) :: zrhsj1

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

    !BH080303    real(c_float) etime,tm1,tm(2) !Dclrtns for lf95 ETIME(). See man etime.
    !BH080303    Using f95 intrinsic subroutine, cpu_time
    real(c_float) :: tm1,tm(2)    !Dclrtns for lf95 ETIME(). See man etime.

    character(len=8) ::  alloc, factor, dalloc
    real(c_double) :: zavarj1(iy+4)
    integer :: ijaj1(iy+4)
    real(c_double) :: zmatcont(12) !XXXX should this be init 0? YuP: done below.
    integer :: janelt(12)
    character(len=1) :: transpose

    integer, save :: icount_imp=0
    integer, save :: inewmax, mlmax, mumax, icoeff, ieq, &
         ieq_tot, icoeff_est, icoeff_est_tot

    real(c_double) :: twopi2
    twopi2 = twopi*twopi
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
    !write(*,'(a,4i5)')'impavnc0: n,icount_imp,ifirst_lr,ilast_lr', &
    !n,icount_imp,ifirst_lr,ilast_lr

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
       inew=iy/2 + itl - 1 ! for itl.lt.iyh
       ipassbnd=iyh-itl+1
    else  !i.e., itl=iyh
       inew=iy ! for itl=iyh
       ipassbnd=0  !i.e., no trapped region for iyh=itl,
       !as in following case.
    endif
    if (symtrap .ne. "enabled") then
       inew=iy ! for symtrap.ne."enabled"
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

!MPIINSERT_IF_RANK_EQ_0      
    if(inewjmax.gt.inewjx)then !YuP[2019-04-24] added printout
      WRITE(*,*)'impavnc0: inewjmax>inewjx', inewjmax,inewjx
    endif 
!MPIINSERT_ENDIF_RANK

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

    !write(*,*)'impavnc0: iy,iyh,itl,itu,ml,mu,inew,jx,inewjx', &
    !iy,iyh,itl,itu,ml,mu,inew,jx,inewjx

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
    !write(*,*)'impavnc0:size(abd_lapack,1),size(abd_lapack,2)=',
    !size(abd_lapack,1),size(abd_lapack,2)
    !write(*,*)'impavnc0:size of a_csr,ja_csr,ia_csr,alu,jlu,ju',
    !size(a_csr),size(ja_csr),size(ia_csr),size(alu),
    !size(jlu),size(ju)
    !write(*,*)'impavnc0:size of jw_ilu,w_ilu_rhs0,sol,vv',
    !size(jw_ilu),size(w_ilu),size(rhs0),size(sol),size(vv)



    !     So, no confusion:
    if (soln_method.eq.'itsol1')  &
         call bcast(abd_lapack(1:icsrij,1),zero,icsrij)

    !write(*,*)'impavnc0: size(ja_csr),icsrip=',size(ja_csr),icsrij
    !write(*,*)'impavnc0: size(ia_csr),icsrip=',size(ia_csr),icsrip

    if (soln_method.eq.'itsol' .or. soln_method.eq.'itsol1') then
       call bcast(a_csr(1:icsrij),zero,icsrij)
       call ibcast(ja_csr(1:icsrij),0,icsrij)
       call ibcast(ia_csr(1:icsrip),0,icsrip)
    elseif (soln_method.eq.'it3dv') then
       if (ifirst_lr.eq.1) then
          call bcast(a_csr(1:icsrij),zero,icsrij)
          call ibcast(ja_csr(1:icsrij),0,icsrij)
          call ibcast(ia_csr(1:icsrip),0,icsrip)
       endif
    elseif (soln_method.eq.'it3drv') then
       if (ifirst_lr.eq.1) then
          !BH080429   Problem:  ipofi only calculated once in code, in tdtranspn.
          !BH080429            call ibcast(ipofi(1,1),0,iy*lrz)
          call bcast(a_csr(1:icsrij),zero,icsrij)
          call ibcast(ja_csr(1:icsrij),0,icsrij)
          call ibcast(ia_csr(1:icsrip),0,icsrip)
          call bcast(ar_csr(1:icsrijr),zero,icsrijr)
          call ibcast(jar_csr(1:icsrijr),0,icsrijr)
          call ibcast(iar_csr(1:icsrip),0,icsrip)
          call bcast(ac_csr(1:icsrijc),zero,icsrijc)
          call ibcast(jac_csr(1:icsrijc),0,icsrijc)
          call ibcast(iac_csr(1:icsrip),0,icsrip)
       endif
    endif

    !     Copy current distribution into f_..   YuP: moved to tdchief:
    !cc      call dcopy(iyjx2*ngen,f(0,0,1,l_),1,f_(0,0,1,l_),1)

    !write(*,*) 'impavnc0, l_,elecfld:', l_ ,elecfld(l_)  
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
          abd=zero !YuP[2019-04-24] was absent
       else ! Not allocated yet
          !write(*,*)'impavnc0 Before allocate; Size req. for abd:',lenabd
          !YuP[2019-04-24] allocate(abd(md1abd,inewjmax),STAT=istat) !Big size ~3*jx*iy^2
          ![2019-04-24] abd() is used in subr.dgbtrf and dgbtrs with dims abd(md1abd,inewjx)
          !So changing allocation to
          allocate(abd(md1abd,inewjx),STAT=istat) !YuP[2019-04-24] !Big size ~3*jx*iy^2
          ![2019-04-24] The problem is - inewjx may change from one flux surface
          ! to another, and that is what inewjmax was trying to account for
          !(it is the maximum of inew*jx = (iyh+itl-1)*jx over all flux surfaces).
          !
          !write(*,*)'impavnc0 abd: Allocate/istat=',istat
          !XXXcall bcast(abd,0,lenabd)   YuP: ok, below
          abd=zero !YuP[2019-04-24] was 0.
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

             !write(*,*)'impavnc0:  factor.eq."disabled", OK?'
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
                           .and. i.eq.itu) rhs(ieq)=zero !YuP[2019-04-24] was 0.

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
                if (istat.ne.0) WRITE(*,*)'impanvc0: ampfda prblm'
                allocate(ampfdd(0:iy,jx),STAT=istat)
                if (istat.ne.0) WRITE(*,*)'impanvc0: ampfdd prblm'
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
                      !YuP[2019-05-30] Note that fh() and fg() are allocated as 
                      ! fh(0:iy+1,0:jx+1,1,1:lrz), i.e. k=1 index only, for now.
                      
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
                      !write(*,'(a,3i4,4e12.3)')'l_,j,i,rhs=',
                      !l_,j,i,rhs_old,rhs(ieq),ampfda(i,j),ampfdd(i,j)

                   endif        ! On if j.eq.1
                enddo           ! On indexsww
             enddo              ! On j
          endif                 ! On ampfmod/kopt

       else  ! factor.ne."disabled", i.e. FACTORIZE

          !.......................................................................
          !
          !     reinitialize to zero matrix and rhs
          !
          !write(*,*)'impavnc0:size(rhs) ',size(rhs)
          if (soln_method.eq.'direct' .or. soln_method.eq.'itsol1') then
             !XXX call bcast(abd,0,lenabd) YuP: ok, below
             abd=zero !YuP[2019-04-24] was 0.
             do i=1,inewjx_(l_)
                rhs(i)=zero !YuP[2019-04-24] was 0.
                ipivot(i)=0
             enddo
          elseif (soln_method.eq.'itsol') then
             do i=1,inewjx_(l_)
                rhs(i)=zero !YuP[2019-04-24] was 0.
             enddo
          elseif (ifirst_lr.eq.1) then  ! for it3dv and it3drv
             do i=1,ieq_tot
                rhs(i)=zero !YuP[2019-04-24] was 0.
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
          !write(*,*)'impavnc0, before j and i loops: ieq,icoeff=',
          !ieq,icoeff

          !%OS
          !          zzzto=0.0
          !          zzz12=second()
          !%OS

          !..................................................................
          !     do loop 12 is the major loop over j (momentum).
          !..................................................................

          zrhsj1=zero !YuP[2019-04-24] was 0.
          !BH040719          call bcast(zavarj1,zero,mu+1)
          !BH040719          call ibcast(ijaj1,0,mu+1)
          call bcast(zavarj1,zero,iy+4)
          zmatcont=zero !YuP[2019-04-24] added
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
                !     for example xmm(i,j,k)*f(i-1,j-1,l_)+ x0m(i,j-1)*f(i,j-1,l_) +..
                !     + seven other terms on left side = z00(i,j,k)
                !     xmm, xpm, t0m, z00 etc are statement functions defined in
                !
                !     abd stores the transposed of the band matrix in the LAPACK sgbtrf format
                !.................................................................

                zmatcont(1)=xmm(i,j,k)
                zmatcont(2)=x0m(i,j,k)
                zmatcont(3)=xpm(i,j,k)
                zmatcont(4)=xm0(i,j,k)
                zmatcont(5)=x00(i,j,k)
                zmatcont(6)=xp0(i,j,k)
                zmatcont(7)=xmp(i,j,k)
                zmatcont(8)=x0p(i,j,k)
                zmatcont(9)=xpp(i,j,k)

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

                zmatcont(1)=xmm(i,j,k)
                zmatcont(2)=x0m(i,j,k)+xpm(i,j,k)
                zmatcont(3)=xm0(i,j,k)
                zmatcont(4)=x00(i,j,k)+xp0(i,j,k)
                zmatcont(5)=xmp(i,j,k)
                zmatcont(6)=x0p(i,j,k)+xpp(i,j,k)

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
                   zmatcont(1)=one !YuP[2019-04-24] was 1.
                   rhs(ieq)=f_(iyh,1,k,l_) ! f_ at previous t-step
                   xnorm=one !YuP[2019-04-24] was 1.
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
                   zmatcont(1)=one !YuP[2019-04-24] was 1.
                   xnorm=one !YuP[2019-04-24] was 1.
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
                   zmatcont(1)=x00(i,j,k)
                   zmatcont(2)=xp0(i,j,k)
                   zmatcont(3)=x0p(i,j,k)
                   zmatcont(4)=xpp(i,j,k)

                   call impnorm(xnorm,zmatcont,rhs(ieq),nvar)
                   !if (ieq.le.840)write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !write(*,*)'HERE0.1: i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                   zmatcont(1)=xm0(i,j,k)
                   zmatcont(2)=x00(i,j,k)
                   zmatcont(3)=xmp(i,j,k)
                   zmatcont(4)=x0p(i,j,k)

                   call impnorm(xnorm,zmatcont,rhs(ieq),nvar)
                   !if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !write(*,*)'HERE0.2:i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                   zmatcont(1)=tm0(j,k)
                   zmatcont(2)=tmp(j,k)
                   zmatcont(3)=t00(j,k)
                   zmatcont(4)=t0p(j,k)
                   zmatcont(5)=tp0(j,k)
                   zmatcont(6)=tpp(j,k)
                   zmatcont(7)=tu0(j,k)
                   zmatcont(8)=tup(j,k)

                   call impnorm(xnorm,zmatcont,rhs(ieq),nvar)
                   !if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !write(*,*)'HERE0.3: i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                   zmatcont(1)=xm0(i,j,k)
                   zmatcont(2)=x00(i,j,k)+xp0(i,j,k)
                   zmatcont(3)=xmp(i,j,k)
                   zmatcont(4)=x0p(i,j,k)+xpp(i,j,k)

                   call impnorm(xnorm,zmatcont,rhs(ieq),nvar)
                   !if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !write(*,*)'HERE0.4: i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                   zmatcont(1)=xm0(i,j,k)
                   zmatcont(2)=x00(i,j,k)
                   zmatcont(3)=xp0(i,j,k)
                   zmatcont(4)=xmp(i,j,k)
                   zmatcont(5)=x0p(i,j,k)
                   zmatcont(6)=xpp(i,j,k)
                   call impnorm(xnorm,zmatcont,rhs(ieq),nvar)
                   !if (ieq.le.840) write(*,*)'ieq,rhs(ieq)',ieq,rhs(ieq)
                   !write(*,*)'HERE0.5: i,j,icntsww,xnorm',i,j,icntsww,xnorm
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
                zmatcont(1)=x0m(i,j,k)
                zmatcont(2)=xpm(i,j,k)
                zmatcont(3)=x00(i,j,k)
                zmatcont(4)=xp0(i,j,k)
                zmatcont(5)=x0p(i,j,k)
                zmatcont(6)=xpp(i,j,k)

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
                zmatcont(1)=xmm(i,j,k)
                zmatcont(2)=x0m(i,j,k)
                zmatcont(3)=xm0(i,j,k)
                zmatcont(4)=x00(i,j,k)
                zmatcont(5)=xmp(i,j,k)
                zmatcont(6)=x0p(i,j,k)

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
                zmatcont(1)=tmm(j,k)
                zmatcont(2)=tm0(j,k)
                zmatcont(3)=tmp(j,k)
                zmatcont(4)=t0m(j,k)
                zmatcont(5)=t00(j,k)
                zmatcont(6)=t0p(j,k)
                zmatcont(7)=tpm(j,k)
                zmatcont(8)=tp0(j,k)
                zmatcont(9)=tpp(j,k)
                zmatcont(10)=tum(j,k)
                zmatcont(11)=tu0(j,k)
                zmatcont(12)=tup(j,k)

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
                zmatcont(1)=xmm(i,j,k)
                zmatcont(2)=xm0(i,j,k)
                zmatcont(3)=x0m(i,j,k)
                zmatcont(4)=x00(i,j,k)
                zmatcont(5)=xpm(i,j,k)
                zmatcont(6)=xp0(i,j,k)

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
                zmatcont(1)=xmm(i,j,k)
                zmatcont(2)=xm0(i,j,k)
                zmatcont(3)=x0m(i,j,k)+xpm(i,j,k)
                zmatcont(4)=x00(i,j,k)+xp0(i,j,k)

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
                zmatcont(1)=tmm(jx,k)
                zmatcont(2)=tm0(jx,k)
                zmatcont(3)=t0m(jx,k)
                zmatcont(4)=t00(jx,k)
                zmatcont(5)=tpm(jx,k)
                zmatcont(6)=tp0(jx,k)
                zmatcont(7)=tum(jx,k)
                zmatcont(8)=tu0(jx,k)

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
                zmatcont(1)=x00(i,j,k)
                zmatcont(2)=xm0(i,j,k)
                zmatcont(3)=x0m(i,j,k)
                zmatcont(4)=xmm(i,j,k)

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
                zmatcont(1)=x00(i,j,k)
                zmatcont(2)=xp0(i,j,k)
                zmatcont(3)=xpm(i,j,k)
                zmatcont(4)=x0m(i,j,k)

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
                   zmatcont(1)=one !YuP[2019-04-24] was 1.0
                   rhs(ieq)=zero !YuP[2019-04-24] was 0.

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
                   zmatcont(1)=tp0(j,k)
                   zmatcont(2)=tpp(j,k)
                   zmatcont(3)=t00(j,k)
                   zmatcont(4)=t0p(j,k)
                   zmatcont(5)=tu0(j,k)
                   zmatcont(6)=tup(j,k)
                   zmatcont(7)=tm0(j,k)
                   zmatcont(8)=tmp(j,k)
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
                   zmatcont(1)=tpm(jx,k)
                   zmatcont(2)=tp0(jx,k)
                   zmatcont(3)=t0m(jx,k)
                   zmatcont(4)=t00(jx,k)
                   zmatcont(5)=tum(jx,k)
                   zmatcont(6)=tu0(jx,k)
                   zmatcont(7)=tmm(jx,k)
                   zmatcont(8)=tm0(jx,k)

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
                   zmatcont(1)=tpm(j,k)
                   zmatcont(2)=tp0(j,k)
                   zmatcont(3)=tpp(j,k)
                   zmatcont(4)=t0m(j,k)
                   zmatcont(5)=t00(j,k)
                   zmatcont(6)=t0p(j,k)
                   zmatcont(7)=tum(j,k)
                   zmatcont(8)=tu0(j,k)
                   zmatcont(9)=tup(j,k)
                   zmatcont(10)=tmm(j,k)
                   zmatcont(11)=tm0(j,k)
                   zmatcont(12)=tmp(j,k)

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
                            zmatcont(1)=half !YuP[2019-04-24] was 0.5
                            zmatcont(2)=-half !YuP[2019-04-24] was -0.5
                            janelt(1)=ieq
                            janelt(2)=inew
                            rhs(ieq)=zero !YuP[2019-04-24] was 0.
                            xnorm=2.d0 !YuP[2019-04-24] was 2.0
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
                               abd(ibandpieq-ijaj1(jcont),ijaj1(jcont))= zavarj1(jcont)
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
                            zmatcont(1)=half !YuP[2019-04-24] was 0.5
                            zmatcont(2)=-half !YuP[2019-04-24] was -0.5
                            janelt(1)=ieq
                            janelt(2)=inew
                            rhs(ieq)=zero !YuP[2019-04-24] was 0.
                            xnorm=2.d0 !YuP[2019-04-24] was 2.0
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
                            zmatcont(1)=half !YuP[2019-04-24] was 0.5
                            zmatcont(2)=-half !YuP[2019-04-24] was -0.5
                            janelt(1)=ieq
                            janelt(2)=inew+(ieq_(l_)-1)
                            rhs(ieq)=zero !YuP[2019-04-24] was 0.
                            xnorm=2.d0 !YuP[2019-04-24] was 2.0
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
                scal(kku,l_)=one/xnorm !!YuP[2019-04-24] was 1./xnorm
                !if(kku.le.0 .or. kku.gt.iyjx*ngen) then
                !write(*,*)'impavnc0: k,l_,ieq,ieq_(l_),kku,ieq_tot=',
                !k,l_,ieq,ieq_(l_),kku,ieq_tot
                !                pause
                !endif
                
                !..................................................................
                !     End of major j and I loops..
                !..................................................................
                !write(*,*)'ij-loop-end: ieq,rhs(ieq) ',ieq,rhs(ieq)
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
                   !write(*,*)'impanvc0: Number of CSR coeffs: icoeff=',icoeff
                   !         icoeff_est=(iy+2)*3+(jx-1)*9*((iyh-itl)+2*itl)  Calc'd above
                   !write(*,*)'impavnc0: Estimate of of icoeff=', icoeff_est
                   !         Check coeff storage:
                   if (icoeff.gt.size(a_csr)) then
!MPIINSERT_IF_RANK_EQ_0
                      WRITE(*,*)'impavnc0:icoeff.gt.size(a_csr)'
!MPIINSERT_ENDIF_RANK
                      STOP
                   endif
                endif  ! on soln_method.ne.'direct'

                !     Check number of equations:
                !        if (soln_method.eq.'itsol' .or. soln_method.eq.'itsol1') then
                !write(*,*)'impavnc0: ieq should equal inewjx: ',ieq,inewjx
                !        elseif(soln_method.eq.'it3dv' .or. soln_method.eq.'it3drv') then
                !write(*,*)'impavnc0: ieq should equal ieq_tot: ',ieq,ieq_tot
                !           pause
                !        endif


                if (soln_method.eq.'itsol') then
                   iunit=39
                   i1=1
                   i2=icoeff
                endif  ! on solution_meth=itsol

                !..................................................................
                !     Next endif is the end for "if(factor.eq."disabled") then".
                !     Thus matrix and rhs are defined
                !..................................................................

             endif     ! on factor

             !write(*,*)'impavnc0:  factor =',factor



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
                ! nzmax = integer. The  length of the arrays c and jc.
                !         aplb will stop if the result matrix C  has a number
                !         of elements that exceeds nzmax. See ierr.
                ! ierr  = integer. serving as error message.
                !         ierr = 0 means normal return,
                !         ierr .gt. 0 means that aplb stopped while computing the
                !         i-th row  of C with i=ierr, because the number
                !         of elements in C exceeds nzmax.
                ! work arrays:
                ! iw    = integer work array of length equal to the number of
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
                ! XXX , changed  a(1),ja(1),ia(1) ~~> a,ja,ia  ! YuP:ok
                call bndcsr(n_rows_A,abd_lapack,nabd,lowd,ml_a,mu_a, &
                     a_csr,ja_csr,ia_csr,icsrij,ierr_csr)
                if (ierr_csr.ne.0) then
!MPIINSERT_IF_RANK_EQ_0
                   WRITE(*,*)'impavnc0/bndcsr: STOP ierr_csr=',ierr_csr
!MPIINSERT_ENDIF_RANK
                   stop
                endif

             endif  !  on soln_method.eq.'itsol1'


             if (soln_method.eq.'itsol' .or. soln_method.eq.'itsol1') then
110             format(10(1pe10.2))
111             format(10i10)
             endif  ! on soln_method itsol or itsol1


             if ( (soln_method.eq.'it3dv')  .and.  &
                  ilast_lr.eq.1 .and. n.eq.1  ) then
                !write out coeffs in a format which can be compared directly 
                !with soln_method=itsol output.
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
                   !XXX again, now passing by array, was passing elem.  YuP:ok
                   call ilut (n_rows_A,a_csr,ja_csr,ia_csr,lfil0, &
                        droptol,alu,jlu,ju,iwk_ilu,w_ilu, &
                        jw_ilu,ierr)
                   !Setup up above, since vv size dependency
                   do i=1,n_rows_A
                      rhs0(i)=rhs(i)    !Copy rhs, for input to pgmres since
                      !pgmres modifies this input.
                   enddo
                elseif(soln_method.eq.'it3dv') then
                   ! perform soln only at last flux surface (l_=lrz)
                   if ( ilast_lr.eq.1 ) then
                      ! XXX passing array.  YuP:ok
                      call ilut (n_rows_A,a_csr,ja_csr,ia_csr,lfil0, &
                           droptol,alu,jlu,ju,iwk_ilu,w_ilu, &
                           jw_ilu,ierr)
                      !              call ilutp(n_rows_A,a_csr(1),ja_csr(1),ia_csr(1),
                      !     +             lfil0,droptol,permtol,mbloc,
                      !     +             alu(1),jlu(1),ju(1),iwk_ilu,w_ilu(1),jw_ilu(1),
                      !     +             iperm_ilu(1),ierr)
!MPIINSERT_IF_RANK_EQ_0
                      !WRITE(*,*)'impavnc0 aft.ilut:  l_,ierr',l_,ierr,soln_method
!MPIINSERT_ENDIF_RANK
                      do i=1,n_rows_A
                         rhs0(i)=rhs(i) !Copy rhs, for input to pgmres since
                         !this input is modified.
                      enddo
                   endif
                elseif(soln_method.eq.'it3drv') then
                   ! perform soln only at last flux surface (l_=lrz)
                   if ( ilast_lr.eq.1 ) then
                      call ilut (n_rows_A,ac_csr,jac_csr,iac_csr,lfil0, &
                           droptol,alu,jlu,ju,iwk_ilu,w_ilu, &
                           jw_ilu,ierr)
                      !              call ilutp(n_rows_A,ac_csr(1),jac_csr(1),iac_csr(1),
                      !     +             lfil0,droptol,permtol,mbloc,
                      !     +             alu(1),jlu(1),ju(1),iwk_ilu,w_ilu(1),jw_ilu(1),
                      !     +             iperm_ilu(1),ierr)
!MPIINSERT_IF_RANK_EQ_0
                      !WRITE(*,*)'impavnc0 aft.ilut:  l_,ierr',l_,ierr,soln_method
!MPIINSERT_ENDIF_RANK
                      do i=1,n_rows_A
                         rhs0(i)=rhs(i) !Copy rhs, for input to pgmres since
                         !this input is modified.
                      enddo
                   endif
                endif

                call cpu_time(tm(1))
                tm1= tm(1) - tm1
                !write(*,*)'impavnc0: icount_imp, time for ilut tm1=',
                !icount_imp,tm1
                if (ierr.ne.0) then
!MPIINSERT_IF_RANK_EQ_0
                   WRITE(*,*)'impavnc0 after ilut: ierr=',ierr
!MPIINSERT_ENDIF_RANK
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
                   !XXXX  YuP: ok
                   call pgmres(n_rows_A,krylov1,rhs0,sol,vv,epsilon, &
                        maxits,iout,a_csr,ja_csr,ia_csr,alu, &
                        jlu,ju,ierr)
                   rhs(1:n_rows_A)=sol(1:n_rows_A)
                elseif(soln_method.eq.'it3dv') then
                   ! perform soln
                   ! only at last flux surface (l_=lrz).
                   if ( ilast_lr.eq.1 ) then
!MPIINSERT_IF_RANK_EQ_0
                      !WRITE(*,*)'impavnc0 before pgmres:  ierr',ierr,soln_method
!MPIINSERT_ENDIF_RANK
                      !XXX  YuP: ok
                      call pgmres(n_rows_A,krylov1,rhs0,sol,vv,epsilon, &
                           maxits,iout,a_csr,ja_csr,ia_csr,alu, &
                           jlu,ju,ierr)
!MPIINSERT_IF_RANK_EQ_0
                      !WRITE(*,*)'impavnc0  after pgmres:  ierr',ierr
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
                      !WRITE(*,*)'impavnc0 before pgmres:  ierr',ierr,soln_method
!MPIINSERT_ENDIF_RANK
                      !XXX YuP:ok
                      call pgmres(n_rows_A,krylov1,rhs0,sol,vv,epsilon, &
                           maxits,iout,ac_csr,jac_csr,iac_csr,alu, &
                           jlu,ju,ierr)
!MPIINSERT_IF_RANK_EQ_0
                      !WRITE(*,*)'impavnc0  after pgmres:  ierr',ierr
!MPIINSERT_ENDIF_RANK

                      do inr = 1,n_rows_A
                         rhs(inr)=sol(inr)
                         !-> rhs(1:n_rows_A)=sol(1:n_rows_A) results in stack overflow
                      enddo
                   endif
                endif

                call cpu_time(tm(1))
                tm1= tm(1) - tm1
                !write(*,*)'impavnc0: time for pgmres tm1=',tm1
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

!MPIINSERT_IF_RANK_EQ_0      
      WRITE(*,'(a,i4,2e19.11)')  &
       'impavnc rhs before dgbtrf l_,MIN(rhs),MAX(rhs)', &
                l_,MINVAL(rhs),MAXVAL(rhs)
      WRITE(*,'(a,i4,2e19.11)')  &
       'impavnc ABD before dgbtrf l_,MIN(abd),SUM(abd)', &
                l_,MINVAL(abd),SUM(abd)
!MPIINSERT_ENDIF_RANK
                
                !     factorize matrix
                !_cray  64-bit-compiler uses sgbtrf (from lapack library).
                !_pc    32-bit-compiler uses dgbtrf (from lapack library).
                if (factor .ne. "disabled") then
                   !          call sgbtrf(inewjx,inewjx,ml,mu,abd,md1abd,ipivot,info)
                   !write(*,*)'impavnc0 before dgbtrf, l_,inewjx,ml=',l_,inewjx,ml
                   ! XXX  YuP:ok
                   call dgbtrf(inewjx,inewjx,ml,mu,abd,md1abd,ipivot,info)
                   if (info .ne. 0) then
                      print *,' warning after sgbtrf in impavnc0: info = ',info
                      stop 'impavnc0 1'
                   endif
                endif
                !       tm1 = etime(tm) - tm1
                !write(*,*)'impavnc0: time for dgbtrf tm1=',tm1

                !     solve system
                !_cray  64-bit-compiler uses sgbtrs (from lapack library).
                !_pc    32-bit-compiler uses dgbtrs (from lapack library).

                !BH080303        tm1 = etime(tm)
                call cpu_time(tm1)

                inbrhs = 1
                transpose = 'n'
                !        call sgbtrs(transpose,inewjx,ml,mu,inbrhs,abd,md1abd,ipivot
                ! XXXX  YuP: ok
                call dgbtrs(transpose,inewjx,ml,mu,inbrhs,abd,md1abd, &
                     ipivot,rhs,inewjx,info)

!MPIINSERT_IF_RANK_EQ_0      
      WRITE(*,'(a,i4,2e19.11)')  &
       'impavnc rhs after  dgbtrs l_,MIN(rhs),MAX(rhs)', &
                l_,MINVAL(rhs),MAXVAL(rhs)
      WRITE(*,'(a,i4,2e19.11)')  &
       'impavnc ABD after  dgbtrs l_,MIN(abd),SUM(abd)', &
                l_,MINVAL(abd),SUM(abd)
!MPIINSERT_ENDIF_RANK
                
                if (info .ne. 0) then
                   print *,' warning after sgbtrs in impavnc0: info = ',info
                   stop 'impavnc0 2'
                endif
                !       tm1 = etime(tm) - tm1
                !write(*,*)'impavnc0: time for dgbtrs tm1=',tm1

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
                !write(*,*)'impavnc0: n,ilast_lr or rhs=>f',n,ilast_lr
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

!MPIINSERT_IF_RANK_EQ_0
     !WRITE(*,'(a,2i6)') 'impavnc[ZOW]: sol.found, f is updated. k,l_=',k,l_
      WRITE(*,'(a,2i4,3e15.7)')  &
       'impavnc[ZOW] f fnd,updatd n,l_,MIN(f),MAX(f),SUM(f)=', &
                n,l_,MINVAL(f),MAXVAL(f),SUM(f)
!MPIINSERT_ENDIF_RANK
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
                   call dcopy(iyjx2*ngen,f(0:iy+1,0:jx+1,1:ngen,l_),1, &
                                        fh(0:iy+1,0:jx+1,1:ngen,l_),1)
                   !YuP[2019-05-30] Note that in case of cqlpmod.eq."enabled", 
                   ! fh() is allocated as fh(0:iy+1,0:jx+1,1:ngen,0:ls+1) in wpalloc
                   factor="disabled"
                   icount_ese=2
                   go to 10

                elseif (icount_ese.eq.2) then

                   !             g function is in f.
                   !             Calculate fluxes, restore f with f_, and return.

                   flux1(l_)=fluxpar( &
                        1,x,coss(1:iy,l_),cynt2(1:iy,l_),cint2,temp1,iy,jx)
                   flux2(l_)=fluxpar( &
                        1,x,coss(1:iy,l_),cynt2(1:iy,l_),cint2,temp2,iy,jx)
                   call dcopy(iyjx2*ngen,f_(0:iy+1,0:jx+1,1:ngen,l_),1, &
                                          f(0:iy+1,0:jx+1,1:ngen,l_),1)

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
                   call dcopy(iyjx2,f(0:iy+1,0:jx+1,kelec,l_),1, &
                                   fh(0:iy+1,0:jx+1,1,l_),1)
                   !YuP[05-2017] corrected dcopy in the above line:
                   !was iyjx2*ngen, but we only copy one species (kelec)
                   factor="disabled"
                   icount_ampf=2 !next, re-define rhs() for the 'g' function equations
                   go to 10

                elseif (icount_ampf.eq.2) then

                   !             g function is in f, store into fg
                   call dcopy(iyjx2,f(0:iy+1,0:jx+1,kelec,l_),1, &
                                   fg(0:iy+1,0:jx+1,1,l_),1)
                   !YuP[05-2017] corrected dcopy in the above line:
                   !was iyjx2*ngen, but we only copy one species (kelec)
                   !BH170312:  Checking effect of zeroing correction part of distn:
                   !BH170312:              call bcast(fg(0,0,1,l_),zero,iyjx2*ngen)
                   !             Restore f with f_, and return.
                   call dcopy(iyjx2,f_(0:iy+1,0:jx+1,kelec,l_),1, &
                                     f(0:iy+1,0:jx+1,kelec,l_),1)
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
  real(c_double) function z00(i,j,k)
    use param_mod
    use cqlcomm_mod
    use advnce_mod !here: in z00(). To get xmm(), etc. (many functions)
    implicit none
    integer :: i
    integer :: j
    integer :: k !XXX example why k is terrible global var name
    ! k is species index, everywhere in the code, except subr. frsubs, ilut, r8subs
    integer :: istat0
    integer :: iboot
    real(c_double) :: z00f
    real(c_double) :: z00t
    real(c_double) :: z00itl
    real(c_double) :: z00itl1
    real(c_double) :: z00itu1

    !cc         save  ! YuP: Not really needed


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
    !YuP[2019-04-24] Changed statement functions into scalars below.
    real(c_double) :: t0ml_
    real(c_double) :: t00l_
    real(c_double) :: t0pl_
    real(c_double) :: t0mu_
    real(c_double) :: t00u_
    real(c_double) :: t0pu_


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
          z00itl1=z00f -xpm(i,j,k)*bsl(j-1,k,l_)-xp0(i,j,k)*bsl(j,k,l_) &
               -xpp(i,j,k)*bsl(j+1,k,l_)
          z00t=z00itl1
       else
          !              itu+1 case:
          if (i.eq.(itu+1))then
             z00itu1=z00f -xmm(i,j,k)*bsu(j-1,k,l_)-xm0(i,j,k)*bsu(j,k,l_) &
                  -xmp(i,j,k)*bsu(j+1,k,l_)
             z00t=z00itu1
             !              itl (or itu)case:
          else
          
         t0ml_=qz(j)*( &
         cl(itl-1,j-1,l_)*dj(itl,j-1,k,l_)*eym5(itl,l_)) &
         +r2y(j,l_)*(-de(itl-1,j)*(1.-di(itl-1,j-1,k,l_))) &
         /(2.*dx(j))

         t00l_= &
         +qz(j)*( &
         -cl(itl-1,j,l_)*dj(itl,j,k,l_)*eym5(itl,l_) &
         +cl(itl-1,j-1,l_)*(1.-dj(itl,j-1,k,l_))*eym5(itl,l_)) &
         +r2y(j,l_)*(dd(itl-1,j)*(1.-di(itl-1,j,k,l_)) &
         +df(itl-1,j)*eym5(itl,l_))

         t0pl_=qz(j)*( &
         -cl(itl-1,j,l_)*eym5(itl,l_)*(1.-dj(itl,j,k,l_))) &
         +r2y(j,l_)*de(itl-1,j)/(2.*dx(j))*(1.-di(itl-1,j+1,k,l_))

         t0mu_=qz(j)*( &
         -cl(itu+1,j-1,l_)*dj(itu,j-1,k,l_)*eyp5(itu,l_)) &
         +r2y(j,l_)*( &
         +de(itu,j)*di(itu,j-1,k,l_))/(2.*dx(j))

         t00u_= &
         +qz(j)*( &
         +cl(itu+1,j,l_)*dj(itu,j,k,l_)*eyp5(itu,l_) &
         -cl(itu+1,j-1,l_)*(1.-dj(itu,j-1,k,l_))*eyp5(itu,l_)) &
         +r2y(j,l_)*( &
         -dd(itu,j) &
         *di(itu,j,k,l_) &
         +df(itu,j)*eyp5(itu,l_))

         t0pu_=qz(j)*( &
         +cl(itu+1,j,l_)*(1.-dj(itu,j,k,l_))*eyp5(itu,l_)) &
         +r2y(j,l_)*( &
         -de(itu,j)*di(itu,j+1,k,l_)/ &
         (2.*dx(j)))
          
             z00itl= z00f -(  &
                   t0ml_*bsl(j-1,k,l_)+t00l_*bsl(j,k,l_)   &
                  +t0pl_*bsl(j+1,k,l_)+t0mu_*bsu(j-1,k,l_) &
                  +t00u_*bsu(j,k,l_)  +t0pu_*bsu(j+1,k,l_)      )
             z00t=z00itl
          endif
       endif
       z00=z00t

    else ! iboot=0
       z00=z00f
    endif

  end function z00

  !      end of file impavnc0.f

  subroutine it3dalloc
      use param_mod
      use cqlcomm_mod
      implicit none
      integer :: i
      integer :: i2
      integer :: istat(18)
      integer :: istat0
!MPIINSERT_INCLUDE

!dir$ nobounds

!.................................................................
!     Setup CSR storage for iterative, sparse matrix solve,
!     coefficient matrix and rhs/solution matrix.
!     Storage is allocated depending on the settings of soln_method.
!     Storage is provided for coefficient matrix for up to the complete
!     set of equations:
!     soln_method='it3dv'  ==> number of eqns = sum l_=1:lrz {inew*jx},
!         where inew=iyh_(l_) + itl_(l_) - 1
!         The number of columns = sum l_=1:lrz {inew*jx*9}
!         (Usually there are 9 coeffs per eqn; the few 12 coeff cases
!          are more than offset the few 2 and 3 coeff cases.)
!     soln_method='it3drv' ==>
!
!
!.................................................................



!     For soln_method='itsol' or 'itsol1', allocate sufficient
!     storage for flux-surface soln with largest number or eqns (l_=1).
!     Specific shape of abd_lapack will be set in impanvc0 for
!     soln_method='itsol1' case.

!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)'it3dalloc: entering...  soln_method=', &
                 soln_method
!MPIINSERT_ENDIF_RANK



      if (   soln_method.eq.'itsol' .or. soln_method.eq.'it3dv' &
        .or. soln_method.eq.'it3drv') then

!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'icsrij,icsrip,icsrikry,lfil,iwk_ilu = ', &
                    icsrij,icsrip,icsrikry,lfil,iwk_ilu
!MPIINSERT_ENDIF_RANK

         istat(1:18)=0
         allocate( a_csr(icsrij),stat=istat(1))
         allocate(ja_csr(icsrij),stat=istat(2))
         allocate(ia_csr(icsrip),stat=istat(3))
         allocate(alu(iwk_ilu),stat=istat(4))
         allocate(jlu(iwk_ilu),stat=istat(5))
         allocate(ju(icsrip),stat=istat(6))
         allocate(jw_ilu(icsri2),stat=istat(7))
         allocate(w_ilu(icsrip),stat=istat(8))
         allocate(rhs0(icsrip),stat=istat(9))
         allocate( sol(icsrip),stat=istat(10))
         allocate(vv(icsrikry),stat=istat(11))
         i2=11
         if (soln_method .eq. 'it3drv') then
            write(*,*)'icsrijr,icsrijc = ',icsrijr,icsrijc
            allocate(ipofi(iy,lrz),stat=istat(12))
            allocate( ar_csr(icsrijr),stat=istat(13))
            allocate(jar_csr(icsrijr),stat=istat(14))
            allocate(iar_csr(icsrip),stat=istat(15))
            allocate(iac_csr(icsrip),stat=istat(18))
            allocate( ac_csr(icsrijc),stat=istat(16))
            allocate(jac_csr(icsrijc),stat=istat(17))
            i2=18
         endif
         istat0=0
         do i=1,i2 !!12
            if (istat(i).ne.0) then
               write(*,*)'it3dalloc:  i,istat(i)=',i,istat(i)
               istat0=max(istat0,1)
            endif
         enddo
         if (istat0.ne.0) stop 'Allocation problem in it3dalloc'
!$$$

      elseif ( soln_method.eq.'itsol1' ) then

!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'lapacki,lapackj,icsrij,icsrip,iwk_ilu = ', &
                    lapacki,lapackj,icsrij,icsrip,iwk_ilu
!MPIINSERT_ENDIF_RANK

         allocate(abd_lapack(lapacki,lapackj),stat=istat(1))
         allocate(a_csr(icsrij),stat=istat(1))
         allocate(ja_csr(icsrij),stat=istat(2))
         allocate(ia_csr(icsrip),stat=istat(3))
         allocate(alu(iwk_ilu),stat=istat(4))
         allocate(jlu(iwk_ilu),stat=istat(5))
         allocate(ju(icsrip),stat=istat(6))
         allocate(jw_ilu(icsri2),stat=istat(7))
         allocate(w_ilu(icsrip),stat=istat(8))
         allocate(rhs0(icsrip),stat=istat(9))
         allocate(sol(icsrip),stat=istat(10))
         allocate(vv(icsrikry),stat=istat(11))
         istat0=0
         do i=1,11
            if (istat(i).ne.0) then
               write(*,*)'it3dalloc:  i,istat(i)=',i,istat(i)
               istat0=max(istat0,1)
            endif
         enddo
         if (istat0.ne.0) stop 'Allocation problem in it3dalloc'

      endif

!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)'it3dalloc: exiting...'
!MPIINSERT_ENDIF_RANK

      return
    end subroutine it3dalloc

!=======================================================================

   subroutine it3ddalloc
      use param_mod
      use cqlcomm_mod
      implicit none
      integer :: istat
!MPIINSERT_INCLUDE

!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)'it3ddalloc: deallocating...'
!MPIINSERT_ENDIF_RANK

      if ( soln_method.eq.'itsol' .or. soln_method.eq.'it3dv' &
           .or. soln_method.eq.'it3drv') then
         deallocate (a_csr,ja_csr,ia_csr,alu,jlu,ju, &
                    jw_ilu,w_ilu,rhs0,sol,vv)
         if (soln_method.eq.'it3drv') then
            deallocate (ar_csr,jar_csr,iar_csr,ac_csr,jac_csr,iac_csr)
         endif

      elseif ( soln_method.eq.'itsol1' ) then
         deallocate (abd_lapack,a_csr,ja_csr,ia_csr,alu,jlu,ju, &
              jw_ilu,w_ilu,rhs0,sol,vv)

      elseif ( soln_method.eq.'it3drv' ) then

      endif

      return
    end subroutine it3ddalloc


!=======================================================================

   subroutine de_alloc ! YuP[11-2017] more deallocation
      use param_mod
      use cqlcomm_mod
      implicit none
      integer :: istat
!MPIINSERT_INCLUDE

!  The purpose of this subroutine is to ensure deallocation of variables
!  at the end of a run.  If running cql3d as a stand-alone code,
!  deallocation would occur automatically.
!  With present TRANSP coupling, invoking cql3d through
!  a subroutine call, allocated memory would build up with
!  each call to cql3d.

!  Generally, all pointers that are defined in comm.h
!  can be deallocated here.
!  Not all of them are allocated during run.
!  It depends on particular setup/cqlinput.
!  So, always check - add if(ASSOCIATED(array_name))
!  in front of deallocate(array_name) statement.
!  Other arrays like abd(),...ampfda() [few lines below]
!  are NOT in comm.h, but rather in local common blocks.
!  So, these common blocks must be added here, too.


!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)'de_alloc: START deallocating...'
!MPIINSERT_ENDIF_RANK

      if(ASSOCIATED(rhs0)) then
        ! rhs0 and other arrays are allocated => deallocate
        deallocate(rhs0)
      endif
      if(ASSOCIATED(abd)) then
        deallocate(abd)
      endif
      if(ASSOCIATED(ipivot)) then
        deallocate(ipivot)
      endif
      if(ASSOCIATED(ampfda)) then
        deallocate(ampfda,ampfdd,ampfaa)
      endif
      if(ASSOCIATED(fh)) then
        deallocate(fh,fg)
      endif
      if(ASSOCIATED(urfb)) then
        deallocate(ilim1,ilim2,ifct1,ifct2,jminray,jmaxray,lloc,llray)
        deallocate(urfb,urfc,cosmz,urftmp,urfpwr,urfpwrc,urfpwrl)
        deallocate(g_,alfag,argmnt,ilim1d,ilim2d,ilim1dd,ilim2dd)
        deallocate(psiloc,scalurf,cwexde,cweyde,cwezde,delpwr,fluxn)
        deallocate(seikon,spsi,sdpwr,sbtot,sene,salphac,salphal)
        deallocate(ws,wr,wz,wnpar,wdnpar,wnper,wphi)
        deallocate(ilowp,iupp,ifct1_,ifct2_)
      endif

      if(ASSOCIATED(f)) then
        deallocate(f)
      endif
      if(ASSOCIATED(favg)) then
        deallocate(favg)
      endif
      if(ASSOCIATED(fxsp)) then
        deallocate(fxsp)
      endif
      if(ASSOCIATED(f_)) then
        deallocate(f_)
      endif
      if(ASSOCIATED(spasou)) then
        deallocate(spasou)
      endif
      if(ASSOCIATED(velsou)) then
        deallocate(velsou)
      endif
      if(ASSOCIATED(velsou2)) then
        deallocate(velsou2)
      endif
      if(ASSOCIATED(source)) then
        deallocate(source)
      endif
      if(ASSOCIATED(gone)) then
        deallocate(gone)
      endif
      if(ASSOCIATED(egylosa)) then
        deallocate(egylosa)
      endif
      if(ASSOCIATED(i0tran)) then
        deallocate(i0tran)
      endif
      if(ASSOCIATED(cal)) then
        deallocate(cal,cbl,ccl,cdl,cel,cfl,eal,ebl,scal,cet,cex)
        deallocate(synca,syncd,taulos,psi0bar)
        deallocate(delecfld0,elecfldn,delecfld0n,elecn,di,dj,ss,dcofleg)
        deallocate(dpcosz,ssy,ssyy,ssyi,ssyyy,pcurr,pcurrm,pdens,pdenm)
        deallocate(pengy,pengym)
        deallocate(wtfl0,wtflm,jflbin)
        deallocate(currv,pwrrf,pwrrfs,wflux,feta,fetb)
        deallocate(sgaint,entr)
      endif
      if(ASSOCIATED(f_aveth)) then
        deallocate(f_aveth)
      endif
      if(ASSOCIATED(densz)) then
        deallocate(densz,waa,wbb,cosz,dtau,sinz,tanz,yz,tot)
        deallocate(vflux,sincosba)
      endif
      if(ASSOCIATED(sovt)) then
        deallocate(sovt)
      endif
      if(ASSOCIATED(sigsxr)) then
        deallocate(sigsxr)
      endif

      if(ASSOCIATED(pentr)) then
        deallocate(pentr)
      endif
      if(ASSOCIATED(sounor)) then
        deallocate(sounor)
      endif
      if(ASSOCIATED(tamt1)) then
        deallocate(tamt1,tamt2)
      endif
      if(ASSOCIATED(cqlb)) then
        deallocate(cqlb,cqlc,cqle,cqlf)
      endif
      if(ASSOCIATED(frn_2)) then
        deallocate(frn_1,frn_2)
      endif
      if(ASSOCIATED(frn)) then
        deallocate(frn)
      endif
      if(ASSOCIATED(fvn)) then
        deallocate(fvn)
      endif
      if(ASSOCIATED(fvn_1)) then
        deallocate(fvn_1)
      endif
      if(ASSOCIATED(dl)) then
        deallocate(dl)
      endif
      if(ASSOCIATED(d_rr)) then
        deallocate(d_rr)
      endif
      if(ASSOCIATED(d_r)) then
        deallocate(d_r)
      endif
      if(ASSOCIATED(f_vtor)) then
        deallocate(f_vtor)
      endif
      if(ASSOCIATED(fnhalf)) then
        deallocate(fnhalf)
      endif
      if(ASSOCIATED(fnp0)) then
        deallocate(fnp0)
      endif
      if(ASSOCIATED(fnp1)) then
        deallocate(fnp1)
      endif
      if(ASSOCIATED(dls)) then
        deallocate(dls)
      endif
      if(ASSOCIATED(fedge)) then
        deallocate(fedge)
      endif
      if(ASSOCIATED(bndmats)) then
        deallocate(bndmats)
      endif
      if(ASSOCIATED(wcqlb)) then
        deallocate(wcqlb,wcqlc,wcqle,wcqlf)
      endif
      if(ASSOCIATED(rdcb)) then
        deallocate(rdcb,rdcc,rdce,rdcf)
      endif
      if(ASSOCIATED(ilpm1ef)) then
        deallocate(ilpm1ef)
      endif
      if(ASSOCIATED(rhspar)) then
        deallocate(rhspar)
      endif
      if(ASSOCIATED(fg_)) then
        deallocate(f_lm,f_lp,f_up,eg_,fg_)
      endif
      if(ASSOCIATED(csv)) then
        deallocate(csv)
      endif
      if(ASSOCIATED(deltarz)) then
        deallocate(deltarho,deltarhop,deltarz)
      endif
      if(ASSOCIATED(fgg)) then
        deallocate(fgg,egg,temp1,temp2,temp3,temp4,temp5,temp6)
        deallocate(xllji,xppji)
      endif
      if(ASSOCIATED(jbm1)) then
        deallocate(jbm1,jb0,jbp1)
      endif
      if(ASSOCIATED(nrayelt)) then
        deallocate(nrayelt,jslofas,nurefls,keiks,jpes,jpis,istarts)
        deallocate(iprmt5,jhlfs,sxxrt,skpsi,skth,skphi)
        deallocate(lrayelt,delpwr0,nrayelt0)
      endif
      if(ASSOCIATED(eqdell)) then
        deallocate(eqdell,eqbpol,solr,solz,drpmconz)
      endif
      if(ASSOCIATED(cynt2_)) then
        deallocate(cynt2_,vpint_,vptb_,cosovb,bovcos,adv)
      endif

!MG added 11/13/2017
      if(ASSOCIATED(dentarget)) then
         deallocate(dentarget)
      endif
      if(ASSOCIATED(sx)) then
         deallocate(sx)
      endif
      if(ASSOCIATED(xmdx)) then
         deallocate(xmdx)
      endif
      if(ASSOCIATED(cosz1)) then
         deallocate(cosz1)
      endif
      if(ASSOCIATED(sinz1)) then
         deallocate(sinz1)
      endif
      if(ASSOCIATED(sinz2)) then
         deallocate(sinz2)
      endif
      if(ASSOCIATED(thtf1)) then
         deallocate(thtf1)
      endif
      if(ASSOCIATED(thtf2)) then
         deallocate(thtf2)
      endif
      if(ASSOCIATED(alfi)) then
         deallocate(alfi)
      endif
      if(ASSOCIATED(alfa)) then
         deallocate(alfa)
      endif
      if(ASSOCIATED(truncd)) then
         deallocate(truncd)
      endif
      if(ASSOCIATED(w_ilu)) then
         deallocate(w_ilu)
      endif
      if(ASSOCIATED(sol)) then
         deallocate(sol)
      endif
      if(ASSOCIATED(vv)) then
         deallocate(vv)
      endif
      if(ASSOCIATED(jw_ilu)) then
         deallocate(jw_ilu)
      endif
      if(ASSOCIATED(ampfln)) then
         deallocate(ampfln)
      endif
      if(ASSOCIATED(ampflg)) then
         deallocate(ampflg)
      endif
      if(ASSOCIATED(ampfa)) then
         deallocate(ampfa)
      endif
      if(ASSOCIATED(ampfb)) then
         deallocate(ampfb)
      endif
      if(ASSOCIATED(ampfc)) then
         deallocate(ampfc)
      endif
      if(ASSOCIATED(ampf2ebar)) then
         deallocate(ampf2ebar)
      endif
              write(*,*)'it3dalloc-1.0'

      if(ASSOCIATED(dym5)) deallocate(dym5)
      if(ASSOCIATED(dyp5)) deallocate(dyp5)
      if(ASSOCIATED(eym5)) deallocate(eym5)
      if(ASSOCIATED(eyp5)) deallocate(eyp5)
      if(ASSOCIATED(y)) deallocate(y)
      if(ASSOCIATED(dy)) deallocate(dy)
      if(ASSOCIATED(yptb)) deallocate(yptb)
      if(ASSOCIATED(coss)) deallocate(coss)
      if(ASSOCIATED(cynt2)) deallocate(cynt2)
      if(ASSOCIATED(batot)) deallocate(batot)
      if(ASSOCIATED(lmax)) deallocate(lmax)
      if(ASSOCIATED(vpint)) deallocate(vpint)
      if(ASSOCIATED(psiiv)) deallocate(psiiv)
      if(ASSOCIATED(psiba)) deallocate(psiba)
      if(ASSOCIATED(psisq)) deallocate(psisq)
      if(ASSOCIATED(psicu)) deallocate(psicu)
      if(ASSOCIATED(psiqu)) deallocate(psiqu)
      if(ASSOCIATED(bavpd)) deallocate(bavpd)
      if(ASSOCIATED(bavdn)) deallocate(bavdn)
      if(ASSOCIATED(psiir)) deallocate(psiir)
      if(ASSOCIATED(vderb)) deallocate(vderb)
      if(ASSOCIATED(sinn)) deallocate(sinn)
      if(ASSOCIATED(tann)) deallocate(tann)
      if(ASSOCIATED(ymid)) deallocate(ymid)
              write(*,*)'it3dalloc-1.1'
      if(ASSOCIATED(tau)) deallocate(tau)
      if(ASSOCIATED(vptb)) deallocate(vptb)
      if(ASSOCIATED(zboun)) deallocate(zboun)
      if(ASSOCIATED(idx)) deallocate(idx)
      if(ASSOCIATED(imax)) deallocate(imax)
      if(ASSOCIATED(dz)) deallocate(dz)
      if(ASSOCIATED(pol)) deallocate(pol)
      if(ASSOCIATED(solrz)) deallocate(solrz)
      if(ASSOCIATED(solzz)) deallocate(solzz)
      if(ASSOCIATED(thtab)) deallocate(thtab)
      if(ASSOCIATED(z)) deallocate(z)
      if(ASSOCIATED(zmid)) deallocate(zmid)
      if(ASSOCIATED(bbpsi)) deallocate(bbpsi)
      if(ASSOCIATED(bpolz)) deallocate(bpolz)
      if(ASSOCIATED(btorz)) deallocate(btorz)
      if(ASSOCIATED(consnp)) deallocate(consnp)
      if(ASSOCIATED(ptime)) deallocate(ptime)
      if(ASSOCIATED(pefld)) deallocate(pefld)
      if(ASSOCIATED(rovsp)) deallocate(rovsp)
      if(ASSOCIATED(restp)) deallocate(restp)
      if(ASSOCIATED(restnp)) deallocate(restnp)
      if(ASSOCIATED(vpov)) deallocate(vpov)
      if(ASSOCIATED(es)) deallocate(es)
      if(ASSOCIATED(bpsi)) deallocate(bpsi)
      if(ASSOCIATED(d2bpsi)) deallocate(d2bpsi)
              write(*,*)'it3dalloc-1.2'
      if(ASSOCIATED(d2solrz)) deallocate(d2solrz)
      if(ASSOCIATED(d2solzz)) deallocate(d2solzz)
      if(ASSOCIATED(d2bpolz)) deallocate(d2bpolz)
      if(ASSOCIATED(d2btorz)) deallocate(d2btorz)
      if(ASSOCIATED(d2thtpol)) deallocate(d2thtpol)
      if(ASSOCIATED(d2es)) deallocate(d2es)
      if(ASSOCIATED(thtpol)) deallocate(thtpol)
      if(ASSOCIATED(esfi)) deallocate(esfi)
      if(ASSOCIATED(psiesfi)) deallocate(psiesfi)
      if(ASSOCIATED(psifi)) deallocate(psifi)
      if(ASSOCIATED(espsifi)) deallocate(espsifi)
      if(ASSOCIATED(soupp)) deallocate(soupp)
      if(ASSOCIATED(pcurra)) deallocate(pcurra)
      if(ASSOCIATED(pdenra)) deallocate(pdenra)
      if(ASSOCIATED(pfdenra)) deallocate(pfdenra)
      if(ASSOCIATED(pfcurra)) deallocate(pfcurra)
      if(ASSOCIATED(pucrit)) deallocate(pucrit)
      if(ASSOCIATED(peoe0)) deallocate(peoe0)
      if(ASSOCIATED(psrc)) deallocate(psrc)
      if(ASSOCIATED(peoed)) deallocate(peoed)
      if(ASSOCIATED(cint2)) deallocate(cint2)
      if(ASSOCIATED(dx)) deallocate(dx)
      if(ASSOCIATED(dxi)) deallocate(dxi)
      if(ASSOCIATED(ifp)) deallocate(ifp)
      if(ASSOCIATED(sg)) deallocate(sg)
      if(ASSOCIATED(sgx)) deallocate(sgx)
      if(ASSOCIATED(sgxx)) deallocate(sgxx)
      if(ASSOCIATED(sh)) deallocate(sh)
      if(ASSOCIATED(shx)) deallocate(shx)
      if(ASSOCIATED(shxx)) deallocate(shxx)
      if(ASSOCIATED(shxxx)) deallocate(shxxx)
              write(*,*)'it3dalloc-1.3'
      if(ASSOCIATED(tam1)) deallocate(tam1)
      if(ASSOCIATED(tam2)) deallocate(tam2)
      if(ASSOCIATED(tam3)) deallocate(tam3)
      if(ASSOCIATED(tam4)) deallocate(tam4)
      if(ASSOCIATED(tam5)) deallocate(tam5)
      if(ASSOCIATED(tam6)) deallocate(tam6)
      if(ASSOCIATED(tam7)) deallocate(tam7)
      if(ASSOCIATED(tam8)) deallocate(tam8)
      if(ASSOCIATED(tam9)) deallocate(tam9)
      if(ASSOCIATED(tam10)) deallocate(tam10)
      if(ASSOCIATED(tam11)) deallocate(tam11)
      if(ASSOCIATED(tam12)) deallocate(tam12)
      if(ASSOCIATED(tam13)) deallocate(tam13)
      if(ASSOCIATED(tam14)) deallocate(tam14)
      if(ASSOCIATED(tam15)) deallocate(tam15)
      if(ASSOCIATED(tam16)) deallocate(tam16)
      if(ASSOCIATED(tam17)) deallocate(tam17)
      if(ASSOCIATED(tam18)) deallocate(tam18)
      if(ASSOCIATED(tam19)) deallocate(tam19)
      if(ASSOCIATED(tam20)) deallocate(tam20)
      if(ASSOCIATED(tam21)) deallocate(tam21)
      if(ASSOCIATED(tam22)) deallocate(tam22)
      if(ASSOCIATED(tam23)) deallocate(tam23)
      if(ASSOCIATED(tam24)) deallocate(tam24)
      if(ASSOCIATED(tam25)) deallocate(tam25)
      if(ASSOCIATED(tam26)) deallocate(tam26)
      if(ASSOCIATED(tam27)) deallocate(tam27)
      if(ASSOCIATED(tam28)) deallocate(tam28)
      if(ASSOCIATED(tam29)) deallocate(tam29)
      if(ASSOCIATED(tam30)) deallocate(tam30)
      if(ASSOCIATED(x)) deallocate(x)
              write(*,*)'it3dalloc-1.4'
      if(ASSOCIATED(xmidpt)) deallocate(xmidpt)
      if(ASSOCIATED(xi)) deallocate(xi)
      if(ASSOCIATED(xsq)) deallocate(xsq)
      if(ASSOCIATED(x3i)) deallocate(x3i)
      if(ASSOCIATED(x2i)) deallocate(x2i)
      if(ASSOCIATED(xcu)) deallocate(xcu)
      if(ASSOCIATED(xcenter)) deallocate(xcenter)
      if(ASSOCIATED(xcensq)) deallocate(xcensq)
      if(ASSOCIATED(xcent3)) deallocate(xcent3)
      if(ASSOCIATED(uoc)) deallocate(uoc)
      if(ASSOCIATED(enerkev)) deallocate(enerkev)
      if(ASSOCIATED(gamma)) deallocate(gamma)
      if(ASSOCIATED(gamsqr)) deallocate(gamsqr)
      if(ASSOCIATED(gamcub)) deallocate(gamcub)
      if(ASSOCIATED(gammi)) deallocate(gammi)
      if(ASSOCIATED(gamm2i)) deallocate(gamm2i)
      if(ASSOCIATED(gamm1)) deallocate(gamm1)
      if(ASSOCIATED(tcsgm1)) deallocate(tcsgm1)
      if(ASSOCIATED(gamefac)) deallocate(gamefac)
      if(ASSOCIATED(ident)) deallocate(ident)
      if(ASSOCIATED(temc1)) deallocate(temc1)
      if(ASSOCIATED(temc2)) deallocate(temc2)
      if(ASSOCIATED(temc3)) deallocate(temc3)
      if(ASSOCIATED(temc4)) deallocate(temc4)
      if(ASSOCIATED(itemc1)) deallocate(itemc1)
      if(ASSOCIATED(itemc2)) deallocate(itemc2)
              write(*,*)'it3dalloc-1.5'
      if(ASSOCIATED(l_lower)) deallocate(l_lower)
      if(ASSOCIATED(lpt)) deallocate(lpt)
      if(ASSOCIATED(mun)) deallocate(mun)
      if(ASSOCIATED(fll)) deallocate(fll,STAT=istat)
              write(*,*)'it3dalloc fll', istat
      if(ASSOCIATED(xpar)) deallocate(xpar,STAT=istat)
              write(*,*)'it3dalloc xpar', istat
      if(ASSOCIATED(rheads)) deallocate(rheads,STAT=istat)
              write(*,*)'it3dalloc rheads', istat
      if(ASSOCIATED(dfvlle)) deallocate(dfvlle,STAT=istat)
              write(*,*)'it3dalloc dfvlle', istat
      if(ASSOCIATED(dfvlli)) deallocate(dfvlli,STAT=istat)
              write(*,*)'it3dalloc dfvlli', istat
!BH180430: Next statement causing Seg fault. Don't see why??
!      if(ASSOCIATED(xperp)) deallocate(xperp,STAT=istat)
!              write(*,*)'it3dalloc xperp', istat
!              write(*,*)xl,jmaxxl
!      if(ASSOCIATED(xl)) deallocate(xl,STAT=istat)
!              write(*,*)'it3dalloc xl', istat
!!      if(ASSOCIATED(jmaxxl)) deallocate(jmaxxl,STAT=istat)
!!              write(*,*)'it3dalloc jmaxxl', istat

!      if(ASSOCIATED(xlm)) deallocate(xlm,STAT=istat)
!              write(*,*)'it3dalloc xlm', istat
!      if(ASSOCIATED(dxl)) deallocate(dxl,STAT=istat)
!              write(*,*)'it3dalloc dxl', istat
!      if(ASSOCIATED(fl)) deallocate(fl,STAT=istat)
!              write(*,*)'it3dalloc fl', istat
!      if(ASSOCIATED(fl1)) deallocate(fl1,STAT=istat)
!              write(*,*)'it3dalloc fl1', istat
! YuP: A problem with deallocation: sometimes ok, sometimes
! crashes (using the same exe file).
!      if(ASSOCIATED(fl2)) deallocate(fl2,STAT=istat)
!          write(*,*)'it3dalloc fl2', istat

      if(ASSOCIATED(pparea)) deallocate(pparea)
      if(ASSOCIATED(faci)) deallocate(faci)
      if(ASSOCIATED(pprps)) deallocate(pprps)
      if(ASSOCIATED(ppars)) deallocate(ppars)
      if(ASSOCIATED(dff)) deallocate(dff)
      if(ASSOCIATED(cthta)) deallocate(cthta)
      if(ASSOCIATED(cah)) deallocate(cah)
      if(ASSOCIATED(xm)) deallocate(xm)
      if(ASSOCIATED(so)) deallocate(so)
      if(ASSOCIATED(gon)) deallocate(gon)
      if(ASSOCIATED(dbb)) deallocate(dbb)
      if(ASSOCIATED(item1)) deallocate(item1)
      if(ASSOCIATED(item2)) deallocate(item2)
      if(ASSOCIATED(item3)) deallocate(item3)
      if(ASSOCIATED(item4)) deallocate(item4)
      if(ASSOCIATED(item5)) deallocate(item5)
      if(ASSOCIATED(item6)) deallocate(item6)
      if(ASSOCIATED(rhs)) deallocate(rhs)
              write(*,*)'it3dalloc-1.7'

      if(ASSOCIATED(da)) deallocate(da)
      if(ASSOCIATED(db)) deallocate(db)
      if(ASSOCIATED(dc)) deallocate(dc)
      if(ASSOCIATED(dd)) deallocate(dd)
      if(ASSOCIATED(de)) deallocate(de)
      if(ASSOCIATED(df)) deallocate(df)
      if(ASSOCIATED(ca)) deallocate(ca)
      if(ASSOCIATED(cb)) deallocate(cb)
      if(ASSOCIATED(cc)) deallocate(cc)
      if(ASSOCIATED(ce)) deallocate(ce)
      if(ASSOCIATED(cd)) deallocate(cd)
      if(ASSOCIATED(cf)) deallocate(cf)
              write(*,*)'it3dalloc-1.8'

      if(ASSOCIATED(bqlm)) deallocate(bqlm)
      if(ASSOCIATED(tem1)) deallocate(tem1)
      if(ASSOCIATED(tem2)) deallocate(tem2)
      if(ASSOCIATED(tem3)) deallocate(tem3)
      if(ASSOCIATED(tem4)) deallocate(tem4)
      if(ASSOCIATED(tem5)) deallocate(tem5)
      if(ASSOCIATED(tem6)) deallocate(tem6)
      if(ASSOCIATED(xhead)) deallocate(xhead)
      if(ASSOCIATED(xtail)) deallocate(xtail)
      if(ASSOCIATED(ytail)) deallocate(ytail)
      if(ASSOCIATED(yhead)) deallocate(yhead)
      if(ASSOCIATED(fpn)) deallocate(fpn)
        write(*,*)'it3dalloc-1.9'
      if(ASSOCIATED(dyi)) deallocate(dyi)
      if(ASSOCIATED(pleg)) deallocate(pleg)
      if(ASSOCIATED(tfl)) deallocate(tfl)
      if(ASSOCIATED(tbl)) deallocate(tbl)
      if(ASSOCIATED(tal)) deallocate(tal)
      if(ASSOCIATED(currvs)) deallocate(currvs)
      if(ASSOCIATED(dxm5)) deallocate(dxm5)
      if(ASSOCIATED(exm5)) deallocate(exm5)
      if(ASSOCIATED(dxp5)) deallocate(dxp5)
      if(ASSOCIATED(exp5)) deallocate(exp5)
      if(ASSOCIATED(pm)) deallocate(pm)
      if(ASSOCIATED(cog)) deallocate(cog)
      if(ASSOCIATED(choose)) deallocate(choose)
      if(ASSOCIATED(constp)) deallocate(constp)
      if(ASSOCIATED(sigmtt)) deallocate(sigmtt)
      if(ASSOCIATED(sigftt)) deallocate(sigftt)
      if(ASSOCIATED(xlndnz)) deallocate(xlndnz)
      if(ASSOCIATED(ix1)) deallocate(ix1)
      if(ASSOCIATED(ix2)) deallocate(ix2)
      if(ASSOCIATED(ix3)) deallocate(ix3)
      if(ASSOCIATED(ix4)) deallocate(ix4)
      if(ASSOCIATED(ix5)) deallocate(ix5)
      if(ASSOCIATED(ix6)) deallocate(ix6)
      if(ASSOCIATED(ix7)) deallocate(ix7)
      if(ASSOCIATED(ix8)) deallocate(ix8)
      if(ASSOCIATED(tom1)) deallocate(tom1)
      if(ASSOCIATED(tom2)) deallocate(tom2)
      if(ASSOCIATED(tom3)) deallocate(tom3)
      if(ASSOCIATED(tom4)) deallocate(tom4)
      if(ASSOCIATED(fctrl)) deallocate(fctrl)

!MG end added by  11/13/2017

      if(ASSOCIATED(delta_bdb0)) then
        deallocate(delta_bdb0)
      endif


!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)'de_alloc:  DONE deallocating'
!MPIINSERT_ENDIF_RANK

      return
    end subroutine de_alloc


end module impavnc0_mod
