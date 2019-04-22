module tdtravct_mod

!
!

contains

      subroutine tdtravct(f1,kprofile,ktransp)
      use param_mod
      use comm_mod
      use r8subs_mod, only : cvmgt
      implicit integer (i-n), real*8 (a-h,o-z)

!..............................................................
!     This routine calculates the advective radial transport
!     coefficient, d_r(i,j,ktransp,l), i.e., radial pinch
!     velocity (positive in dirn of increasing rho).
!     For pinch="simple", the advective coeff (and the radial
!        diff coeff) are independent of velocity.
!     For pinch="case1", the advective coeff and the
!        radial diffusion coeff are velocity independent.
!        Include geometry+relaxation, not included in "simple".
!     For pinch="case2", the advective coeff is indep of
!        velocity but the radial diffusion coeff can be
!        velocity dependent.
!     Else:  we assume that the advective coeff has the same
!        velocity dependence as the radial diffusion coeff, i.e.,
!     d_r=adv(l)*d_rr. The coefficient is computed at
!     mid-mesh radial points, i.e., at radial bin boundaries,
!     and is determined so that the density profile of
!     species ktransp is forced to equal the profile
!     reden(kprofile,l). The procedure for attaining the desired
!     density is necessarily iterative (over time step n) and
!     a relaxation coefficient relaxden (0<relaxden.le.1) may
!     be utilized. The "old" distribution is f1 (defined on the
!     transport velocity mesh).
!..............................................................


      dimension f1(0:iyp1,0:jxp1,ngen,0:*)


      dimension xx(lrza)
      logical check
      character*8 ifirst
      data ifirst /"first"/

      f2(i,j,k,l)=f1(i,j,k,l)*vptb_(i,l)/zmaxpsi(l)

!     Following ytr and xtr3 statement functions are consistent with
!     trans.h, except xtr3 uses d_rr in place of d_r (as in xtr()),
!     for case3nn below.
!     ztra(i,j,l) is ztrp(i,l)/zmaxpsi*d^3u [cynt2_=cynt2 for ipacktp=0].
      ztra(i,j,l)=cosovb(idx(i,l),l)/dvol(lrindx(l))*4.*pi**2*radmaj* &
        cynt2_(i,l)*cint2(j)
      ytr(i,l)=h_r(lrindx(l))*drrt(k)*d_rr(idx(i,l),j,k,indxlr(l))/ &
        drp5(lrindx(l))*bovcos(idx(i,l),l)
      xtr3(i,l)=-h_r(lrindx(l))*drrt(k)*d_rr(idx(i,l),j,k,indxlr(l))* &
        bovcos(idx(i,l),l)
      if (lrz.ne.lrzmax) stop 'tdtravct: Fixes reqd for lrz.ne.lrzmax'

!..............................................................
!     Compute two integrals that are defined on mid-mesh points.
!..............................................................

      if (pinch.eq.'disabled')  return
      if (ngen.gt.1) stop 'STOP in tdtravct: ngen.gt.1'

      if (ifirst.eq."first") then
         call bcast(d_r,zero,iyjx2*ngen*(lrz+1))
         ifirst="notfirst"
      endif

      dtfactor=1.e0

!     Set target density at the beginning of transport.
!     If colmodl=0, set it equal to initial density of the
!     general species.
!     (This can get generalized to a time dependent target density,
!      or to multiple species.)
      kp=kprofile
      if (n.eq.nontran .and. colmodl.eq.0) then
         kp=1
         do l=1,lrz
            dentarget(l)=reden(kp,l)
         enddo
      else
         do l=1,lrz
            dentarget(l)=reden(kp,l)
         enddo
      endif

!$$$c     Reducing target density by error on previous step
!$$$      do l=1,lrz
!$$$         dentarget(l)=dentarget(l)-
!$$$     +        relaxden*(reden(ktransp,l)-dentarget(l))
!$$$      enddo


      k=ktransp

!      write(*,*)'tdtravct:dtr,dtrtr,dtr_eff',dtr,dtrtr,dtreff
!      write(*,*)'tdtravctt:kp,k,dentarget(1),reden(k,1)',
!     +                    kp,k,dentarget(1),reden(k,1)

!-----------------------------------------------------------------------
!    Smoothing of density and distn in radius introduced, to avoid a
!    numerical instability which occurs with non-iterated pinch term,
!    pinch=case1,case2,case3.   YuP, 091209.
!-----------------------------------------------------------------------

      if (pinch.eq.'case1' .or. pinch.eq.'case2' &
                           .or. pinch.eq.'case3') then

      do i=1,iytr(lrz)
      do j=1,jx
         do l=2,lrz-1 !-YuP: smoothening the f1 radial profile:
            tr1(l)=(f1(i,j,k,l-1) + f1(i,j,k,l) + f1(i,j,k,l+1))/3.d0
         enddo !-YuP: tr1 here is used as a working array
         tr1(1)=f1(i,j,k,1)
         tr1(lrz)=f1(i,j,k,lrz)
         do l=1,lrz
!-YuP            f1(i,j,k,l)=tr1(l) !-YuP: f1 now has a smooth radial profile
!-YuP-100326: In some cases, smoothening is not good: MST tests.
         enddo
      enddo
      enddo

      do l=2,lrz-1 !-YuP: smoothening the density profile:
         tr1(l)=(reden(k,l-1) + reden(k,l) + reden(k,l+1))/3.d0
      enddo !-YuP: tr1 here is used as a working array
      tr1(1)=reden(k,1)
      tr1(lrz)=reden(k,lrz)
      do l=1,lrz
!-YuP         reden(k,l)=tr1(l) !-YuP: reden now is the smooth density profile
!-YuP-100326: In some cases, smoothening is not good: MST tests.
      enddo

      endif  ! On pinch


!=======================================================================
      if (pinch.eq."simple".or. pinch.eq."simplen") then
!=======================================================================

!       Check diff coeff is independent of vel.
        do ii=1,4
           if (difus_vshape(ii).ne.0.) stop 'tdtravct:pinch=simple'
        enddo
        do 160 l=1,lrz-1
          tr1(l)= reden(1,l)                !-YuP: definition changed
          tr2(l)= (reden(1,l+1)-reden(1,l-1)) / &
                      ( drp5(l)+drp5(l-1) ) !-YuP: definition of dn/dr changed
          adv(k,l)=tr2(l)/tr1(l) *advectr   !-YuP: sign reversed
          do 170 i=1,iy_(l)
            do 180 j=1,jx
              d_r(i,j,k,l)=adv(k,l)*drrt(k)*d_rr(i,j,k,l)
 180        continue ! j
 170      continue ! i
 160    continue ! l

!=======================================================================
      elseif (pinch.eq."case1" .or. pinch.eq."case1n") then
!=======================================================================

!       Check diff coeff is independent of vel.
        do ii=1,4
           if (difus_vshape(ii).ne.0.) &
                stop 'tdtravct:pinch=case1, chx difus_vshape'
        enddo
        do l=1,lrz-1
            tr3(l)=(dentarget(l)-reden(k,l))*relaxden*dvol(l)/ &
                 (4.*pi**2*radmaj)/dtreff*advectr  !-YuP: *advectr
            tr1(l)=h_r(l)*drrt(k)*d_rr(1,1,k,l)* &
                 0.5*(reden(k,l)+reden(k,l+1))
            tr2(l)=h_r(l)*drrt(k)*d_rr(1,1,k,l)* &
                   (reden(k,l+1)-reden(k,l))/drp5(l)*advectr !-YuP: *advectr
!           Other version (cancelling-out some terms. But:~same result):
!           tr3(l)=(dentarget(l)-reden(k,l))*relaxden*drp5(l)*advectr !-YuP: *advectr
!           tr1(l)=d_rr(1,1,k,l)*0.5*(reden(k,l)+reden(k,l+1))
!           tr2(l)=d_rr(1,1,k,l)*
!     +            (reden(k,l+1)-reden(k,l))/drp5(l)*advectr !-YuP: *advectr
         enddo

!..............................................................
!     Now recursively compute the advection scale factor.
!..............................................................

         adv(k,1)=(-tr3(1)+tr2(1))/tr1(1)   !-YuP: signs corrected
         do l=2,lrz-1                       !-YuP: signs corrected:
           adv(k,l)=(-tr3(l)+tr2(l)+adv(k,l-1)*tr1(l-1)-tr2(l-1))/tr1(l)
           !!! adv(k,l)=(-tr3(l)+tr2(l))/tr1(l) ! ~same result
         enddo
!         write(*,'(a,8x,a,7x,a,5x,a,5x,a,7x,a,7x,a)')
!     +   'tdtravct: tr1','tr2','tr3','dentarget','reden','h_r','d_r'
         do l=1,lrz-1
            do i=1,iy_(l)
               do j=1,jx
                  d_r(i,j,k,l)=adv(k,l)*drrt(k)*d_rr(i,j,k,l)
               enddo ! j
            enddo ! i
!          write(*,'(i4,7e11.3)')
!     +    l, tr1(l),tr2(l),tr3(l), dentarget(l),reden(k,l),
!     +    h_r(l), d_r(iy/4,jx/2,k,l)
         enddo ! l

!=======================================================================
      elseif (pinch.eq."case2" .or. pinch.eq."case2n") then
!=======================================================================

!***NEED TO FIX THIS CASE TO TAKE D_RR OUT OF ADV.***
!        if (pinch.eq."case2n") stop 'tdtravct: Dont forget more work'
         do l=1,lrz-1
            lp=l+1
            tr1(l)=0.
            tr2(l)=0.
            call bcast(tam1,zero,jx)
            call bcast(tam2,zero,jx)
            do  i=1,iytr(lrz)
               id=idx(i,l)
               ie=idx(i,l+1)
               if (id.eq.0 .or. ie.eq.0) go to 21
               do j=1,jx
                  tam1(j)=tam1(j)+(dl(id,j,k,l)*f2(id,j,k,l) &
                       *cynt2_(id,l) &
                       *cosovb(id,l)+cosovb(ie,lp)* &
                       (1.-dl(id,j,k,l))*f2(ie,j,k,lp)*cynt2_(ie,lp)) &
                       *bovcos(id,l)
                  tam2(j)=tam2(j)+(cosovb(ie,lp)*f2(ie,j,k,lp) &
                       *cynt2_(ie,lp) &
                       -cosovb(id,l)*f2(id,j,k,l)*cynt2_(id,l)) &
                       *bovcos(id,l) &
                       /drp5(l)*drrt(k)*d_rr(id,j,k,l)
               enddo ! j
 21            continue
            enddo ! i
            do  j=1,jx
               tr1(l)=tr1(l)+tam1(j)*cint2(j)
               tr2(l)=tr2(l)+tam2(j)*cint2(j)*advectr  !-YuP: *advectr
            enddo ! j
         enddo ! l

!..............................................................
!     Compute the number of particles (modulo 4*pi**2*radmaj)
!     which need to be gained this time step at mesh point "l".
!     Relaxation coefficient slows down the "correction" term.
!..............................................................

         do l=1,lrz-1
            tr3(l)=(dentarget(l)-reden(k,l))*relaxden*dvol(l)/ &
                 (4.*pi**2*radmaj)/dtreff*advectr    !-YuP: *advectr
         enddo

!..............................................................
!     Now recursively compute the advection scale factor.
!..............................................................

         adv(k,1)=(h_r(1)*tr2(1)-tr3(1))/(h_r(1)*tr1(1))
         do l=2,lrz-1
          adv(k,l)=(-tr3(l)+h_r(l)*tr2(l)-h_r(l-1)*(tr2(l-1) &
            -adv(k,l-1)*tr1(l-1)))/(h_r(l)*tr1(l))
         enddo
!         write(*,'(a,8x,a,7x,a,5x,a,5x,a,7x,a,7x,a)')
!     +   'tdtravct: tr1','tr2','tr3','dentarget','reden','h_r','d_r'
         do l=1,lrz-1
            do i=1,iy_(l)
               do j=1,jx
                  d_r(i,j,k,l)=adv(k,l)
               enddo ! j
            enddo ! i
!          write(*,'(i4,7e11.3)')
!     +    l, tr1(l),tr2(l),tr3(l), dentarget(l),reden(k,l),
!     +    h_r(l), d_r(iy/4,jx/2,k,l)
         enddo ! l

!=======================================================================
      elseif (pinch.eq."case3" .or. pinch.eq."case3n") then
!=======================================================================

        do 10 l=1,lrz-1
          lp=l+1
          tr1(l)=0.
          tr2(l)=0.
          call bcast(tam1,zero,jx)
          call bcast(tam2,zero,jx)
          do 20 i=1,iytr(lrz)
            id=idx(i,l)
            ie=idx(i,l+1)
            if (id.eq.0 .or. ie.eq.0) go to 20
            do 30 j=1,jx
              tam1(j)=tam1(j)+(dl(id,j,k,l)*f2(id,j,k,l)*cynt2_(id,l) &
                *cosovb(id,l)+cosovb(ie,lp)* &
                (1.-dl(id,j,k,l))*f2(ie,j,k,lp)*cynt2_(ie,lp)) &
                *drrt(k)*d_rr(id,j,k,l)*bovcos(id,l)
              tam2(j)=tam2(j)+(cosovb(ie,lp)*f2(ie,j,k,lp)*cynt2_(ie,lp) &
                -cosovb(id,l)*f2(id,j,k,l)*cynt2_(id,l))*bovcos(id,l) &
                /drp5(l)*drrt(k)*d_rr(id,j,k,l)
 30         continue !on j
 20       continue !on i
          do 40 j=1,jx
            tr1(l)=tr1(l)+tam1(j)*cint2(j)
            tr2(l)=tr2(l)+tam2(j)*cint2(j)*advectr  !-YuP: *advectr
 40       continue !on j
 10     continue !on l

!..............................................................
!     Compute the number of particles (modulo 4*pi**2*radmaj)
!     which need to be gained this time step at mesh point "l".
!     Relaxation coefficient slows down the "correction" term.
!..............................................................

        do 50 l=1,lrz-1
!BH010423  Didn't like relaxden coeff being a divisor.
!BH010423  (Not a big problem as most runs now using relaxden=1.)
!BH010423          tr3(l)=(reden(kp,l)-reden(k,l))/relaxden*dvol(l)/
          tr3(l)=(dentarget(l)-reden(k,l))*relaxden*dvol(l)/ &
            (4.*pi**2*radmaj)/dtreff*dtfactor*advectr  !-YuP: *advectr
!        write(*,*)'tdtravct: l,dentarget(l),reden(k,l),tr3(l)',
!     +                       l,dentarget(l),reden(k,l),tr3(l)
 50     continue

!..............................................................
!     Now recursively compute the advection scale factor.
!..............................................................

       adv(k,1)=(h_r(1)*tr2(1)-tr3(1))/(h_r(1)*tr1(1))
       do l=2,lrz-1
          adv(k,l)=(-tr3(l)+h_r(l)*tr2(l)-h_r(l-1)*(tr2(l-1) &
               -adv(k,l-1)*tr1(l-1)))/(h_r(l)*tr1(l))
       enddo
        do 60 l=1,lrz-1
          do 70 i=1,iy_(l)
            do 80 j=1,jx
              d_r(i,j,k,l)=adv(k,l)*drrt(k)*d_rr(i,j,k,l)
 80         continue
 70       continue
 60     continue


!=======================================================================
      elseif (soln_method.eq.'it3drv' .and. pinch.eq."case3nn") then
!=======================================================================

!BH080503:  For soln_method=it3drv,  idx(i,l_)=i.   Using idx
!BH080503:  below for legacy purposes, although maybe not needed.
       do l=1,lrz-1
          lp=l+1
          tr1(l)=0.
          tr2(l)=0.
          call bcast(tam1,zero,jx)
          call bcast(tam2,zero,jx)
          do i=1,iytr(lrz)
             id=idx(i,l)
             ie=idx(i,l+1)
             if (id.eq.0 .or. ie.eq.0) go to 223
             do j=1,jx
                tam1(j)=tam1(j)+ztra(id,j,l)*xtr3(i,l)*(dl(id,j,k,l)* &
                     f2(id,j,k,l)+(1.-dl(id,j,k,l))*f2(ie,j,k,lp))
                tam2(j)=tam2(j)+ztra(id,j,l)*ytr(i,l)* &
                     (f2(ie,j,k,lp)-f2(id,j,k,l))
             enddo              !on j
 223         continue
          enddo                 !on i
          do j=1,jx
             tr1(l)=tr1(l)+tam1(j)
             tr2(l)=tr2(l)+tam2(j)
          enddo                 !on j
       enddo                    !on l

!..............................................................
!     Compute the particle density which is desired
!     to be gained this time step at mesh point "l".
!     Relaxation coefficient slows down the "correction" term.
!..............................................................

       do l=1,lrz-1
          tr3(l)=(dentarget(l)-reden(k,l))*relaxden/dtreff
!          write(*,*)'tdtravct: l,dentarget(l),reden(k,l),tr3(l)',
!     +         l,dentarget(l),reden(k,l),tr3(l)
       enddo

!..............................................................
!     Now recursively compute the advection scale factor.
!..............................................................

       adv(k,1)=(tr3(1)-tr2(1))/tr1(1)
       do l=2,lrz-1
          adv(k,l)=(tr3(l)-tr2(l) +(tr2(l-1)+adv(k,l-1)*tr1(l-1))) &
               /(h_r(l)*tr1(l))
       enddo
!         write(*,'(a,8x,a,7x,a,5x,a,5x,a,7x,a,7x,a)')
!     +   'tdtravct: tr1','tr2','tr3','dentarget','reden','h_r','d_r'
       do l=1,lrz-1
          do i=1,iy_(l)
             do j=1,jx
                d_r(i,j,k,l)=adv(k,l)*drrt(k)*d_rr(i,j,k,l)
             enddo ! j
          enddo ! i
          write(*,'(i4,7e11.3)') &
          l, tr1(l),tr2(l),tr3(l), dentarget(l),reden(k,l), &
          h_r(l), d_r(iy/4,jx/2,k,l)
       enddo ! l

!=======================================================================
      endif                     ! on pinch
!=======================================================================


      do l=1,lrz-1
!        write(*,*)'tdtravct: l,dentarget(l),reden(k,l),tr1(l),tr2(l)'//
!     +    ',tr3(l),adv(k,l)',
!     +    l,dentarget(l),reden(k,l),tr1(l),tr2(l),tr3(l),adv(k,l)
      enddo

      do l=1,lrz
         do j=1,jx
            do i=1,iy
               fxsp(i,j,1,l)=f2(i,j,1,l)
            enddo
         enddo
      enddo


!.......................................................................
!     Test if Newton Iteration for d_r to be performed.
!.......................................................................

      if (pinch.ne."simplen" .and. &
          pinch.ne."case1n" .and. &
          pinch.ne."case2n" .and. &
          pinch.ne."case3n" ) go to 999

!.......................................................................
!     Newton Iteration section, to improve estimate of radial
!     velocities required to maintain target density profile.
!     Subroutine newt is used to find the approximate solution for
!     adv(k,l) such that (dentarget(l)-reden(k,l)=0.).
!     The calculation of adv at the beginning of this subroutine
!     provides a first guess for the
!     lrz values of adv(k,l).
!     Normalization is applied to the input vector, and
!     the density solution vector.
!     See Numerical Recipes in Fortran, Sec.9.7.
!.......................................................................

      write(*,*) 'tdtravct: [l,adv,d_r,drp5] BEFORE NEWTON ITERATIONS:'

      do l=1,lrz-1
      write(*,'(i4,3e15.5)') l,adv(k,l),d_r(iy/4,jx/2,k,l),drp5(l)
         reden_norm(l)=dentarget(l)
         adv_norm(l)=adv(k,l)
         xx(l)=1.0
      enddo

      reden_norm(lrz)=dentarget(l)
      nn=lrz-1
      call newt(xx,nn,iters,check)

      write(*,*)'tdtravct:[l,adv,d_r,xx] AFTER NEWTON ITERATIONS:',iters

      if (check) then
         stop 'tdtravct: Newton Iteration failed'
      elseif (pinch.eq."simplen" .or. pinch.eq."case1n" .or. &
              pinch.eq."case3n") then
         do l=1,lrz-1
            adv(k,l)=xx(l)*adv_norm(l)
            do i=1,iy_(l)
               do j=1,jx
                  d_r(i,j,k,l)=adv(k,l)*drrt(k)*d_rr(i,j,k,l)
               enddo
            enddo ! i
            write(*,'(i4,3e15.5)') l,adv(k,l),d_r(iy/4,jx/2,k,l),xx(l)
         enddo ! l
      elseif (pinch.eq."case2n") then
         do l=1,lrz-1
            adv(k,l)=xx(l)*adv_norm(l)
            do i=1,iy_(l)
               do j=1,jx
                  d_r(i,j,k,l)=adv(k,l)
               enddo
            enddo ! i
            write(*,'(i4,3e15.5)') l,adv(k,l),d_r(iy/4,jx/2,k,l),xx(l)
         enddo ! l
      endif


 999  return
      end



!=======================================================================
!=======================================================================
      subroutine funcv(nn,xx,ffvec)
      use param_mod
      use comm_mod, only: lrz, ngen, pinch
      use comm_mod, only: adv, adv_norm, d_r, drrt, d_rr
      use comm_mod, only: lrindx,tam2, frn, cynt2, coss, tau
      use comm_mod, only: cint2, zmaxpsi, iy_
      use comm_mod, only: dentarget, reden, reden_norm
      implicit integer (i-n), real*8 (a-h,o-z)
      real*8 :: ffvec(nn)

!.......................................................................
!     This user input subroutine is used by newt, giving the function
!     values fvec(1:nn) resulting for input vector xx(1:nn).
!     In our case:   nn=lrz-1
!                    xx(1:nn)=adv(k,1:lrz-1)
!                    fvec=(dentarget(l)-reden(ktransp,l)),l=1,lrz-1
!.......................................................................

      dimension xx(nn)

      if (nn.ne.(lrz-1)) stop 'funcv: nn.ne.(lrz-1)'

!     If want to transport multiple species, will need to
!     pass species index to here.
      if (ngen.ne.1) stop 'funcv: Have to fix for ngen.ne.1, BH010426'
      k=1

      do l=1,lrz-1

         if (pinch.eq."simplen" .or. pinch.eq."case1n" .or. &
              pinch.eq."case3n") then
            adv(k,l)=xx(l)*adv_norm(l)
            do i=1,iy_(l)
               do j=1,jx
                  d_r(i,j,k,l)=adv(k,l)*drrt(k)*d_rr(i,j,k,l)
               enddo
            enddo

         elseif (pinch.eq."case2n") then
            adv(k,l)=xx(l)*adv_norm(l)
            do i=1,iy_(l)
               do j=1,jx
                  d_r(i,j,k,l)=adv(k,l)
               enddo
            enddo

         endif

      enddo

!      write(*,*)'funcv: xx,l=1,nn:',(xx(l),l=1,nn)
!      write(*,*)'funcv: d_r(1,4,1,l),l=1,nn:',(d_r(1,4,1,l),l=1,nn)

      call tdtransp1              ! Gets new distn in frn on vel mesh

!     Compute density from frn, and form fvec:

      do l=1,nn
         call bcast(tam2,zero,jx)
         lr=lrindx(l)
         do j=1,jx
            do i=1,iy_(l)
            tam2(j)=tam2(j)+frn(i,j,k,l)*cynt2(i,l)* &
              abs(coss(i,l))*tau(i,lr)
            enddo
         enddo

         hn=0.
         do j=1,jx
            hn=hn+tam2(j)*cint2(j)
         enddo

         reden(k,lr)=hn/zmaxpsi(lr)

         ffvec(l)=(dentarget(l)-reden(k,lr))/reden_norm(lr)

!      write(*,*)'funcv: lr,dentarget(l),reden(k,lr),xx(lr),fvec(lr)',
!     +                 lr,dentarget(l),reden(k,lr),xx(lr),fvec(lr)


      enddo

      return
      end


!=======================================================================
!=======================================================================
      subroutine tdtransp1
      use param_mod
      use comm_mod
      use r8subs_mod, only : cvmgt, dcopy
      implicit integer (i-n), real*8 (a-h,o-z)

!..............................................................
!     Time advancement splitting step for transport.
!     This is slightly modified version of subroutine tdtransp,
!     prepared to Newton iteration for velocities to maintain
!     the target densities.
!..............................................................

      include 'trans.h90'
!      data nobind /"disabled"/    !Need BLOCK DATA or following stmt
      nobind="disabled"

!.......................................................................

      call dcopy(iyjx2*ngen*lrors,f(0:iyjx2*ngen*lrors-1,0,1,1),1, &
           fvn(0:iyjx2*ngen*lrors-1,0,1,1),1)

!..............................................................
!     The velocity split for this time step is done and we assume
!     that the results of this split resides in fvn(i,j,k,l) and
!     that this is defined on the velocity or complete mesh. Since
!     the tranport equation is not advanced in the vicinity of the
!     pass/trapped boundary we first interpolate onto a new mesh
!     that that does not include these mesh points, and this is
!     accomplished in a density conserving manner.
!..............................................................

!%OS
!      write(*,*) 'f(i,3,1,1),i=1,iy',(f(i,3,1,1),i=1,iy)
      if (nonadi .eq. 6) then
        call tdtrvtor2(fvn(0,0,1,1),frn(0,0,1,1),vpint,vpint_,1)
      else
        call tdtrvtor(fvn,frn)
      endif
!%OS
!      write(*,*) 'frn(i,3,1,1),i=1,iy',(frn(i,3,1,1),i=1,iy)

!..............................................................
!     Generate advection coefficients - then determine the
!     Chang-Cooper weights (radial direction).
!..............................................................

!BH010426      call tdtravct(frn,kelecm,kelecg)
      call tdtrwtl

!..............................................................
!     Keep a copy for accuracy check later...
!..............................................................

!BH010426      call dcopy(iyjx2*ngen*(lrors+1),frn,1,frn_2,1)

!..............................................................
!     Loop over species index.
!..............................................................

      do 100 k=1,ngen

!..............................................................
!     Loop over momentum
!..............................................................

        do 200 ic=1,iytr(lrors)/2
          if (l_lower(ic).eq.lrors) go to 200
          i_=iytr(lrors)+1-ic
          do 201 ii=1,2
            i=ic
            if (ii.eq.2) i=i_
            lu=l_lower(i)

!..............................................................
!     Initialize the calculation (zero flux) at the innermost
!     radial mesh point seen by a particle with a given adiabatic
!     invariant if meshy="fixed_mu" or at l=1 if meshy="fixed_y".
!     In the event that there is a binding equation ((lpt(i).ne.
!     lrors) meaning a particle with a given "i" is passing at the center
!     (l=1) and at some point (lpt(i_)) becomes trapped) we will
!     also initialize at l=lrors. There will then be two initial
!     sweeps, instead on just one. One sweep will be up in "l" the
!     other down, and both terminate at lpt(i).
!..............................................................


            dtfactor=1.e0     !delr is propto 1/dttr
                              !reduces time step by 1/dtfactor
            do 222 j=1,jx
              fg_(j,ii,lu)=delr(i,lu)*dtfactor/betr(i,lu)
              eg_(j,ii,lu)=alpr(i,lu)/betr(i,lu)
              frn(idx(i,lrors),j,k,lrors)=vptb_(idx(i,lrors), &
                lrindx(lrors))/zmaxpsi(lrindx(lrors)) &
                *frn(idx(i,lrors),j,k,lrors)
              eg_(j,ii,lrors)=0.
              fg_(j,ii,lrors)=frn(idx(i,lrors),j,k,lrors)
 222        continue

!..................................................................
!     Now complete the initial sweep over 1 (or 2) regions.
!..................................................................

            do 240 l=l_lower(i)+1,lpt(i)-1
              do 230 j=1,jx
                eg_(j,ii,l)=alpr(i,l)/(betr(i,l)-gamr(i,l) &
                  *eg_(j,ii,l-1))
                fg_(j,ii,l)=(delr(i,l)*dtfactor+gamr(i,l)*fg_(j,ii,l-1)) &
                  /(betr(i,l)-gamr(i,l)*eg_(j,ii,l-1))
 230          continue
 240        continue
            if (lpt(i).eq.lrors) go to 256

!..................................................................
!     second region, if necessary...
!..................................................................

            do 250 l=lrors-1,lpt(i)+1,-1
              do 255 j=1,jx
                eg_(j,ii,l)=gamr(i,l)/(betr(i,l)-alpr(i,l) &
                  *eg_(j,ii,l+1))
                fg_(j,ii,l)=(delr(i,l)*dtfactor+alpr(i,l) &
                  *fg_(j,ii,l+1)) &
                  /(betr(i,l)-alpr(i,l)*eg_(j,ii,l+1))
 255          continue
 250        continue
 256        continue
 201      continue
!..................................................................
!     The binding equation (at lpt(i)) is next..
!..................................................................


          i=ic
          if (lpt(i).eq.lrors) go to 260
          lv=lpt(i)
          if (nobind.eq."disabled") then

            do 265 j=1,jx
              tam1(j)=-(gamr(i,lv)*eg_(j,1,lv-1)+gamr(i_,lv) &
                *eg_(j,2,lv-1))*.5 &
                +(betr(i,lv)+betr(i_,lv))*.5 - alpr(i,lv)*eg_(j,1,lv+1)
              tam2(j)=(gamr(i,lv)*fg_(j,1,lv-1)+gamr(i_,lv) &
                *fg_(j,2,lv-1))*.5 &
                +alpr(i,lv)*fg_(j,1,lv+1)+delr(i,lv)*dtfactor
              frn(idx(i,lpt(i)),j,k,lpt(i))=tam2(j)/tam1(j)
              frn(idx(i_,lpt(i_)),j,k,lpt(i_))=tam2(j)/tam1(j)
 265        continue
          else
            do 266 j=1,jx
              tam1(j)=-(gamr(i,lv)*eg_(j,1,lv-1)) &
                +(betr(i,lv)) - alpr(i,lv)*eg_(j,1,lv+1)
              tam2(j)=(gamr(i,lv)*fg_(j,1,lv-1)) &
                +alpr(i,lv)*fg_(j,1,lv+1)+delr(i,lv)*dtfactor
              frn(idx(i,lpt(i)),j,k,lpt(i))=tam2(j)/tam1(j)
              tam3(j)=-(gamr(i_,lv)*eg_(j,2,lv-1)) &
                +(betr(i_,lv)) - alpr(i_,lv)*eg_(j,2,lv+1)
              tam4(j)=(gamr(i_,lv)*fg_(j,2,lv-1)) &
                +alpr(i_,lv)*fg_(j,2,lv+1)+delr(i_,lv)*dtfactor
              frn(idx(i_,lpt(i_)),j,k,lpt(i_))=tam4(j)/tam3(j)
 266        continue
          endif
 260      continue

!..................................................................
!     solve for the new distribution...
!..................................................................
          do 285 ij=1,2
            i=ic
            if (ij.eq.2) i=i_
            do 280 l=lpt(i)-1,l_lower(i),-1
              do 270 j=1,jx
                frn(idx(i,l),j,k,l)=eg_(j,ij,l)*frn(idx(i,l+1),j,k,l+1) &
                  +fg_(j,ij,l)
 270          continue
 280        continue
            do 295 l=lpt(i)+1,lrors
              do 290 j=1,jx
                frn(idx(i,l),j,k,l)=eg_(j,ij,l)*frn(idx(i,l-1),j,k,l-1) &
                  +fg_(j,ij,l)
 290          continue
 295        continue
            do 310 l=l_lower(i),lrors
              ii=idx(i,l)
              do 312 j=1,jx
                frn(ii,j,k,l)=frn(ii,j,k,l)/vptb_(ii,lrindx(l))* &
                  zmaxpsi(lrindx(l))
 312          continue
 310        continue
 285      continue
 200    continue
 100  continue

!BH010426      call tdtrchk

      if (nobind.eq."enabled") then
        call tdtrsym
      endif

!.................................................................
!     Redefine frn at v=0 so it is unique.
!.................................................................

!%OS
      if (nonadi .ne. 3) then

        do 400 k=1,ngen
          do 410 l=1,lrors
            zs=0.
            zt=0.
            do 420 i=1,iytr(l)
              zs=zs+vpint_(idx(i,l),lrindx(l))
              zt=zt+vpint_(idx(i,l),lrindx(l))*frn(idx(i,l),1,k,l)
 420        continue
            do 430 i=1,iytr(l)
              frn(idx(i,l),1,k,l)=zt/zs
 430        continue
 410      continue
 400    continue

      endif

!......................................................................
!     The distribution frn is currently defined on the transport
!     mesh - interpolate onto the velocity mesh, returning it in frn.
!......................................................................

      if (nonadi .eq. 6) then
        call tdtrrtov2(frn(0,0,1,1),frn(0,0,1,1),vpint,vpint_,1)
      else
        call tdtrrtov(frn)
      endif

      return
      end subroutine


!=======================================================================
!=======================================================================
      subroutine newt(x,n,iters,check)
      use comm_mod, only : NP, fvec, newtv_nn
!
!..................................................................
!     A Newton-Raphson iteration, here used with subroutine funcv
!     to obtain radial velocities (transp='enabled') which keep
!     plasma density close to a target profile.
!     This algorithm is adapted (slightly) from Numerical Recipes.
!     See section 9.7.
!     From NR:
!     Given an initial guess x(1:n) for a root in n dimensions, find
!     the root by a GLOBALLY CONVERGENT Newton's method.  The vector
!     of functions to be zeroed, fvec(1:n), is returned by a user-
!     supplied subroutine funcv(n,x,fvec).
!     newt output variable  check is "true" for converged to local
!     minimum rather than zero, else "false" for normal return.
!
!..................................................................
!
      INTEGER n,MAXITS, iters
      LOGICAL check
      REAL*8 x(n),TOLF,TOLMIN,TOLX,STPMX
!     Choose following TOLX to match value in lnsrch.
!     Choose following NP to match value in fffmin.
!BH091214      PARAMETER (NP=100,MAXITS=25,TOLF=1.e-2,TOLMIN=1.e-5,TOLX=1.e-3,
      PARAMETER (MAXITS=100,TOLF=1.e-6,TOLMIN=1.e-6,TOLX=1.e-10, &
                 STPMX=100.)

!     USES fdjac,fffmin,lnsrch,lubksb,ludcmp


      INTEGER i,its,j,indx(NP)
      REAL*8 d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP), &
      		g(NP),p(NP),xold(NP),fffmin
      REAL*8 fn,fnp5,one

      one=1.d0
      iters=0
      if (n.gt.NP) stop 'newt: Increase NP.'

      newtv_nn=n
      f=fffmin(x)            ! Calls funcv
!      write(*,*)'newt: fffmin',f
      test=0.
      do   i=1,n
         if(abs(fvec(i)).gt.test)test=abs(fvec(i))
      enddo
!      write(*,*)'newt: fvec',(fvec(i),i=1,n)
      if(test.lt..01*TOLF)return
      sum=0.
      do   i=1,n
         sum=sum+x(i)**2
      enddo
      fn=n
      stpmax=STPMX*max(sqrt(sum),fn)
      do   its=1,MAXITS
         iters=its
         call fdjac(n,x,fvec,NP,fjac)
         do i=1,n
            sum=0.
            do   j=1,n
               sum=sum+fjac(j,i)*fvec(j)
            enddo
            g(i)=sum
         enddo
         do   i=1,n
            xold(i)=x(i)
         enddo
         fold=f
         do   i=1,n
            p(i)=-fvec(i)
         enddo
         call ludcmp(fjac,n,NP,indx,d)
         call lubksb(fjac,n,NP,indx,p)
!         write(*,*)'newt, Before lnsrch: f,x',f,(x(i),i=1,n)
         call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fffmin)
!         write(*,*)'newt, After lnsrch: f,x',f,(x(i),i=1,n)
         test=0.
         do   i=1,n
            if(abs(fvec(i)).gt.test)test=abs(fvec(i))
         enddo
         if(test.lt.TOLF)then
            check=.false.
            return
         endif
         if(check)then
            test=0.
            fnp5=.5*n
            den=max(f,fnp5)
            do   i=1,n
               temp=abs(g(i))*max(abs(x(i)),one)/den
               if(temp.gt.test)test=temp
            enddo
            if(test.lt.TOLMIN)then
               check=.true.
            else
               check=.false.
            endif
            return
         endif
         test=0.
         do   i=1,n
            temp=(abs(x(i)-xold(i)))/max(abs(x(i)),one)
            if(temp.gt.test)test=temp
         enddo
         if(test.lt.TOLX)return
      enddo
      write(*,*) 'MAXITS exceeded in newt'
      END


!=======================================================================
!=======================================================================
      subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
!                                                                      !
!    Description:                                                      !
!                                                                      !
!    Adopted from Numerical Recipes in FORTRAN, Chapter 9.7, 2nd ed.   !
!                                                                      !
!    Given an n-dimensional point XOLD(1:N), the value of the function !
!    and gradient there, FOLD and G(1:N), and a direction P(1:N),      !
!    finds a new point X(1:N) along the direction P from XOLD where    !
!    the function FUNC has decreased 'sufficiently'. The new function  !
!    value is returned in F. STPMAX is an input quantity that limits   !
!    the length of the steps so that you do not try to evaluate the    !
!    function in regions where it is undefined or subject to overflow. !
!    P is usually the Newton direction. The output quantity CHECK is   !
!    false on a normal; exit. It is true when X is too close to XOLD.  !
!    In a minimization algorithm, this usually signals convergence and !
!    can be ignored. However, in a zero-finding algorithm the calling  !
!    program should check whether the convergence is spurious.         !
!                                                                      !
!    Called by:       NEWT                                             !
!                                                                      !
!    Calls:           FUNC                                             !
!                                                                      !
      INTEGER n
      LOGICAL check
      REAL*8 f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
!BH091214      PARAMETER (ALF=1.e-4,TOLX=1.e-3)
      PARAMETER (ALF=1.e-5,TOLX=1.e-10)
      EXTERNAL func
!     USES func

      INTEGER i
      REAL*8 a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope, &
      	       sum,temp,test,tmplam
      REAL*8 zero,one
      zero=0.d0
      one=1.d0
      check=.false.
      sum=zero
      do    i=1,n
         sum=sum+p(i)*p(i)
      enddo
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
         do   i=1,n
            p(i)=p(i)*stpmax/sum
         enddo
      endif
      slope=0.
      do   i=1,n
         slope=slope+g(i)*p(i)
      enddo
      test=0.
      do   i=1,n
         temp=abs(p(i))/max(abs(xold(i)),one)
         if(temp.gt.test)test=temp
      enddo
      alamin=TOLX/test
      alam=1.
 1    continue
      do   i=1,n
         x(i)=xold(i)+alam*p(i)
      enddo
      f=func(x)
      if(alam.lt.alamin)then
         do   i=1,n
            x(i)=xold(i)
         enddo
         check=.true.
         return
      else if(f.le.fold+ALF*alam*slope)then
         return
      else
         if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
         else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/ &
                  (alam-alam2)
            if(a.eq.zero)then
               tmplam=-slope/(2.*b)
            else
               disc=b*b-3.*a*slope
               tmplam=(-b+sqrt(disc))/(3.*a)
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
         endif
      endif
      alam2=alam
      f2=f
      fold2=fold
      alam=max(tmplam,.1*alam)

      goto 1

      END


!=======================================================================
!=======================================================================
      subroutine fdjac(n,x,fvec,NP,df)
      INTEGER n,NP
      REAL*8 df(NP,NP), fvec(n),x(n),EPS
      PARAMETER (EPS=1.e-7)
!     USES funcv
      INTEGER i,j
      REAL*8 h,temp,f(NP)
      do   j=1,n
         temp=x(j)
         h=EPS*abs(temp)
         if(h.eq.0.)h=EPS
         x(j)=temp+h
         h=x(j)-temp
         call funcv(n,x,f)
         x(j)=temp
         do   i=1,n
            df(i,j)=(f(i)-fvec(i))/h
         enddo
      enddo
      return
      END

!=======================================================================
!=======================================================================
      REAL*8 function fffmin(x)
      use comm_mod, only : NP, fvec, newtv_nn
      REAL*8 x(*)
!     USES funcv
      INTEGER i
      REAL*8 sum
      call funcv (newtv_nn,x,fvec)
      sum=0.
      do   i=1,newtv_nn
         sum=sum+fvec(i)**2
      enddo
      fffmin=0.5*sum
      return
      END

!=======================================================================
!=======================================================================
      subroutine ludcmp(a,n,NP,indx,d)
      INTEGER n,NP,indx(n)
      REAL*8 d,a(np,np),TINY
      PARAMETER (TINY=1.0e-20)

      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NP),zero
      save icall_count !-YuP Counting calls to subroutine


      zero=0.d0
      d=1.
      imax=1 ! YuP-2011: added. Will be re-defined.
      do  i=1,n
         aamax=0.
         do   j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         enddo
         if (aamax.eq.zero) then
            icall_count=icall_count+1 !-YuP Counting calls
            write(*,*) 'singular matrix in ludcmp',icall_count
            !pause
            aamax=TINY ! YuP-2011: added
         endif
         vv(i)=1./aamax
      enddo
      do  j=1,n
         do   i=1,j-1
            sum=a(i,j)
            do   k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
         enddo
         aamax=0.
         do i=j,n
            sum=a(i,j)
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo
         if (j.ne.imax)then
            do   k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            enddo
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.zero)a(j,j)=TINY
         if(j.ne.n)then
            dum=1./a(j,j)
            do   i=j+1,n
               a(i,j)=a(i,j)*dum
            enddo
         endif
      enddo
      return
      END

!=======================================================================
!=======================================================================
      subroutine lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum,zero
      zero=0.d0
      ii=0
      do   i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do   j=ii,i-1
               sum=sum-a(i,j)*b(j)
            enddo
         else if (sum.ne.zero) then
            ii=i
         endif
         b(i)=sum
      enddo
      do   i=n,1,-1
         sum=b(i)
         do   j=i+1,n
            sum=sum-a(i,j)*b(j)
         enddo
         b(i)=sum/a(i,i)
      enddo
      return
      END

end module tdtravct_mod
