module tdtranspn_mod

  !---BEGIN USE

  use r8subs_mod, only : dcopy
  use tdtravct_mod, only : tdtravct
  use tdtrvtor_mod, only : tdtrvtor
  use tdtrwtl_mod, only : tdtrwtl

  !---END USE


contains

      subroutine tdtranspn
      use param_mod
      use comm_mod
      use r8subs_mod, only : cvmgt, dcopy

      implicit integer (i-n), real*8 (a-h,o-z)
      dimension  zmatcont(12), janelt(12)
      include 'trans.h90'
      save ifirst
      data ifirst /1/
!..............................................................
!     Radial transport matrix coefficients are generated in
!     CSR format, to be added to the velocity space coefficients
!     for full 2DVel-1Drho matrix prepartory to solution by iterative
!     sparse-matrix methods:    soln_method="it3drv".
!     "l" is the index for the radial mesh bins.
!     Presently, the set up for transport handles only one general
!     species (k) at a time.
!     tdtranspn provides an alternative (new) soln method to the
!     splitting soln implemented in tdtransp.
!     BH071105
!.......................................................................
      call dcopy(iyjx2*ngen*lrors,f(0:iyjx2*ngen*lrors-1,0,1,1),1, &
           fvn(0:iyjx2*ngen*lrors-1,0,1,1),1)
!..............................................................
!     The distribution from the previous time step resides in
!     fvn(i,j,k,l).  This is now (071111) defined on the complete
!     velocity mesh (ipacktp=0,  see tdtrmuy.f).
!     Finite difference equations have been constructed for a radial
!     diffusion + advection equation as given in the cql3d manual.
!     At each time step, the radial advection term is chosen to
!     keep the density of the diffusing particles near to a specified
!     target radial profile [see cqlinput_help and tdtravct.f].
!     For the j=1 (v=0) points, only the theta=pi (i=iy) point
!     is explicity transported, as the other theta points are
!     set equal to f(i=iy,j=1, l) in the velocity-space setting
!     of coefficients (in impanvnc0.f). [lbdry()="conserv"]
!     For each theta point (i) on the mesh at the outermost
!     radial bin (l=lrz), there is a smallest minor radius
!     l_lower(i) with corresponding i mesh point.
!     [Presently, only use ymesh="fixed_y" (not "fixed_mu) is fully
!     implemented, which gives l_lower(i)=1.]
!     Boundary conditions are zero-radial-flux at the inner
!     edge of the rya(l=l_lower(i)) radial bin, and a specified
!     Maxwellian at the outer radius rya(l=lrz).
!     Better BCs might be future mod [BH071107]:
!        Boundary conditions presently are zero-radial-flux at the
!        inner edge of the rya(l=l_lower(i)) radial bin, and a
!        specified Maxwellian at the outer radius rya(l=lrz+1/2)=1.0
!        [The outer bc is slightly different from the splitting method
!        setup, in which the distribution at rya(l=lrz) is taken
!        to be Maxwellian for the radial split (but general for the
!        velocity split).]
!
!     At the outer radius, l=lrz, the trapping fraction is
!     greater than towards the magnetic axis.  Therefore, some of the
!     i-points will transform from trapped particles (for which
!     the distribution is symetric about theta=pi/2) to passing
!     particles.  At these radii, given by l=lpt(i), the trapped
!     particles at l.ge.lpt,  with given i, are connected (bound) to
!     corresponding transitting particle distributions at l=lpt-1
!     in the co-  and counter-dirns, and there is a special equation for
!     f(i,j,k,l=lpt(i)).  At other points in i,j,l space,
!     a simple 3-point difference equation holds connecting distrns
!     at l-1,l,l+1, for each i,j. [However, the passing particles
!     at lpt(i)-1 with i.gt.iyh are connected to trapped particles
!     at lpt(i) with i.lt.iyh.]
!..............................................................
!.......................................................................
!  Copy fvn to frn as in tdtransp, although they are now on same grid
!.......................................................................

      call tdtrvtor(fvn,frn)

!..............................................................
!     Generate advection coefficients - then determine the
!     Chang-Cooper weights (radial direction).
!     Presently only set up for 1 general species,
!       except for pinch.eq."disabled" (adv()=0.).
!..............................................................
      if (ngen.gt.1) then
         write(*,*)'Presently set up for only 1 gen species: STOP'
         stop
      endif
      kprofile=kelecm
      if (niong.ne.0) kprofile=kionm(1)
      if (colmodl.eq.0) kprofile=1
      call tdtravct(frn,kprofile,1)  !Obtain d_r advection coeffs
                                     !No-op if pinch="disabled"
      call tdtrwtl              !  Sets all dl(,,,)=0.5
!..............................................................
!     Keep a copy for accuracy check later... [following tdtransp].
!..............................................................
      call dcopy(iyjx2*ngen*(lrors+1),frn,1,frn_2,1)
!.......................................................................
!.......................................................................
!  Following do loop setup from impavnc0 might be expanded for
!  other cases.  Use same eqn order as impavnc0 in full radius
!  calcs, for compatibility of the resulting CSR matrix.
!.......................................................................
!.......................................................................
!.......................................................................
!  ipofi(i,l)  (Iprime, i.e., indexsww) of i; set up this array.
!  That is, ipofi(i,l) gives the indexsww value corresponding to
!  to each i, i=1,iy.  For i within the neg v-par trapped region,
!  it uses the the indexsww for theta at symmetic points about y=pi/2.
!  The ipofi(i,l) variable facilitates calculation of corresponding
!  (i.e., given i,j)equation numbers displaced by one radial
!  mesh point.
!  It also simplifies radial differencing for each i up to the
!  flux surface where l=lpt(i) [if such condition
!  exists].
!  Recall that i counts counter-clockwise in the vpar-vperp plane,
!    i.e., in the direction of increasing pitch angle theta.
!    ipofi(1,l) will be = inew_(l)-1, ipofi(inew_(l),l)=inew_(l).
!  The indexsww values are as in subroutine impavnc0.
!.......................................................................
      if (ifirst.eq.1) then

         do l=1,lrz
            itl=itl_(l)
            itu=itu_(l)
            iyh=iyh_(l)
            iyy=iy_(l) !-YuP-101215: Don't use iy=; it's in common /params/
                       ! Don't let overwrite the cqlinput value!
            if (itl.lt.iyh ) then
               inew=iyh + itl - 1
               ipassbnd=iyh-itl+1
            else                !i.e., itl=iyh
               inew=iyy
               ipassbnd=0       !i.e., no trapped region for iyh=itl.
                                !      however, this case not checked out.
            endif

            do i = 1,iymax
            ipofi(i,l)=0
            enddo
            ihalf=0
            icntsww=1
            if (iyh.eq.itl) icntsww=0
            do  indexsww=1,inew
               if(indexsww.le.ipassbnd) then
                  i=iyh+1-indexsww
                  ipofi(i,l)=indexsww
                  ipofi(iyy-i+1,l)=indexsww !connects counter-i to co-i
                                           !in the trapped region
               elseif(ihalf.eq.0) then
                  i = itl - icntsww
                  ihalf = 1
                  ipofi(i,l)=indexsww
               else
                  i = itu + icntsww
                  ihalf = 0
                  icntsww = icntsww + 1
                  ipofi(i,l)=indexsww
               endif
            enddo               ! on indexsww
         enddo                  ! on l

         ifirst=0
      endif                     ! on ifirst
!.......................................................................
!  Main loops setting up radial transport coefficients
!.......................................................................
      k=1  ! only one gen species for time being
      ieq=0     !equation counter
      icoeff=0  !CSR coefficient counter

      do l=1,lrz   !loop over flux surfaces

!         write(*,*)'tdtranspn: l,ieq+1,ieq_(l),icoeff',
!     +                         l,ieq+1,ieq_(l),icoeff
         itl=itl_(l)
         itu=itu_(l)
         iyh=iyh_(l)
         iyy=iy_(l)
!.......................................................................
!     define effective number of unknowns in the theta direction and
!     resulting band widths
!     see below for comments on inew, etc.
!.......................................................................
!BH071029:  During work on 3d fully-implicit solve, found a problem
!BH071029:  with iyh=itl cases.  Fixing this case, and moving ipassbnd
!BH071029:  setting to here [See further notes below for ipassbnd]:
!BH071029:  [But subsequently found had problems with iyh=itl elsewhere
!BH071029:  in the code, so suggest bigger rya(1) or iy for time being.]

      if (itl.lt.iyh ) then
         inew=iyh + itl - 1
         ipassbnd=iyh-itl+1
      else  !i.e., itl=iyh
         inew=iyy
         ipassbnd=0  !i.e., no trapped region for iyh=itl.
                     !      however, this case not checked out.
      endif

      do j=1,jx
            ihalf=0
            icntsww=1
            if (iyh.eq.itl) icntsww=0
!..................................................................
!     indexsww is the major loop over I (theta).
!..................................................................
            do indexsww=1,inew
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

              ieq=ieq+1
              iar_csr(ieq)=icoeff+1   !Coeff counter at beginning of row

!             index for scaling of coeffs (saved in impavnc0):
              if (soln_method.eq."itsol" .or. soln_method.eq."itsol1" &
                 .or. soln_method.eq."direct") then
                 kku=ieq+(k-1)*iyjx  !-YuP-> version from impavnc
              else
                 kku=ieq-(ieq_(l)-1)+(k-1)*iyjx
                 ! k.gt.1 not enabled for soln_method it3dv or it3drv
              endif

              if (l.eq.1) then

                 if (j.eq.1) then
                    if (ieq.lt.inew .or. lbdry(k).ne."conserv") then  !0. diagonal element for 1st
                                           !inew-1 equations
                       icoeff=icoeff+1
                       ar_csr(icoeff)=0.
                       jar_csr(icoeff)=ieq
                    else        !i.e., ieq.eq.inew (l=1,i=iyy),
                                !link to 2nd flux surface, i-iyy
                       icoeff=icoeff+1
                       ar_csr(icoeff)=scal(kku,l)*betrp(i,l)
                       jar_csr(icoeff)=ieq
                       icoeff=icoeff+1
                       ar_csr(icoeff)=-scal(kku,l)*alprp(i,l)
                       jar_csr(icoeff)=(ieq_(2)-1)+ipofi(i,2)
                   endif
                 else  ! i.e., j.gt.1
!                   No need to distinguish passing or trapped for l=1.
                    icoeff=icoeff+1
                    ar_csr(icoeff)=scal(kku,l)*betrp(i,l)
                    jar_csr(icoeff)=ieq
                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*alprp(i,l)
                    jar_csr(icoeff)=ieq_(2)-1+(j-1)*inew_(2)+ipofi(i,2)
                 endif


              elseif (l.lt.lrz) then  !i.e., 1<l<lrz
                 if (j.eq.1 .and. ((ieq-(ieq_(l)-1)).lt.inew &
                                .or. lbdry(k).ne."conserv")) then
                       icoeff=icoeff+1
                       ar_csr(icoeff)=0.
                       jar_csr(icoeff)=ieq


                 elseif (l.eq.lpt(i)) then !at inner radius for trapping
                                           !at the given i: couple to passing
                                           !at smaller l.  CHECK OUT MORE!!!
                                           !-YuP: corrected. Was:  elseif (i.eq.lpt(i))
                    ii=iyy+1-i

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*0.5*gamrp(i,l)
                    jar_csr(icoeff)=ieq_(l-1)-1+ &
                                    (j-1)*inew_(l-1)+ipofi(i,l-1)

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*0.5*gamrp(ii,l)
                    jar_csr(icoeff)=ieq_(l-1)-1+ &
                                    (j-1)*inew_(l-1)+ipofi(ii,l-1)

                    icoeff=icoeff+1
                    ar_csr(icoeff)= &
                              scal(kku,l)*0.5*(betrp(i,l)+betrp(ii,l))
                    jar_csr(icoeff)=ieq

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*alprp(i,l)
                    jar_csr(icoeff)=ieq_(l+1)-1+ &
                                    (j-1)*inew_(l+1)+ipofi(i,l+1)
                 else  ! i.e., regular points

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*gamrp(i,l)
                    jar_csr(icoeff)=ieq_(l-1)-1+ &
                                    (j-1)*inew_(l-1)+ipofi(i,l-1)

                    icoeff=icoeff+1
                    ar_csr(icoeff)=scal(kku,l)*betrp(i,l)
                    jar_csr(icoeff)=ieq

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*alprp(i,l)
                    jar_csr(icoeff)=ieq_(l+1)-1+ &
                                    (j-1)*inew_(l+1)+ipofi(i,l+1)

                 endif

              elseif (l.eq.lrz) then  !no transp for now, just vel distn bc.
                                      !Put placeholders in the CSR matrix

                 icoeff=icoeff+1
                 ar_csr(icoeff)=0.0
                 jar_csr(icoeff)=ieq

              endif

!       if(l.eq.10 .and. j.eq.jx/2) then
!       write(*,'(3i9,e13.5)')
!     ~  i,iar_csr(ieq),jar_csr(icoeff),ar_csr(icoeff)
!       endif


           enddo    ! on indexsww

        enddo       ! on j

      enddo         ! on l

!  Final element for CSR storage
      iar_csr(ieq+1)=icoeff+1
!      pause

      return
      end
end module tdtranspn_mod
