module baviorbt_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast
  use comm_mod
  use lookup_mod, only : lookup
  use micgnbnd_mod, only : micgnbnd
  use param_mod
  use psif_mod, only : psif
  use r8subs_mod, only : lug, luf
  use zcunix_mod, only : terp1
  use zcunix_mod, only : terp2

  !---END USE

contains

  subroutine baviorbt
    implicit integer (i-n), real(c_double) (a-h,o-z)
    real(c_double) ::wk_lug(nstps) ! YuP[2019-04-22] working array for lug()
    save

!  *** NOTE: Bunch of write(*,*) statements related to checking
!            out effect of taunew.  For now, might use taunew="disabled".
!            See comments in a_change.h.

!...................................................................
!     New fully numerical calculation of tau and dtau. (BH: 970211).
!     (Consider special integration for orbit intervals with
!      turning points.  See, for example, Numerical Recipes for
!      Fortran, p. 135ff).
!
!     This routine computes the quantity equal to the product of the
!     bounce time with the particle velocity. This depends only on the
!     midplane pitch angle and is stored in tau(i,lr_). The bounce
!     time is an average over dy(i,lr_) for each pitch angle i.
!     [BH091031: But see below, nii=1, i.e., only calc at y(i,l_).]
!     The local average time weights (bounce time increment
!     times in the distance increment dz(l,lr_) evaluated in the nbhd.
!     of z(l,lr_) and y(i,l_)) are also computed here and are stored
!     in dtau(i,l,lr_).  The bounce average of a function g(i,l)
!     is the sum over l from 1 to lz of g(i,l)*dtau(i,l,lr_)/tau(i,lr_).
!     Caution: for trapped particles we must have for some l
!     z(l,lr_) .lt. zboun(i,lr_) .lt. z(l+1,lr_). There will be a
!     nonzero value in dtau(i,l+1,lr_) in general. Thus a reasonable
!     value of g must be stored in g(i,l+1). In this code we use the
!     value of g at the bounce point if possible. Otherwise we set
!     g(i,l+1)=g(i,l).
!     Caution: The above introduces nonzero elements into dtau(l,i,lr_)
!     at i=imax(l,lr_)+1, in some cases (Bobh, 970126).
!     BH,090901:  Might not agree very well with the original
!                 baviorbto below, where intervals considered are
!                 [z_l,z_l+1],l=1,lmax where z_1=0, z_lmax+1=zboun.
!                 In the bounce-average of a function g, g is to be
!                 tabulated on z_l points.
!                 On the other hand, baviorbto does not properly
!                 fill in dtau for i=itl pitch angle (tp-bndry).
!                 Need to check this. (start with mnemonic.nc op).
!...................................................................

!cc      common/temp_imax_old/ imax_old(lza)

!cc      dimension lmax_old(iy) ! for tests/debug only

!cc      real(c_double), dimension(iy):: prnt1,prnt2,prnt3,prnt4 !For gdb printing
!cc      real(c_double), dimension(iy):: prnt5,prnt6             !For gdb printing

!...................................................................
!     ERROR CONTROL - change nstps if bigger than iy*jx, because
!       of use of tem1,tem2
!...................................................................

!     nii is number of pitch angle subintervals used in the
!       calculation of dtau and tau.  nstps is the number of poloidal
!       mesh subintervals (in namelist setup. default=100).
!      nii=25   !VALUE used until 090907
      nii=25 ! nii=1 gives bad profile of j_bs(lr) in bootcalc test
      if (nii*nstps .gt. iyjx-1) then
         nstps=iyjx/nii
         write(*,*) 'WARNING: nstps in baviorbt reduced to',nstps
      endif

      if (nii.gt.iy) then
         nii=iy
         write(*,*) 'WARNING: nii in baviorbt reduced to iy',nii
      endif

!...................................................................
!     Limit the range of integration at the pass/trapped boundary
!...................................................................

!      call bcast (tau(1,lr_),zero,iyh)
      call bcast(tau(1:iy,lr_),zero,iy)
      call bcast(dtau(1:iy*lz,1,lr_),zero,iy*lz)
      call bcast(tem2,zero,iyjx)
!      call ibcast(itemc1,1,iy)    ! Test
      call ibcast(itemc1,0,iy)

!...................................................................
!     Begin loop over the particle orbit, doing upper "half" of
!     the flux surface.
!     If eqsym.eq."none", then lower in second section done below.
!...................................................................

      if (eqsym.ne."none") then
         ilzhfs=lz
         if (numclas .eq. 1) ilzhfs=lz/2+1
      else  !  i.e., non-up-down symmetric case
         ilzhfs=lz_bmax(lr_) ! index in pol.mesh corresponding to Bmax
      endif

      do 70 l=1,ilzhfs  ! index in poloidal mesh.

!...................................................................
!     Evaluate the bbpsi function where it will normally be required for
!     integration. Doing this outside the "60" loop saves time.
!     BH,090901:  This does not square very well with the original
!                 baviorbto below, where intervals considered are
!                 [z_l,z_l+1],l=1,lmax where z_1=0, z_lmax+1=zboun.
!                 In the bounce-average of a function g, g is to be
!                 tabulated on z_l points.
!...................................................................

            nsteps=nstps !-YuP
         if (l.eq.1) then
            z2=0.
            z3=zmid(1,lr_)
            !nsteps=nstps/2 !-YuP: because of half-interval z3-z2
         else
            z2=zmid(l-1,lr_)
            z3=zmid(l,lr_)
         endif
         if (l.eq.ilzhfs) then
            z2=zmid(ilzhfs-1,lr_)
            z3=z_bmax(lr_)   !-YuP Was: z(ilzhfs,lr_)  (should be same)
         endif

         step2=(z3-z2)/nsteps
         pt=z2-step2*.5
         do il=1,nsteps
            pt=pt+step2
            tem1(il)=psif(pt) ! B/Bmin
            !if(l.eq.1)write(*,*)il,pt,tem1(il)
         enddo

!c            if(lr_.eq.lrz) then
!c            write(*,*)'baviorbt  B/Bmin',l,lr_,bpsi(l,lr_)
!c            endif
         !write(*,*)'bav_up: l,step2',l,step2

!...................................................................
!     Begin loop over theta
!     For imax: see subroutine micxinil
!...................................................................

!..................................................................
!     Uses library routines:
!     luf(x,table,n) (MATHLIB) which is a function returning the index
!     of the first element in the table that is greater than x.
!     Table elements must be strictly increasing.  x.gt.table(n)==>n+1.
!     lug(x,table,n,iguess) (MATHLIB) same as luf, but with guess index iguess.
!..................................................................


         do 60 i=1,iyh ! index in pitch-angle grid; scan 0 to pi/2

            wp=0.

!     Maybe need to refine theta0 only for trapped or nearly
!     trapped particles. For now, do all theta0:

            yyp=y(i,l_)+0.5*dyp5(i,l_)
            yym=y(i,l_)-0.5*dym5(i,l_)
            uup=cos(yyp)
            uum=cos(yym)
            duu=(uup-uum)/nii
            uu=uum-0.5*duu

            call bcast(tem2,zero,nii*nsteps)
            do ii=1,nii
               uu=uu+duu
               !if(nii.gt.1) sinyy2=1.-uu*uu
               !if(nii.eq.1) sinyy2=sinn(i,l_)**2 !-YuP Much more accurate
               sinyy2=1.-uu*uu
               do il=1,nsteps
                  xs=tem1(il) ! B/Bmin
                  tem2(il+(ii-1)*nsteps)=-1.+sinyy2*xs  ! should be <= 0
               enddo
            enddo

!           Determine last il-index of passing particles for each ii
!           (If there are several zeros in a row (a roundoff effect)
!            step back to the first nonzero tem2).
            do ii=1,nii
              iguess=1  ! YuP: ??? what for? doesn't work in lug()
              ii0= 1+(ii-1)*nsteps
              ii2= ii0+nsteps-1 !YuP[2019-04-22]should be nsteps elements
              wk_lug(1:nsteps)= tem2(ii0:ii2)
              itemc1(ii)=lug(zero,wk_lug,nsteps,iguess)-1
               !write(*,*)'b1 after LUG itemc1=',lr_,l,itemc1(ii)
               if (itemc1(ii).ge.1) then
                  if (tem2(itemc1(ii)+(ii-1)*nsteps).eq.zero) then
 59                  itemc1(ii)=itemc1(ii)-1
                     if (tem2(itemc1(ii)+(ii-1)*nsteps).eq.zero &
                          .and. itemc1(ii).gt.0)  goto 59
                  endif
               endif
            enddo

            uu=uum-0.5*duu
            do ii=1,nii
               uu=uu+duu
               do il=1,itemc1(ii)
                  costh=sqrt(-tem2(il+(ii-1)*nsteps))
                  !if(nii.gt.1) dtau(i,l,lr_)=dtau(i,l,lr_)+uu/costh
                  !if(nii.eq.1) dtau(i,l,lr_)=dtau(i,l,lr_)+step2/costh
                   dtau(i,l,lr_)=dtau(i,l,lr_)+uu/costh
               enddo
            enddo

            !if(nii.gt.1) then
               duum=-duu
               dtau(i,l,lr_)=dtau(i,l,lr_)*step2*duum*twopi/ &
                 (coss(i,l_)*cynt2(i,l_))
            !endif

            !-YuP Added:
            if(l.eq.1) dtau(i,l,lr_)=dtau(i,l,lr_)*2.0
            !-YuP Should be *2 because first node is over half-interval?

            tau(i,lr_)=tau(i,lr_)+dtau(i,l,lr_)

 60      continue ! i=1,iyh
 70   continue ! l=1,ilzhfs


!...................................................................
!     Determine imax and lmax by searching for nonzero dtau
!     (defined in micxinil.f).
!...................................................................

      do l=1,ilzhfs
!cc         imax_old(l)=imax(l,lr_)  ! Can check with debugger
         ii=0
         do i=1,iyh
            if (dtau(i,l,lr_).le.0.) goto 80
            ii=i
         enddo
 80      imax(l,lr_)=ii
!BH091031         write(*,*)'baviorbt:l,imax,imax_old',l,imax(l,lr_),imax_old(l)
      enddo


      do i=1,iyh
!cc         lmax_old(i)=lmax(i,lr_)  ! Can check with debugger
         ll=0
         do l=1,ilzhfs
            if (dtau(i,l,lr_).le.0.) goto 90
            ll=l
         enddo
 90      lmax(i,lr_)=ll
!BH091031         write(*,*)'baviorbt:i,lmax,lmax_old',i,lmax(i,lr_),lmax_old(i)
      enddo



!.......................................................................
!     Do lower half of equilibrium, if eqsym.eq."none"
!.......................................................................

      if (eqsym.eq."none") then

      do l=lz,ilzhfs,-1 ! start at Bmin, scan towards Bmax
!     Check, need special treatment at ilzhfs, to avoid conflict
!     with "upper" integration.

!......................................................................
!     Evaluate the bbpsi function where it will normally be required for
!     integration. Doing this outside the "60" loop saves time.
!     BH,090901:  This does not square very well with the original
!                 baviorbto below, where intervals considered are
!                 [z_l,z_l+1],l=1,lmax where z_1=0, z_lmax+1=zboun.
!                 In the bounce-average of a function g, g is to b
!                 tabulated on z_l points.
!......................................................................


!        Order z2,z3 so (z3-z2) is positive, as for l=1,ilzhrs intgratn.
!        Half-steps at ends of z-integration.
            nsteps=nstps !-YuP
         if (l.eq.lz) then
            z2=zmid(lz-1,lr_) !-YuP: this point is between z(lz-1) and zmax
            z3=zmax(lr_)  !-YuP Was: zmid(lz,lr_) which is defined as zmax(lr_)
            !nsteps=nstps/2 !-YuP: because of half-interval
         else
            z2=zmid(l-1,lr_)
            z3=zmid(l,lr_)
         endif
         if (l.eq.ilzhfs) then
            z2=z_bmax(lr_)   !-YuP Was: z(ilzhfs,lr_)  (should be same)
            z3=zmid(ilzhfs,lr_)
         endif

         step2=(z3-z2)/nsteps
         !-YuP was: pt=z2-step2*.5
         pt=z3+step2*0.5 !-YuP
         do il=1,nsteps
            !-YuP was: pt=pt+step2
            pt=pt-step2 !-YuP
            !-YuP: decreasing pt corresponds to a scan from Bmin to Bmax
            tem1(il)=psif(pt)
            !-YuP: tem1 is increasing now, as required by lug() function
         enddo
!...................................................................
!     Begin loop over theta
!     Tempoarily store tau and dtau data at i+iyh
!     For imax: see subroutine micxinil
!...................................................................

!..................................................................
!     Uses library routines:
!     luf(x,table,n) (MATHLIB) which is a function returning the index
!     of the first element in the table that is greater than x.
!     Table elements must be strictly increasing.  x.gt.table(n)==>n+1.
!     lug(x,table,n,iguess) (MATHLIB) same as luf, but with guess index iguess.
!..................................................................


         do i=1,iyh

            wp=0.

!     Maybe need to refine theta0 only for trapped or nearly
!     trapped particles. For now, do all theta0:

            yyp=y(i,l_)+0.5*dyp5(i,l_)
            yym=y(i,l_)-0.5*dym5(i,l_)
            uup=cos(yyp)
            uum=cos(yym)
            duu=(uup-uum)/nii
            uu=uum-0.5*duu

            call bcast(tem2,zero,nii*nsteps)
            do ii=1,nii
               uu=uu+duu
               !if(nii.gt.1) sinyy2=1.-uu*uu
               !if(nii.eq.1) sinyy2=sinn(i,l_)**2 !-YuP Much more accurate
               sinyy2=1.-uu*uu
               do il=1,nsteps
                  xs=tem1(il)
                  tem2(il+(ii-1)*nsteps)=-1.+sinyy2*xs
               enddo
            enddo

!           Determine last il-index of passing particles for each ii
!           (If there are several zeros in a row (a roundoff effect
!            step back to the first nonzero tem2).
            do ii=1,nii
              iguess=1
              ii2= ii0+nsteps-1 !YuP[2019-04-22]should be nsteps elements
              wk_lug(1:nsteps)= tem2(ii0:ii2)
              itemc1(ii)=lug(zero,wk_lug,nsteps,iguess)-1
               !write(*,*)'b2 after LUG itemc1=',lr_,l,itemc1(ii)
               if (itemc1(ii).ge.1) then
                  if (tem2(itemc1(ii)+(ii-1)*nsteps).eq.zero) then
 58                  itemc1(ii)=itemc1(ii)-1
                     if (tem2(itemc1(ii)+(ii-1)*nsteps).eq.zero &
                          .and. itemc1(ii).gt.0)  goto 58
                  endif
               endif
            enddo


            uu=uum-0.5*duu
            do ii=1,nii
               uu=uu+duu
               do il=1,itemc1(ii)
                  costh=sqrt(-tem2(il+(ii-1)*nsteps))
            !if(nii.gt.1) dtau(i+iyh,l,lr_)=dtau(i+iyh,l,lr_)+uu/costh
            !if(nii.eq.1) dtau(i+iyh,l,lr_)=dtau(i+iyh,l,lr_)+step2/costh
                  dtau(i+iyh,l,lr_)=dtau(i+iyh,l,lr_)+uu/costh
               enddo
            enddo

!           Put tau,dtau  results in second quadrant of theta0
!           NOTE:  for l=ilzhfs, following dtau(,,) will add in half
!                  z-step contribution from upper half flux surface.

            !if(nii.gt.1) then
              duum=-duu
              dtau(i+iyh,l,lr_)=dtau(i+iyh,l,lr_)*step2*duum*twopi/ &
                 (coss(i,l_)*cynt2(i,l_))
            !endif

            !-YuP Added:
            if(l.eq.lz) dtau(i+iyh,l,lr_)=dtau(i+iyh,l,lr_)*2.0
            !-YuP Should be *2 because last node is over half-interval?

            tau(i+iyh,lr_)=tau(i+iyh,lr_)+dtau(i+iyh,l,lr_)

         enddo  ! on i=1,iyh
      enddo  ! on l=lz,ilzhfs


!     Put dtau(i+iyh,l,lr_) into dtau(i,l,lr_), adding upper and
!     lower contributions at l=ilzhfs
      do i=1,iyh
         do l=lz,ilzhfs+1,-1
            dtau(i,l,lr_)=dtau(i+iyh,l,lr_)
         enddo
         l=ilzhfs
         dtau(i,l,lr_)=dtau(i,l,lr_)+dtau(i+iyh,l,lr_)
      enddo


!...................................................................
!     Determine imax and lmax by searching for nonzero dtau
!...................................................................

      do l=lz,ilzhfs,-1
!cc         imax_old(l)=imax(l,lr_)  ! Can check with debugger
         ii=0
         do i=1,iyh
            if (dtau(i,l,lr_).le.0.) goto 81
            ii=i
         enddo
 81      imax(l,lr_)=ii
      enddo


      do i=1,iyh
!cc         lmax_old(i)=lmax(i,lr_)  ! Can check with debugger
         ll=0
         do l=lz,ilzhfs,-1
            if (dtau(i,l,lr_).le.0.) goto 91
            ll=l
         enddo
 91      lmax(i+iyh,lr_)=ll
!BH091031         write(*,*)'baviorbt:i,lmax,lmax_old',i,lmax(i,lr_),lmax_old(i)
      enddo


      endif  ! On eqsym.eq."none"



!......................................................................
!     Reassign the bounce point for particles at pass/trapped
!     boundary to zmax(lr_)+abit
!......................................................................

      zboun(itl,lr_)=z_bmax(lr_)+1.e-10
      if (eqsym.ne."none") then
         zboun(itu,lr_)=zboun(itl,lr_)
      else
         zboun(itu,lr_)=z_bmax(lr_)-1.e-10
      endif


!......................................................................
!     Symmetrize the bounce times
!     For eqsym.eq."none" case, first add up the tau contributions
!     from the "upper" and "lower" integration.
!     Pitch angle symmetries are unchanged from eqsym.ne."none" case.
!......................................................................

      if (eqsym.eq."none") then
         do i=1,iyh
            tau(i,lr_)=tau(i,lr_)+tau(i+iyh,lr_)
         enddo
      endif


!DIR$ IVDEP

      do 120 i=1,iyh
         tau(iy+1-i,lr_)=tau(i,lr_)
 120  continue

      do 121 i=1,iy
         vptb(i,lr_)=abs(coss(i,l_))*tau(i,lr_)
 121  continue

      do 131 l=1,lz             !  For all l
         do i=1,iyh
            dtau(iy+1-i,l,lr_)=dtau(i,l,lr_)
         enddo
!cc         do i=1,iy
!cc            prnt5(i)=dtau(i,l,lr_)
!cc         enddo
!cc         lll=l  !added as stop point for gdb   ???? not used here?
 131  continue  ! l=1,lz


!cc      do i=1,iy
!cc         prnt1(i)=tau(i,1)
!cc         prnt2(i)=tau(i,2)
!cc         prnt3(i)=tau(i,3)
!cc         prnt4(i)=tau(i,4)
!cc      enddo
!BH091031      do i=1,iy
!BH091031         do ll=1,lr_
!BH091031            prnt5(i,ll)=tau(i,ll)
!BH091031         enddo
!BH091031      enddo


      return
      end subroutine baviorbt
!
!
!=======================================================================
!=======================================================================
!
!
      subroutine baviorbto
      use param_mod
      use comm_mod
      use r8subs_mod, only : luf
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!...................................................................
!     This routine computes the quantity equal to the product of the
!     bounce time with the particle velocity. This depends only on the
!     midplane pitch angle and is stored in tau(i,lr_). However at the
!     pass/trapped boundary only the partial integral, integrated
!     from 0 to zstar is computed. The remainder of the calculation
!     is done in subroutine micgnbnd. The local time weights associated
!     with z(l,lr_) and y(i,l_) are also computed here and are stored in
!     dtau(i,l,lr_). We have assumed in the calculation that the
!     function about which the bounce average is to be computed is
!     linear in bbpsi between z(l,lr_) and z(l+1,lr_). The bounce average
!     is the sum from 1 to lz of g(i,l)*dtau(i,l,lr_)/tau(i,lr_).
!     Caution: for trapped particles we must have for some l
!     z(l,lr_) .lt. zboun(i,lr_) .lt. z(l+1,lr_). There will be a nonzero
!     value in dtau(i,l+1,lr_) in general. Thus a reasonable value
!     must be stored in g(i,l+1). In this code we use the value
!     of g at the bounce point if possible. Otherwise we set
!     g(i,l+1)=g(i,l).
!     Refer to p. 96-7, Killeen, Kerbel, Mccoy, Mirin book.
!...................................................................

!cc      real(c_double), dimension(iy):: prnt1,prnt2,prnt3,prnt4  !For gdb printing
!cc      real(c_double), dimension(iy,lrzmax):: prnt5   !For gdb printing


      character*8 bnc

      if(eqsym.eq."none") STOP 'baviorbto not revised for eqsym.eq.none'

!...................................................................
!     ERROR CONTROL - change nstps if bigger than iyjx
!...................................................................

      if (nstps .gt. iyjx) nstps=iyjx

!...................................................................
!     Limit the range of integration at the pass/trapped boundary
!...................................................................

      zboun(itl,lr_)=zstar
      call bcast(tau(1:iyh,lr_),zero,iyh)
      !tau = 0
      call bcast(dtau(1:iy*lz,1,lr_),zero,iy*lz)
      !dtau = 0

!...................................................................
!     Begin loop over the particle orbit
!...................................................................

      ilzhfs=lz
      if (numclas .eq. 1) ilzhfs=lz/2+1
      do 70 l=1,ilzhfs-1

!...................................................................
!     Evaluate the bbpsi function where it will normally be required for
!     integration. Doing this outside the "60" loop saves time.
!...................................................................

        z2=z(l,lr_)
        x2=bbpsi(l,lr_)
        z3=z(l+1,lr_)
        step2=(z3-z2)/nstps
        pt=z2+step2*.5
        do 20 il=1,nstps
          tem1(il)=psif(pt)
          pt=pt+step2
 20     continue

!...................................................................
!     Begin loop over theta
!...................................................................

        do 60 i=1,imax(l,lr_)
!...................................................................
!     The integration for bounce time associated with z(l,lr_) involves
!     two integral contributions: the first over the interval
!     [z(l-1,lr_),z(l,lr_] is stored in wm; the second over
!     [z(l,lr_),zb] is stored in wp.
!     Here zb= minimum(z(l,lr_),zboun(i,lr_)).
!...................................................................

          wm=0.
          wp=0.

!.................................................................
!     Proceed with calculation of wm
!     wm contribution is void for l=1
!.................................................................

          if(l.ne.1) wm=temc1(i)
          temc1(i)=0.0

!...................................................................
!     Reassign z3 if the particle bounces in the current region
!...................................................................

          if (zboun(i,lr_) .lt. z(l+1,lr_)) then
            if (zboun(i,lr_) .lt. z(l,lr_)) goto 60
            bnc="enabled"
            z3=zboun(i,lr_)
            x3=psif(z3)
          else
            z3=z(l+1,lr_)
            x3=bbpsi(l+1,lr_)
            bnc="disabled"
          endif
!...................................................................
!     Proceed with the calculation of wp
!...................................................................

          if(bnc.eq."enabled") then
            step2sww=(z3-z2)/nstps
            pt=z2+step2sww*.5
            do 30 il=1,nstps
              xs=psif(pt)
              wp=wp+step2sww*(1.-(xs-x2)/(x3-x2)) &
                /sqrt(1.-sinn(i,l_)**2*xs)
              pt=pt+step2sww
 30         continue
          else
            step2sww=step2
            do 31 il=1,nstps
              xs=tem1(il)
              swwtemp2=1.0/sqrt(1.-sinn(i,l_)**2*xs)
              swwtemp3=step2*(xs-x2)/(x3-x2)*swwtemp2
              temc1(i)=temc1(i)+swwtemp3
              wp=wp+step2*(1.-(xs-x2)/(x3-x2))*swwtemp2
 31         continue
          endif
          dtau(i,l,lr_)=wp+wm
          tau(i,lr_)=tau(i,lr_)+dtau(i,l,lr_)

!...................................................................
!     Include the last contribution for passing particles and for
!     particles such that z(l+1,lr_) .gt. zboun(i,lr_) .gt. z(l,lr_)
!...................................................................

          if (l .eq. ilzhfs-1 .and. lmax(i,lr_) .eq. ilzhfs) goto 40
          if (lmax(i,lr_) .eq. l) goto 40
          goto 60
 40       continue
          pt=z2+.5*step2sww
          wm=0.
          do 50 il=1,nstps
            xs=psif(pt)
            wm=wm+step2sww*(xs-x2)/(x3-x2)/sqrt(1.-sinn(i,l_)**2*xs)
            pt=pt+step2sww
 50       continue
          dtau(i,l+1,lr_)=wm
          tau(i,lr_)=tau(i,lr_)+dtau(i,l+1,lr_)
 60     continue
 70   continue

!...................................................................
!     Reassign the bounce point for for particles at pass/trapped
!     boundary to .ge. zmax(lr_)
!...................................................................

      zboun(itl,lr_)=zmax(lr_)+1.e-10
      zboun(itu,lr_)=zboun(itl,lr_)

!...................................................................
!     Symmetrize the bounce times
!...................................................................

!DIR$ IVDEP
      do 120 i=1,iyh
        tau(iy+1-i,lr_)=tau(i,lr_)
 120  continue

      do 131 l=1,ilzhfs
        do 130 i=1,iyh
          dtau(iy+1-i,l,lr_)=dtau(i,l,lr_)
 130    continue
 131  continue

!..................................................................
!     Determine bounce time in trapped/passing sliver region.
!..................................................................

      call micgnbnd

!cc      do i=1,iy
!cc         prnt1(i)=tau(i,1)
!cc         prnt2(i)=tau(i,2)
!cc         prnt3(i)=tau(i,3)
!cc         prnt4(i)=tau(i,4)
!cc      enddo
!cc      do i=1,iy
!cc         do ll=1,lr_
!cc            prnt5(i,ll)=tau(i,ll)
!cc         enddo
!cc      enddo

      return
      end subroutine baviorbto


!
!
      subroutine deltar
      use param_mod
      use comm_mod
      use r8subs_mod, only : luf
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save


!     Calculate the first order orbit-width shift, deltarho, of
!     guiding center orbits from their associated zero-orbit-width,
!     bounce-averaged flux surface, divided by particle speed.
!     That is, for a particle at given radius rya(lr) and distance
!     along the magnetic surface, z(l,lr) [i.e., at given cyl coords
!     R,Z], and given midplane velocity coordinates, v0,theta0,
!     obtain the shift FROM the bounce-averaged flux surface for
!     the particle (thus particles in the co-current electric drift
!     direction will be shifted outwards from their BA flux surface).
!     Velocity is divided out of the orbit shift, reducing the shift,
!     deltarho(theta0,z,rya), to depend only on three variables.
!
!     deltarhop(theta,z,rya), where theta is local pitch, is derived
!       from deltarho.
!     deltarz(R,Z,theta), orbit shift as function of R,Z and local pitch
!       angle is then interpolated from deltarhop.
!       Adjustments are made to deltarz to ensure that it does not
!       point to a flux surface with negative rho.
!
!     This subroutine is called after above calc of tau/dtau (taunew=
!     "enabled"), and call to tdmshst.
!     Is is assumed (has been checked in ainsetva.f) that lrzdiff
!     .eq."disabled", giving lrz=lrzmax.
!
!     Trapped and transiting particles are treated separately.
!     See BH notes: 090821.


      write(*,*)'deltar: Beginning of deltarho orbit-shift calc'
      k=1  !ONLY setup for one species, ngen=1

!......................................................................
      do lr=1,lrzmax
!......................................................................

      iyy=iy_(lr) !-YuP-101215: Don't use iy=; it's in common /params/
      itl=itl_(lr)
      iyh=iyh_(lr)
      call bcast(deltarho(1:iy*lrzmax,1,lr),zero,iy*lrzmax) ! YuP: is iy*lrzmax correct?.. umm you tell me XXX


!.......................................................................
!     First evaluate absolute values of partial(epsi)/partial(Z)
!     (=tz1()), and 1./(R*omegac) (=tz2() along z().
!.......................................................................


      do l=1,lz  !lz is independent of lr
         tz1(l)=terp2(solrz(l,lr),solzz(l,lr),nnr,er,nnz,ez,epsi, &
              epsirr,epsizz,epsirz,nnra,0,1)
         tz1(l)=abs(tz1(l))
         tz2(l)=1./(solrz(l,lr)*bnumb(k)*charge*bbpsi(l,lr)*bmod0(lr) &
              /(fmass(k)*clight))
         tz2(l)=abs(tz2(l))
      enddo

!      write(*,*)'deltar:dPsi/dZ (tz1) =',(tz1(l),l=1,lz)
!      write(*,*)'deltar:1/R*omegac (tz2) =',(tz2(l),l=1,lz)


!     Trapped particles:  orbit shift given by vertical drift off
!     the flux surface from initial point to banana tip.
!


!     For transiting particles, setup partial integrals from
!     z=0 to z=z(l) of dtau:
!     Use temp1 storage, so check lz.le.jx+1
      if (lz.gt.jx+1) STOP 'In baviorbt: Check lz versus jx+1'
      call bcast(temp1(0:iyjx2-1,0),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)
      !temp1 = 0
      do i=1,itl
         temp1(i,1)=dtau(i,1,lr)
         do l=2,lz
            temp1(i,l)=temp1(i,l-1)+dtau(i,l,lr)
         enddo
      enddo

      call bcast(temp2(0:iyjx2-1,0),zero,iyjx2)
      call bcast(temp3(0:iyjx2-1,0),zero,iyjx2)
      !temp3 = 0
      do i=1,iyh

         if (i.le.itl) then !transiting
!           First compute partial integrals z(1) to z(lz) in temp2
!           and temp3, then form orbit shift.
!           Note:  temp2(i,l=0)=temp3(i,l=0)=0.
            do l=1,lz
              temp2(i,l)=temp2(i,l-1)+dtau(i,l,lr)*tz1(l)*tz2(l) &
                   *(1.d0-0.5*sinz(i,l,lr)**2)*(temp1(i,lz)-temp1(i,l))
              temp3(i,l)=temp3(i,l-1)+dtau(i,l,lr)*tz1(l)*tz2(l) &
                   *(1.d0-0.5*sinz(i,l,lr)**2)*temp1(i,l)
            enddo
            do l=1,lz
               deltarho(i,l,lr)=(1/tau(i,lr))* &
                    (temp2(i,lz)-temp2(i,l-1)-temp3(i,l))
            enddo


         else                   !trapped

!           lmax(i,lr) will be less than lz for trapped particles
!BH160605: Try this change.  Not justified yet.
!BH160605:            if (lmax(i,lr).ge.lz) STOP 'deltarho calc problem'
            if (lmax(i,lr).gt.lz) STOP 'deltarho calc problem'
            do l=lmax(i,lr),1,-1
               deltarho(i,l,lr)= deltarho(i,l+1,lr)+dtau(i,l,lr) &
                    *tz1(l)*tz2(l)*(1.d0-0.5*sinz(i,l,lr)**2)
            enddo

         endif

      enddo  ! on i
!$$$      write(*,*)'baviorbt: Part dtau intrl:temp1=',(temp1(itl,l),l=1,lz)
!$$$      write(*,*)'baviorbt: Part dtau intrl:temp2=',(temp2(itl,l),l=1,lz)
!$$$      write(*,*)'baviorbt: Part dtau intrl:temp3=',(temp3(itl,l),l=1,lz)
!$$$      write(*,*)'baviorbt: deltarho(transiting)=',
!$$$     +     (deltarho(itl,l,lr),l=1,lz)
!$$$      write(*,*)'baviorbt: deltarho(trapped)=',
!$$$     +     (deltarho(itl+1,l,lr),l=1,lz)

!     Changing deltarho from shift in poloidal flux psi to shift
!     to rho (as specified by radcoord). [eqrho is in cms, for
!     radcoord=sqtorflx, eqpsi in cgs.]
!     Splines setup with eqpsi sign temporarily reversed,
!     to give monotonically increasing function.
!      call coeff1(nconteqn,eqpsi,eqrho,d2eqrho,i1p,1,workk)


      itab(1)=1 ! just to check rho value in tab(1) with debugger.
      itab(2)=1
      itab(3)=0
      do j=1,nconteqn
         eqpsi(j)=-eqpsi(j)
      enddo
      call terp1(nconteqn,eqpsi,eqrho,d2eqrho,-epsicon(lr),1,tab,itab)
!      write(*,*)'deltar: eqpsi(1:nconteqn)=',(eqpsi(l),l=1,nconteqn)
!      write(*,*)'deltar: eqrho(1:nconteqn)=',(eqrho(l),l=1,nconteqn)
      do j=1,nconteqn
         eqpsi(j)=-eqpsi(j)
      enddo
      drhodpsi=tab(2)
!$$$      write(*,*)'deltar: HERE1'
      do l=1,lz
         do i=1,iyh
            deltarho(i,l,lr)=abs(drhodpsi)*deltarho(i,l,lr)/rhomax
         enddo
      enddo
!$$$      write(*,*)'deltar: HERE2'

!     Adjusting sign of deltarho for co-/counter-initial velocity,
!     and for electrons or ions.  sgn(vpar)=+1 for i=1,iyh

      do l=1,lz
         do i=1,iyh
            isgn_deltarho=-cursign*sign(one,bnumb(k))
            deltarho(i,l,lr)=deltarho(i,l,lr)*isgn_deltarho
            deltarho(iyy-(i-1),l,lr)=-deltarho(i,l,lr)
         enddo
      enddo



!$$$      write(*,*)"baviorbt:deltarho(tp_bndry) Norm=rhomax*x*vnorm",
!$$$     +     (deltarho(itl,l,lr),l=1,lz)
!$$$      write(*,*)"baviorbt:deltarho(tp_bndry+1) Norm=rhomax*x*vnorm",
!$$$     +     (deltarho(itl+1,l,lr),l=1,lz)
!$$$      write(*,*)"baviorbt:deltarho(tp_bndry-1) Norm=rhomax*x*vnorm",
!$$$     +     (deltarho(itl-1,l,lr),l=1,lz)
!$$$      write(*,*)"baviorbt:deltarho(i=1-passing) Norm=rhomax*x*vnorm",
!$$$     +     (deltarho(1,l,lr),l=1,lz)
!$$$      write(*,*)"baviorbt:deltarho(i=iy-passing) Norm=rhomax*x*vnorm",
!$$$     +     (deltarho(1,l,lr),l=1,lz)


!......................................................................
      enddo  ! on lr
!......................................................................



!......................................................................

!  Interpolate above deltarho data to orbit shift deltrz(ir=1,nr_delta,
!  iz=1,nz_delta,it=1,nt_delta) on a regular R,Z,theta grid. theta is
!  local pitch angle.  nr_delta,nz_delta,nt_delta are namelist input
!  variables.
!  Only do this after filling in all above lr=1,lrzmax values.
!  (Subroutine calls are with lr in reverse order from lrzmax to 1.)
!......................................................................

!     First, interpolate deltarho data to deltarhop(i,l,lr) where
!     i is a equispaced local pitch angle grid: sin**2(theta)=
!     (B/B0)*sin**2(theta0), and theta0 is the corresponding
!     equatorial plane value. (This is 1st order in small banana
!     width calculation, so don't calculate B/B0 change to shifted
!     flux surface.)
!     Then, can interpolate deltarhop to R,Z grid at same pitch
!     angles i.

      dt_delta=pi/(nt_delta-1)
      do i=1,nt_delta
         t_delta(i)=(i-1)*dt_delta
      enddo
!      write(*,*)'deltar:t_delta(1:nt_delta)=',t_delta(1:nt_delta)

      do lr=1,lrzmax   !Coding setup for lrzdiff.ne."enabled"
         iyy=iy_(lr) !YuP-101215: Don't use iy=; it's in common /params/
!         deltarhop(1,1,lr)=deltarho(1,1,lr)
!         deltarhop(nt_delta,1,lr)=deltarho(iyy,1,lr)
         do l=1,lz
            do i=1,nt_delta/2
               sin_theta02=sin(t_delta(i))**2/bbpsi(l,lr)
               theta0=asin(sqrt(sin_theta02))
               call lookup(theta0,y(1:iyy,lr),iyy,wtu,wtl,ii)
               deltarhop(i,l,lr)= &
                    wtl*deltarho(ii-1,l,lr)+wtu*deltarho(ii,l,lr)
               deltarhop(nt_delta+1-i,l,lr)=-deltarhop(i,l,lr)
            enddo
         enddo
!      write(*,*)'deltar: lr,lz,iyy=',lr,lz,iyy
!      write(*,*)'deltar: deltarhop(1:nt_delta,1,lr)',
!     +                   deltarhop(1:nt_delta,1,lr)
!      write(*,*)'deltar: deltarho(1:iyy,1,lr)',deltarho(1:iyy,1,lr)
!      write(*,*)
!      write(*,*)'deltar: deltarhop(1:nt_delta,10,lr)',
!     +                   deltarhop(1:nt_delta,10,lr)
!      write(*,*)'deltar: deltarho(1:iyy,10,lr)',deltarho(1:iyy,10,lr)
!      write(*,*)
!      write(*,*)'deltar: deltarhop(1:nt_delta,lz-1,lr)',
!     +                   deltarhop(1:nt_delta,lz-1,lr)
!      write(*,*)'deltar: deltarho(1:iyy,lz-1,)',deltarho(1:iyy,lz-1,lr)
!      write(*,*)
!      write(*,*)'deltar: deltarhop(1:nt_delta,lz,lr)',
!     +                   deltarhop(1:nt_delta,lz,lr)
!      write(*,*)'deltar: deltarho(1:iyy,lz,lr)',deltarho(1:iyy,lz,lr)
!      write(*,*)

      enddo  !On lr


!......................................................................

!     Need to apply following to deltarho and deltarz, after
!     multiplication by velocity.


!     Modify negative values of deltarhop so that inward radial shift
!     is sensible: not greater than the distance to the magnetic axis.
!     [Choose that inwards shift is not further than to center
!     of the 1st radial bin.]
!
!$$$      do lr=1,lrzmax            !Coding setup for lrzdiff.ne."enabled"
!$$$         do l=1,lz
!$$$            do i=1,nt_delta
!$$$               deltarhop(i,l,lr)=min(rya(1)-rya(lr),deltarhop(i,l,lr))
!$$$            enddo
!$$$         enddo
!$$$      enddo  !On lr
!......................................................................
!  OR
!......................................................................
!
!     Simplest to modify orbit shift after interpolation for
!     orbit shift:
!
!     deltarz_interp=min(rya(1)-rho,deltarz_int*v)
!     where
!     rya(1) in innermost cql3d flux surface (as the .nc file)
!     rho is radial flux surface coord at R,Z location in question
!     deltarz_int is interpolate to given R,Z,theta, of deltarz array.
!     v is particle velocity.
!
!......................................................................

!     Setup R,Z grid.  Grid covers approx 10 percent outside LCFS.
      dr_delta=(1.1*rgeom2-0.9*rgeom1)/(nr_delta-1)
      dz_delta=(1.1*zgeomp-(-1.1*zgeomp))/(nz_delta-1) !symmetry assumed
      do i=1,nr_delta
         r_delta(i)=0.9*rgeom1+(i-1)*dr_delta
      enddo
      do i=1,nz_delta
         z_delta(i)=-1.1*zgeomp+(i-1)*dz_delta
      enddo

!     Simple loop over R,Z points, determining nearest l,lr points in
!     the z(l=1,lz,lr=1,lrzmax) grid.
!     NB:  up-down symmetry assumed
!     If outside LCFS, extrapolate rho shift at constant poloidal
!     angle from outermost computational flux surface (rya(lrzmax)).

!     tr2() below contains bin boundaries, in ascending order (for luf)
      do lr=1,lrzmax
         tr2(lr)=psimag-psivalm(lr)
      enddo
!      write(*,*)'deltar: tr2 (bin bndries)',tr2(1:lrzmax)

!      do lr=1,lrzmax
!         write(*,*)'deltar: pol(1:lz,lr)=',pol(1:lz,lr)
!      enddo

      do ir=1,nr_delta
         do iz=1,nz_delta
            rpos=r_delta(ir)
            zpos=z_delta(iz)

!           pol flux value and radial bin lr
            psival=terp2(rpos,zpos,nnr,er, &
              nnz,ez,epsi,epsirr,epsizz,epsirz,nnra,0,0)
            apsi=psimag-psival
            lr=min(luf(apsi,tr2(1),lrzmax),lrzmax)

!           Determine poloidal angle about mag axis:
!           (see eqorbit.f and micxiniz.f, only setup for up-down symmetry)
!$$$            tang=(zmag-zpos)/(rmag-rpos)
!$$$            if (tang.ge.0.) then
!$$$               thetpol=atan(tang)
!$$$            else
!$$$               thetpol=pi-atan(-tang)
!$$$            endif
            if (eqsym.ne."none") then  !up-down symmetrize
               thetpol=abs(atan2(zpos-zmag,rpos-rmag)) ! in [0,pi]
            else  !i.e., no up-down symmetrization
               thetpol=atan2(zpos-zmag,rpos-rmag)
               if (thetpol.lt.zero) thetpol=thetpol+twopi ! in [0,pi2)
            endif

!           Determine pol angle bin
            ll=min(luf(thetpol,pol(1:lz,lr),lz),lz-1)
!            write(*,*)'deltar:rpos,rmag,zpos,zmag,tang,thetpol=',
!     +                        rpos,rmag,zpos,zmag,tang,thetpol

!           Bilinear interpolate between rya(lr) and rya(lr+/-1) for
!           psival outside/inside ray(lr).
!           NB:  For poloidal position, the normalized z-values are
!                same (to roundoff) for each flux surface, keeping
!                the poloidal points close to each other for each ll.

!           ione quantity enables interpolation to nearest
!           lr,lr+ione (ione=+/- 1) flux surfaces with one coding:
            ione=-1
            if (apsi.gt.(psimag-equilpsi(lr))) ione=+1
!            write(*,*)'deltar:ir,iz,ll,lr,apsi,psimag,psim-equilpsi(lr)'
!     +                 ,ir,iz,ll,lr,apsi,psimag,psimag-equilpsi(lr)
!            write(*,*)'deltar: ione',ione

            if ((lr.ne.1 .and. ione.lt.0) &
                 .or. (lr.ne.lrzmax .and. ione.gt.0)) then
               area1=abs((rpos-solrz(ll,lr))*(zpos-solzz(ll,lr)))
               area2=abs((solrz(ll,lr+ione)-rpos)* &
                    (zpos-solzz(ll,lr+ione)))
               area3=abs((solrz(ll+1,lr+ione)-rpos)* &
                    (solzz(ll+1,lr+ione)-zpos))
               area4=abs((rpos-solrz(ll+1,lr))*(solzz(ll+1,lr)-zpos))
               areatot=area1+area2+area3+area4

!              Use same interpolation for delta_bdb0(R,Z)=B/B0
               delta_bdb0(ir,iz)=(area3*bbpsi(ll,lr) &
                    +area4*bbpsi(ll,lr+ione) &
                    +area1*bbpsi(ll+1,lr+ione) &
                    +area2*bbpsi(ll+1,lr))/areatot

               do i=1,nt_delta
                  deltarz(ir,iz,i)=(area3*deltarhop(i,ll,lr) &
                       +area4*deltarhop(i,ll,lr+ione) &
                       +area1*deltarhop(i,ll+1,lr+ione) &
                       +area2*deltarhop(i,ll+1,lr))/areatot
               enddo
!               write(*,*)'deltar: area1,area2,area3,area4,areatot',
!     +                            area1,area2,area3,area4,areatot
!               write(*,*)'deltar: deltarz(ir,iz,1:4),deltarz(ir,iz,10)',
!     +                           deltarz(ir,iz,1:4),deltarz(ir,iz,10)

            elseif (lr.eq.1 .and. ione.lt.0 ) then
!               write(*,*)'deltar: lr,lrzmax,ione =',lr,lrzmax,ione

!              Interpolate deltarhop at same poloidal angle,
!              at lr=1 flux surface
               call lookup(thetpol,pol(1:lz,lr),lz,weightu,weightl,lll)
               delta_bdb0(ir,iz)=weightl*bbpsi(lll-1,lr) &
                                  +weightu*bbpsi(lll,lr)
               do i=1,nt_delta
                  deltarz(ir,iz,i)=weightl*deltarhop(i,lll-1,lr) &
                                  +weightu*deltarhop(i,lll,lr)
               enddo
!               write(*,*)'deltar: deltarz(ir,iz,1:4),deltarz(ir,iz,10)',
!     +                           deltarz(ir,iz,1:4),deltarz(ir,iz,10)

            elseif (lr.eq.lrzmax .and. ione.gt.0 ) then
!               write(*,*)'deltar: lr,lrzmax,ione =',lr,lrzmax,ione


!              Interpolate deltarhop at same poloidal angle,
!              at lr=lrzmax flux surface
               call lookup(thetpol,pol(1:lz,lr),lz,weightu,weightl,lll)
!               write(*,*)'deltar: tang,weightl,weightu,lll =',
!     +                            tang,weightl,weightu,lll
               delta_bdb0(ir,iz)=weightl*bbpsi(lll-1,lr) &
                                  +weightu*bbpsi(lll,lr)
               do i=1,nt_delta
                  deltarz(ir,iz,i)=weightl*deltarhop(i,lll-1,lr) &
                                  +weightu*deltarhop(i,lll,lr)
               enddo
               write(*,*)'deltar: deltarz(ir,iz,1:4),deltarz(ir,iz,10)', &
                                 deltarz(ir,iz,1:4),deltarz(ir,iz,10)

            else
               write(*,*)"deltar:  Shouldn't happen!"
               STOP
            endif


         enddo ! On iz
      enddo    ! On ir

!      write(*,*)'delt1r: deltarz(1:nr_delta,nz_delta/2+1,1)=',
!     +                   deltarz(1:nr_delta,nz_delta/2+1,1)
!      write(*,*)'delt1r: deltarz(1:nr_delta,nz_delta/2+2,1)=',
!     +                   deltarz(1:nr_delta,nz_delta/2+2,1)


 999  return
      end subroutine deltar


      real(c_double) function deltarz_interp(rr,zz,thet,v)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save


!     Interpolate deltarz array to given r,z,thet point,
!     multiply by velocity v, and adjust calc if necessary
!     so deltarz_interp distance from rya(1) flux surface does
!     not exceed distance from rya(1), i.e., orbit shift does
!     not point to a flux surface smaller than rya(1).
!     r,z (cms), thet (radians) v (cm/sec).

!     Tri-linear interpolation.

      ir=(rr-r_delta(1))/dr_delta +1
!      ir=max(ir,1)     ! Assume on grid, as it goes 10% outside LCFS
!      ir=min(ir,nr_delta-1)
      iz=(zz-z_delta(1))/dz_delta +1
!      iz=max(iz,1)     ! Assume on grid, as it goes 10% outside LCFS
!      iz=min(iz,nz_delta-1)
      it=(thet-t_delta(1))/dt_delta +1
!      it=max(it,1)
      it=min(it,nt_delta-1)  ! to treat thet=pi case.

      rweightl=(rr-r_delta(ir))/dr_delta
      rweightr=1.d0-rweightl
      zweightl=(zz-z_delta(iz))/dz_delta
      zweightr=1.d0-zweightl
      tweightl=(thet-t_delta(it))/dt_delta
      tweightr=1.d0-tweightl

      deltai=tweightl*(rweightl*zweightl*deltarz(ir,iz,it) &
                      +rweightr*zweightl*deltarz(ir+1,iz,it) &
                      +rweightl*zweightr*deltarz(ir,iz+1,it) &
                      +rweightr*zweightr*deltarz(ir+1,iz+1,it)) &
            +tweightr*(rweightl*zweightl*deltarz(ir,iz,it+1) &
                      +rweightr*zweightl*deltarz(ir+1,iz,it+1) &
                      +rweightl*zweightr*deltarz(ir,iz+1,it+1) &
                      +rweightr*zweightr*deltarz(ir+1,iz+1,it+1))

!     Probably better to apply following step outside of the
!     interpolation of deltarz(,,) since rho not known here.
!      deltarz_interp=min(rya(1)-rho,deltai*v)

      deltarz_interp=deltai

      return
      end





end module baviorbt_mod
