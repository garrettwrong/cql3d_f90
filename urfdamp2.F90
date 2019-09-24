module urfdamp2_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast
  use tdnflxs_mod, only : tdnflxs
  use urfmidv_mod, only : urfmidv_db
  use urfmidv_mod, only : urfmidv_dc
  !XXXXXXXXXXXXXXX these require disabling TKR checks!
  !use urfpackm_mod, only : unpack
  !use urfpackm_mod, only : unpack16
  external unpack
  external unpack16

  !---END USE

!
!

contains

      subroutine urfdamp2(krf)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!yup      save

!..................................................................
!     This routine computes the damping for all
!     ray elements on all rays for UNIT power delpwr. This is put into
!     array urfpwr(is,iray,krf). The algorithm follows that of urfb0.
!
!     Power absorption for each ray element is obtained by integration
!     over the resonance region of velocity space of energy*(df/dt)_QL.
!     Local values of the distribution function are obtained from
!     midplane values.  The algorithm follows that of urfb0.
!.......................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

      complex*16 cwz,cwxyp,cwxym,cei

      integer*1 ipwr(nrayelts) ! YuP: For MPI:  indicator  (0 or 1)
                               ! that a ray element contributed any power

      dimension gftp(0:jx) !, fppj(0:iy+1,0:jx+1) ! local working arrays

      fppj(i,j)=g_(i,j+1,k,l_)  + &
               (g_(i,j,k,l_)-g_(i,j+1,k,l_))*dj(i,j,k,l_)

      data nray0 /1/

      cei=(0.,1.)

!..................................................................
!     fill out increasing psi array of bin boundaries, as in urffflx.
!..................................................................
         do l=1,setup0%lrzmax
            tr2(l)=psimag-psivalm(l)
         enddo

!     k is general species which this krf mode damps on:
      k=nrfspecies(krfn(krf))
      qmc= charge*bnumb(k)/(fmass(k)*clight) !store in common/orb0/
      qmca= abs(qmc)
      o65535= 1.d0/65535.d0
      vnorm4i= 1.d0/vnorm4
!..................................................................
!     Express the Chang-Cooper weighted average g_(i,j+1/2,l_)
!..................................................................
!      call bcast(fppj,zero,iyjx2)
!      do i=1,iy
!      do j=1,jx
!         fppj(i,j)=g_(i,j+1,k,l_)  +
!     1            (g_(i,j,k,l_)-g_(i,j+1,k,l_))*dj(i,j,k,l_)
         !YuP: above is same as below, but only one '*':
!        fppj(i,j)=g_(i,j+1,k,l_)*(1.-dj(i,j,k,l_)) +
!     1            g_(i,j,k,l_)*dj(i,j,k,l_)
!      enddo
!      enddo
      !YuP-101228: Setting fppj as an array rather than a function
      !reduces cpu time - fppj does not need to be calc. along each ray

!.......................................................................
!     Up-down symmetry factor symm= is defined in aingeom
!**bh050820:
!**bh050820:  Trapped-particle bounce-time factor, trapfac=1.,
!**bh050820:  accounting for bounce time twice transitting bounce time.
!**bh091031:  Bounce-averages are divided by bounce times, so this
!**bh091031:  trapfrac factor does not enter into the problem.
!**bh091031:  symm factor though involves allocating half the ray power
!**bh091031:  above the midplane (and half mirrored below), for
!**bh091031:  up-down-symmetric case.
!..................................................................

      trapfac=one

      alf=0.5*pi*qmca**2/symm

!.......................................................................
!     Set indicators of electrons or ions
!.......................................................................

!BH080920      if (kionn.eq.1) then
      if (kiongg(k).eq.k) then
        signi=1.0
      else
        signi=0.0
      endif

!BH080920      if (kelecg.eq.1) then
      if (kelecg.eq.k) then
        signe=1.0
      else
        signe=0.0
      endif

!..................................................................
!     Loop over rays
!..................................................................

      do 10 iray=nray0,nray(krf)
        ipwr=0 ! For MPI: initialize
#ifdef __MPI
         iraykrf= iray + nrayn*(krf-1)
         if(mpisize.gt.1) then
            mpiworker= MOD(iraykrf-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
      !if(mpirank.eq.mpiworker) then
      !PRINT *,'n,iray,krf,mpiworker=',n,iray,krf,mpiworker
      !endif
#endif
#ifdef __MPI
      if(mpirank.eq.mpiworker) then

#endif

!..................................................................
!     Loop of ray elements - jump out if ray element does not contribute
!     to current flux surface.
!..................................................................

        do 20 is=1,nrayelt(iray,krf)

!cc          if(lr_.ne.lloc(is,iray,krf)) go to 20
          lr_=lloc(is,iray,krf)
          l_=indxlr(lr_)
          call tdnflxs(lmdpln(l_)) !-> Determine itl,itu, etc.

          jmin=jminray(is,iray,krf)
          if(jmin.gt.jx)  go to 20

          jmax=jmaxray(is,iray,krf)
          jmax0=min0(jmax,jx-1)
          jmaxp=min0(jmax+1,jx)
          jminm=max0(jmin-1,1)

          locatn=(jjx*(is-1)+jjx*nrayelts*(iray-1))/ibytes+1
          locatn16=(jjx*(is-1)+jjx*nrayelts*(iray-1))/ibytes16+1
          call bcast(db(1:iy,0:jx),zero,iyjxp1)
          call bcast(dc(1:iy,0:jx),zero,iyjxp1)
          call ibcast(ilim1d,0,jx)
          call ibcast(ilim2d,0,jx)

!..................................................................
!     Unpack stored data - see subroutine urfpack
!..................................................................
          !if(urfb_version.eq.1)then ! 2 is the new version developed by YuP
             ! if 1, it will use the original version
             call unpack(ilowp(locatn,krf),8,ilim1(1),jjx)
             call unpack(iupp(locatn,krf),8,ilim2(1),jjx)
             call unpack16(ifct1_(locatn16,krf),8,ifct1(1),jjx)
             call unpack16(ifct2_(locatn16,krf),8,ifct2(1),jjx)
          !endif

!..................................................................
!     rr is the major radius R location of the ray element.
!     ll if the index of the array pol (poloidal angle array).
!..................................................................
          ! rr=wr(is,iray,krf) ! YuP: Not used?
          ll=llray(is,iray,krf)

!..................................................................
!     Determine the length of the ray element, slngth
!..................................................................
          slngth=0.0
          if (is.eq.1) then
            slngth=ws(1,iray,krf)
          elseif (is.lt.nrayelt(iray,krf)) then
            slngth=0.5*(ws(is+1,iray,krf)-ws(is-1,iray,krf))
          elseif (is.eq.nrayelt(iray,krf)) then
            slngth=0.5*(ws(is,iray,krf)-ws(is-1,iray,krf))
          endif
               ![YuP 08-27-2015] Added skipping condition:
               if(slngth.le.0.d0) goto 20 !-> next ray element
               ! From print-out: found that sometimes slngth=0.

!..................................................................
!     Obtain the polarizations..
!..................................................................
          cwz=cwezde(is,iray,krf)
          cwxyp=0.5*(cwexde(is,iray,krf)+cei*cweyde(is,iray,krf))
          cwxym=0.5*(cwexde(is,iray,krf)-cei*cweyde(is,iray,krf))

!..................................................................
!     Compute the cyclotron frequency.
!..................................................................
!**bh930729wce=charge*sbtot(is,iray,krf)/fmass(kelecg)/clight
          wce=qmca*sbtot(is,iray,krf)
          wcen_omega= nharm(krf)*wce/omega(krf)
!..................................................................
!     alfa(i) will have in it all the leading factors in the contribution
!     to the B coefficient, including ray power and geometry, see notes.
!     Here we have set what was delpwr(is,iray,krf) in subroutine urfb0
!     to 1, in order to normalize the power absorption. Later in
!     urfdamp0 we will trace the ray and compute a new
!     delpwr(is,iray,krf) and this will be used to scale what
!     is computed HERE.
!..................................................................

          alfa_work= urfmult*wnpar(is,iray,krf)**2 &
            *slngth*sbtot(is,iray,krf)/(dpsi(lr_) &
            *fluxn(is,iray,krf)*twopi*wdnpar(is,iray,krf)*omega(krf)) &
            *clight

!..................................................................
!     Determine most of the argument of the Bessel functions, also
!     some leading coefficient factors.
!..................................................................

          do 40 j=jmin,jmax
            i2=ilim2(j)
            i1=min0(iy,ilim1(j))
            rrr= (1.d0-wcen_omega*gammi(j))**2
            div=1.
            if (j.gt.1) then
              div=div+1.
              rrr=rrr+(1.d0-wcen_omega*gammi(j-1))**2
            endif
            if (j.lt.jx) then
              div=div+1.
              rrr=rrr+(1.d0-wcen_omega*gammi(j+1))**2
            endif
            rrr=rrr/div

            alfaj=alfa_work*truncd(j)*gammi(j)*xcu(j)*vnorm3/rrr
            arg_bssl=wnper(is,iray,krf)*x(j)* &
                     vnorm*omega(krf)/(wce*clight*bsslstp(krf))
!..................................................................
!     Add in the contribution to urfb (B_0,indxlr_)  (innermost loop)
!..................................................................
            do i=i2,i1
              temc2(i)=1.
            enddo
            temc2(i2)=DBLE(ifct2(j))*o65535
            temc2(i1)=DBLE(ifct1(j))*o65535
            tem=wcen_omega/(gamma(j)*bbpsi(ll,lr_))
            do 50 i=i2,i1
              thtf1i=1.
              if (i.ge.itl .and. i.le.itu) then
                thtf1i=.5
              endif
              indx=sinz(i,ll,lr_)*arg_bssl+1
              db(i,j)= alfaj*alf*thtf1i*abs(coss(i,l_)) &
                   *abs( cosz(i,ll,lr_)*cwz*jb0(indx,krf) &
                        +sinz(i,ll,lr_)* &
             ( cwxyp*(signe*jbp1(indx,krf)+signi*jbm1(indx,krf)) &
              +cwxym*(signe*jbm1(indx,krf)+signi*jbp1(indx,krf)) ) &
                         )**2 &
                    * temc2(i)
              sinn_=sinn(i,l_)
              if (sinn_.ne.zero) then
              dc(i,j)=db(i,j)*(tem-sinn_**2)/(x(j)*coss(i,l_)*sinn_)
              endif
 50         continue

 40       continue ! j

!..................................................................
!     Need a slightly larger integration range...
!     This is to accomodate the interpolation that is done in
!     the local damping calculation, below.
!..................................................................

!     Account for extension of j-values by -1.
          do 4 j=jminm,jmax0
            ilim2d(j)=min0(ilim2(j),ilim2(j+1))
            ilim1d(j)=max0(ilim1(j),ilim1(j+1))
 4        continue
!     Extend  i-bounds by 1
          do 5 j=jminm,jmax
            ilim2d(j)=max0(ilim2d(j)-1,1)
            ilim1d(j)=min0(ilim1d(j)+1,iy)
 5        continue
          ilim2d(jmaxp)=ilim2d(jmax)
          ilim1d(jmaxp)=ilim1d(jmax)
!     Account for extension of j-values by +1.
          do 6 j=jmin,jmaxp
            ilim2dd(j)=min0(ilim2d(j),ilim2d(j-1))
            ilim1dd(j)=max0(ilim1d(j),ilim1d(j-1))
 6        continue
          ilim2dd(jminm)=ilim2d(jminm)
          ilim1dd(jminm)=ilim1d(jminm)
          do 7 j=jminm,jmaxp
            ilim2dd(j)=max0(ilim2dd(j),1)
            ilim1dd(j)=min0(ilim1dd(j),iy)
 7        continue

          call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)
          call urfmidv_db(jmin,jmax)
!cc          call bcast(temp1(0,0),zero,iyjx2) ! Not really necessary
          call urfmidv_dc(jmin,jmax)

          do j=jminm-1,jmaxp
         !Setting gftp as an array rather than function reduces cpu time
         ! advnce.f90:  cl(i,j)=.25*vptb(itl,lr_)/vptb(i,lr_)*dc(i,j)
             cl_itlm1= .25*vptb(itl,lr_)/vptb(itl-1,lr_)*dc(itl-1,j)
             cl_itlp1= .25*vptb(itl,lr_)/vptb(itl+1,lr_)*dc(itl+1,j)
             cl_itup1= .25*vptb(itl,lr_)/vptb(itu+1,lr_)*dc(itu+1,j)
             gftp(j)= &
                 +cl_itlm1*eyp5(itl-1,l_)*(fppj(itl,j)-fppj(itl-1,j)) &
              +2.*cl_itlp1*eyp5(itl,l_) * (fppj(itl+1,j)-fppj(itl,j)) &
                 +cl_itup1*eyp5(itu,l_) * (fppj(itu+1,j)-fppj(itu,j))
             ! Note: cl(i,j) is proportional to dc(i,j)
          enddo

!..................................................................
!     Express the velocity flux at (i,j+1/2)
!..................................................................
!
          prf_rayel=0.d0 !sum-up contribution from a given ray element.
          do 80 j=jminm,jmaxp
          sum_i=0.d0
          do 90 i=ilim2dd(j),ilim1dd(j)
              dgfp= db(i,j)*exp5(j)*(g_(i,j+1,k,l_)-g_(i,j,k,l_)) &
                 -db(i,j-1)*exp5(j-1)*(g_(i,j,k,l_)-g_(i,j-1,k,l_))
              if(i .ne. itl .and. i .ne. itu) then
                dgfp= dgfp &
                   +dc(i,j)*0.5*dyi(i,l_)*(fppj(i+1,j)-fppj(i-1,j)) &
                 -dc(i,j-1)*0.5*dyi(i,l_)*(fppj(i+1,j-1)-fppj(i-1,j-1))
              else
                dgfp= dgfp + (gftp(j)-gftp(j-1))
              endif
              thtf2i=1.
              if (i.ge.itl .and. i.le.itu) then
                 thtf2i=2.
              endif
              zfact=.5*zmaxpsii(lr_)*dvol(lr_)*fmass(k)/vnorm2
              sum_i=sum_i + thtf2i*dgfp*tcsgm1(j)*cynt2(i,l_)*zfact
 90       continue
          prf_rayel= prf_rayel +sum_i
 80       continue

          urfpwr(is,iray,krf)=prf_rayel

          if(prf_rayel.ne.0.d0) ipwr(is)=1 ! For MPI

!         See scaleurf in cqlinput_help for explanation of following
!         procedure.
          urfpt=urfpwr(is,iray,krf)+urfpwrc(is,iray,krf) &
                +urfpwrl(is,iray,krf)
          if (urfpt.gt.one .and. scaleurf.eq."enabled")  then
            scalurf(is,iray,krf)=1.d0/urfpt
            urfp1=urfpwr(is,iray,krf)
            urfp2=urfpwrc(is,iray,krf)
            urfp3=urfpwrl(is,iray,krf)
            urfpwr(is,iray,krf)=urfp1/urfpt
            urfpwrc(is,iray,krf)=urfp2/urfpt
            urfpwrl(is,iray,krf)=urfp3/urfpt
          endif

!.............................................................
!     Overwrite salphac if iurfcoll.eq."damp_out".  This is added
!      for code diagnostic purposes and gives o/p to the
!      netcdf rf file.
!.............................................................
          if (iurfcoll(krfn(krf)).eq."damp_out") then
             if (slngth.eq.zero) then
                salphac(is,iray,krf)=zero
             else
                salphac(is,iray,krf)=urfpwr(is,iray,krf)/slngth
             endif
          endif

 20     continue ! do is=1,nrayelt(iray,krf)

#ifdef __MPI
         !PRINT*,'SEND_urfpwr: mpirank,krf,iray=',mpirank,krf,iray
         ! Count number of elements for this ray
         irayis=0
         do is=1,nrayelt(iray,krf)! Loop over ray elements
            if(ipwr(is).eq.1) irayis=irayis+1 ! incremented up to nrayis
         enddo
         nrayis=irayis
         mpisz=nrayis ! number of elements
         irayis=0
         do is=1,nrayelt(iray,krf)! Loop over ray elements
           if(ipwr(is).eq.1)  then
             irayis=irayis+1 ! incremented up to nrayis
             urftmp(0*mpisz+irayis)= urfpwr(is,iray,krf)
             urftmp(1*mpisz+irayis)= urfpwrc(is,iray,krf)
             urftmp(2*mpisz+irayis)= urfpwrl(is,iray,krf)
             urftmp(3*mpisz+irayis)= scalurf(is,iray,krf)
             urftmp(4*mpisz+irayis)= salphac(is,iray,krf)
           endif
         enddo
         mpitag= iray + nrayn*(krf-1) ! combined: ray-number + wave-mode
         call MPI_SEND(ipwr,nrayelts, MPI_INTEGER1,0,mpitag,MPI_COMM_WORLD,mpiierr)
         call MPI_SEND(urftmp,5*mpisz,MPI_DOUBLE_PRECISION,0,mpitag,MPI_COMM_WORLD,mpiierr)
         !PRINT*,'SEND_urfpwr: krf,iray,mpirank=',krf,iray,mpirank
      !-----------------------------------------------------------

#endif
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif
#ifdef __MPI
      if(mpirank.eq.0) then !-------------------------------------------
         !PRINT*,'recv_urfpwr: mpirank,krf,iray=',mpirank,krf,iray
         call MPI_RECV(ipwr,nrayelts,MPI_INTEGER1,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,mpiierr)
         mpitag=mpistatus(MPI_TAG)
         call MPI_RECV(urftmp,nrayelts*5,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,mpitag,MPI_COMM_WORLD,mpistatus,mpiierr)
         mpiray=MOD(mpitag-1,nrayn)+1  ! determine which ray
         mpikrf=(mpitag-mpiray)/nrayn +1 ! determine which krf wave mode
         ! Get mpisz5 (Number of elements received)
         call MPI_GET_COUNT(mpistatus,MPI_DOUBLE_PRECISION,mpisz5,mpiierr)
         mpisz= mpisz5/5 ! urftmp contains 5 arrays
         irayis=0
         do is=1,nrayelt(mpiray,mpikrf)! Loop over ray elements
           if(ipwr(is).eq.1)  then
             irayis=irayis+1 ! incremented up to nrayis==mpisz
             urfpwr(is,mpiray,mpikrf)=  urftmp(0*mpisz+irayis)
             urfpwrc(is,mpiray,mpikrf)= urftmp(1*mpisz+irayis)
             urfpwrl(is,mpiray,mpikrf)= urftmp(2*mpisz+irayis)
             scalurf(is,mpiray,mpikrf)= urftmp(3*mpisz+irayis)
             salphac(is,mpiray,mpikrf)= urftmp(4*mpisz+irayis)
           endif
         enddo
         !PRINT*,'recv_urfpwr: mpikrf,mpiray,mpisz=', &
         !                    mpikrf,mpiray,mpisz
      endif !-----------------------------------------------------------

#endif
 10   continue ! do iray=nray0,nray(krf)


      return
      end subroutine urfdamp2


end module urfdamp2_mod
