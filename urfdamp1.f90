module urfdamp1_mod

  !---BEGIN USE

  use tdnflxs_mod, only : tdnflxs
  !XXXXXXXXX these require disbaling TKR checks !!
  !use urfpackm_mod, only : unpack
  !use urfpackm_mod, only : unpack16
  external unpack
  external unpack16

  !---END USE

!
!

contains

      subroutine urfdamp1(krf)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
!yup      save

!.......................................................................
!     This routine computes the damping for all
!     ray elements on all rays for UNIT power delpwr. This is put into
!     array urfpwr(is,iray,krf). The algorithm follows that of urfb0.
!
!     Power absorption for each ray element is obtained by integration
!     over the resonance region of velocity space of energy*(df/dt)_QL.
!     Local values of the distribution function are obtained from
!     midplane values.  The algorithm follows that of urfb0.
!.......................................................................

!MPIINSERT_INCLUDE

      complex*16 cwz,cwxyp,cwxym,cei

      integer*1 ipwr(nrayelts) ! YuP: For MPI:  indicator  (0 or 1)
                               ! that a ray element contributed any power

      data nray0 /1/

      cei=(0.,1.)

!..................................................................
!     fill out increasing psi array of bin boundaries, as in urffflx.
!..................................................................
         do l=1,lrzmax
            tr2(l)=psimag-psivalm(l) ! store in comm.h
         enddo
!.......................................................................
!     Up-down symmetry factor symm= is defined in aingeom
!**bh050820:
!**bh050820:  Trapped-particle bounce-time factor, trapfac, here set =1.,
!**bh050820:  accounts for bounce time being twice the transitting bounce
!**bh091031:  time.  However, bounce-averages are divided by bounce times,
!**bh091031:  so this trapfac factor cancels out.
!**bh091031:  The symm factor though involves allocating half the ray power
!**bh091031:  above the midplane (and half mirrored below), for
!**bh091031:  up-down-symmetric case.
!..................................................................

      trapfac=one
      o65535= 1.d0/65535.d0
      vnorm2i= 1.d0/vnorm2

!     k is general species to which this krf mode is to be applied:
      k=nrfspecies(krfn(krf))
      qmc= charge*bnumb(k)/(fmass(k)*clight) !store in common/orb0/
      qmca= abs(qmc)
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
!MPIINSERT_MPIWORKER_IRAYKRF
!MPIINSERT_IF_RANK_EQ_MPIWORKER
        prf_ray_sum=0.d0 !For print-out only: sum of prf_rayel over rayel

!..................................................................
!     Loop of ray elements - jump out if ray element does not contribute
!     to current flux surface.
!..................................................................

        do 20 is=1,nrayelt(iray,krf)

!cc          if(lr_.ne.lloc(is,iray,krf)) go to 20

!..................................................................
!     rr is the major radius R location of the ray element.
!     ll if the index of the array pol (poloidal angle array).
!..................................................................

          rr=wr(is,iray,krf)
          zz=wz(is,iray,krf)
          bray=sbtot(is,iray,krf) ! B-field at the ray element
          !-> For FOW (some values at the ray-element coord):
          !if (fow.ne.'disabled') then
             R_loc= rr   !-> store in common/orb_loc/
             Z_loc= zz   !-> store in common/orb_loc/
             B_loc=bray  !-> store in common/orb_loc/
          !endif
          !->
          ll=llray(is,iray,krf)
          lrloc=lloc(is,iray,krf) !local to ray element point.
          jmin=jminray(is,iray,krf) !jmin,jmax, etc.,
                                    !are local to ray element point.
          jmax=jmaxray(is,iray,krf)
          jmax0=min(jmax,jx-1)
          anpar0=wnpar(is,iray,krf) ! n_par0
          anpar02=anpar0**2
          dnpar=wdnpar(is,iray,krf) ! delta(n_par)
          anper0=wnper(is,iray,krf) ! n_prp0
          ! Determine the length of the ray element, slngth:
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
          ! Obtain the polarizations:
          cwz=cwezde(is,iray,krf)
          cwxyp=0.5*(cwexde(is,iray,krf)+cei*cweyde(is,iray,krf))
          cwxym=0.5*(cwexde(is,iray,krf)-cei*cweyde(is,iray,krf))
          ! Compute the (cyclotron frequency)/omega  :
          wc_w=qmca*bray/omega(krf)
          wcn_w= nharm(krf)*wc_w
          abssl=anper0/(wc_w*clight*bsslstp(krf))
          del_pwr=delpwr(is,iray,krf)
          scal_urf=scalurf(is,iray,krf)
          alf_rayelt=alf*urfmult*slngth*bray*clight*clight/ &
                    abs(fluxn(is,iray,krf)*twopi*omega(krf)*anpar0)
          alf_rayel=alf_rayelt*(anpar02/dnpar)*abs(anpar0/clight)
!..................................................................
!     alf_rayel has in it all the leading factors in the contribution
!     to B coefficient, including the ray power and geometry, see notes.
!     For examination of possible non-quasilinear effects,
!     urfmult is a multiplier for the QL diffusion coefficients and
!     damping, defaulted to 1.0.
!..................................................................
          ! Define vpar_res (resonance Vpar corresponding to npar0)
          !vpar_res= clight*(1.d0-wcn_w)/anpar0
          ! Check that |vpar_res| is within u-grid
          ! ( |vpar_res| < vnorm  ):
          if( abs(clight*(1.d0-wcn_w)).ge.abs(vnorm*anpar0) &
              .or. jmin.gt.jx) then
                ! Resonance is outside of u-grid.
                ! The above condition includes the case of Npar=0.
                jminray(is,iray,krf)=jx+1 !To indicate: out of grid
                !WRITE(*,*)'iray,is,jmin=', iray,is,jmin
                goto 20 !-> next ray element
          endif
          vpar_res= clight*(1.d0-wcn_w)/anpar0

!          if(urfb_version.eq.2)then ! new version developed by YuP
!              ! if 1, it will continue with the original version
!              lr_=lrloc !default value for ZOW; to be found for FOW.
!              !!!! For now, Only setup for ICRH/Non-Relativistic,
!              !!!! which is the case if(kiongg(k).eq.k)
!              iurfb=1 !1 for calc. prf_rayel
!              call urfb_add(anpar0,vpar_res,
!     +             abssl,alf_rayelt,
!     +             bray,wcn_w,cwz,cwxyp,cwxym,krf,
!     +             iurfb,prf_rayel)
!              goto 200 !-> Define urfpwr(is,iray,krf)=prf_rayel;
!                       ! Then - next ray-element.
!          endif

!..................................................................
!     Unpack stored data - see subroutine urfpack
!..................................................................

          locatn=(jjx*(is-1)+jjx*nrayelts*(iray-1))/ibytes+1
          locatn16=(jjx*(is-1)+jjx*nrayelts*(iray-1))/ibytes16+1
          ! locatn16 specifies the storage location
          ! in the compressed arrays
          !if(urfb_version.eq.1)then ! 2 is the new version developed by YuP
            ! if 1, it will use the original version
            call unpack(ilowp(locatn,krf),8,ilim1(1),jjx)
            call unpack(iupp(locatn,krf),8,ilim2(1),jjx)
            call unpack16(ifct1_(locatn16,krf),8,ifct1(1),jjx)
            call unpack16(ifct2_(locatn16,krf),8,ifct2(1),jjx)
          !endif

          prf_rayel=0.d0 !sum-up contribution from a given ray element.

          do 40 j=jmin,jmax0 !YuP: was jmax !ggfij not def. for j=jx
            i2=ilim2(j)
            i1=min0(iy,ilim1(j))
            ulocg=vnorm*x(j)*gammi(j) ! V (speed, not momentum/mass)
            ! Smooth out the diffusion coefficient in j to avoid
            ! dividing by a near zero value in the evaluation of alfaj.
            !YuP: This "smoothening" is useless in low-T
            !     non-relativ. plasma
            !     where gamma(j)=1 with good accuracy.
            !     As a result, we can get a nearly 1/0 from 1/rrr term
            !     at region of crossing the resonance wcn_w=1.
            rrr= (1.d0-wcn_w*gammi(j))**2
            div=1.
            if (j.gt.1) then
              div=div+1.
              rrr=rrr+(1.d0-wcn_w*gammi(j-1))**2
            endif
            if (j.lt.jx) then
              div=div+1.
              rrr=rrr+(1.d0-wcn_w*gammi(j+1))**2
            endif
            rrr=rrr/div

!$$$ BH: 170101: Reverting back to older method, keeps backwards compatibility.
!$$$                  !YuP[04-2016] Added averaging of rrr over [i2;i1] points.
!$$$                  !The value of sqrt(rrr) is |vpar_res*npar/c|.
!$$$                  !Instead of using |vpar_res|, we can use |v*cos(theta_loc)|
!$$$                  !and average it over few points in theta_local.
!$$$                  lr_=lrloc
!$$$                  l_=indxlr(lr_)
!$$$                  call tdnflxs(lmdpln(l_)) !Determine itl,itu,indxlr_, etc.
!$$$                  upar_loc_sum=0.d0 ! initialize
!$$$                  theta_loc_sum=0.d0
!$$$                  !Adjust the range in i, to be at least 3 pts.
!$$$                  !(from printout: sometimes i2=i1, and it can be vpar~0)
!$$$                  i22=i2
!$$$                  i11=i1
!$$$                  if(i2-i1 .le. 1) then
!$$$                    ! One or two points only. Add two more points:
!$$$                    if(i2.ge.2 .and. i1.le.iy-1) then
!$$$                      i22=i2-1
!$$$                      i11=i1+1
!$$$                    endif
!$$$                    if(i2.eq.1) then
!$$$                      i22=i2
!$$$                      i11=i1+2
!$$$                    endif
!$$$                    if(i1.eq.iy) then
!$$$                      i22=i2-2
!$$$                      i11=i1
!$$$                    endif
!$$$                  endif
!$$$                  do i=i22,i11 ! loop in theta_local, for a given u-level
!$$$                     im1=max(i-1,1)  ! i-1, but never a 0
!$$$                     ip1=min(i+1,iy) ! i+1, but never iy+1
!$$$                     dtheta_loc= 0.5*(yz(ip1,ll,lr_)-yz(im1,ll,lr_))
!$$$                     cos_loc= cosz(i,ll,lr_) ! can be 0 ! cos(theta_local)
!$$$                     upar_loc=ulocg*cos_loc !ulocg is speed vnorm*x(j)/gamma
!$$$                     upar_loc_sum= upar_loc_sum+abs(upar_loc)*dtheta_loc
!$$$                     theta_loc_sum= theta_loc_sum+dtheta_loc
!$$$                  enddo ! iloc
!$$$                  upar_loc_avg= upar_loc_sum/theta_loc_sum
!$$$                  !The range theta_loc_sum in theta_local
!$$$                  !is supposed to be small.
!$$$                  !Now define "averaged" rrr:
!$$$                  rrr= (upar_loc_avg*anpar0/clight)**2
!$$$                  !YuP[04-2016] Added averaging of |upar_loc| and rrr: done
!$$$                  !This part over-writes the original definition of rrr.
!$$$                  !Same procedure could be added to urfdamp2 and vlf.

!.................................................................
!     Determine most of the argument of the Bessel functions, also
!     some leading coefficient factors.
!.................................................................
            arg_bssl= abssl*x(j)*vnorm
            !-------------------------------------------ZOW
            lr_=lrloc
            l_=indxlr(lr_)
            call tdnflxs(lmdpln(l_)) !-> Determine itl,itu, etc.
            do i=i2,i1
               temc2(i)=1.
            enddo
            temc2(i2)=DBLE(ifct2(j))*o65535
            temc2(i1)=DBLE(ifct1(j))*o65535
            uloc3=xcu(j)*vnorm3
            uloc=x(j)*vnorm
            tem=wcn_w/(gamma(j)*bbpsi(ll,lr_))
            alfaj=alf_rayel*truncd(j)*gammi(j)*uloc3/rrr
            sum_i=0.d0
            do 50 i=i2,i1
              thtf1i=1.
              thtf2i=1.
              if (i.ge.itl .and. i.le.itu) then
                thtf1i=.5/trapfac
                thtf2i=2.
              endif
              indx=sinz(i,ll,lr_)*arg_bssl+1
              ! Make certain indx is within table bounds
              indx=min(indx,nbssltbl)
              !!! if (indx.gt.nbssltbl) call urfwrong(1)
              sinn_=sinn(i,l_)
              coss_=coss(i,l_)
              dbij= alfaj*thtf1i/dpsi(lr_) &
                   *abs(coss_) &
                   *abs( cosz(i,ll,lr_)*cwz*jb0(indx,krf) &
                        +sinz(i,ll,lr_)* &
             ( cwxyp*(signe*jbp1(indx,krf)+signi*jbm1(indx,krf)) &
              +cwxym*(signe*jbm1(indx,krf)+signi*jbp1(indx,krf)) ) &
                         )**2 &
                    *temc2(i)

              coss_sinn= coss_*sinn_
              dcij=0.d0
              if(coss_sinn.ne.zero) then
                dcij= dbij*(tem-sinn_**2)/(x(j)*coss_sinn)
              endif
              !YuP[2017-11-22] Added: check proximity to loss cone.
              jp1=min(j+1,jx) !Not exceeding jx (j+/-1 will be used for df/dv)
              jm1=max(j-1,1)  !Not smaller than 1
              ip1=min(i+1,iy)
              im1=max(i-1,1)
              !Check: if at least one of points across differencing
              !is in the loss cone, skip this point:
              if(gone(i,jp1,k,l_)+gone(i,jm1,k,l_).lt.zero) goto 50 !-> next i
              if(gone(ip1,j,k,l_)+gone(im1,j,k,l_).lt.zero) goto 50 !-> next i
              !YuP[2017-11-22] Added: check proximity to loss cone.
              ggfij= &
                   +dbij*0.5*dxi(j)*(f(i,jp1,k,l_)-f(i,jm1,k,l_)) &
                   +dcij*0.5*dyi(i,l_)*(f(ip1,j,k,l_)-f(im1,j,k,l_))
              zfact=zmaxpsii(lr_)*dvol(lr_)*fmass(k)*vnorm2i
!bug          sum_i= sum_i -4.*thtf1i*ggfij*cynt2(i,l_)*zfact
              sum_i= sum_i    -thtf2i*ggfij*cynt2(i,l_)*zfact !ZOW
!         write(*,*) dcij,(f(i+1,j,k,l_)-f(i-1,j,k,l_)),ggfij
!         if(sum_i.lt.0.d0)then
!           write(*,*)'urfdamp1: sum_i<0; ggfij=',ggfij
!           pause
!         endif
 50         continue ! i=i2,i1
            !-------------------------------------------ZOW


            prf_rayel= prf_rayel+ x(j)*gammi(j)*dx(j)*sum_i
 40       continue ! j

200       urfpwr(is,iray,krf)=prf_rayel

          prf_ray_sum= prf_ray_sum+prf_rayel ! For print-out only.

          if(prf_rayel.ne.0.d0) ipwr(is)=1 ! For MPI


!         See scaleurf in cqlinput_help for explanation of following
!         procedure.
          urfpt=urfpwr(is,iray,krf)+urfpwrc(is,iray,krf) &
                +urfpwrl(is,iray,krf)
          if (urfpt.gt.one .and. scaleurf.eq."enabled")  then
            scal_urf=1.d0/urfpt
            scalurf(is,iray,krf)=scal_urf !=1.d0/urfpt
            urfp1=urfpwr(is,iray,krf)
            urfp2=urfpwrc(is,iray,krf)
            urfp3=urfpwrl(is,iray,krf)
            urfpwr(is,iray,krf)= urfp1*scal_urf
            urfpwrc(is,iray,krf)=urfp2*scal_urf
            urfpwrl(is,iray,krf)=urfp3*scal_urf
            ipwr(is)=1 ! For MPI
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
             ipwr(is)=1 ! For MPI
          endif

 20     continue ! do is=1,nrayelt(iray,krf)

!MPIINSERT_IF_RANK_EQ_0
      if(prf_ray_sum.gt.0.1)then
      WRITE(*,'(a,3i4,e12.4)')'URFDAMP1: n,iray,krf, prf_ray_sum=', &
                           n,iray,krf, prf_ray_sum
      endif
!MPIINSERT_ENDIF_RANK


!MPIINSERT_SEND_URFPWR
!MPIINSERT_ENDIF_RANK
!MPIINSERT_RECV_URFPWR
 10   continue ! do iray=nray0,nray(krf)




      return
      end
end module urfdamp1_mod
