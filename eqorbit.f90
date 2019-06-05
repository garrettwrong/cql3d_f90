module eqorbit_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use bcast_mod, only : bcast
  use eqfpsi_mod, only : eqfpsi
  use eqjac_mod, only : eqjac
  use eqrhs_mod, only : eqrhs
  use eqwrng_mod, only : eqwrng
  use exlin_mod, only : exlin
  use zcunix_mod, only : coeff1
  use zcunix_mod, only : terp1
  use zcunix_mod, only : terp2

  !---END USE

!
!

contains

      subroutine eqorbit(epsicon_)
      use param_mod
      use cqlcomm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

      dimension slrz(2),rwork(80),iwork(30),epsierr(lfielda)

!     Temp storage
      dimension d2bpsi_(:),d2solr_(:),d2solz_(:)
      dimension es_new(:),bpsi_new(:),solr_new(:),solz_new(:)
      pointer d2bpsi_,d2solr_,d2solz_
      pointer es_new,bpsi_new,solr_new,solz_new


      data iwarn/0/

!................................................................
!     This routine determines field line data, namely B(s)/B(0)
!     where s is the field line arc length by integrating out
!     the orbit equations from the equilibrium solution for the
!     poloidal flux coordinate epsi.
!     "Orbit" equations are following along the B-field, and the
!     R,Z values of the flux surface are traced out.
!..................................................................


!     Setup temp storage for non-updown-symm case
      if (eqsym.eq."none") then
         allocate( d2bpsi_(lfield),STAT=istat)
         if (istat.ne.0)  STOP 'eqorbit allocate problem'
         call bcast(d2bpsi_,zero,SIZE(d2bpsi_))
         allocate( d2solr_(lfield),STAT=istat)
         if (istat.ne.0)  STOP 'eqorbit allocate problem'
         call bcast(d2solr_,zero,SIZE(d2solr_))
         allocate( d2solz_(lfield),STAT=istat)
         if (istat.ne.0)  STOP 'eqorbit allocate problem'
         call bcast(d2solz_,zero,SIZE(d2solz_))

         allocate( es_new(lfield),STAT=istat)
         if (istat.ne.0)  STOP 'eqorbit allocate problem'
         call bcast(es_new,zero,SIZE(es_new))
         allocate( bpsi_new(lfield),STAT=istat)
         if (istat.ne.0)  STOP 'eqorbit allocate problem'
         call bcast(bpsi_new,zero,SIZE(bpsi_new))
         allocate( solr_new(lfield),STAT=istat)
         if (istat.ne.0)  STOP 'eqorbit allocate problem'
         call bcast(solr_new,zero,SIZE(solr_new))
         allocate( solz_new(lfield),STAT=istat)
         if (istat.ne.0)  STOP 'eqorbit allocate problem'
         call bcast(solz_new,zero,SIZE(solz_new))
      endif

      if(epsicon_.le.psilim) then ! YuP[2015/05/03] Just in case:
         write(*,*)'eqorbit: WARNING: epsicon_.le.psilim'
         write(*,*)'eqorbit: epsicon_, psilim =', epsicon_, psilim
         write(*,*)'eqorbit: Resetting epsicon_ to be inside LCFS'
         epsicon_= psilim + 0.001*(psimag-psilim)
         write(*,*)'eqorbit: epsicon_, psilim =', epsicon_, psilim
      endif

!..................................................................
!     Determine f(psi) and df(psi)/dpsi.
!..................................................................

      call eqfpsi(epsicon_,fpsi_,fppsi_)

!..................................................................
!     Determine the first radial mesh point iless such that the
!     psi(iless,zmag) is less than epsicon_ (the psi associated with
!     the contour of interest), that is, point is outside epsicon_
!     contour.  The search is performed at the magnetic axis height,
!     zmag, at the er() radial grid points.
!..................................................................

      do 10 i=imag+1,nnr
         psi2=terp2(er(i),zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
              epsirz,nnra,0,0)

         if (psi2.lt. epsicon_) then
            iless=i
            go to 11
         elseif (i.eq.nnr) then
            call eqwrng(1)
         endif
 10   continue
 11   continue

!..................................................................
!     Do a Newton's iteration to find the R=rcon where the contour
!     intersects the Z=zmag plane.
!..................................................................

      iter=0
      imore=iless-1
      ivalue=imore
      if (imore.eq.imag) ivalue=iless
      rcon=er(ivalue)
      epsitst=terp2(rcon,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
           epsirz,nnra,0,0)
      if( ivalue.eq.imore) then
         psi1=epsitst
      else
         psi1=terp2(er(iless),zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
           epsirz,nnra,0,0)
      endif

 15   continue
      dpsidr=terp2(rcon,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,1,0)
      a=(epsicon_-epsitst)/dpsidr

!..................................................................
!     New value for rcon..
!..................................................................

      rcon=rcon+a
      epsitst=terp2(rcon,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,0,0)
      err=abs((epsicon_-epsitst)/(psi1+psi2))
      if (err .lt. 1.e-6) go to 20
      iter=iter+1
      if (iter.gt.100 .and. err.lt.0005) go to 20
      if (iter.gt.150) call eqwrng(2)
      go to 15
 20   continue

!..................................................................
!     (rcon,zmag) is the coordinate of the initial value for the orbit
!     integrator. Before doing the integration make a guess as to the
!     length of the field line.
!     Guess B-pol for this flux surface; Guess field line length.
!     From these determine a step size for the integrator which
!     will hopefully finish integration after lfield steps or so.
!..................................................................

      rad=abs(rcon-rmag)
      dpsidr=terp2(rmag,rad,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,1,0)
      dpsidz=terp2(rmag,rad,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,0,1)
      bpoltst=sqrt(dpsidr**2+dpsidz**2)/rmag
      qsafty=(rad*fpsi_/rmag)/(rmag*bpoltst)
      if (eqsym.ne."none") then
         zmaxtst=pi*qsafty*rmag
      else
         zmaxtst=2.*pi*qsafty*rmag
      endif
      zstep=zmaxtst/(lfield*0.9) !-YuP Was zmaxtst/(lfield-5)

!..................................................................
!     solr_(nn) and solz_(nn) are R(s) and Z(s), the solutions to
!     the orbit integrator at the nn'th step. At this step
!     s=es_(nn). We will use namelist input variables to test for
!     accuracy. rtol is the relative tolerance parameter, atol
!     is the absolute tolerance parameter, with
!     rtol*abs(solr_(nn) or solz_(nn)) + atol .gt. local error.
!     neq is the number of equations which is 2.
!     methflag=10 for nonstiff Adams method (no Jacobian used).
!     =21 for stiff (BDF) method, user supplied Jacobian.
!     =22 for stiff method, internally generated Jacobian.
!..................................................................

      es_(1)=0.
      neq=2
      solr_(1)=rcon
      solz_(1)=zmag
      s=0.
      slrz(1)=rcon
      slrz(2)=zmag
      sout=zstep
      itol=1
      itask=1
      istate=1
      iopt=0
      lrw=80
      liw=30
      iter=0
      mf=methflag
      nn=2
      nmag=0
!..................................................................
!     Call the orbit solver for the next zstep...
!..................................................................

 30   continue ! handle for iterations

!BH160310
      iwork(6)=2500  !Optional input, max steps during lsode call.
                     !Default is 500.
      iopt=1   ! to use optional input.

!     cursign not used in eqrhs poloidal mag fld components, so
!     also integrate in starting at outer equatorial plane in pos
!     Z-dirn.
      call lsode(eqrhs,neq,slrz,s,sout,itol,rtol,atol,itask, &
        istate,iopt,rwork,lrw,iwork,liw,eqjac,mf)
! YuP This is a standard subroutine from
!     http://www.netlib.no/netlib/odepack/
!     See r8lsode.f - Contains lsode and related subroutines.
!     If such subroutine (and related subroutines) are already
!     present in your library suite, remove r8lsode.f file from
!     your makefile and use your own lsode.

!      if(istate.lt.0)then
!      write(*,*)'eqorbit: after lsode: iter,istate=',iter,istate
!      endif

      solr_(nn)=slrz(1)
      solz_(nn)=slrz(2)
      if (solr_(nn).lt.rmag .and. nmag.eq.0) then
        nmag=nn  ! Gives step where mag axis rmag is passed.
!BH091011: Following variable is z-coord along the epsicon_ field line.
!BH091011: Don't think it is used elsewhere, as such, so changing name
!BH091011: so it doesn't conflict with z-cooord of the magnetic axis.
!BH091011: (Been this way since pre-early 90's).
!BH091011:        zmag=(solz_(nmag)+solz_(nmag-1))*.5
        zmag_nn=(solz_(nmag)+solz_(nmag-1))*.5 !Z-location, passing rmag
      endif
      es_(nn)=sout
      if (nn.eq.2) then
        zfrst=solz_(nn)-zmag
      endif
      zrat=(solz_(nn)-zmag)/zfrst

!..................................................................
!     Test if integration is complete by checking if sign of Z has
!     changed.
!..................................................................

      if (nn.eq.2) istop=0
      if (eqsym.ne."none")  then
         if (zrat.lt.zero) istop=2
      else   ! eqsym.eq."none"
         if (zrat.gt.zero .and. istop.eq.1) istop=2
         if (zrat.lt.zero) istop=1  !Passing zmag in inboard side
      endif
!..................................................................
!     If there are more than lfield/2 orbit points, quit..
!..................................................................

      if (istop.eq.2) then  ! Test on integration achieved end-point
         if (nn.gt.lfield/2) then
            go to 40
         else
            xrat=DBLE(nn)/DBLE(lfield)
            zstep=xrat*zstep
            slrz(1)=rcon
            slrz(2)=zmag
            s=0.
            istate=1
            itask=1
            sout=zstep
            istate=1
            nn=2
            iter=iter+1
            nmag=0

!..................................................................
!     If zstep has been changed too often, call eqwrng.
!..................................................................

!BH160310            if (iter.gt.10) call eqwrng(3)
            if (iter.gt.25) call eqwrng(3)
            go to 30
         endif

!...................................................................
!     If integrator has needed too many orbit points, increase zstep.
!     (Allow for addition of one point in eqsym.eq."none" case below)
!..................................................................

      elseif (nn+2.gt.lfield) then
        zstep=zstep*10.
        sout=zstep
        s=0.
        itask=1
        istate=1
        nn=2
        slrz(1)=rcon
        slrz(2)=zmag
        iter=iter+1
        nmag=0
        if (iter.gt.10) call eqwrng(3)
        go to 30


!..................................................................
!     If calculation is not finished, do the next step..
!..................................................................

      else
        sout=sout+zstep
        nn=nn+1
        go to 30
      endif
 40   continue

!..................................................................
!     The orbit integrator has overshot zmax_, the maximum orbit
!     point. Backtrack by linear interpolation.
!..................................................................

      delz=solz_(nn-1)-solz_(nn)
      zratio=(solz_(nn-1)-zmag)/delz
!     similar triangles:
      es_(nn)=es_(nn-1)+zratio*(-es_(nn-1)+es_(nn))
      solz_(nn)=zmag
!     "l"/"u" suffices refer to lower/upper in es_-traj.
      drdsl=(solr_(nn-1)-solr_(nn-2))/(es_(nn-1)-es_(nn-2))
      sl=(es_(nn-2)+es_(nn-1))*.5
      if (eqsym.eq."none") then
         drdsu=(solr_(nn)-solr_(nn-1))/(es_(nn)-es_(nn-1))
      else
         drdsu=zero
      endif
      su=es_(nn)
      sval=(es_(nn)+es_(nn-1))*.5
!     project dr/ds derivative to sval
      drdsval=exlin(drdsl,drdsu,sl,su,sval)
      solr_(nn)=solr_(nn-1)+drdsval*(es_(nn)-es_(nn-1))

      if (eqsym.eq."none") then
!        First, print out how well the orbit integrator did at closure:
!BH091031         write(*,*)'eqorbit: solr_(1),solr_(nn)=',solr_(1),solr_(nn)

!        Adjust endpoint values to be max of the two (so don't create
!        lower min-|B| point on the FS).
!        [Could smooth our the join abit, if necessary.
!        Alternatively, start orbit integration at different location,
!        say, at 45 or 90 deg poloidal angle rays.]
!YuP         solr_(1)=max(solr_(1),solr_(nn))
!YuP         solr_(nn)=solr_(1)

         ! YuP[2015/05/03] Modified the adjustment of 1st/last point:
         ! At small rho, the "orbit" usually "drifts" towards magn.axis,
         ! so that solr_(nn)<solr_(1).
         ! But at outer surfaces, the "orbit" diverges to the outside,
         ! so that solr_(nn)>solr_(1).
         ! In such a case, the orbit can get outside of LCFS,
         ! if the starting point solr_(1) is close to the LCFS.
         ! So, setting solr_(1)=max(solr_(1),solr_(nn)) would have
         ! a bad result. Simply move the last point to the starting point:
         solr_(nn)=solr_(1)
      endif

      lorbit_=nn

!      write(*,*)
!      write(*,*)'eqorbit: lorbit_=',lorbit_
!      do j=1,lorbit_
!         write(*,*)'j,solr_(j),solz_(j)=',
!     +        j,solr_(j),solz_(j)
!      enddo

      !YuP-120402: combine two last steps if the last step is too small
      if(      es_(lorbit_)-es_(lorbit_-1)   .lt. &
          0.2*(es_(lorbit_-1)-es_(lorbit_-2)) ) then
          ! The last step is < 0.2 of previous step.
          es_(lorbit_-1)= es_(lorbit_) ! Merge two last points together.
          solr_(lorbit_-1)= solr_(lorbit_)
          solz_(lorbit_-1)= solz_(lorbit_)
          lorbit_= lorbit_-1 ! Less one point.
      endif


!..................................................................
!     Compute bpsi_(nn)=B(es_(nn))/B(0.) and thtpol_(nn), the poloidal
!     angle. Also determine dl=sqrt(dr**2+dz**2) and B-pol (eqbpol_).
!..................................................................

      ! l=1 point:
      l=1
      bpsi_(1)=1.d0
      thtpol_(1)=0.d0
      eqdell_(1)=0.d0
      eqdells_=0.d0
      dpsidr=terp2(solr_(1),solz_(1),nnr,er,nnz,ez,epsi, &
          epsirr,epsizz,epsirz,nnra,1,0)
      dpsidz=terp2(solr_(1),solz_(1),nnr,er,nnz,ez,epsi, &
          epsirr,epsizz,epsirz,nnra,0,1)
      az=dpsidr**2+dpsidz**2
      eqbpol_(1)=sqrt(az)/solr_(1)
      bval=sqrt(az+(fpsi_)**2)/solr_(1)
      bmidplne_=bval

      lbad=0

      do 50 l=2,lorbit_
        dpsidr=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi, &
          epsirr,epsizz,epsirz,nnra,1,0)
        dpsidz=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi, &
          epsirr,epsizz,epsirz,nnra,0,1)
!        epsierr(l)=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi,
!     1    epsirr,epsizz,epsirz,nnra,0,0)  !NOT USED
        az=dpsidr**2+dpsidz**2
        eqbpol_(l)=sqrt(az)/solr_(l)
        bval=sqrt(az+(fpsi_)**2)/solr_(l)
        bpsi_(l)=bval/bmidplne_

        if(bpsi_(l).le.1.d0 .and. lbad.eq.0)then  !Added YuP[2015/05/03]
           lbad= l ! remember this point number
        endif
        ! Sometimes, at ending point, in case of eqsym=none,
        ! an orbit/surface cannot get to the starting point
        ! because of accuracy of the orbit integrator.
        ! Orbit gets to a smaller R, where the value of bval/bmidplne_
        ! is a bit smaller than at the starting point, say 0.9999
        ! instead of 1.d0. Although the ending point is forced
        ! to be the same as the starting point (later in this subr),
        ! the previous point (lorbit_-1) may still have the value of
        ! bpsi_(lorbit_-1) = 0.9999 or so.
        ! In any case, bpsi_(l) cannot be less than 1.d0 -
        ! the code cannot handle such equilibria.
        if(lbad.gt.2 .and. thtpol_(l-1).gt.pi*1.5) then
          ! This adjustment only works for eqsym=none
          ! (during 2nd half of surface).
          ! Adjust (R,Z) points by a straight line connecting
          ! the first point where bpsi_<1 is detected
          ! to the starting point
          ![to be exact, connecting point(lbad-1) to the point(lorbit_)]
          solr_(l)= solr_(lbad-1) &
                  -(l-lbad+1)*(solr_(lbad-1)-solr_(1))/(lorbit_-lbad+1)
          solz_(l)= solz_(lbad-1) &
                  -(l-lbad+1)*(solz_(lbad-1)-solz_(1))/(lorbit_-lbad+1)
          dpsidr=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi, &
                 epsirr,epsizz,epsirz,nnra,1,0)
          dpsidz=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi, &
                 epsirr,epsizz,epsirz,nnra,0,1)
          az=dpsidr**2+dpsidz**2
          eqbpol_(l)=sqrt(az)/solr_(l)
          bval=sqrt(az+(fpsi_)**2)/solr_(l)
          bpsi_(l)=bval/bmidplne_
        endif


        if(eqsym.ne.'none')then ! only for half-surface tracing:
!.......................................................................
!     CQL3D cannot handle multiple wells. If this happens - such as
!     would in a PBX run - create a simple fix that will keep the
!     code from failing.
!     980919:  Can't specifically remember what the failure was.
!              But this simple fix causes runs with radial coord.
!              proportional to sqrt(area)... to fail when there is
!              minor non-monoticity in bpsi_.  The reduction of lorbit_
!              reduces the area of the flux surface causing
!              non-monoticity of the radial coordinate.  For now,
!              simply print out warnings of non-monoticity in bpsi_
!              (BobH).
!     YuP[2015/05/03]: adjusting bpsi_(l) to maintain monoticity.
!.......................................................................
            if (bpsi_(l-1).gt.bpsi_(l)) then
               iwarn=iwarn+1
               if (iwarn.eq.1) then
                 write(*,1000) l,rmag,rcon ! XXX had , bug  YuP:agreed
                 !print*,bpsi_(1:l)
               end if
               bpsi_(l)=bpsi_(l-1)+em40 !YuP[2015/05/03] redefine: increasing
 1000 format(//,1x,'eqorbit/WARNING: Non-Monoton B(s)/B(0)',i6,2e17.10)
            endif
        endif ! eqsym

        eqdell_(l)=sqrt((solr_(l)-solr_(l-1))**2+ &
            (solz_(l)-solz_(l-1))**2)
        eqdells_=eqdells_+eqdell_(l)
        if (eqsym.ne."none") then  !up-down symmetrize
             thtpol_(l)=abs(atan2(solz_(l)-zmag,solr_(l)-rmag)) ! in [0,pi]
        else                  !i.e., no up-down symmetrization
             thtpol_(l)=atan2(solz_(l)-zmag,solr_(l)-rmag)
             if (thtpol_(l).lt.zero) thtpol_(l)=thtpol_(l)+twopi !in [0,pi2)
        endif
        if (thtpol_(l).gt.pi*.5 .and. thtpol_(l-1).lt.pi*.5) lpi2=l

 50   continue ! l


      if (eqsym.ne."none") then
         thtpol_(lorbit_)=pi
      else
         thtpol_(lorbit_)=2*pi
      endif
! 60   continue

!......................................................................
!     zmax_ is field line length from 0 to pi pol angle(eqsym.ne."none")
!                                     0 to 2*pi        (eqsym.eq."none")
!     es_bmax_ is fld line len from min to max |B| point, in either case
!......................................................................

      zmax_=es_(lorbit_)
      if (eqsym.ne."none") es_bmax_=abs(zmax_)

!..................................................................
!     set es_ and zmax_ explicitly to positive values (BobH,970620)
!..................................................................

      zmax_=abs(zmax_)
      do l=1,lorbit_
         es_(l)=abs(es_(l))
      enddo
!.......................................................................
!     If eqsym.ne."none", then we can assume the minimum |B| is on the
!     equatorial plane flux surface at zmag=0. where solz_(1)=0.
!     The maximum and minimum R points on the flux surface are at
!     solr_(1) and solr_(nn), respectively.
!
!     For eqsym.eq."none", the field line trace was started at rcon,zmag,
!     in general not at the minimum magnetic field/max major radius
!     point on the flux surface.  We find the minimum |B| point on
!     the flux surface, and then re-grid the above data onto distance grid,
!     es_(nn), shifted to start at the minimum-B point on the flux
!     surface.
!.......................................................................

      if (eqsym.ne."none") then
         rmcon_=solr_(lorbit_)
         rpcon_=solr_(1)
         zmcon_=solz_(lorbit_)
         zpcom_=solz_(1)
         lbpsi_min_= 1         ! YuP added [140710]
         lbpsi_max_= lorbit_   ! YuP added [140710]
         bpsi_min_=  bpsi_(lbpsi_min_) ! YuP added [140710] ==1.0
         bpsi_max_=  bpsi_(lbpsi_max_) ! YuP added [140710] ==Bmax/Bmin

      else                      !eqsym.eq."none"

!        Spline bpsi(l) vs es(l), and find minimum by Newton interation.
!        Use periodic bc's (i1p(1:2)=3)
         i1p(1)=3
         i1p(2)=3
         call coeff1(lorbit_,es_,bpsi_,d2bpsi_,i1p,1,work)
!        Find index of min bpsi_:
         call aminmx(bpsi_,1,lorbit_,1,bpsi_min_,bpsi_max_, &
                     lbpsi_min_,lbpsi_max_)

!         do l=1,lorbit_
!          write(*,'(a,i5,4e15.8)')'eqorbit: l,R,Z,thtpol,bpsi',
!     +     l,solr_(l)-rmag,solz_(l)-zmag,thtpol_(l),bpsi_(l)
!         enddo
!         write(*,*)'lbpsi_min_,lbpsi_max_,bpsi_min_,bpsi_max_',
!     +              lbpsi_min_,lbpsi_max_,bpsi_min_,bpsi_max_
!         pause

!        Iterate to find location of min spline value of bpsi_ value:
!        (Taylor expand in distance s around initial point, take deriv,
!        and set it = 0.)
!        Remember, es_ can jump near this point:
!          If es_min_.lt.zmax_/2, add zmax_, so es_ is continuous in vicinity.
         zmaxd2=zmax_/two
         es_min_=es_(lbpsi_min_)
!         if (es_min_.lt.zmaxd2) es_min_=es_min_+zmax_
         itab(1)=1
         itab(2)=1
         itab(3)=1
         do iter=1,8
            call terp1(lorbit_,es_,bpsi_,d2bpsi_,es_min_,1,tab,itab)
            es_min_=es_min_-tab(2)/tab(3)
            if (es_min_.lt.zero) es_min_=es_min_+zmax_
            if (es_min_.ge.zmax_) es_min_=es_min_-zmax_
!BH091031            write(*,*)'eqorbit: iter,es_min_,bpsi_min_',es_min_,tab(1)
         enddo

!        Reset zero location and re-grid es_ array:
         des_new=zmax_/(lorbit_-1)
         es_new(1)=zero
         do l=2,lorbit_
            es_new(l)=(l-1)*des_new
         enddo

!        Re-grid orbit data:
         call coeff1(lorbit_,es_,solr_,d2solr_,i1p,1,work)
         call coeff1(lorbit_,es_,solz_,d2solz_,i1p,1,work)

         itab(1)=1
         itab(2)=0
         itab(3)=0
         do l=1,lorbit_
!MAYBE SHOULD RECALC BPSI_ FROM R,Z, SO GET BPSI_(1)=1.0??
!$$$            call terp1(lorbit_,es_,bpsi_,d2bpsi_,es_new(l),1,tab,itab)
!$$$            bpsi_new(l)=tab(1)
            call terp1(lorbit_,es_,solr_,d2solr_,es_new(l),1,tab,itab)
            solr_new(l)=tab(1)
            call terp1(lorbit_,es_,solz_,d2solz_,es_new(l),1,tab,itab)
            solz_new(l)=tab(1)
!$$$            if (l.eq.1) then
!$$$               thtpol_(1)=0.
!$$$               eqdell_(1)=0.
!$$$               eqdells_=0.
!$$$            else
!$$$               eqdell_(l)=sqrt((solr_(l)-solr_(l-1))**2+
!$$$     1              (solz_(l)-solz_(l-1))**2)
!$$$               eqdells_=eqdells_+eqdell_(l)
!$$$               thtpol_(l)=atan2(solz_(l)-zmag,solr_(l)-rmag)
!$$$               if(thtpol_(l).lt.zero) thtpol_(l)=thtpol_(l)+twopi !in [0,pi2)
!$$$ 1             if(thtpol_(l).gt.pi*.5.and.thtpol_(l-1).lt.pi*.5) lpi2=l
!$$$            endif
!$$$         enddo

         enddo  !On l

         do l=1,lorbit_
!            bpsi_(l)=bpsi_new(l)
            solr_(l)=solr_new(l)
            solz_(l)=solz_new(l)
         enddo

!.......................................................................
!        Re-calc bpsi_ for shifted es_,solr_,solz_ grids.
!.......................................................................
         do  l=1,lorbit_
            dpsidr=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi, &
                 epsirr,epsizz,epsirz,nnra,1,0)
            dpsidz=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi, &
                 epsirr,epsizz,epsirz,nnra,0,1)
!            epsierr(l)=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi,
!     1           epsirr,epsizz,epsirz,nnra,0,0)  !NOT USED
            az=dpsidr**2+dpsidz**2
            eqbpol_(l)=sqrt(az)/solr_(l)
            bval=sqrt(az+(fpsi_)**2)/solr_(l)
            if (l.eq.1) then
               bpsi_(1)=1.
               thtpol_(1)=0.
               bmidplne_=bval
               eqdell_(1)=0.
               eqdells_=0.
            else
               bpsi_(l)=bval/bmidplne_

!.......................................................................
!     CQL3D cannot handle multiple wells. If this happens - such as
!     would in a PBX run - create a simple fix that will keep the
!     code from failing.
!     980919:  Can't specifically remember what the failure was.
!              But this simple fix causes runs with radial coord.
!              proportional to sqrt(area)... to fail when there is
!              minor non-monoticity in bpsi_.  The reduction of lorbit_
!              reduces the area of the flux surface causing
!              non-monoticity of the radial coordinate.  For now,
!              simply print out warnings of non-monoticity in bpsi_
!              (BobH).
!.......................................................................

               if(eqsym.ne.'none')then ! only valid for half-surface case
               if (bpsi_(l-1).gt.bpsi_(l)) then
                 bpsi_(l)=bpsi_(l-1)+em40 !YuP[2015/05/03] redefine: increasing
               endif
               endif

               bpsi_(l)= max(bpsi_(l), 1.d0) ! YuP[2015/05/03]

               eqdell_(l)=sqrt((solr_(l)-solr_(l-1))**2+ &
                    (solz_(l)-solz_(l-1))**2)
               eqdells_=eqdells_+eqdell_(l)
               if (eqsym.ne."none") then !up-down symmetrize
                  thtpol_(l)=abs(atan2(solz_(l)-zmag,solr_(l)-rmag)) ! [0,pi]
               else             !i.e., no up-down symmetrization
                  thtpol_(l)=atan2(solz_(l)-zmag,solr_(l)-rmag)
                  if (thtpol_(l).lt.zero) thtpol_(l)=thtpol_(l)+twopi ! [0,pi2)
               endif
               if (thtpol_(l).gt.pi*.5.and.thtpol_(l-1).lt.pi*.5) lpi2=l
            endif               ! On l

         enddo                  ! On l
         thtpol_(lorbit_)=2*pi

!        Re-spline bpsi_
         call coeff1(lorbit_,es_,bpsi_,d2bpsi_,i1p,1,work)
         call aminmx(bpsi_,1,lorbit_,1,bpsi_min_,bpsi_max_, &
                     lbpsi_min_,lbpsi_max_)

!        Determine major radius at inboard and outboard side of flux
!        surface.  We will define it by the min and max |B|-points
         rpcon_=solr_(1)
         zpcon_=solz_(1)
!        Iterate for max |B| point.
         es_bmax_=es_(lbpsi_max_)
         itab(1)=1
         itab(2)=1
         itab(3)=1
         do iter=1,8
            call terp1(lorbit_,es_,bpsi_,d2bpsi_,es_bmax_,1,tab,itab)
            es_bmax_=es_bmax_-tab(2)/tab(3)
!BH091031            write(*,*)'eqorbit: iter,es_bmax_,bpsi_max_',es_bmax_,tab(1)
         enddo
         bpsi_max_=tab(1)

!        Re-spline solr_,solz_ to obtain rmcon_,zmcon_
         call coeff1(lorbit_,es_,solr_,d2solr_,i1p,1,work)
         call coeff1(lorbit_,es_,solz_,d2solz_,i1p,1,work)
         itab(1)=1
         itab(2)=0
         itab(3)=0
         call terp1(lorbit_,es_,solr_,d2solr_,es_bmax_,1,tab,itab)
         rmcon_=tab(1)
         call terp1(lorbit_,es_,solz_,d2solz_,es_bmax_,1,tab,itab)
         zmcon_=tab(1)

!......................................................................
!     The maximum |B| point is critical in the bounce-averaged
!     calculation, along with minimum-|B| point, as they define
!     the trapped-passing boundary.
!     Add extra point to es_() array at es_bmax_, and save array location
!     and related data.
!......................................................................

         if (es_bmax_.gt.es_(lbpsi_max_))  lbpsi_max_=lbpsi_max_+1

         do l=lorbit_,lbpsi_max_,-1
            es_(l+1)=es_(l)
            bpsi_(l+1)=bpsi_(l)
            solr_(l+1)=solr_(l)
            solz_(l+1)=solz_(l)
            eqdell_(l+1)=eqdell_(l)
            thtpol_(l+1)=thtpol_(l)
            eqbpol_(l+1)=eqbpol_(l)
         enddo
         l=lbpsi_max_
         es_(l)=es_bmax_
         bpsi_(l)=bpsi_max_
         solr_(l)=rmcon_
         solz_(l)=zmcon_
         thtpol_(l)=atan2(solz_(l)-zmag,solr_(l)-rmag)
         if (thtpol_(l).lt.zero) thtpol_(l)=thtpol_(l)+twopi !in [0,pi2)

         dpsidr=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi, &
              epsirr,epsizz,epsirz,nnra,1,0)
         dpsidz=terp2(solr_(l),solz_(l),nnr,er,nnz,ez,epsi, &
              epsirr,epsizz,epsirz,nnra,0,1)
         az=dpsidr**2+dpsidz**2
         eqbpol_(l)=sqrt(az)/solr_(l)
!     bval=sqrt(az+(fpsi_)**2)/solr_(l)  !Use above value
!     bpsi_(l)=bval/bmidplne_            !Use above value

!        Minus two previous eqdell_(), plus three new eqdell_()
         eqdells_=eqdells_-eqdell_(l-1)-eqdell_(l+1)
         do l=lbpsi_max_-1,lbpsi_max_+1
            eqdell_(l)=sqrt((solr_(l)-solr_(l-1))**2+ &
                 (solz_(l)-solz_(l-1))**2)
            eqdells_=eqdells_+eqdell_(l)
         enddo
         lorbit_=lorbit_ + 1
         zmax_=es_(lorbit_)


         !YuP[2015/05/03] Another adjustment:
         ! a very persistent bad point (lorbit_-1)
         if( bpsi_(lorbit_-1) .le. bpsi_(lorbit_) ) then
            bpsi_(lorbit_-1)= 0.5*(bpsi_(lorbit_-2)+bpsi_(lorbit_))
            ! bpsi_ array should be strictly monotonic
            ! in each half of the surface.
            ! Two points with same value of bpsi_() will trigger
            ! a warning message from luf()
         endif

!         do l=1,lorbit_
!          write(*,'(a,2i5,3e15.8)')'eqorbit:',
!     +     l_,l,solr_(l),thtpol_(l),bpsi_(l)
!         enddo
!         write(*,*)'lbpsi_max_,bpsi_min_,bpsi_max_',
!     +              lbpsi_max_,bpsi_min_,bpsi_max_
!      write(*,*)'eqorbit: epsicon_,solr_(1)',epsicon_,solr_(1)


      endif  ! On eqsym


!..................................................................
!     Determine bpol and btor at poloidal angle=90 deg (R=rmag)
!..................................................................

      bthr_=(eqbpol_(lpi2)+eqbpol_(lpi2-1))*.5
      btoru_=fpsi_/rmag

!..................................................................
!     Determine bpol and btor at the outer midplane.
!..................................................................

      bthr0_=eqbpol_(1)
      btor0_=fpsi_/rpcon_

      if (eqsym.eq."none") then
         deallocate(d2bpsi_,d2solr_,d2solz_,STAT=istat)
         deallocate(es_new,bpsi_new,solr_new,solz_new,STAT=istat)
      endif

      return
      end subroutine eqorbit


end module eqorbit_mod
