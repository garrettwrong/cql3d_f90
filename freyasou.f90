!
!
      subroutine freyasou(qx,qy,qz,qr,vx,vy,vz,iqts,curdep, &
        bmsprd,multiply,multiplyn)
      use bcast_mod, only : bcast
      use bcast_mod, only : ibcast
      use cqlcomm_mod
      use cqlconf_mod, only : setup0
      use eqfpsi_mod, only : eqppsi
      use eqfpsi_mod, only : eqfpsi
      use param_mod
      use r8subs_mod, only : dscal
      use sourcpwr_mod, only :sourcpwr
      use tdnflxs_mod, only : tdnflxs
      use tdtoaray_mod, only : tdtoaray
      use urfb0_mod, only : luf_bin
      use zcunix_mod, only : coeff1
      use zcunix_mod, only : terp1
      use zcunix_mod, only : terp2
      implicit none

!MPIINSERT_INCLUDE

      character*8 multiply,method

      real(c_double) :: qx(*),qy(*),qz(*),qr(*),vx(*),vy(*),vz(*)
      integer llbrth(0:lrza)
      real(c_double) :: curbrth(lrza)
      real(c_double) :: d2vphipl(lrza),workk(3*lrza+1)

      real(c_double) :: ranvc(:)
      pointer ranvc
      real(c_double) :: px(:)
      pointer px
      real(c_double) :: py(:)
      pointer py
      real(c_double) :: pz(:)
      pointer pz
      integer lbrth(:)
      pointer lbrth
      real(c_double) :: psipt(:)
      pointer psipt
      real(c_double) :: pr(:)
      pointer pr
      
      integer i,j,l,ll,k,ipts,istat,mrans,ipar,mpar,jpar !local
      integer mran,ma,multiplyn !local
      integer iless,imore,iter !local
      integer iqts ! in arg. list
      real(c_double) :: twopi2,qmc,qmca,rbrth,rbrth2 ! local
      real(c_double) :: dpsidr,dpsidz,bxcomp,bycomp,bzcomp,bmag2,one_wcb !local
      real(c_double) :: vmag,bmag, bmra,bmsprd,apsi,psibrth,Bpol_factor !local
      real(c_double) :: rhox_gyro,rhoy_gyro,rhoz_gyro ! local
      real(c_double) :: vphi,vphip,vrp,vphiplas,dotprd,a,aa,s,scalfact !local
      real(c_double) :: cosp,sinp,curnorm,curdep,cospitch,rcon
      real(c_double) :: epsitst,err
      real(c_double) :: bzero,psb,bratio,tmdplne

      ! YuP[2019-06-12] Getting lrz and lrzmax from setup0 type:
      integer :: lrz
      integer :: lrzmax
      lrz = setup0%lrz
      lrzmax =  setup0%lrzmax

!
!..................................................................
!
!     This routine computes the bounce averaged source suitable for
!     time advancement from the ion birth points as determined by the
!     NFREYA code as it exists in the transport code ONETWO. For the
!     moment the value of psi at which the particles are born is
!     assumed to be the same as the bounce average of psi. qx qy and
!     qz are the locations of birth in a cartesian system, and vx,vy and
!     and vz are the velocity coordinates. kfrsou is the species index.
!
!..................................................................
!     psivalm(l) is the value of psi that delimits the EXACT OUTER EDGE
!     of flux volume l. The volume, area and incremental volume and
!     area arrays, vol, area, dvol and darea (functions of l) are
!     computed with this boundary in mind. Thus there is a PRECISE
!     correspondence between the volume between psival(l-1) and
!     psival(l) and dvol(l).
!..................................................................

      itab(1)=1
      itab(2)=0
      itab(3)=0
      method="meth1"

      twopi2=twopi*twopi

      k=kfrsou

      ! q/mc for source particles:
      qmc= charge*bnumb(k)/(fmass(k)*clight)
      qmca=abs(qmc)
      do 10 l=1,lrz
        tr(l)=psimag-psivalm(l)
 10   continue
      tr(lrz)=psimag-psilim

!..................................................................
!     If variable multiply.ne."disabled" and multiplyn.ge.1,
!     we assume that each particle birth
!     point represents the mean value of a normally distributed
!     beamlet with standard deviation equal to beamsprd. If multiply
!     is "disabled" then the code goes with the actual NFREYA deposition
!     with no alterations.
!     In the case that multiply is activated, multiply particles are
!     born distributed as described above. Generate space for these new
!     particles and for the random number vector.
!..................................................................

      ma=multiplyn
      if (multiply.eq."disabled") ma=1
      ipts=iqts*ma
      mrans=ipts*3

      write(*,*)'freyasou: allocating ipts 6 times, ipts=',ipts
      allocate(px(ipts),STAT=istat)
      call bcast(px,zero,SIZE(px))
      allocate(py(ipts),STAT=istat)
      call bcast(py,zero,SIZE(py))
      allocate(pz(ipts),STAT=istat)
      call bcast(pz,zero,SIZE(pz))
      allocate(lbrth(ipts),STAT=istat)
      call ibcast(lbrth,0,SIZE(lbrth)) !-YuP: was 'zero', which is 0.d0
      allocate(psipt(ipts),STAT=istat)
      call bcast(psipt,zero,SIZE(psipt))
      allocate(pr(ipts),STAT=istat)
      call bcast(pr,zero,SIZE(pr))
      if (multiply.ne."disabled") then
         allocate(ranvc(mrans),STAT=istat)
         call bcast(ranvc,zero,SIZE(ranvc))


!..................................................................
!     Call random number generator - scale to obtain desired
!     standard deviation
!..................................................................

!     iseed=123457        ***old**bh
!     call rnset(iseed)   ***old**bh
!     call rnnoa(mrans,ranvc)   ***old**bh
        call frnnoa(mrans,ranvc)
        bmra=bmsprd*radmin
        call dscal(mrans,bmra,ranvc,1)

!..................................................................
!     Determine the birth points....
!..................................................................

        mpar=0
        mran=0
        do 30 ipar=1,ipts
          jpar=(ipar-1)/multiplyn+1
          px(ipar)=qx(jpar)+ranvc((ipar-1)*3+1)
          py(ipar)=qy(jpar)+ranvc((ipar-1)*3+2)
          pz(ipar)=qz(jpar)+ranvc((ipar-1)*3+3)
 30     continue
      else
        do 45 ipar=1,ipts
          px(ipar)=qx(ipar)
          py(ipar)=qy(ipar)
          pz(ipar)=qz(ipar)
 45     continue
      endif  !on multiply
      write(*,*)'freyasou: After allocate'

!      write(*,*)
!      write(*,*)'freyasou:i,qx(i),qy(i),qz(i),qr(i),vx(i),vy(i),vz(i):'
!      do i=1,5
!         write(*,*) i,qx(i),qy(i),qz(i),qr(i),vx(i),vy(i),vz(i)
!      enddo

!..................................................................
!     Mark birth particles according to their birth flux surface bin.
!     lbrth(ipar) is the label of the relevant flux surface.
!..................................................................

      do 60 ipar=1,ipts

        rbrth2= px(ipar)**2 +py(ipar)**2
        rbrth=sqrt(rbrth2)
        pr(ipar)=rbrth
        psibrth=terp2(rbrth,pz(ipar),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
          epsirz,nnra,0,0)
        psipt(ipar)=psibrth
        apsi=psimag-psibrth
        if (psibrth.lt.psilim) then  !i.e., outside LCFS
          lbrth(ipar)=0
        else
          lbrth(ipar)=luf_bin(apsi,tr(1:lrz))
        endif

        if(fr_gyrop.eq.'enabled')then ! YuP[03/13/2015]
           ! Note: fr_gyrop=fr_gyro through frnfreya call.
           ! Adjust coordinate from local position to g.c. position
           ! First, Determine components of B:
           dpsidr=terp2(rbrth,pz(ipar),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
             epsirz,nnra,1,0)
           dpsidz=terp2(rbrth,pz(ipar),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
             epsirz,nnra,0,1)
           call eqfpsi(psibrth,fpsi_,fppsi_)
!           Bpol_factor=0d0   !Test: zero Bpol effect, leaving only btor.
           Bpol_factor=1.d0
           bxcomp=(-cursign*Bpol_factor*dpsidz*px(ipar)-bsign*fpsi_* &
                py(ipar))/rbrth2
           bycomp=(-cursign*Bpol_factor*dpsidz*py(ipar)+bsign*fpsi_* &
                px(ipar))/rbrth2
           bzcomp=+cursign*Bpol_factor*dpsidr/rbrth
           if (multiply.eq."disabled") then
              jpar=ipar
           else
              jpar=(ipar-1)/multiplyn+1
           endif
           bmag2= bxcomp**2+bycomp**2+bzcomp**2 ! B^2
           !===> This part is to add gyro-correction
           !    (to find g.c. position from a local position of an ion)
           one_wcb= 1.d0/(qmc*bmag2) ! = 1/(omega_c*B)
           ! Components of rho_gyro are found
           ! from  rho_gyro= (1/omega_c)[v x B]/|B|
           rhox_gyro= one_wcb*(vy(jpar)*bzcomp - vz(jpar)*bycomp)
           rhoy_gyro= one_wcb*(vz(jpar)*bxcomp - vx(jpar)*bzcomp)
           rhoz_gyro= one_wcb*(vx(jpar)*bycomp - vy(jpar)*bxcomp)
           ! The vector rho_gyro is directed from ion position
           ! towards the gyro-center (for a positive charge q)
           !-> Gyro-radius correction:
           px(ipar)= px(ipar) +rhox_gyro ! Now px is at g.c.
           py(ipar)= py(ipar) +rhoy_gyro ! Now py is at g.c.
           pz(ipar)= pz(ipar) +rhoz_gyro ! Now pz is at g.c.
           ! Update other values. Now - at g.c. position:
           rbrth2= px(ipar)**2 +py(ipar)**2
           rbrth=sqrt(rbrth2)
           pr(ipar)=rbrth
           psibrth= terp2(rbrth,pz(ipar),nnr,er,nnz,ez,epsi, &
                          epsirr,epsizz,epsirz,nnra,0,0)
           psipt(ipar)=psibrth ! at g.c. now
           apsi=psimag-psibrth
           if (psibrth.lt.psilim) then  !i.e., g.c. outside LCFS
              !Or maybe still use the actual particle position for lbrth?
              lbrth(ipar)=0
           else
              lbrth(ipar)=luf_bin(apsi,tr(1:lrz))
           endif
           ! Also, update gx,qy,qz - it will update xpts,ypts,zpts
           qx(ipar)=px(ipar)
           qy(ipar)=py(ipar)
           qz(ipar)=pz(ipar)  ! Need to check multiply.eq."enabled" case
        endif ! fr_gyro

 60   continue


!      do ipar=1,10
!         write(*,*)'freyasou: ipar,psipt(ipar),lbrth(ipar)'//
!     +        'vx(ipar),vy(ipar),vz(ipar)=,',
!     +        ipar,psipt(ipar),lbrth(ipar),vx(ipar),vy(ipar),vz(ipar)
!      enddo




!..................................................................
!     Subtract rotation velocity vphipl from the toroidal
!     velocity of the particles, if n.ge.nonvphi .and. n.lt.noffvphi.
!..................................................................

      if ((n+1).ge.nonvphi .and. (n+1).lt.noffvphi) then
        if (multiply.ne."disabled") stop 'rerun with multiply=disabled'
        i1p(1)=4
        i1p(2)=4
        itab(1)=1
        itab(2)=0
        itab(3)=0
        call  coeff1(lrzmax,equilpsp(1),vphipl,d2vphipl,i1p,1,workk)

        do 70 ipar=1,ipts
          psibrth=psipt(ipar)

!         prevent possibility of unhealthy extrapolation
          if (psibrth.lt.psilim) psibrth=psilim
          if(psibrth.gt.(-equilpsp(1))) psibrth=-equilpsp(1)

          psibrth=-psibrth


          call terp1(lrzmax,equilpsp(1),vphipl,d2vphipl,psibrth, &
            1,tab,itab)
          if(rmag.gt.1.d-8)then
            vphiplas=tab(1)*(pr(ipar)/rmag) ! rmag=0 in mirror machine
            cosp=px(ipar)/pr(ipar)
            sinp=py(ipar)/pr(ipar)
            vphi=-sinp*vx(ipar)+cosp*vy(ipar)
            vphip=vphi-vphiplas
            vrp=cosp*vx(ipar)+sinp*vy(ipar)
            vx(ipar)=vrp*cosp-vphip*sinp
            vy(ipar)=vrp*sinp+vphip*cosp
!           vz(ipar)=vz(ipar)
          else
            goto 70 ! skip, for now
          endif

!          if (ipar.lt.10) then
!             write(*,*)'freyasou: ipar,psibrth,vphiplas,'//
!     &            'vx(ipar),vy(ipar),vz(ipar)=,',
!     &            ipar,psibrth,vphiplas,vx(ipar),vy(ipar),vz(ipar)
!          endif

 70     continue ! ipar=1,ipts
      endif ! ((n+1).ge.nonvphi .and. (n+1).lt.noffvphi)

!.......................................................................
!BH110614:  Loop structure is changed.  Main loop is now over
!BH110614:  particles, and data is accumulated for each flux
!BH110614:  surface as advance through the particle (ipar) loop.
!
!xx     Begin loop over the flux surfaces: Recall that information
!xx     pertaining to individual flux surfaces (including S) is
!xx     stored on disk. This introduces some awkwardness in this
!xx     routine (such as the existence of do loop 60).
!.......................................................................

!.......................................................................
!     curnorm will be the integral over the source over velocity and
!     configuration space volume. To get the actual current, the
!     source arrays which are being computed below will have to
!     be divided by curnorm (giving total current=1.) and then
!     multiplied by curdep, curdep being the total number of
!     particles deposited per second in the plasma.
!     asor will (after rescaling at the end of the routine) have in it
!     the local current (part/sec) for the flux surface in question.
!.......................................................................

      curnorm=0.
      do 200 ll=1,lrz
        call tdnflxs(ll)
        llbrth(ll)=0
!BH171014        call bcast(source(0,0,kfrsou,ll),zero,iyjx2)
        call bcast(source(0:iy+1,0:jx+1,k,ll),zero,iyjx2)

!..................................................................
!     Begin loop over the birth particles.
!     Determine the value of psi and R at birth
!..................................................................

        k=kfrsou
        do 100 ipar=1,ipts

          if (lbrth(ipar).ne.ll) go to 100
          llbrth(ll)=llbrth(ll)+1
          rbrth=pr(ipar)
          rbrth2= rbrth*rbrth
          psibrth=psipt(ipar)
!          write(*,*)'freyasou, do 100, ll,ipar,ipts=',ll,ipar,ipts
!          write(*,*)'freyasou, do 100, rbrth,psibrth=',rbrth,psibrth

!..................................................................
!     Determine the pitch angle (cos(theta)) at birth. Compute
!     B dot v / mag B/ mag v
!..................................................................

          dpsidr=terp2(rbrth,pz(ipar),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
            epsirz,nnra,1,0)
          dpsidz=terp2(rbrth,pz(ipar),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
            epsirz,nnra,0,1)

!..................................................................
!     Determine f(psi)
!..................................................................

          call eqfpsi(psibrth,fpsi_,fppsi_)

!..................................................................
!     Determine components of B.
!..................................................................
          if(rbrth.gt.1.d-8)then
!BHandYuP110715: Account for bsign/cursign

          bxcomp=(-cursign*dpsidz*px(ipar)-bsign*fpsi_*py(ipar))/rbrth2
          bycomp=(-cursign*dpsidz*py(ipar)+bsign*fpsi_*px(ipar))/rbrth2
          bzcomp=+cursign*dpsidr/rbrth
          else ! in a mirror machine, can get to R=0
             go to 100 ! skip, for now
          endif

!.......................................................................
!     Compute the dot product and the cosine
!     BH110614:  But, pitch angle is measure from B-field,
!     BH110614:  however, it is pos. in pos. tor. dirn, by convention.
!     BH110614:  In present NSTX study, bsign=-1., but get pos cospitch
!     BH110614:  for part with pos vphi.  Maybe ok, by logic needs
!     BH110614:  verifying.  Also, effect of cursign?
!.......................................................................

          if (multiply.eq."disabled") then
            jpar=ipar
          else
            jpar=(ipar-1)/multiplyn+1
          endif
          dotprd=bxcomp*vx(jpar)+bycomp*vy(jpar)+bzcomp*vz(jpar)
          vmag=sqrt(vx(jpar)**2+vy(jpar)**2+vz(jpar)**2)
          bmag2= bxcomp**2+bycomp**2+bzcomp**2 ! B^2
          bmag=sqrt(bmag2)
!BHandYuP110715: cospitch measured with vpar pos in pos vphi
!BHandYuP110715: (CCW dirn viewed from top).
          cospitch=bsign*dotprd/bmag/vmag ! local cos(pitch)
!          write(*,*)'freyasou: ll,ipar,vmag,bmag,cospitch=',
!     +                         ll,ipar,vmag,bmag,cospitch

!..................................................................
!     Translate to midplane coordinates. This requires knowledge
!     of B-ratio = B(at birth)/B(at midplane). We require rcon,
!     the value of R at which the flux surface of birth pierces
!     the midplane.
!     Determine the first mesh point iless such the epsi(iless,nzc)
!     is less than psibrth (the epsi associated with the contour of
!     interest).
!     Two methods can be used; default 1-d splines (meth1)
!     Other method is probably much slower (meth2).
!..................................................................

          if (method.eq."meth2") then   !NOT TRUE, method hardwired.
            do 110 i=imag+1,nnr
              if (epsi(i,nzc).lt. psibrth) then
                iless=i
                go to 111
              elseif (i.eq.nnr) then
                call frwrong(3)
              endif
 110        continue
 111        continue

!..................................................................
!     Do a Newton's iteration to find the R=rcon where the contour
!     intersects the Z=0. axis.
!..................................................................

            iter=0
            imore=iless-1
            rcon=er(iless)
            epsitst=epsi(iless,nzc)
 115        continue
            dpsidr=terp2(rcon,zero,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
              epsirz,nnra,1,0)
            a=(psibrth-epsitst)/dpsidr

!..................................................................
!     New value for rcon..
!..................................................................

            rcon=rcon+a
            epsitst=terp2(rcon,zero,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
              epsirz,nnra,0,0)
            err=abs((psibrth-epsitst)/(epsi(imore,nzc)+epsi(iless,nzc)))
            if (err .lt. 1.e-3) go to 120
            iter=iter+1
            if (iter.gt.5 .and. err.gt..01) call frwrong(4)
            if (iter.gt.5) go to 120
            go to 115
 120        continue


!..................................................................
!     Determine B0, i.e., B at midplane
!..................................................................

!BH990429            dpsidr=terp2(rcon,rbrth,0.,nnr,er,
!BH990429                         nnz,ez,epsi,epsirr,epsizz,
            dpsidr=terp2(rcon,zero,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
              epsirz,nnra,1,0)
            bzero=sqrt(fpsi_**2+dpsidr**2)/rcon

          elseif (method.eq."meth1") then   !TRUE

            psb=-psibrth
            call terp1(lrzmax,equilpsp(1),bmdplne,d2bmdpl,psb, &
              1,tab,itab)
            bzero=tab(1)

          endif

          bratio=bmag/bzero
!         bratio, physically, is .ge.1.  Due to inaccuracies, it can
!         numerically get .lt.1.
!         Allowing for some tolerance, these cases are adjusted.
          if (bratio.lt.1.) then
            bratio=1.0
            if (bratio.lt.0.99) &
            write(*,*)'freyasou: bratio=[should be=1, else stop]',bratio
!BH041215            if (bratio.le.0.9) stop 'freyasou: bratio.le.0.9'
          endif
          tmdplne=asin(sqrt((1.-cospitch**2)/bratio))
          if (cospitch.lt.0.) then
            tmdplne=pi-tmdplne
          endif

!..................................................................
!     Determine the indices of the midplane meshes to which the
!     velocity coordinates at the birth point transform.
!     i for theta, j for velocity.
!     Locate the flux surface to which this particle will be
!     assigned.
!..................................................................

          j=luf_bin(vmag/vnorm,xmidpt(1:jx)) !jx
          if (j.gt.jx) go to 100
          i=luf_bin(tmdplne,ymid(1:iy,l_)) !iy
          if (i.gt.iy) then
             write(*,*)'freyasou: ll,i,j,iy,ipar,tmdplne=', &
                  ll,i,j,iy,ipar,tmdplne
             call frwrong(5)
          endif

!..................................................................
!     Increment the source by the birth particle weight.
!     2. in denominator added because zmaxpsi(lr_), vptb are half orbit
!     quantities in CQL
!..................................................................

!VickieLynch060824 found this should be before do 100:          k=kfrsou
          aa=bmdplne(lr_)/(cynt2(i,l_)*cint2(j) &
            *dpsi(lr_)*twopi*vptb(i,lr_)*one_*symm)
          !YuP[03/27/2015] Changed 2.0->symm,
          !which is symm=2 in case of up-dn symmetrical equilibrium.
          !In this case,
          !for passing particles vptb*symm corresponds to the whole orbit;
          !for trapped: vptb*symm corresponds to half orbit (from tip to tip)

!..................................................................
!     Symmetrize in trapped region...
!..................................................................

          if (i.ge.itl_(l_) .and. i.le.itu_(l_)) then   !trapped
            source(i,j,k,ll)=source(i,j,k,ll)+aa/2. ! ZOW
            ! Factor 1/2 because vptb for trapped - over 1/4 physical orbit
            !(in case of up-dn symmetrical equil)
            source(iy+1-i,j,k,ll)=source(i,j,k,ll)  ! ZOW
          else
            source(i,j,k,ll)=source(i,j,k,ll)+aa    ! ZOW
          endif
 100   continue ! ipar=1,ipts
!=========================================== Done LOOP IN PARTICLES


        s=0.
        do 150 j=1,jx
          do 160 i=1,iy
            s=s+source(i,j,k,ll)*cynt2(i,l_)*cint2(j)*vptb(i,lr_)
 160      continue
 150    continue

!..................................................................
!     asor(k,1,lr_) will temporarily hold the non-normalized current
!..................................................................

        asor(k,1,lr_)=s/zmaxpsi(lr_)*one_
        asorz(k,1,lr_)=asor(k,1,lr_)
        curnorm=curnorm+dvol(lr_)*asor(k,1,lr_)
        call tdtoaray
 200   continue ! ll

      deallocate(px,STAT=istat)
      deallocate(py,STAT=istat)
      deallocate(pz,STAT=istat)
      deallocate(pr,STAT=istat)
      deallocate(lbrth,STAT=istat)
      deallocate(psipt,STAT=istat)
      if (multiply.ne."disabled") then
         deallocate(ranvc,STAT=istat)
      endif
      tr3=one !YuP[2019-06-08]was call bcast(tr3(1),one,lrz)

!..................................................................
!     Smooth the source if input variable smooth .ge. .001
!..................................................................

      call frsmooth(k,curnorm)

!..................................................................
!     Compute the factor necessary to scale asor so that the source
!     current is equal to curdep (a quantity determined by NFREYA).
!     If frsmooth has been utilized, tr3 has been altered as has
!     curnorm.
!..................................................................

      scalfact=curdep/curnorm

!..................................................................
!     Scale all the sources (at each flux surface) by scalfact*tr3.
!     This forces the source current equal to curdep.
!..................................................................

      do 300 ll=1,lrz
        call tdnflxs(ll)
        curbrth(ll)=DBLE(llbrth(ll))/dvol(lr_)*scalfact
        call dscal(iyjx2,scalfact*tr3(ll),source(0:iy+1,0:jx+1,k,ll),1)
        asor(k,1,lr_)=asorz(k,1,lr_)*scalfact*tr3(ll)
        xlncur(k,lr_)=asor(k,1,lr_)*zmaxpsi(lr_)
        xlncurt(lr_)=xlncur(k,lr_)
        call sourcpwr(k)
        call tdtoaray
 300  continue

      return
      end subroutine freyasou


!======================================================================
!======================================================================


      subroutine frnnoa(mrans,ranvc)
      use iso_c_binding, only : c_double
      implicit none
!     Generate a vector of mrans normal (0,1) pseudo-random numbers
!     using a method based on the central limit theorem
!     and power residue method.  Output becomes truly normal
!     as k (below) goes to infinity.
!     General formula for the deviate is
!     y=((sum of x(i),i=1,k)-k/2.)/sqrt(k/12.),
!     where x(i) are uniformally distributed on 0,1.
!     Method is borrowed through ibm.  They use k=12.
!
!     parameter (k=12,rtkd12=1.00)
!
      integer mrans ! arg. in subr.
      real(c_double) :: ranvc(mrans) ! arg. in subr.
      real(c_double) :: RANDOM_my ! external function
      integer i,j,k ! local
      real(c_double) :: seed,rtkd12,a ! local

!YuP110316      integer iflag

!YuP110316      iflag=0
      seed=0.d0  ! I.E., continue RN sequence in RANDOM_my

!990131      iseed=123457
!990131      call  ranset(iseed)
!     drand() is a unix library function (see Absoft f77 manual).
!     drand(0) returns next real*8 random number of the
!     sequence.  The input parameter must be integer.
!     Presumably each realization of the
!     sequence is the same.  Rand(iflag.ne.0) could be
!     used to start a different sequence.

!BH060821:
!     random_number(harvest) is [0.,1.] r.n., an f90 intrinsic sub
!     call random_seed(size,put,get) see f90 intrinsics
!     call random_seed() initializes the rn generator.

!      call random_seed()

!BH060821:  BUT, instead of random_number will use RANDOM_my,
!           which is from the GA portlib.f (see function in zfreya.f)
      k=12
      rtkd12=1.0
      do 100 i=1,mrans
        a=0.0
        do 10 j=1,12
!990131          a=a+ranf()
!BH060821          a=a+drand(iflag)
!BH060821           call random_number(harvest)
!BH060821           a=a+harvest

!YuP110316           a=a+RANDOM_my(iflag)
           a=a+RANDOM_my(seed)  !input should be real*8
!-YuP           a=a+drand(iflag) ! YuP: unresolved drand
 10     continue
        ranvc(i)=(a-0.5*k)/rtkd12
 100  continue
      return
      end subroutine frnnoa
