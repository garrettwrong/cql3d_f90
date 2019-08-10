module tdsxr0_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use eqfpsi_mod, only : eqfpsi
  use r8subs_mod, only: dcopy
  use tdsetsxr_mod, only : tdsetsxr
  use tdsxr_mod, only : tdsxr
  use tdsxrplt_mod, only : tdsxrplt
  use tdsxrplt_mod, only : tdsxrvw
  use zcunix_mod, only : terp2

  !---END USE

!
!

contains

  subroutine tdsxr0(rb,ene,icall,iplotsxr)
    use cqlconf_mod, only : setup0
    use param_mod
    use cqlcomm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)


!.......................................................................
!     This subroutine calculates and plots a SXR energy spectrum.
!
!     Toroidal geometry with circular or non-circular flux surfaces
!     is assumed (depending on eqmod namelist variable).
!     The detector is assumed to lie in the x-z plane at y=0.

!     The radial flux surface bins are traversed by the viewing cord up
!     to four times as in the case where the viewing cord is incident
!     from the outer equatorial plane, grazes the inside wall of the
!     tokamak, and continues on to the outside equatorial plane.
!
!     For efficiency, contributions to the SXR spectrum are computed
!     only once for each radial bin (defined by input rb) crossed
!     by the viewing cord.
!     The first task in the computation is obtaining the lengths of
!     the viewing chord in each radial bin, the weighted direction
!     of the viewing cord w.r.t. the B-field, and the weighted poloidal
!     angle on the flux surface.
!     These quantities are indexed by (bin-index,crossing-index 1:4)
!     [ibin1,ibin2 as below], and are then used to obtain the
!     contributions to the SXR flux.
!
!
!
!     rb=array of minor radius bin boundaries from 0. to radmin, in cms.
!        It is only used in the circular plasma case, and/or for
!        determination of the step size along the sight line for
!        both circ and eqdsk equilibria.
!     Dimension is 0:setup0%lrzmax.
!     zeff=array of effective charge, dimension setup0%lrzmax.
!     nv = number of viewing chords.
!     x_sxr(1:nv) = x-coordinate of detector
!     z_sxr(1:nv) = z-coordinate of detector
!
!     if x_sxr(1) and z_sxr(1) are both zero, then rd,thetd are used instead
!
!     rd(1:nv) = minor radius of detector (cms.),  .ge.radmin.
!     thetd(1:nv) = poloidal location of detector (radians).
!     thet1(1:nv) =angle in the poloidal plane including the
!          antenna, measured from the z-axis (degs).
!     thet2(1:nv) = toroidal angle of the viewing coord, measured
!           from the x-z plane (y=0.) in right-hand-sense with
!           respect to the z-axis direction (degs).
!     We have detector direction vector x,y,z-components:
!           detec_x=sin(thet1)*cos(thet2)
!           detec_y=sin(thet1)*sin(thet2)
!           detec_z=cos(thet2)
!     enmin,enmax= are minimum and maximum detector energies (keV).
!       (same for all chords).
!     fds= specifies the step size along the viewing cord, as a fraction
!     of the average radial width of each radial bin (0.2 is a
!     reasonable value).
!
!     icall=  "first" means set up the time-independent and viewing
!     chord-independent piece of the sxr diagnostic, otherwise skip it.
!     icall will be reset to "notfrst", upon exit from this subroutine.
!
!     Ionic charge in included in the cross-section calculation using
!     the input variable zeff(lr_).
!.......................................................................


#ifdef __MPI
      include 'mpilib.h'
#endif

      character*8 icall,iplotsxr,icalls
      dimension rb(0:setup0%lrzmax),ene(setup0%lrzmax)
      dimension efluxwk(nen,setup0%lrzmax) ! local working array
      real(c_double), allocatable, dimension(:):: tempp4,tempp5,tempp6

      if(nv.gt.nva) stop 'nv.gt.nva in tdsxr0'

!.......................................................................
!     Allocate temporary storage, to be deallocated at end of subroutine
!.......................................................................
!BH121115      itempp=int(25*setup0%lrzmax*nv/fds)
      itempp=int(100*setup0%lrzmax*nv/fds)
      allocate(tempp4(itempp),STAT=istat4)
      allocate(tempp5(itempp),STAT=istat5)
      allocate(tempp6(itempp),STAT=istat6)
      if (istat4.ne.0 .or. istat5.ne.0 .or. istat6.ne.0) then
         write(*,*)'tdsxr0, tempp1: allocation problem'
         STOP
      endif
!.......................................................................



      icalls=icall
      iistep=0

!.......................................................................
!*********loop over viewing angles      **************************
!.......................................................................
      do 420  nn=1,nv

!.......................................................................
!     Start source at detector position and calculate viewing cosines
!.......................................................................

       xs_(2)=0.0
       if (x_sxr(1).eq.zero .and. z_sxr(1).eq.zero) then
          xs_(1)=radmaj+rd(nn)*cos(thetd(nn)*pio180)
          xs_(3)=rd(nn)*sin(thetd(nn)*pio180)
       else
          xs_(1)=x_sxr(nn)
          xs_(3)=z_sxr(nn)
       endif
!.......................................................................
!     Calculate viewing cosines.
!.......................................................................

        alphad(1)=sin(thet1(nn)*pio180)*cos(thet2(nn)*pio180)
        alphad(2)=sin(thet1(nn)*pio180)*sin(thet2(nn)*pio180)
        alphad(3)=cos(thet1(nn)*pio180)

!.......................................................................
!     Stepping along the sightline obtaining distance in each radial
!     bin, the weighted angle with respect to the b-field, and the
!     average polidal position in the radial bin....
!.......................................................................

        if(fds.ge.1.)  stop 'fds step size in sxr too large'
        ibin1=setup0%lrzmax
        ibin2=1

!.......................................................................
!     ibin1= the radial bin of the source point along the viewing cord,
!     numbered from 1 for the bin with inner edge at the magnetic
!     axis, to setup0%lrzmax for the bin with outer limit equal to radmin.
!     (Bin ibin occupies the region from .ge.rb(ibin-1) to
!     .lt.rb(ibin)).
!     Emission will only be calculated for each crossing of a radial
!     bin.  Need to keep track of several possible crossings per bin.
!     ibin2= refers to successive classes of radial motion along the
!     viewing cord.
!     = 1, then this is first pass from larger to smaller minor radius.
!       2, minor radius is increasing as the viewing cord is transitted
!       3, viewing cord has not passed out of the plasma, and the
!          minor radius is decreasing.
!       4, minor radius increasing again.
!.......................................................................

        s=0.0
        istart=1
        istep=0
        do 99  j=1,4
          do 991  i=1,setup0%lrzmax
            ibin(i,j)=0
            sxry(i,j)=0.0
            sang(i,j)=0.0
            spol(i,j)=0.0
 991      continue
 99     continue

!.......................................................................
!     main loop for determination of view chord increments
!.......................................................................

 100    continue

!990131        ds=fds*amin1(rb(ibin1)-rb(ibin1-1),radmin/setup0%lrzmax)
!BH081106        ds=fds*min(rb(ibin1)-rb(ibin1-1),radmin/setup0%lrzmax)
        ds=fds*min(rb(ibin1)-rb(ibin1-1),(rb(ibin1)-rb(0))/setup0%lrzmax)
        s=s+ds
        istep=istep+1  ! counts for each view chord
        do 102  j=1,3
 102    xs_(j)=xs_(j)+alphad(j)*ds

!.......................................................................
!     Save view cords for plotting at end of subroutine
!.......................................................................

        iistep=iistep+1  ! counts up over all view chords

!        write(*,*)'tdsxr0: nn,istep,iistep,ibin1,ds=',
!     &       nn,istep,iistep,ibin1,ds

        if (iistep.le.itempp) then
           lensxr(nn)=istep
           tempp4(iistep)=xs_(1)
           tempp5(iistep)=xs_(2)
           tempp6(iistep)=xs_(3)
        else
           stop 'tdsxr0: TOO MANY STEPS to save in SXR cords'
        endif

!.......................................................................
!     Toroidal angle of source point
!.......................................................................

        phis=atan2(xs_(2),xs_(1))



!.......................................................................
!     Beginning of if-statement for circular/noncircular
!     Non-circular (eqmod.eq."enabled") case starts at line 349.
!.......................................................................

        if (eqmod.ne."enabled")  then

!.......................................................................
!     Minor radius of source point
!.......................................................................

          rs=sqrt((xs_(1)-radmaj*cos(phis))**2 &
            +(xs_(2)-radmaj*sin(phis))**2 +xs_(3)**2)

!.......................................................................
!     Continue stepping inwards until the plasma is entered,
!     or more than 1000/fds steps are taken.
!.......................................................................

!BH180207          if(istep.gt.(1000/fds).and.istart.eq.1)  then
!BH180207: Needed a few more steps for a difficult case.
          if(istep.gt.(2000/fds).and.istart.eq.1)  then
             write(*,*) 'Warning: SXR sightline',nn, &
                        " missed or didn't reach plasma"
             go to 420
          endif
          if(rs.ge.radmin.and.istart.eq.1)  go to 100
          istart=0

!.......................................................................
!     Exit from this loop if the source point leaves the plasma.
!.......................................................................

          if(rs.ge.radmin)  go to 200

!.......................................................................
!     Decide which radial bin the source point
!     is in:
!.......................................................................

          if(ibin2.eq.2.or.ibin2.eq.4)  go to 130

!.......................................................................
!     ibin2 is equal to 1 or 3. Thus at previous step the bin disignator
!     was in a decreasing mode.
!     Check if in the same bin as previous step.
!.......................................................................

          rs1=rs-rb(ibin1)
          rs2=rs-rb(ibin1-1)
          if(rs1*rs2.le.0..and.rs1.ne.zero)  go to 150

!.......................................................................
!     Have changed bins.
!.......................................................................

          if(rs1.ge.0.)  go to 115
!.......................................................................
!     Continuing in decreasing mode.
!.......................................................................
          ibin1=ibin1-1
          if(ibin1.lt.1)  stop 'stop 1 in tdsxr0'
          go to 150
!..................................................................
!     Have changed directions in minor radius
!..................................................................

 115      ibin2=ibin2+1
          ibin1=ibin1+1
          go to 150

!..................................................................
!     ibin2 is equal to 2 or 4.  At previous step bin designator was in
!     increasing mode.
!     Check if source is in same bin as previous step.
!..................................................................

 130      rs1=rs-rb(ibin1-1)
          rs2=rs-rb(ibin1)
          if(rs1*rs2.le.0..and.rs2.ne.zero)  go to 150

!.......................................................................
!     Have changed bins
!.......................................................................

          if(rs1.lt.0.0)  go to 135

!.......................................................................
!     Continuing on is same direction
!.......................................................................

          ibin1=ibin1+1
          if(ibin1.gt.setup0%lrzmax)  stop 'stop 2 in tdsxr0'
          go to 150

!.......................................................................
!     Have changed directions
!.......................................................................

 135      ibin2=ibin2+1
          ibin1=ibin1-1

!.......................................................................
!     Add in ds, and increment weighted angles.
!.......................................................................
 150      if(ibin2.gt.4)  stop 'stop 3 in tdsxr0'
          sxry(ibin1,ibin2)=sxry(ibin1,ibin2)+ds
!.......................................................................
!     Angle between viewing dirn.and toroidal (approx. b-field) direction
!     (Detector views sxr's emitted towards it).
!.......................................................................

          angle=acos(-sin(phis)*alphad(1)+cos(phis)*alphad(2))
          sang(ibin1,ibin2)=sang(ibin1,ibin2)+ds*angle

!.......................................................................
!     Poloidal angle
!.......................................................................

          arg1=xs_(3)
          arg2=sqrt(xs_(1)**2+xs_(2)**2)-radmaj
          if(arg1.eq.zero .and. arg2.eq.zero)  then
            polang=0.
          else
            polang=atan2(arg1,arg2)
          endif

!.......................................................................
!     assume up-down symmetry
!.......................................................................

          polang=abs(polang)
          if(polang.gt.pi)  polang=pi
          spol(ibin1,ibin2)=spol(ibin1,ibin2)+ds*polang

!.......................................................................








!.......................................................................
!     If non-circular:
!.......................................................................

        else


          if ((psimag-psilim).le.0.0) stop 'psimag.lt.psilim in tdsxr0'
          tr(0)=zero
          do 160  l=1,setup0%lrzmax
 160      tr(l)=psimag-psivalm(l)  !Increases from 0. to psimag-psilim at edge.
          xx=xs_(1)
          yy=xs_(2)
          rr=sqrt(xx**2+yy**2)
          zz=xs_(3)
!.......................................................................
!     Check if on eqdsk grid
!.......................................................................

          iongrid=0
          if ((rr.lt.er(nnr).and. rr.gt.er(1)) .and. &
               (zz.lt.ez(nnz).and.zz.gt.ez(1))) iongrid=1
          if (iongrid.eq.0) go to 165
!.......................................................................
!     Minor radius of source point
!.......................................................................

          ppsi=terp2(rr,zz,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
            epsirz,nnra,0,0)
          apsi=psimag-ppsi     !Increases from 0. to psimag at edge.
          apsi=max(apsi,zero)  !Sometimes there might be a small
                               !error in psimag causing a problem
                               !near the magnetic axis [BH, 121115].

!.......................................................................
!     Continue stepping inwards until the plasma is entered,
!     or more than 1000/fds steps are taken.
!.......................................................................

!BH180207     165   if(istep.gt.(1000/fds).and.istart.eq.1)  then
!BH180207: Needed a few more steps for a difficult case.
 165      if(istep.gt.(2000/fds).and.istart.eq.1)  then
             write(*,*) 'Warning: SXR sightline',nn, &
                        " missed or didn't reach plasma"

             go to 420
          endif
          if(((iongrid.eq.0).or.(ppsi.le.psilim)).and. &
                                              istart.eq.1)  go to 100

          istart=0  ! in the plasma

!.......................................................................
!     Exit from this loop if the source point leaves the plasma.
!.......................................................................

          if(ppsi.le.psilim)  go to 200

!.......................................................................
!     Decide which radial bin the source point
!     is in:
!.......................................................................

          if(ibin2.eq.2.or.ibin2.eq.4)  go to 180

!.......................................................................
!     ibin2 is equal to 1 or 3. Thus at previous step the bin designator
!     was in a decreasing mode.
!     Check if in the same bin as previous step.
!.......................................................................

          rs1=apsi-tr(ibin1)     !  ibin1 is initialized = setup0%lrzmax
          rs2=apsi-tr(ibin1-1)
          if(rs1*rs2.le.zero.and.rs1.ne.zero) go to 190 !in [ibin1,ibin1+1)

!.......................................................................
!     Have changed bins.
!.......................................................................

          if(rs1.ge.zero)  go to 170  ! Changed to increasing mode
!.......................................................................
!     Continuing in decreasing mode.
!.......................................................................
          ibin1=ibin1-1
          if(ibin1.lt.1) then
             write(*,*) 'stop 4 in tdsxr0'
             stop
          endif
          go to 190
!.......................................................................
!     Have changed directions in minor radius
!.......................................................................

 170      ibin2=ibin2+1
          ibin1=ibin1+1
          go to 190

!.......................................................................
!     ibin2 is equal to 2 or 4.  At previous step bin designator was in
!     increasing mode.
!     Check if source is in same bin as previous step.
!.......................................................................

 180      rs1=apsi-tr(ibin1-1)
          rs2=apsi-tr(ibin1)
          if(rs1*rs2.le.zero .and. rs2.ne.zero)  go to 190

!.......................................................................
!     Have changed bins
!.......................................................................

          if(rs1.lt.zero)  go to 185

!.......................................................................
!     Continuing on is same direction
!.......................................................................

          ibin1=ibin1+1
          if(ibin1.gt.setup0%lrzmax)  stop 'stop 5 in tdsxr0'
          go to 190

!.......................................................................
!     Have changed directions
!.......................................................................

 185      ibin2=ibin2+1
          ibin1=ibin1-1

!.......................................................................
!     Add in ds, and increment weighted angles.
!.......................................................................
 190      if(ibin2.gt.4)  stop 'stop 6 in tdsxr0'
          sxry(ibin1,ibin2)=sxry(ibin1,ibin2)+ds
!.......................................................................
!     Angle between (negative) viewing dirn, and b-field direction,
!     (considering b-field oriented with positive toroidal component).
!     (Detector views sxr's emitted towards it).
!.......................................................................

!.......................................................................
!     Determine the pitch angle (cos(theta)) at birth. Compute
!     B dot v / mag B/ mag v
!.......................................................................

          dpsidr=terp2(rr,zz,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
            epsirz,nnra,1,0)
          dpsidz=terp2(rr,zz,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
            epsirz,nnra,0,1)

!.......................................................................
!     Determine f(psi)
!.......................................................................

          call eqfpsi(ppsi,fpsi_,fppsi_)
          fpsi_=abs(fpsi_)

!.......................................................................
!     Determine components of B.
!.......................................................................

          xcomp=(-cursign*dpsidz*xx-bsign*fpsi_*yy)/rr**2
          ycomp=(-cursign*dpsidz*yy+bsign*fpsi_*xx)/rr**2
          zcomp=+cursign*dpsidr/rr

!.......................................................................
!     Compute the dot product and the cosine
!.......................................................................

          dotprd=-(xcomp*alphad(1)+ycomp*alphad(2)+zcomp*alphad(3))
          bbmag=sqrt(xcomp**2+ycomp**2+zcomp**2)
!BHandYuP110715: cospitch measured with vpar pos in pos vphi
!BHandYuP110715: (CCW dirn viewed from top).
          cospitch=bsign*dotprd/bbmag
          angle=acos(cospitch)
          sang(ibin1,ibin2)=sang(ibin1,ibin2)+ds*angle

!.......................................................................
!     Poloidal angle   (rmag is major radius of magnetic axis).
!.......................................................................

          arg1=zz
          arg2=rr-rmag
          if(arg1.eq.zero .and. arg2.eq.zero)  then
            polang=0.
          else
            polang=atan2(arg1,arg2)
          endif

!.......................................................................
!     assume up-down symmetry
!.......................................................................

          polang=abs(polang)
          if(polang.gt.pi)  polang=pi
          spol(ibin1,ibin2)=spol(ibin1,ibin2)+ds*polang

!.......................................................................
!     end of circular/non-circular if-statement
!.......................................................................

        endif

!.......................................................................
!     Radial bin.
!.......................................................................

        ibin(ibin1,ibin2)=ibin1
        go to 100

!.......................................................................
!     *******     Have finished sightline determination   ******
!.......................................................................

 200    continue

!.......................................................................
!     Obtain weighted angle between viewing and toroidal directions,
!     and poloidal angle, in each or the setup0%lrzmax radial bins.
!.......................................................................

        do 210  j=1,ibin2
          do 211  i=1,setup0%lrzmax
            if(sxry(i,j).eq.zero)  go to 211
            sang(i,j)=sang(i,j)/sxry(i,j)
            spol(i,j)=spol(i,j)/sxry(i,j)
            if (spol(i,j) .gt. pi) spol(i,j)=pi
 211      continue
 210    continue

!.......................................................................
!     Integrations over velocity space for several detector energies.
!     Set up energies linearly from enmin to enmax.
!     (Check nen.le.nena)
!.......................................................................

        if(nen.le.nena)  go to 299
        write (*,1000) nen, nena
 1000   format("nen = ",i3, "is too large. Reset automatically to ",i5)
        nen=nena
 299    den=(enmax-enmin)/(nen-1.)

        do 300  ien=1,nen
          en_(ien)=enmin+(ien-1)*den

!.......................................................................
!     kev
!     Units of m0*c**2
!.......................................................................

          enk(ien)=en_(ien)/restmkev  !YuP[01-2011] was 512.

!.......................................................................
!     Legendre decomposition of XR cross sections due to ei  and ee
!     collisions.  Done only once for each energy.
!     Results stored in sigsxr(1:jx,1:msxr,1:nen,2) by sub tdsetsxr.
!.......................................................................

          if (icall .eq. "first") call tdsetsxr(ien)
 300    continue ! ien=1,nen

!.......................................................................
!     Obtain energy flux for each energy and sightline
!.......................................................................

        efluxwk=zero !YuP[2019-06-08]was call bcast(efluxwk(1,1),zero,nen*lrzmax)
        mpisz= nen ! number of elements in efluxwk(nen,i)

        do i=1,setup0%lrzmax
#ifdef __MPI
         if(mpisize.gt.1) then
            mpiworker= MOD(i-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
#endif
#ifdef __MPI
      if(mpirank.eq.mpiworker) then

#endif
        do j=1,ibin2
        if(sxry(i,j).ne.zero)  then
          do ien=1,nen !should be the innermost: setting plegg at ien=1
              call tdsxr(ien,sang(i,j),spol(i,j),ibin(i,j), &
                         dedotei,dedotee) !-> out
              efluxwk(ien,i)= efluxwk(ien,i)+ &
                1.6e-12*enk(ien)*sxry(i,j)*ene(ibin(i,j))* &
                (zeff(ibin(i,j))*dedotei +dedotee)
              ! efluxwk(ien,i): summation over j=1,ibin2 only
          enddo ! ien=1,nen
        endif
        enddo ! j=1,ibin2
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif
#ifdef __MPI
      if(mpirank.eq.0) then !-------------------------------------------
        call MPI_RECV(tem2, mpisz,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,mpiierr)
        mpitag=mpistatus(MPI_TAG)
        mpil_=mpitag ! determine which radial surface sent the data
        call dcopy(mpisz,tem2(1:mpisz),1,efluxwk(1:mpisz,mpil_),1)
      endif !-----------------------------------------------------------

#endif
#ifdef __MPI
      if(mpirank.eq.mpiworker) then !-----------------------------------
        call dcopy(mpisz,efluxwk(1:mpisz,i),1,tem2(1:mpisz),1)
        mpitag= i ! i=1,setup0%lrzmax
        call MPI_SEND(tem2, mpisz,MPI_DOUBLE_PRECISION,0,mpitag,MPI_COMM_WORLD,mpiierr)
      endif !-----------------------------------------------------------
#endif
        enddo ! i=1,setup0%lrzmax

        call bcast(eflux(1:nen,nn),zero,nen) !YuP[2019-06-08] range 1:nen

        do ien=1,nen  ! for each energy bin
        do i=1,setup0%lrzmax ! sum-up contributions from all flux surfaces
          eflux(ien,nn)= eflux(ien,nn) + efluxwk(ien,i)
        enddo ! i=1,setup0%lrzmax
        enddo ! ien=1,nen

#ifdef __MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
#endif
#ifdef __MPI
      call MPI_BCAST(eflux(1:nena,nn),nena,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
#endif

!.......................................................................
!     Integrate energy flux spectra from enmin to enmax.
!     Assume the flux varies as exp(-e/t) between flux values,
!     to improve accuracy.
!.......................................................................

        inegsxr(nn)=nen
        efluxt(nn)=0.0

        do 400  ien=1,nen
          if (eflux(ien,nn).le.0.) then
            inegsxr(nn)=ien-1
            go to 402
          endif
          if (ien.eq.nen) go to 402

!         Avoid some overflows:
          if (eflux(ien+1,nn).lt.0.) go to 402
          if (eflux(ien,nn).lt.em100*em100) go to 400
!BH160507     if(eflux(ien,nn).eq.eflux(ien+1,nn))  go to 401
!BH160507 Allowing for possible increase in eflux_npa with en.
!BH160507 Could modify to fork with eflux_npa(i:i+1) are too close
          if(eflux(ien,nn).le.eflux(ien+1,nn))  go to 401
          if(eflux(ien+1,nn).eq.zero)  go to 401

!990131          t=(en_(ien+1)-en_(ien))/alog(eflux(ien,nn)/eflux(ien+1,nn))
          t=(en_(ien+1)-en_(ien))/log(eflux(ien,nn)/eflux(ien+1,nn))
!BH160507          c=eflux(ien,nn)*exp(en_(ien)/t)
!BH160507          write(*,*)'tdnpa0.f: i,t,c,eflux_npa(i,nn)=',
!BH160507     +                         i,t,c,eflux_npa(i,nn)
!BH160507          efluxt(nn)=efluxt(nn)+c*t*
!BH160507     +      (exp(-en_(ien)/t)-exp(-en_(ien+1)/t))*1.e3
          efluxt(nn)=efluxt(nn)+t*(eflux(ien,nn)-eflux(ien+1,nn))*1.e3
          go to 400
 401      efluxt(nn)=efluxt(nn)+eflux(ien,nn)*(en_(ien+1)-en_(ien))*1.e3
 400    continue !   ien=1,nen
 402    continue
!
        icall="notfirst"

!.......................................................................
!*********END loop over viewing angles      **************************
!.......................................................................

 420  continue ! nn=1,nv

!.......................................................................
!     Plot eflux........
!.......................................................................

      nenaa=nena
      if (iplotsxr.eq.'yes') &
          call tdsxrplt(en_,eflux,nen,nenaa, &
                  efluxt,nv,inegsxr,softxry,setup0%lnwidth)

!.......................................................................
!     Plot SXR view cords in poloidal cross-section, if eqmod="enabled"
!.......................................................................

      if (eqmod.eq."enabled" .and. icalls.eq."first" .and. &
           iplotsxr.eq.'yes') then
         call tdsxrvw(tempp4,tempp5,tempp6)
      endif

      deallocate(tempp4,tempp5,tempp6,STAT=istat1)

      return
      end subroutine tdsxr0


end module tdsxr0_mod
