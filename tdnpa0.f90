module tdnpa0_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use eqfpsi_mod, only : eqfpsi
  use tdnpa_mod, only : tdnpa
  use tdsetnpa_mod, only : tdsetnpa
  use tdsxrplt_mod, only : tdsxrplt
  use tdsxrplt_mod, only : tdsxrvw

  !---END USE

!
!

contains

      subroutine tdnpa0(rb,ene,icall,iplotnbi)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save
!VT....V.Tang, NPA diagnostic 9-25-05.............................
!VT    Currently, I have only made changes to the non-circular section
!VT    Checks out 5/12/06, Version 1.0, NPA Mod
!VT...............................................................
!BH100425: Modified and simplified treatment:  see BH100423 below.
!BH100425: Results changed substantially.  CX cross-section updated
!BH100425: from IAEA ALLADIN data base, nominally good from 0.12 eV
!BH100425: to 630 keV.
!.......................................................................
!     [Earlier SXR routines modified for NPA. Comments from tdsxr0, etc.]
!     This subroutine calculates and plots a SXR///NPA energy spectrum.
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
!     [BH100423:  The CX cross-section is simple, compared to that
!                 for XRays:  forward scattering only, with no energy
!                 change, i.e., scattering cross-section per energy
!                 and solid angle has delta-functions in energy and
!                 solid angle.
!                 Hence, NPA is a relatively short calculation compared
!                 to SXR, so  the calculation can be performed
!                 as a simple integration at each step along the
!                 viewing cord.]
!
!     rb=array of minor radius bin boundaries from 0. to radmin, in cms.
!     Dimension is 0:lrzmax.
!     zeff=array of effective charge, dimension lrzmax.
!     rd_npa= minor radius of detector (cms.),  .ge.radmin.
!     thetd_npa= poloidal location of detector (radians).
!     nv_npa= number of viewing chords.
!     thet1_npa(1:nv_npa)=angle in the poloidal plane including the
!          antenna, measured from the z-axis (degs).
!     thet2_npa(1:nv_npa)=toroidal angle of the viewing coord, measured
!           from the x-z plane (y=0.) in right-hand-sense with
!           respect to the z-axis direction (degs).
!     We have detector direction vector x,y,z-components:
!           detec_x=sin(thet1_npa)*cos(thet2_npa)
!           detec_y=sin(thet1_npa)*sin(thet2_npa)
!           detec_z=cos(thet1_npa)
!     enmin_npa,enmax_npa= are minimum and maximum detector energies (keV).
!       (same for all chords).
!     fds_npa= specifies the step size along the viewing cord, as a fraction
!     of the average radial width of each radial bin (0.2 is a
!     reasonable value).
!
!     icall=  "first" means set up the time-independent and viewing
!     chord-independent piece of the sxr diagnostic, otherwise skip it.
!     icall will be reset to "notfrst", upon exit from this subroutine.
!     NPA:  only presently called with icall="notfrst" [BH100425].
!     NPA:  References to icall can be removed. [BH100425].
!
!     Ionic charge in included in the cross-section calculation using
!     the input variable zeff(lr_).
!.......................................................................


!MPIINSERT_INCLUDE

      character*8 icall,icalls,iplotnbi

      real*8,dimension(lrzmax)::ene
      real*8,dimension(0:lrzmax)::rb
      real*8,dimension(lrzmax,4,nen_npa)::atten
      real*8,dimension(lrzmax,nen_npa)::stoplamda
      real*8,dimension(nen_npa)::vel
      real*8,dimension(nen_npa,npaproc)::sigmanpa,emitnpa
      real*8,dimension(nen_npa,nv_npa)::eflux_npa
      dimension efluxwk(nen_npa,lrzmax) ! local working array

      integer::jj
      real*8, allocatable, dimension(:):: tempp4,tempp5,tempp6

      print*, 'tdnpa0::Doing NPA diagnostic, CQL at step', n
      eflux_npa(:,:)=0.

      if(nv_npa.gt.nva) stop 'nv_npa.gt.nva in tdnpa0'

      if(nen_npa.gt.nena) then
         write (*,1000) nen_npa, nena
 1000    format("nen_npa = ",i3, "is too large. Reset to ",i5)
         nen_npa=nena
      endif

!.......................................................................
!     Allocate temporary storage, to be deallocated at end of subroutine
!.......................................................................
!BH121115      itempp=int(15*lrzmax*nv_npa/fds)
      itempp=int(100*lrzmax*nv_npa/fds)
      allocate(tempp4(itempp),STAT=istat4)
      allocate(tempp5(itempp),STAT=istat5)
      allocate(tempp6(itempp),STAT=istat6)
      if (istat4.ne.0 .or. istat5.ne.0 .or. istat6.ne.0) then
         write(*,*)'tdsxr0, tempp1: allocation problem'
         STOP
      endif
!.......................................................................


      den=(enmax_npa-enmin_npa)/(nen_npa-1.) !energy bin width(keV)

!BH   Take en_() to be the center of each energy bin
      do 5  ien=1,nen_npa
         en_(ien)=enmin_npa+(ien-1)*den  ! (keV)
 5    continue

      do 223 ien=1,nen_npa
         vel(ien)=(2.*en_(ien)*ergtkev/fmass(1))**0.5  ! For ion#1 !!!
 223  continue

!.......................................................................
!     Obtain CX cross-section(s) versus energy
!.......................................................................

!     if (icall .eq. "first")
      call tdsetnpa(sigmanpa)


!VT..............................................
!VT....call lam subroutine to get the needed attenuation cross
!VT    sections for NPA diagnostic
!VT.................................................................

      call lam(en_,ene,stoplamda)

      icalls=icall
      iistep=0


!.......................................................................
!*********main loop over viewing angles      ***************************
!.......................................................................

      do 420  nn=1,nv_npa

!.......................................................................
!     Start source at detector position and calculate viewing cosines
!.......................................................................

!     if x_npa(1) and z_npa(1) are both zero, then rd_npa,thetd_npa
!     are used instead.
       xs_(2)=0.0d0
       if (x_npa(1).eq.zero .and. z_npa(1).eq.zero) then
          xs_(1)= rmag + rd_npa(nn)*cos(thetd_npa(nn)*pio180)
          xs_(3)= zmag + rd_npa(nn)*sin(thetd_npa(nn)*pio180)
          ! YuP June 2013: In above, changed radmaj -> rmag, added zmag
       else
          xs_(1)=x_npa(nn)
          xs_(3)=z_npa(nn)
       endif
       print*,'tdnpa0: starting x,y,z',xs_(1),xs_(2),xs_(3)

!.......................................................................
!     Calculate viewing cosines.
!.......................................................................

        alphad(1)=sin(thet1_npa(nn)*pio180)*cos(thet2_npa(nn)*pio180)
        alphad(2)=sin(thet1_npa(nn)*pio180)*sin(thet2_npa(nn)*pio180)
        alphad(3)=cos(thet1_npa(nn)*pio180)

!.......................................................................
!     Stepping along the sightline obtaining distance in each radial
!     bin, the weighted angle with respect to the b-field, and the
!     average polidal position in the radial bin....
!.......................................................................

        if(fds_npa.ge.1.)  stop 'tdfds_npa step size in npa too large'
        ibin1=lrzmax
        ibin2=1

!.......................................................................
!     ibin1= the radial bin of the source point along the viewing cord,
!     numbered from 1 for the bin with inner edge at the magnetic
!     axis, to lrzmax for the bin with outer limit equal to radmin/
!     (rhomax for eqmod.eq.'enabled').
!     (Bin ibin occupies the region from .ge.rb(ibin-1) to
!     .lt.rb(ibin)).
!     Emission will only be calculated for each crossing of a radial
!     bin.  Need to keep track of several possible crossings per bin.
!     ibin2= refers to successive classes of radial motion along the
!     viewing cord.
!     = 1, then this is first pass from larger to smaller minor radius.
!     2, minor radius is increasing as the viewing cord is transitted
!     3, viewing cord has not passed out of the plasma, and the
!     minor radius is decreasing.
!     4, minor radius increasing again.
!.......................................................................

        s=0.0
        istart=1
        istep=0

!VT indicators for when we switch ibin2 and ibin1 so we can properly
!VT account for attenuation;
        iclass12=0
        iclass23=0
        ichangebin13=0
        ichangebin24=0
!VT...............................................

        do 99  j=1,4
          do 991  i=1,lrzmax  ! Initialize:
            ibin(i,j)=0
            sxry(i,j)=0.0
            sang(i,j)=0.0
            spol(i,j)=0.0
!VT attenuation factor
            do 992 kk=1,nen_npa
               atten(i,j,kk)=0.0
 992        continue
!VT................................
 991     continue
 99      continue

!.......................................................................
!     main loop for determination of view chord increments
!.......................................................................

 100    continue

        ds=fds_npa*min(rb(ibin1)-rb(ibin1-1),radmin/lrzmax)
        s=s+ds
        istep=istep+1  ! counts for each view chord
        do 102  j=1,3
 102    xs_(j)=xs_(j)+alphad(j)*ds

!.......................................................................
!     Save view cords for plotting at end of subroutine
!.......................................................................

        iistep=iistep+1  ! counts up over all view chords
        if (iistep.le.itempp) then
           lensxr(nn)=istep
           tempp4(iistep)=xs_(1)
           tempp5(iistep)=xs_(2)
           tempp6(iistep)=xs_(3)
        else
           stop 'tdnpa0: TOO MANY STEPS to save in NPA cords'
        endif

!.......................................................................
!     Toroidal angle of source point
!.......................................................................

        phis=atan2(xs_(2),xs_(1))

!.......................................................................
!     Beginning of if-statement for circular/noncircular
!.......................................................................

        if (eqmod.ne."enabled")  then

           write(*,*)
           write(*,*)'tdnpa0:CAUTION: Circular SXR NOT ADJUSTED FOR NPA'
           write(*,*)

!.......................................................................
!     Minor radius of source point
!.......................................................................

          rs=sqrt((xs_(1)-radmaj*cos(phis))**2 &
            +(xs_(2)-radmaj*sin(phis))**2 +xs_(3)**2)

!.......................................................................
!     Continue stepping inwards until the plasma is entered,
!     or more than 1000/fds_npa steps are taken.
!.......................................................................

          if(istep.gt.(1000/fds_npa).and.istart.eq.1)  then
             write(*,*) 'Warning: NPA sightline',nn, &
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
          if(ibin1.lt.1)  stop 'stop 1 in tdnpa0'
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
          if(ibin1.gt.lrzmax)  stop 'stop 2 in tdnpa0'
          go to 150

!.......................................................................
!     Have changed directions
!.......................................................................

 135      ibin2=ibin2+1
          ibin1=ibin1-1

!.......................................................................
!     Add in ds, and increment weighted angles.
!.......................................................................
 150      if(ibin2.gt.4)  stop 'stop 3 in tdnpa0'
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
!     If non-circular:   (Only checked out implementation for NPA.)
!.......................................................................

        else  ! eqmod.ne.disabled


          if ((psimag-psilim).le.0.0) stop 'psimag.lt.psilim in tdnpa0'
          tr(0)=0.0
          do 160  l=1,lrzmax
 160      tr(l)=psimag-psivalm(l)
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
!     or more than 1000/fds_npa steps are taken.
!.......................................................................

 165      if(istep.gt.(1000/fds_npa).and.istart.eq.1)  then
             write(*,*) 'Warning: npa sightline',nn, &
                        " missed or didn't reach plasma"

             go to 420
          endif
          if(((iongrid.eq.0).or.(ppsi.le.psilim)).and. &
                                              istart.eq.1)  go to 100

          istart=0  ! in the plasma, i.e., no longer starting.

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
!     ibin2 is equal to 1 or 3. Thus at previous step the bin disignator
!     was in a decreasing mode.
!     Check if in the same bin as previous step.
!.......................................................................

          rs1=apsi-tr(ibin1)
          rs2=apsi-tr(ibin1-1)
          if(rs1*rs2.le.0..and.rs1.ne.zero)  go to 190

!.......................................................................
!     Have changed bins.
!.......................................................................

          if(rs1.ge.0.)  go to 170
!.......................................................................
!     Continuing in decreasing mode.
!.......................................................................
          ichangebin13=1 !cVT
          ibin1=ibin1-1
          if(ibin1.lt.1)  stop 'stop 4 in tdnpa0'
          go to 190
!.......................................................................
!     Have changed directions in minor radius
!.......................................................................

 170      ibin2=ibin2+1
          ibin1=ibin1+1
          iclass12=1 !cVT
          go to 190

!.......................................................................
!     ibin2 is equal to 2 or 4.  At previous step bin designator was in
!     increasing mode.
!     Check if source is in same bin as previous step.
!.......................................................................

 180      rs1=apsi-tr(ibin1-1)
          rs2=apsi-tr(ibin1)
          if(rs1*rs2.le.0..and.rs2.ne.zero)  go to 190

!.......................................................................
!     Have changed bins
!.......................................................................

          if(rs1.lt.0.0)  go to 185

!.......................................................................
!     Continuing on is same direction
!.......................................................................
          ichangebin24=1 !cVT
          ibin1=ibin1+1
          if(ibin1.gt.lrzmax)  stop 'stop 5 in tdnpa0'
          go to 190

!.......................................................................
!     Have changed directions
!.......................................................................

 185      ibin2=ibin2+1
          ibin1=ibin1-1
          iclass23=1 !cVT

!.......................................................................
!     Add in ds, and increment weighted angles.
!.......................................................................

 190      if(ibin2.gt.4)  stop 'stop 6 in tdnpa0'
          sxry(ibin1,ibin2)=sxry(ibin1,ibin2)+ds ! cumulative length (cms)
                                                 ! in each radial bin.
          do 195 ien=1,nen_npa

!BH          print*,'tdnpa0::Stepping thru energies for atten'
!VT check bin changes and increment attenuation factor depending on case
!VT stoplamda stands for lamda and it's the mean-free path for a
!VT flux surface as a function of neutral energy (n_en)
          if (iclass12.eq.1) then !changed ibin2, class from 1 to 2  or 3 to 4
          atten(ibin1,ibin2,ien)=ds/stoplamda(ibin1,ien) &
           +atten(ibin1-1,ibin2-1,ien)  !going from into plasma to out...
!VT                                     ibin1 flux sufrace now going up
!          print*,'tdnpa0::ibin2 change from 1 to 2 or 3 to 4 in atten'
!          print*,'tdnpa0::Previous class attenuation factor for
!     1     energy bin:',en_m(ien),atten(ibin1-1,ibin2-1,ien)
          elseif (iclass23.eq.1) then !changed ibin2, class from 2 to 3
          atten(ibin1,ibin2,ien)=ds/stoplamda(ibin1,ien) &
            +atten(ibin1+1,ibin2-1,ien) !going from out of plasma to in
!VT                         again...inbin flux sufraces now going down

!          print*,'tdnpa0::ibin2 change from 2 to 3 in atten'
!          print*,'tdnpa0::Previous class attenuation factor for energy
!     1     bin:',en_m(ien),atten(ibin+1,ibin2-1,ien)

          elseif (ichangebin13.eq.1) then !decreasing ibin1, in class 1 or 3
          atten(ibin1,ibin2,ien)=ds/stoplamda(ibin1,ien) &
             +atten(ibin1+1,ibin2,ien) !get result from last flux surface
!          print*,'tdnpa0::ibin1 has decreased by 1 in atten'

          elseif (ichangebin24.eq.1) then !increasing ibin1, in class 2 or 4
          atten(ibin1,ibin2,ien)=ds/stoplamda(ibin1,ien) &
           +atten(ibin1-1,ibin2,ien) !get result from last flux sufrace
!          print*,'tdnpa0::ibin1 has increased by 1 in atten'

          else !no change in bins or direction
          atten(ibin1,ibin2,ien)=ds/stoplamda(ibin1,ien) &
           +atten(ibin1,ibin2,ien)
          endif

 195      continue

          iclass12=0
          iclass23=0
          ichangebin13=0
          ichangebin24=0

!VT........................................................................

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
!BHandYuP110715: Account for bsign/cursign
          xcomp=(-cursign*dpsidz*xx-bsign*fpsi_*yy)/rr**2
          ycomp=(-cursign*dpsidz*yy+bsign*fpsi_*xx)/rr**2
          zcomp=+cursign*dpsidr/rr

!.......................................................................
!     Compute the dot product and the cosine
!.......................................................................

          dotprd=-(xcomp*alphad(1)+ycomp*alphad(2)+zcomp*alphad(3))
          bbmag=sqrt(xcomp**2+ycomp**2+zcomp**2)
!BHandYuP110715: cospitch measured with vpar pos in pos vphi dirn
!BHandYuP110715: (CCW dirn viewed from top).
          cospitch=bsign*dotprd/bbmag
          cospitch=max(cospitch,-1.d0) ! to be sure cos>=-1
          cospitch=min(cospitch,+1.d0) ! to be sure cos<=+1
          sinpitch=sqrt(1.d0-cospitch**2)
          angle=acos(cospitch)
          sang(ibin1,ibin2)=sang(ibin1,ibin2)+ds*angle

!.......................................................................
!     Poloidal angle   (rmag is major radius of magnetic axis).
!.......................................................................

          arg1=zz-zmag
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
!     *******Have finished sightline and attenuation determination******
!.......................................................................

 200    continue

!        do ien=1,nen_npa
!           do l=1,lrzmax
!            write(*,*)'l,ien,atten(l,j=1:4,ien)=',l,ien,atten(l,1:4,ien)
!           enddo
!        enddo


!.......................................................................
!     Obtain weighted angle between viewing and toroidal directions,
!     and poloidal angle, in each of the lrzmax radial bins.
!.......................................................................

        do 210  j=1,ibin2  ! Radial bin class, ibin2.le.4
          do 211  i=1,lrzmax
            if(sxry(i,j).eq.zero)  go to 211
            sang(i,j)=sang(i,j)/sxry(i,j)
            spol(i,j)=spol(i,j)/sxry(i,j)
            if (spol(i,j) .gt. pi) spol(i,j)=pi
 211      continue
 210    continue

!        do l=1,lrzmax
!           write(*,*)'l,sxry(l,j=1:4)=',l,sxry(l,1:4)
!        enddo
!        do l=1,lrzmax
!           write(*,*)'l,sang(l,j=1:4)=',l,sang(l,1:4)
!        enddo
!        do l=1,lrzmax
!           write(*,*)'l,spol(l,j=1:4)=',l,spol(l,1:4)
!        enddo


!.......................................................................
!     Obtain energy flux for each energy and sightline
!.......................................................................

        call bcast(emitnpa,zero,nen_npa*npaproc)
        call bcast(efluxwk(1,1),zero,nen_npa*lrzmax)
        mpisz= nen_npa ! number of elements in efluxwk(:,i)

!       set atten bins that will be double counted to some
!       very high attenuation

        atten(1,2,:)=100
        atten(lrzmax,3,:)=100
        atten(1,4,:)=100

        do i=1,lrzmax   !lrzmax
!MPIINSERT_MPIWORKER_I
!MPIINSERT_IF_RANK_EQ_MPIWORKER
        do j=1,ibin2  ! Radial bin class, ibin2.le.4
              if(sxry(i,j).eq.zero)  go to 311
!     Obtain neutral emittance emitnpa (#/sec*cm**3*ster*eV),
!     for each npa_process().
!     enn density pertains to each npa_process(), generally
!     a namelist input neutral density.  For npa_process(5)=
!     'radrecom', enn(1:lrzmax,5) is set to the electron density.
              call tdnpa(vel,sang(i,j),ibin(i,j),spol(i,j), &
                   sigmanpa,emitnpa)
              do kk=1,npaproc
              do 320 ien=1,nen_npa
!              if (j.EQ.1 .AND. i.EQ.1) THEN
!     VT specifically prohibit double counting of l=1 surface
!     VT.......... have to interpolate values for emitnpa since the grid
!     VT           for emitnpa is velocity and we want energy....

                 if (npa_process(kk).ne.'notset') then
                 if (atten_npa.eq."enabled") then
                     efluxwk(ien,i)= efluxwk(ien,i)+ &
                         enn(ibin(i,j),kk) &
                         *sxry(i,j)*emitnpa(ien,kk)*exp(-atten(i,j,ien))
                 else
                     efluxwk(ien,i)= efluxwk(ien,i)+ &
                         enn(ibin(i,j),kk) &
                         *sxry(i,j)*emitnpa(ien,kk)
                 endif
                 endif  ! On npa_process(kk)
!     VT................................................................
!                 print*,"eflux_npa(ien,nn)=",eflux_npa(ien,nn)
!                 print*,"enn(ibin(i,j),1)=",enn(ibin(i,j),1)
!                 print*,"enn(ibin(i,j),2)=",enn(ibin(i,j),2)
!                 print*,"sxry(i,j)=",sxry(i,j)
!                 print*,"exp(-attend..)=",exp(-atten(i,j,ien))
!                 if (ien.eq.1)print*,"emitnpa(,1)=",emitnpa(1:nen_npa,1)
!                 print*,"----------------------------"

 320        continue !ien=1,nen_npa
            enddo  ! kk=1,npaproc

 311    continue
        enddo ! j=1,ibin2
!MPIINSERT_ENDIF_RANK
!MPIINSERT_RECV_EFLUXWK
!MPIINSERT_SEND_EFLUXWK
        enddo ! i=1,lrzmax

        call bcast(eflux_npa(1,nn),zero,nen_npa)

        do ien=1,nen_npa  ! for each energy bin
        do i=1,lrzmax ! sum-up contributions from all flux surfaces
           eflux_npa(ien,nn)=eflux_npa(ien,nn)+ efluxwk(ien,i)
        enddo ! i=1,lrzmax
        enddo ! ien=1,nen

!MPIINSERT_BARRIER
!MPIINSERT_BCAST_EFLUX_NPA




!.......................................................................
!     Integrate energy flux spectra from enmin_npa to enmax_npa.
!     Assume the flux varies as exp(-e/t) between flux values,
!     to improve accuracy.
!.......................................................................

        inegsxr(nn)=nen_npa
        efluxt(nn)=0.0
        do 400  i=1,nen_npa
!          write(*,*)'tdnpa0:i,nn,eflux_npa(i,nn),zero=',
!     +          i,nn,eflux_npa(i,nn),zero
          if (eflux_npa(i,nn).le.zero) then
            inegsxr(nn)=i-1
            go to 402
          endif
!BH100811         if (i.eq.nen) go to 402
!BH100811  Added following line.
         if (i.eq.nen_npa) go to 402
!         write(*,*)'tdnpa0: i=',i
!BH100419:  Wonder why I have this (also in tdsxr0.f)??
!BH100811:  An answer:  get out-of-bounds for i=nen_npa.
!BH100419          if (i.eq.nen_npa) go to 402

!         Avoid some overflows:
          if(eflux_npa(i+1,nn).lt.zero) go to 402
          if (eflux_npa(i,nn).lt.em100*em100) go to 400
!BH160507     if(eflux_npa(i,nn).eq.eflux_npa(i+1,nn))  go to 401
!BH160507 Allowing for possible increase in eflux_npa with en.
!BH160507 Could modify to fork with eflux_npa(i:i+1) are too close
          if(eflux_npa(i,nn).le.eflux_npa(i+1,nn))  go to 401
          if(eflux_npa(i+1,nn).eq.zero)  go to 401

          t=(en_(i+1)-en_(i))/log(eflux_npa(i,nn)/eflux_npa(i+1,nn))
!BH160507     c=eflux_npa(i,nn)*exp(en_(i)/t)
!BH160507          write(*,*)'tdnpa0.f: i,t,c,eflux_npa(i,nn)=',
!BH160507     +                         i,t,c,eflux_npa(i,nn)
!BH160507     efluxt(nn)=efluxt(nn)+c*t*(exp(-en_(i)/t)-exp(-en_(i+1)/t))
          efluxt(nn)=efluxt(nn)+t*(eflux_npa(i,nn)-eflux_npa(i+1,nn)) &
            *1.e3
          go to 400
 401      efluxt(nn)=efluxt(nn)+eflux_npa(i,nn)*(en_(i+1)-en_(i))*1.e3
 400    continue
 402    continue

!        write(*,*)
!        write(*,*)'tdnpa0: emitnpa(1:30,1:npaproc)'
!        do i=1,npaproc
!           write(*,*)'tdnpa0: emitnpa(*,i)',(emitnpa(jj,i),jj=1,30)
!        enddo
!        write(*,*)

!        print*,'npa_mesh=[',en_(1:nen_npa),']'
!        print*,'eflux_npa=[',eflux_npa,']'
!        write(*,*)'tdnpa0: efluxt(1:nn)=',efluxt(1:nn)
!
        icall="notfirst"

!.......................................................................
!*********END loop over viewing angles      **************************
!.......................................................................
!        print*,'------------------------------'
!        print*,'Done with sightline',nn

 420  continue

!$$$      iistep=0
!$$$      do nn=1,nv
!$$$         print*,'tempp4 x=[',tempp4(iistep+1:iistep+lensxr(nn)),']'
!$$$         print*,'tempp5 y=[',tempp5(iistep+1:iistep+lensxr(nn)),']'
!$$$         print*,'tempp6 z=[',tempp6(iistep+1:iistep+lensxr(nn)),']'
!$$$         iistep=iistep+lensxr(nn)
!$$$      enddo
!$$$          iiii=int(nen_npa)
!$$$          print*, "tdnpa0::atten for Neutral energy:",en_(iiii)
!$$$          print*,"atten_ibin21=[",atten(1:lrzmax,1,iiii)
!$$$          print*,"atten_ibin22=[",atten(1:lrzmax,2,iiii)
!$$$          print*,"atten_ibin23=[",atten(1:lrzmax,3,iiii)
!$$$          print*,"atten_ibin24=[",atten(1:lrzmax,4,iiii)
!$$$          iiii=int(nen_npa/2)
!$$$          print*, "tdnpa0::atten for Neutral energy:",en_(iiii)
!$$$          print*,"atten_ibin21=[",atten(1:lrzmax,1,iiii)
!$$$          print*,"atten_ibin22=[",atten(1:lrzmax,2,iiii)
!$$$          print*,"atten_ibin23=[",atten(1:lrzmax,3,iiii)
!$$$          print*,"atten_ibin24=[",atten(1:lrzmax,4,iiii)





!.......................................................................
!     Plot eflux_npa........
!.......................................................................

      nenaa=nen_npa
!BH100418:  copied following call here.  Probably should restore
!BH100418:  icall logic for npa....
!BH100418:      call tdsxrvw

!BH100418: Follwing call to tdsxrplt was commented out by VT

!BH100418: Normalize nps flux
!BH100418:      do nn=1,nv_npa
!BH100418:         do i=2,nen_npa
!BH100418:            eflux_npa(i,nn)=eflux_npa(i,nn)/eflux_npa(1,nn)
!BH100418:         enddo
!BH100418:         eflux_npa(1,nn)=1.0
!BH100418:      enddo

!     iplt3d set in tdchief
      if (iplt3d.ne.0 .or. n.eq.0 .or. n.eq.nstop &
           .and. iplotnbi.eq.'yes') then
         call tdsxrplt(en_,eflux_npa,nen_npa,nenaa, &
              efluxt,nv_npa,inegsxr,softxry,lnwidth)

!.......................................................................
!     Plot NPA view cords in poloidal cross-section, if eqmod="enabled"
!.......................................................................

         if (eqmod.eq."enabled") then
            nv=nv_npa
            call tdsxrvw(tempp4,tempp5,tempp6)
         endif
      endif ! On iplt3d, etc.

      deallocate(tempp4,tempp5,tempp6,STAT=istat1)

!.......................................................................
!     Put eflux_npa into eflux in order to pass the data to netcdfrw2.
!     (It would be cleaner to put eflux_npa in common, as done for
!     eflux.  BH100816)
!.......................................................................

      do nn=1,nv_npa
         do ien=1,nen_npa
            eflux(ien,nn)=eflux_npa(ien,nn)
         enddo
      enddo

      return
      end
end module tdnpa0_mod
