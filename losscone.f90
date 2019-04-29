module losscone_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast
  use lossorbm_mod, only : lossorbm
  use zcunix_mod, only : allocate_error

  !---END USE

!
!

contains

      subroutine losscone
      use param_mod
      use comm_mod
      use r8subs_mod, only : luf
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine computes the location of the loss orbits and puts
!     the information into gone(i,j,k,indxlr_).
!     [Subroutine called from ainitial. indxlr_ is set by tdnflxs(ll)
!      in a loop ll=lrors,1,-1 in tdinitl.]
!     Several loss models are available and are invoked by setting
!     lossmode(k) to the chosen "string" in the input:
!     "snk0"= orbits whose energy is less than esink
!     are lost from the system with characteristic time equal to
!     the bounce time.
!     "snk0accl"=same as "snk0" except the loss mechanism is acceler-
!     ated by factor xsinkm (input).
!     "electron" forces a loss cone for general species electrons
!     for perpendicular energies greater that eperc and parallel
!     energies greater than eparc.
!     "mirrorcc"=potential square well model for tandem mirror central
!     cell calculations; orbits passing through a potential jump
!     at the throat of magnitude ephicc(kev) are lost.
!     "mirrsnk"=combination of mirrorcc and snk0.
!     "simplban"=
!     The array gone is utilized in subroutine coefload, which
!     determines suitable Krook operator loss times, depending on the
!     local value of gone.
!..................................................................


!     Pointers for dynamic memory allocation, local usage:
      real(c_double), dimension(:), pointer :: upar,uprp,rho_a
      integer, dimension(:,:,:), pointer :: notlost
      real(c_double), dimension(:,:,:), pointer :: dnotlost

      do 100 k=1,ngen  !Loops down to end of the subroutine

         !call bcast(gone(0,0,k,indxlr_),zero,iyjx2) ! YuP commented-out
         !YuP[09-10-2014] do not initialize gone array here.
         ! it was already done in ainalloc.
         xmul=1.

         !===========================================================
         if (lossmode(k).eq."snk0" .or. lossmode(k).eq."snk0accl" .or. &
              lossmode(k) .eq. "mirrsnk") then
            if (lossmode(k) .eq. "snk0accl") xmul=xsinkm
            if (esink.gt.em90) then
               xsink=sqrt(esink/fions(k))
!            else
!               xsink as specified in namelist
            endif

!..................................................................
!     Find the first mesh point x "outside" the sink hole.
!..................................................................

            do 10 j=2,jx
               if (x(j) .ge. xsink) go to 20
 10         continue
            j=jx
 20         continue
            do 30 jj=1,j
               do 31 i=1,iy
                  gone(i,jj,k,indxlr_)=-xmul
 31            continue
 30         continue


!..................................................................
!     Electron losses defined by eparc and eperc
!..................................................................

         elseif (lossmode(k).eq."electron" .and. k.eq.kelecg) then

            v2=vnorm*vnorm/(clight*clight)
            do 40  i=1,iy
               do 41  j=1,jx
                  gampar=sqrt(1.+(coss(i,l_)*x(j))**2*v2)
                  gamper=sqrt(1.+(sinn(i,l_)*x(j))**2*v2)
                  epar=(gampar-1.)*restmkev !YuP[01-2011] was 512.
                  eper=(gamper-1.)*restmkev !YuP[01-2011] was 512.
                  if(epar.gt.eparc(k,lr_).or.eper.gt.eperc(k,lr_)) then
                     gone(i,j,k,indxlr_)=-1.
                  else
                     gone(i,j,k,indxlr_)=0.
                  endif
 41            continue
 40         continue

!.......................................................................
!     Simple Banana + Gyro-Orbit loss model:
!     -trapped particle bananas plus gyro-orbit width (circ plasma)
!     -greater than dist to plasma edge are lost.
!     -Circulating particles with gyro-orbit width greater than
!      distance to the edge of the plasma are lost.
!     Banana width uses poloidal magnetic field which is an average
!     of Bpol at each radius and of the edge Bpol.
!     bthr0 is at outer midplane.  Might be that bth ---at thetapol
!       =pi/2--- would be more appropriate.  Might want to check this
!       by orbit integrations.  (BobH, 011125).
!.......................................................................

         elseif (     lossmode(k).eq."simplban" &
                 .or. lossmode(k).eq."simplbn1" &
                 .or. lossmode(k).eq."simplbn2") then
!$$$            if (ionce.eq.0) then
!$$$                 ionce=1
!$$$                 write(*,*)'losscone:bthr(ll),ll=1,lrzmax',
!$$$     +                               (bthr(ll),ll=1,lrzmax)
!$$$                 write(*,*)'losscone:bthr0(ll),ll=1,lrzmax',
!$$$     +                               (bthr0(ll),ll=1,lrzmax)
!$$$                 write(*,*)'Losscone:ll,eqbpol(1:lorbit,ll=1,lrzmax)'
!$$$                 do ll=1,lrzmax
!$$$                    write(*,*) ll,(eqbpol(lll,ll),lll=1,lorbit(lr_))
!$$$                 enddo
!$$$            endif

!         Project to (Last Closed Flux Surface)*simpbfac from flux
!         surface radii tabulated on the rya() mesh
            rpconz_max=rpconz(lrzmax)+(simpbfac-rya(lrzmax))* &
                       (rpconz(lrzmax)-rpconz(lrzmax-1))/ &
                       (rya(lrzmax)-rya(lrzmax-1))
            bthr0edge=abs(bthr0(lrzmax)+(simpbfac-rya(lrzmax))* &
                       (bthr0(lrzmax)-bthr0(lrzmax-1))/ &
                       (rya(lrzmax)-rya(lrzmax-1)))
            bmod0edge=abs(bmod0(lrzmax)+(simpbfac-rya(lrzmax))* &
                       (bmod0(lrzmax)-bmod0(lrzmax-1))/ &
                       (rya(lrzmax)-rya(lrzmax-1)))
            bthr0avg=0.5*(abs(bthr0(lr_))+bthr0edge)
            x_loss=abs(bnumb(k)*charge) &
                 *(rpconz_max-rpconz(lr_))/(fmass(k)*vnorm*clight)
            r_loss=rpconz_max-rpconz(lr_)
            omcnst=abs(bnumb(k)*charge)/(fmass(k)*vnorm*clight)
!$$$            write(*,*)'Losscone: lr_,rpconz_max,x_loss',
!$$$     +                           lr_,rpconz_max,x_loss
!$$$            write(*,*)'losscone: bthr0avg,bmod0edge',
!$$$     +                           bthr0avg,bmod0edge

!BH160403: Rethinking gyro-orbit loss.  Should evaluate bmag0 at outside
!BH160403: of banana orbit.  Banana is centered on lr_ flux surface.
!BH160403: Linearly interpolate for bmod0 at outside of the the banana.

            drpconz0=(rpconz(2)-rpconz(1))/(rya(2)-rya(1))
            drpconz1=(rpconz(lrzmax)-rpconz(lrzmax-1))/ &
                     (rya(lrzmax)-rya(lrzmax-1))
            dbmod00=(bmod0(2)-bmod0(1))/(rya(2)-rya(1))
            dbmod01=(bmod0(lrzmax)-bmod0(lrzmax-1))/ &
                     (rya(lrzmax)-rya(lrzmax-1))
            do 50  i=1,iy
               do 51  j=3,jx    !was j=1,jx ! YuP
                  ! YuP[7-22-2014] skip j=1:2 (do not apply a loss cone)
                  !Sometimes a loss cone at v~0 gives unstable solution
                  !and issues for boundary condition at v=0.
                  !Other than that, the effect from skipping j=1:2
                  !is very small - difference in ~6th digit.
                  if (i.gt.itl .and. i.lt.itu) then  !ZOW trapped particles
!BH160529
!BH160529
!BH160529
!BH160529
!BH160529
!BH160529
                  ! BH: We take the abs() for delb1, since the orbit
                  !     shift is from BA radial flux surface, and outside
                  !     scrapeoff will occur on first or second part of
                  !     the bounce.
! YuP[01-2017] Commented with ! by YuP: does not work properly.
! YuP  One of problems: Failure in ngen=2 multiURF tests
! YuP  because subr. deltar is setup/called for one general species only,
! YuP  and so all deltar* arrays are saved for k=1 gen.species only.
!                     delb1=x(j)*abs(coss(i,indxlr_))/(omcnst*bthr0(lr_))
!                     call lookup(rpconz(lr_)+delb1,rpconz,lrzmax,
!     +                           wtu,wtl,linterp)
!                     rpconz1=wtl*rpconz(linterp-1)+wtu*rpconz(linterp)
!                     bmod01=wtl*bmod0(linterp-1)+wtu*bmod0(linterp)
!                     gyrorad1=x(j)*abs(sinn(i,indxlr_))/(omcnst*bmod01)
!                     gone1=0.
!                     if (delb1+gyrorad1 .gt. r_loss) gone1=-1.

!                     if (lossmode(k).eq.'simplbn1')then !YuP[2017-11-21]corr.lossmode(1)
!                        if (delb1+gyrorad1 .gt. r_loss) then
!                           gone(i,j,k,indxlr_)=-1.
!                        else
!                           gone(i,j,k,indxlr_)=0.
!                        endif   ! on delb1
!                     endif  ! on lossmode(k).eq.'simplbn1'

!                     gyrorad2=x(j)*abs(sinn(i,indxlr_))/bmod0edge/omcnst
!                     delb2=x(j)*abs(coss(i,indxlr_))/bthr0avg/omcnst
!                     gone2=0.
!                     if ((gyrorad2 + delb2).gt. x_loss/omcnst) gone2=-1.

!                     if (lossmode(k).eq.'simplbn2') then !YuP[2017-11-21]corr.lossmode(1)
!                        if ((gyrorad2 + delb2)
!     +                       .gt. x_loss/omcnst) then
!                           gone(i,j,k,indxlr_)=-1.
!                        else
!                           gone(i,j,k,indxlr_)=0.
!                        endif
!                     endif  ! on lossmode(k).eq.'simplbn2'

!BH160529 There are problems with application of NB for these models.
!BH160529 Need to first have NB birth points shifted inwards to their
!BH160529 BA radius, if they are going to be scraped off at a banana width.
!BH160529 First, compare delb with deltaRZ calc, which we will subsequently
!BH160529 implement for scrapeoff.
!BH160529 Calc orbit shift with deltarho, to be consistent with calc
!BH160529 of NB birth points (with ndeltarho='freya')
!BH160529 Note: cql3d-FOW convention for deltarho shift.
!BH170708 Note: deltarho calc in subroutine deltar enabled for ngen=1 at this time.
!BH170708       See email YuP: 2017-01-04 12:27
!YuP: rho1 is not used, commenting out
!YuP: (can give a problem in case of lossmode(1)=disabled; lossmode(2)=simplban)
!YuP                     rho1=rya(lr_)+x(j)*vnorm*abs(deltarho(i,1,lr_))
                  ! BH: We take the abs() for deltarho, since the orbit
                  !     shift is from BA radial flux surface, and outside
                  !     scrapeoff will occur on first or second part of
                  !     the bounce. deltarho is anti-symmetric about pi/2.
!YuP                     if (rho1.lt.zero) rho1=zero
!YuP                     if (rho1.gt.1.1d0) rho1=1.1d0
!YuP                     call lookup(rho1,rya(1),lrzmax,wtu,wtl,linterp)
!YuP                    !Interpolate, or extrapolate a bit.
!YuP                     if (linterp.eq. lrzmax .and. wtu.eq.1d0) then
!YuP                        R1=rpconz(linterp)+(rho1-rya(linterp))*drpconz1
!YuP                        b1=bmod0(linterp)+(rho1-rya(linterp))*dbmod01
!YuP                     elseif(linterp.eq.2 .and. wtl.eq.1d0) then
!YuP                        R1=rpconz(linterp)+(rho1-rya(linterp))*drpconz1
!YuP                        b1=bmod0(1)+(rho1-rya(1))*dbmod00
!YuP                        if (R1.lt.rmag) then
!YuP                           R1=rmag
!YuP                           b1=bmag
!YuP                        endif
!YuP                     else
!YuP                        R1=wtl*rpconz(linterp-1)+wtu*rpconz(linterp)
!YuP                        b1=wtl*bmod0(linterp-1)+wtu*bmod0(linterp)
!YuP                     endif
!YuP                     gyrorad=x(j)*abs(sinn(i,indxlr_))/(omcnst*bmod01)
!YuP                     delb_deltarho=R1-rpconz(lr_)
!YuP                     gone0=0.
!YuP                     if (delb_deltarho+gyrorad.gt.r_loss) gone0=-1.
!BH170708  For backward compatibility, reset  gyrorad and delb_deltarho
                    gyrorad=x(j)*abs(sinn(i,indxlr_))/(omcnst*bmod0edge)
                     delb_deltarho=x(j)*abs(coss(i,indxlr_))/(omcnst* &
                                   bthr0avg)
                     if (lossmode(k).eq.'simplban') then
                        if (delb_deltarho+gyrorad .gt. r_loss) then
                           gone(i,j,k,indxlr_)=-1.
                        else
                           gone(i,j,k,indxlr_)=0.
                        endif
                     endif
!BH160529
!BH160529
!BH160529      if (i.eq.(itl+1) .and. j.eq.3) then
!BH160529      write(*,*)'rpconz_max,r_loss',rpconz_max,r_loss
!BH160529      write(*,*)'i,j,lr_,rya,rho1,R,R1,delb_deltarho,'//
!BH160529     +  'delb1,delb2,gyrorad,gyrorad1,gyrorad2,gone0,gone1,gone2'
!BH160529      endif
!BH160529      write(*,*) i,j,lr_,rya(lr_),rho1,rpconz(lr_),R1,delb_deltarho,
!BH160529     +     delb1,delb2,gyrorad,gyrorad1,gyrorad2,gone0,gone1,gone2
!BH160529
!BH160529  Conclusion from comparison of above results, particularly gone0
!BH160529  (simplban, new), gone1(simplbn1, 1st adj), gone2(simplbn2, old):
!BH160529  All about the same, although simplbn2 least severe, simplbn1 is
!BH160529  most, and new method is intermediate. See log.1.simplban_extended.
!BH160529  Application of ndeltarho='freya' (a relevant correction) will
!BH160529  reduce the co-inj losses.
!BH160529

                  else  ! On i

!BH-YuP140807: Should have used bmod0() for transitting particles.
!BH-YuP140807              if (x(j)*abs(sinn(i,indxlr_))/bmod0edge
!BH160529 Could do more for passing particles with deltarho....
                     if (x(j)*abs(sinn(i,indxlr_))/bmod0(indxlr_) &
                        .gt. x_loss) then
                        gone(i,j,k,indxlr_)=-1.
                     else
                        gone(i,j,k,indxlr_)=0.
                     endif
                  endif  ! on i
 51            continue
 50         continue
!..................................................................
!     Now for mirror cases.
!..................................................................

         else if (lossmode(k).eq."mirrorcc") then

            call lossorbm(ephicc,k)

!..................................................................
!     Now for file specification of prompt losses
!..................................................................

         else if (lossmode(k).eq."frm_file") then

!     Allocate local pointered memory for reading lossfile(k)
            allocate (uprp(n_uprp),STAT = istat)
            if(istat .ne. 0) &
                 call allocate_error("uprp, sub rdc_multi",0,istat)
            call bcast(uprp,zero,SIZE(uprp))

            allocate (upar(n_upar),STAT = istat)
            if(istat .ne. 0) &
                 call allocate_error("upar, sub rdc_multi",0,istat)
            call bcast(upar,zero,SIZE(upar))

            allocate (rho_a(n_psi),STAT = istat)
            if(istat .ne. 0) &
                 call allocate_error("rho_a, sub rdc_multi",0,istat)
            call bcast(rho_a,zero,SIZE(rho_a))

            allocate (notlost(n_uprp, n_upar, n_psi),STAT = istat)
            if(istat .ne. 0) &
                 call allocate_error("notlost, sub rdc_multi",0,istat)
            call ibcast(notlost,1,SIZE(notlost))

            allocate (dnotlost(n_uprp, n_upar, n_psi),STAT = istat)
            if(istat .ne. 0) &
                 call allocate_error("dnotlost, sub rdc_multi",0,istat)
            call bcast(dnotlost,one,SIZE(notlost))


! Open file to read lost regions, specified on an equispaced
! uprp,upar grid at each radial flux surface.

            iunit=14
            open(unit=iunit,file='lossfile(k)',status='unknown')
            read(iunit,309) n_uprp !number of uprp points
            read(iunit,309) n_upar !number of upar points
            read(iunit,309) n_R !number of radial (PSI) points
            read(iunit, 3310) vc_cgs !Max vel on the grid
            read(iunit, 3310) upar_min,upar_max !Generally, -1. +1.
                                                ! Max uprp is assumed
                                                ! to = max upar.
                                                ! min uprp=0.
            read(iunit,310) (((notlost(i_prp, i_par, i_psi), &
                  i_prp =0,n_uprp-1), i_par =0,n_upar-1), i_psi =1,n_R)
! notlost=0 for loss cone (if at least one i_phase or i_posn is lost),
! notlost=1 otherwise
            close(iunit)        ! Close file
!-----------------------------------------------------------------
 3310       format(1p6e18.10)
 309        format(10i10)
! BH remove 1p descriptor for gfortran 4..5.1 (20101208)
![Nonsense for i decriptor].
! 310        format(1p6i2)
 310        format(6i2)


            uprp_min=0.d0
            uprp_max=upar_max

!     Check enorm for notlost data:

            vc_cgs2=vc_cgs*vc_cgs
            gammac=sqrt(1.d0+vc_cgs2/clite2)
            if ((gammac-1.d0).lt.1.e-6) then
               enormc=0.5d0*fmass(1)*(vc_cgs)**2/ergtkev
            else
               enormc=(gammac-1.d0)*fmass(1)*clite2/ergtkev
            endif
            write(*,*)'Enorm in losscone, gammac =',enormc,gammac

!.................................................................
!     The rho_a radial mesh and associated coeffc will be reduced
!     by a factor of 2, if lrz.eq.n_psi/2, enabling cql3d to run
!     on half the number of flux surfaces used in the full-wave code.
!     Check if n_psi.eq.lrz.
!     If not, check n_psi/2.eq.lrz.  If so, omit every 2nd radial point
!       of the diffusion coeff grid and coeffs, and reset n_psi.
!       This enables factor of 2 reduction of the cql3d radial mesh.
!       The cql3d and du0u0_input radial meshes are assumed to be
!         the same (or close).
!       Future modification:  Interpolate the du0u0_input radial mesh
!                             to the cql3d radial mesh.
!     Check if n_psi.eq.lrz:  if not, STOP

      if (n_psi.ne.lrz .and. &
          int((n_psi+1.1)/2) .eq. lrz) then

!.................................................................
!     Adjust the radial mesh to use every second one ==> 32 radial
!     points from 64, or 64 radial points from 128.
!.................................................................

         n_psi=n_psi/2
         write(*,*)'losscone: n_psi,lrz=',n_psi,lrz
         do i_psi=1,n_psi
            rho_a(i_psi)=rho_a(2*i_psi)
            do i_upar=1,n_upar
               do i_uprp=1,n_uprp
                  notlost(i_uprp,i_upar,i_psi)= &
                      notlost(i_uprp,i_upar,2*i_psi)
               enddo
            enddo
         enddo

         write(*,*)'losscone: rho_a, no. elements/2:', &
              (rho_a(i),i=1,n_psi)

      endif

      if (n_psi.ne.lrz) then
         WRITE(*,*)'losscone:  n_psi.ne.lrz'
         STOP
      endif


!.................................................................
!  The velocity grids are assumed to be equi-spaced.  Renormalize
!  to account for possible vc_cgs.ne.vnorm

!  Setup arrays of initial upar0,uprp0, momentum-per-rest-mass at the
!  minimum B-field point.
!.................................................................

      dupar=(upar_max-upar_min)/(n_upar-1)*vc_cgs/vnorm
      do i=1,n_upar-1
         upar(i)=upar_min+(i-1)*dupar
      enddo

      duprp=(uprp_max-uprp_min)/(n_uprp-1)*vc_cgs/vnorm
      do i=1,n_uprp-1
         uprp(i)=uprp_min+(i-1)*duprp
      enddo
      dupdup=dupar*duprp

!.................................................................
!  Bilinear interpolate notlost onto the cql3d gone(i,j,k,ll) grid.
!  [Coding modified from rdc_multi.]
!  If any quantity remains on the grid which is less than 1.d0
!  then set it to zero.
!  If vc_cgs.lt.vnorm, and values of gone(i,j=jc,k,ll)=0.d0, then
!  set gone (i,j=jc+1:jx,k,ll)=0.d0.
!.................................................................

      do i_psi=1,n_psi
         do i_upar=1,n_upar
            do i_uprp=1,n_uprp
               if (notlost(i_uprp,i_upar,i_psi).eq.0) then
                  dnotlost(i_uprp,i_upar,i_psi)=-one !i.e.,
                                                     !gone-like lost
               else
                  dnotlost(i_uprp,i_upar,i_psi)=zero !i.e.,
                                                     !gone-like not lost
               endif
            enddo
         enddo
      enddo

      do ll=1,lrz
      do j=1,jx
         do i=1,iy
            xupar=x(j)*coss(i,ll)
            xuprp=x(j)*sinn(i,ll)
            if (xupar.lt.upar_max .and. xupar.gt.upar_min) then
               if (xuprp.lt.uprp_max .and. xuprp.ge.uprp_min) then
                  ipar=(xupar-upar_min)/dupar+1
                  iprp=(xuprp-uprp_min)/duprp+1
                  w00=(upar(ipar+1)-xupar)*(uprp(iprp+1)-xuprp)/dupdup
                  w01=(upar(ipar+1)-xupar)*(xuprp-uprp(iprp))/dupdup
                  w10=(xupar-upar(ipar))*(uprp(iprp+1)-xuprp)/dupdup
                  w11=(xupar-upar(ipar))*(xuprp-uprp(iprp))/dupdup
                  wsum=w00+w01+w10+w11
                  gone(i,j,k,ll)=(w00*dnotlost(iprp,ipar,ll) &
                           +w01*dnotlost(iprp+1,ipar,ll) &
                           +w10*dnotlost(iprp,ipar+1,ll) &
                           +w11*dnotlost(iprp+1,ipar+1,ll))
                  if (gone(i,j,k,ll).ne.zero) gone(i,j,k,ll)=-one
               endif
            endif
         enddo
      enddo
      enddo

!  If vc_cgs.lt.vnorm, haven't covered the whole iy,jx grid.
!  If last point ga vc_cgs indicates prompt loss, then extend
!    this at constant pitch angle to the jx edge of the mesh.
         jc=luf(vc_cgs/vnorm,x,jx)-1 !Will be less than jx if
                                  !vc_cgs/vnorm.lt.x(jx-1)

         write(*,*)'losscone: jx,jc=',jx,jc

         do i=1,iy
         if (jc.lt.jx .and. gone(i,jc,k,ll).eq.-one) then
            do j=jc,jx
               gone(i,j,k,ll)=-one
            enddo
         endif
         enddo




         endif      ! on lossmode(k)
         !===========================================================

         if (lossmode(k) .eq. "mirrsnk") call lossorbm(ephicc,k)

 100  continue ! k  species

!BH160529      if (lr_.eq.1) STOP  !This is for stop after extended o/p
      return
      end


end module losscone_mod
