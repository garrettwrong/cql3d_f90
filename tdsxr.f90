module tdsxr_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use cfpleg_mod, only : cfpleg
  use r8subs_mod, only : dcopy
  use tdnflxs_mod, only : tdnflxs
  use tdwrng_mod, only : tdwrng
  use zcunix_mod, only : terp2

  !---END USE

!
!

contains

      subroutine tdsxr(ien,chi,thpol,ll,edotei,edotee)
      use param_mod
      use comm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
!..................................................................
!     This routine computes the BREMSSTRAHLUNG contribution to the
!     SXR spectrum from e-i and e-e collisions.
!     The integration requires expansion of the
!     electron distribution function and the differential cross
!     section in Legendre polynomials - the latter calculation is
!     done in subroutine TDSETSXR.
!
!     ien is the index of photon energy of interest
!
!     chi is the angle between B and the viewing angle
!
!     thpol is the poloidal angle (representing the point of
!     intersection of the viewing angle and the relevant flux
!     surface.
!
!     ll is the index of the flux surface (rz(ll))
!..................................................................

      character*8 ione
      allocatable :: plegg(:) ! local array

      save  ! need to save plegg at least

      if (.NOT. ALLOCATED(plegg)) then
        allocate( plegg(0:msxr) )
        call bcast(plegg,zero,(msxr+1))
      endif


      if (ien .gt. 1) go to 700 ! skip some of calc., including plegg();
                                ! They are done at first call with ien=1

!..................................................................
!     Set up variables local to flux surface ll
!     Special treatment if f not solved on ll (case lrzdiff=enabled)
!...................................................................

      if (lrzdiff .ne. "enabled") then
        illcal=1  ! No lrzdiff
        go to 1
      else
        illcal=0  ! lrzdiff
      endif

      do 2 il=1,lrz
        if (lrindx(il) .eq. ll) illcal=1  !No lrzdiff for this ll
        if (lrindx(il) .gt. ll) go to 3
 2    continue
 3    if (illcal .eq. 1) go to 1
      ilplus=il-1/(lrz+2-il)
      ilrrigt=lrindx(ilplus)
      ilrleft=lrindx(ilplus-1+1/ilplus)
      if (abs(rz(ilrrigt)-rz(ll)) .le. abs(rz(ilrleft)-rz(ll))) then
        ilrnear=ilrrigt
      else
        ilrnear=ilrleft
      endif

 1    continue
      if (illcal .eq. 1) then
        call tdnflxs(indxlr(ll))
      else
        call tdnflxs(indxlr(ilrnear))
        lr_=ll
      endif

!...................................................................
!     compute the legendre coefficients of f at thpol
!     sax is the field line intersection point with the viewing
!     angle.
!...................................................................

      sax=thpol*zmax(lr_)/pi
      if (abs(sax-zmax(lr_)) .lt. 5.e-11) sax=zmax(lr_)

!....................................................................
!     Copy the electron distribution function into a temporary array
!     used for Legendre coefficient calculation, and obtain poloidal
!     bin index lll.
!     If distribution not computed (ll not in lrindx) then assume a
!     non-relativistic Maxwellian distribution function.
!....................................................................

      if (illcal .eq. 1) then
        call dcopy(iyjx2,f(0:iyjx2-1,0,kelec,l_),1,temp2(0:iyjx2-1,0),1)
      else
        zthta=0.5*fmass(kelec)*vnorm**2/(temp(kelec,lr_)*ergtkev)
        zcoef=vnorm**3*reden(kelec,lr_)/ &
          (2.*pi*temp(kelec,lr_)*ergtkev/fmass(kelec))**1.5
        do 7 m=0,msxr
          call bcast(feta(1:jx,m),zero,jx)
 7      continue
!     Legendre-coefficient of Maxwellian is the Maxwellian for m=0
        do 4 j=jval_(ien),jx
          feta(j,0)=zcoef*exp(-xsq(j)*zthta)
 4      continue
      endif

!.......................................................................
!     Loop over the Legendre index
!.......................................................................
      coschi=cos(chi)
      plegg(0)=one
      if(msxr .eq. 0) go to 600
      plegg(1)=coschi
      do 500 m=2,msxr
        plegg(m)=((2*m-1)*coschi*plegg(m-1)-(m-1)*plegg(m-2))/m
 500  continue
 600  continue

!.......................................................................
!     Compute Legendre coefficient for temp2=f (if lr in lrindx only)
!.......................................................................

      if (illcal .eq. 0) go to 700

      ione="no"
      if (thpol .gt. pi .or. thpol .lt. zero) call tdwrng(4)
      do 10 l=1,lz
        lll=l
        !write(*,*)'tdsxr  z,sax=',lr_,z(l,lr_)-sax
        if (abs(z(l,lr_)-sax) .lt. 5.e-11) go to 11 !-YuP: see line 72
!        if (z(l,lr_) .eq. sax) go to 11 !-YuP: sometimes fails
        if (z(l,lr_) .gt. sax) go to 12
 10   continue
      call tdwrng(5)
 11   ione="yes"
 12   continue

!......................................................................
!     Copy the distribution into temp3 (for cfpleg)
!......................................................................

      call dcopy(iyjx2,temp2(0:iyjx2-1,0),1,temp3(0:iyjx2-1,0),1)

      do 20 m=0,msxr
        call bcast(feta(1:jx,m),zero,jx)
        call bcast(fetb(1:jx,m),zero,jx)

!.......................................................................
!     Obtain Legendre coefficients (in tam1, from cfpleg)  at z(lll,lr_)
!.......................................................................

        jzval=jval_(ien)
        call cfpleg(m,lll,jzval)
        call dcopy(jx,tam1,1,feta(1:jx,m),1)

        if (ione .eq. "yes") go to 30

!.......................................................................
!     Obtain Legendre coefficient at z(lll-1,lr_)
!.......................................................................
        lm1=lll-1
        jzval=jval_(ien)
        call cfpleg(m,lm1,jzval)
        call dcopy(jx,tam1,1,fetb(1:jx,m),1)

!BH110331 Checking feta/b calc in accord with msxr setting.
!        write(*,*)'tdsxr: m,feta(jx/2:jx/2+5,m),fetb(jx/2:jx/2+5,m)=',
!     &       m,feta(jx/2:jx/2+5,m),fetb(jx/2:jx/2+5,m)

!.......................................................................
!     Interpolate to obtain Legendre coefficients at sax
!.......................................................................
        dz1=z(lll,lr_)-sax
        dz2=sax-z(lll-1,lr_)
        dzz=z(lll,lr_)-z(lll-1,lr_)
        dz1=dz1/dzz
        dz2=dz2/dzz
        do 23 j=jval_(ien),jx
          feta(j,m)=feta(j,m)*dz2+fetb(j,m)*dz1
 23     continue
 30     continue


!      write(*,*) 'tdsxr:m,(feta(j,m),j=jval_(ien),jval_(ien)+10)',
!     +     m,(feta(j,m),j=jval_(ien),jval_(ien)+10)

 20   continue

 700  continue ! skipping handle

      edotei=zero
      edotee=zero

!.......................................................................
!     do the velocity integrations
!.......................................................................

!     But first, put elwert factor adjustment for Z.eq.1 into tam6:
!     (beta0 = pre-collision elec momentum/energy in relativistic units,
!      beta1 = post-collision momentum/energy)
      call bcast(tam6,zero,jx)
      do 601 j=jval_(ien),jx
!YuP         beta0=x(j)*cnormi*gammi(j) !cnormi=0 when relativ.eq."disabled"
         beta0=x(j)*gammi(j)/cnorm !YuP[07-2016] cannot have beta0=0, see below
         e1=gamma(j)-enk(ien)
         e12=e1*e1
         p12=e12-1.d0
         if (p12.le.zero) goto 601
         p1=sqrt(p12)
         beta1=p1/e1
         ff0=-2.*pi/(137*beta0)
         ff1=-2.*pi/(137*beta1)
         zz=zeff(ll)
         tam6(j)=(1.-exp(ff0))*(1.-exp(zz*ff1))/ &
                 ((1.-exp(zz*ff0))*(1.-exp(ff1)))
 601  continue

      call bcast(tam5,zero,jx)
      do 701 m=0,msxr
        do 41 j=jval_(ien),jx
          tam5(j)=feta(j,m)*sigsxr(j,m,ien,1)*plegg(m)+tam5(j)
 41     continue
 701  continue

      do 42 j=jval_(ien),jx
        edotei=edotei+tam5(j)*tam6(j)
 42   continue

      call bcast(tam5,zero,jx)
      do 702 m=0,msxr
        do 43 j=jval_(ien),jx
          tam5(j)=feta(j,m)*sigsxr(j,m,ien,2)*plegg(m)+tam5(j)
 43    continue
 702  continue
      do 44 j=jval_(ien),jx
        edotee=edotee+tam5(j)
 44   continue
!MG added 11/13/2017
      if(allocated(plegg) &
       .and. n.eq.nstop .and. ien.eq.nen) deallocate (plegg)
      ! YuP[2017-12-04] Corrected, added n.eq.nstop. .and. ien.eq.nen
      ! plegg is defined once and then re-used many times.
!MG end  added 11/13/2017

      return
      end
end module tdsxr_mod
