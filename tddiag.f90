module tddiag_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdboothi_mod, only : tdboothi
  use tdbootst_mod, only : tdbootst

  !---END USE

!
!

contains

      subroutine tddiag
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This routine computes radial diagnostics (total rf power
!     and current)
!.......................................................................


!.......................................................................
!     Compute the total current in AMPS, (currza).
!.......................................................................

      if (cqlpmod .eq. "enabled") print *," WARNING in tddiag: routine", &
        " not yet ready for CQLP"


!.......................................................................
!     Bootstrap current calculation:
!     If jhirsh=88, use Hirshman '88 banana regime calc,
!     If jhirsh=99, use Sauter et al '99, in banana regime, calc,
!     If jhirsh=0 , use Hinton and Haseltine multi-regime,
!     high-aspect ratio formula.
!     (default=0, but 99 is best for comparison with cql3d).
!.......................................................................

      if (jhirsh.ne.0) then
         call tdboothi
      else
         call tdbootst
      endif

      psynct=0.0
      do 1 k=1,ngen
        currza(k)=0.
        rfpwrt(k)=0.
        gkpwrt(k)=0.
        pegyt(k)=0.0
        pplosst(k)=0.0
        wparzt(k)=0.
        wperpzt(k)=0.

!..................................................................
!     Compute the total RF power absorbed.
!     Also compute the integrals of the perpendicular and
!     parallel energy (units changed from ergs to joules).
!..................................................................

        do 2 ll=1,lrzmax
          currza(k)=currza(k)+darea(ll)*currz(k,ll)
          rfpwrt(k)=rfpwrt(k)+dvol(ll)*rfpwrz(k,ll)
          gkpwrt(k)=gkpwrt(k)+dvol(ll)*gkpwrz(k,ll)
          pegyt(k)=pegyt(k)+dvol(ll)*pegyz(k,ll)
          pplosst(k)=pplosst(k)+dvol(ll)*pplossz(k,ll)
          wparzt(k)=wparzt(k)+dvol(ll)*wpar(k,ll)*1.e-7
          wperpzt(k)=wperpzt(k)+dvol(ll)*wperp(k,ll)*1.e-7
          if(k.eq.kelecg)  psynct=psynct+dvol(ll)*psyncz(ll)
 2      continue
 1    continue

!..................................................................
!     Compute the total current (Amps) - currtza
!     Compute the total current with electron compensating effects
!     Compute the current integral up to rz(l); currtzi(l)
!     Compute the compensated current up to rz(l); currtpzi(l)
!
!     ccurtor(lr_) is the cumulative toroidal current, integrating
!                   curtor in poloidal cross-section (Amps),
!                   accounting for pol variation of tor current.
!     ccurpol(lr_) is the cumulative poloidal current, integrating
!                   currpol over area in toroidal cross-section
!                   at the outer equatorial plane (Amps).
!..................................................................

!     Presently using electron component of bootstrap in
!     total current.  Need to adjust this (BobH, 990821).

      do k=1,2
      do kk=1,2
         bscurmi(0,k,kk)=0.
      enddo
      enddo
      currtzi(0)=0.
      currtpzi(0)=0.
      totcurzi(0)=0.
      ccurtor(0)=0.
      ccurpol(0)=0.
      do 4 ll=1,lrzmax
        currtzi(ll)=currtzi(ll-1)+darea(ll)*currtz(ll)
        currtpzi(ll)=currtpzi(ll-1)+darea(ll)*currtpz(ll)
        do k=1,2
        do kk=1,2
           bscurmi(ll,k,kk)=bscurmi(ll-1,k,kk)+darea(ll)*bscurm(ll,k,kk)
        enddo
        enddo
        totcurz(ll)=bscurm(ll,1,1)+currtpz(ll)
        totcurzi(ll)=totcurzi(ll-1)+darea(ll)*totcurz(ll)
        ccurtor(ll)=ccurtor(ll-1)+darea(ll)*curtor(ll)* &
             rpcon(ll)*onovrp(2,ll)/onovrp(1,ll)
        ccurpol(ll)=ccurpol(ll-1)+twopi*rpconz(ll)* &
             (rpmconz(ll)-rpmconz(ll-1))*curpol(ll)
 4    continue
      do k=1,2
      do kk=1,2
         bscurma(k,kk)=bscurmi(lrzmax,k,kk)
      enddo
      enddo
      currtza=currtzi(lrzmax)
      currtpza=currtpzi(lrzmax)
      totcurza=currtpza+bscurma(1,1)
!
!..................................................................
!     Compute the total source power integrated over space and summed
!     over all beam species
!     Compute the inductance, li. (Needs work for eqmod.ne."enabled":
!        expressions for bpolsqaz and bpolsqlm).
!..................................................................

      do k=1,ngen
        sorpw_rfi(k,0)=0.0
        sorpw_nbii(k,0)=0.0
        do ll=1,lrzmax
           sorpw_rfi(k,ll)=sorpw_rfi(k,ll-1)+sorpw_rf(k,ll)*dvol(ll)
           sorpw_nbii(k,ll)=sorpw_nbii(k,ll-1)+sorpw_nbi(k,ll)*dvol(ll)
        enddo
      enddo

      sorpwti(0)=0.0
      do 11 ll=1,lrzmax
        sorpwti(ll)=sorpwti(ll-1)+sorpwt(ll)*dvol(ll)
 11   continue

      sorpwtza=sorpwti(lrzmax)

      volume=0.
      do 14 ll=1,lrzmax
         volume=volume+dvol(ll)
 14   continue

      li=0.
      if (eqmod.eq."enabled") then
         do 15 ll=1,lrzmax
            li=li+bpolsqaz(ll)*dvol(ll)
 15      continue
         li=li/volume/bpolsqlm
      endif

!..................................................................
!     Compute beam current drive figure of merit, eta, for both
!     cases (with compensating electrons and without).
!     =N(10**14)*radmaj(10**2)*currtza/sorpwtza
!
!     edenlavg=line average density.
!..................................................................

      eden=0.
      edenlavg=0.0
      etemp=0.
      ethtemp=0.
      edntmp=0.
      pden=0.
      pdntmp=0.
      do 30 k=1,ntotal
        if ((k.eq.kelecg.and.kelecm.eq.0) .or. k.eq.kelecm) then
          if (k.eq.kelecg.and.(colmodl.eq.1.or.colmodl.eq.3)) go to 21
          do 20 ll=1,lrzmax
            eden=eden+reden(k,ll)*dvol(ll)
            etemp=etemp+energy(k,ll)*dvol(ll)
            ethtemp=ethtemp+temp(k,ll)*dvol(ll)
            edntmp=edntmp+reden(k,ll)*dvol(ll)*energy(k,ll)
 20       continue
          edenlavg=(rpcon(1)-rmcon(1))*reden(k,0)
          do 22 ll=2,lrzmax
            edenlavg=edenlavg+(rpcon(ll)-rpcon(ll-1)+rmcon(ll-1) &
              -rmcon(ll))*0.5*(reden(k,ll)+reden(k,ll-1))
 22       continue
 21       continue
        elseif (k.gt.ngen) then
          do 40 ll=1,lrzmax
            pden=pden+reden(k,ll)*dvol(ll)
            pdntmp=pdntmp+reden(k,ll)*dvol(ll)*energy(k,ll)
 40       continue
        endif
 30   continue
      edntmp=edntmp*2./3./eden
      if (pden.ne.zero) then
        pdntmp=pdntmp*2./3./pden
        pden=pden/volume
      endif
      eden=eden/volume
      edenlavg=edenlavg/(rpcon(lrzmax)-rmcon(lrzmax))
      etemp=(2./3.)*etemp/volume
      ethtemp=ethtemp/volume
      coef=eden/1.e+14*radmaj*.01/(sorpwtza+em90)

!..................................................................
!     Compute the figures of merit for c.d. efficiency
!..................................................................

      fom=coef*currtza
      fomp=coef*currtpza
      fompla=fomp*edenlavg/eden
      fomtot=coef*(currtpza+bscurma(1,1))

!..................................................................
!     Compute total plasma energy (joules) in each species
!..................................................................

      do 50  k=1,ntotal
        energyt(k)=0.0
        do 51  ll=1,lrzmax
          energyt(k)=energyt(k)+energy(k,ll)*reden(k,ll)*dvol(ll) &
            *1.6e-16
 51     continue
 50   continue

      return
      end
end module tddiag_mod
