module diaggnde2_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use cfpleg_mod, only : cfpleg
  use diagdenz_mod, only : diagdenz
  use diagcfac_mod, only : diagcfac
  use diagwrng_mod, only : diagwrng
  use r8subs_mod, only : dcopy
  use r8subs_mod, only : dscal

  !---END USE

!
!

contains

      subroutine diaggnde2
      use param_mod
      use comm_mod
      use r8subs_mod, only : dscal, dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Same as subroutine diaggnde but use Legendre decomposition for
!     computing the density and parallel current, so that the density
!     is the same as the one used in computing the collision operator.
!     It also allows one to have about a constant startup density
!     for the CQLP case, when the y-mesh is not the same on each
!     l mesh point.
!
!     Main diagnostic routine; computes densities, currents,
!     energies, z=effective, J/P. Also redefies electron density
!     to maintain quasineutrality if necessary.
!..................................................................

!.......................................................................

!..................................................................
!l    1. Compute integrals for 0th and 1st moments using the same
!     Legendre decomposition as in computing the collision term
!     Compute the other terms including bounce-averaged term in the same
!     way as diaggnde.
!     tam11 will contain midplane density moment (using Legendre coeff.)
!     tam13 will contain parallel current - at midplane. (with Leg. coeff.)
!..................................................................

      ilegen=1
      if (cqlpmod .eq. "enabled") ilegen=l_

      currmt(l_)=0.
      currmtp(l_)=0.
      if (l_ .eq. lmdpln_) then
        currt(lr_)=0.
        currtp(lr_)=0.
      endif

      do 90 k=1,ngen
        call bcast(tam2,zero,jx)
        call bcast(tam4,zero,jx)
        call bcast(tam5,zero,jx)

!     compute Legendre decomposition from temp3 to tam1
!     m=0 contribution
        call dcopy(iyjx2,f(0:iyjx2-1,0,k,l_),1,temp3(0:iyjx2-1,0),1)
        call cfpleg(0,ilegen,1)
        do 10 j=1,jx
          tam11(j)=tam1(j)*4.*pi
 10     continue
        if (mx .eq. 0) go to 11
!     m=1 contribution
        call cfpleg(1,ilegen,1)
        do 12 j=1,jx
          tam13(j)=tam1(j)*4.*pi/3.
 12     continue
 11     continue

!     Compute other terms integrating explicitly over theta
        if (l_ .eq. lmdpln_) then
          do 20 i=1,iy
            do 21 j=1,jx
!     include vptb=|cos(th0)| * tau
              tam2(j)=tam2(j)+f(i,j,k,l_)*cynt2(i,l_)* &
                abs(coss(i,lmdpln_))*tau(i,lr_)
              tam4(j)=tam4(j)+f(i,j,k,l_)*cynt2(i,l_)* &
                abs(coss(i,lmdpln_))*tau(i,lr_)* &
                (1.-sinn(i,l_)**2*psiba(i,lr_))
              tam5(j)=tam5(j)+f(i,j,k,l_)*cynt2(i,l_)* &
                abs(coss(i,lmdpln_))*tau(i,lr_)* &
                (sinn(i,l_)**2*psiba(i,lr_))
 21         continue
 20       continue
        endif

        gn=0.
        en=0.
        hn=0.
        sn=0.
        cn=0.
        wpar_=0.
        wperp_=0.
        do 40 j=1,jx

!..................................................................
!     midplane density
!..................................................................

          gn=gn+tam11(j)*cint2(j)

!..................................................................
!     miplane current
!     currvs contains partial sums of the current; cn is the full
!     integral.
!..................................................................


          currv(j,k,l_)=tam13(j)*x(j)*gammi(j)*cint2(j)*dxi(j)
          if (j.eq.1) then
            currvs(j,k)=dx(1)*currv(j,k,l_)
          else
            currvs(j,k)=currvs(j-1,k)+currv(j,k,l_)*dx(j)
          endif
          cn=cn+currv(j,k,l_)*dx(j)

!..................................................................
!     field line density (particles/cm**2)
!..................................................................

          hn=hn+tam2(j)*cint2(j)

!..................................................................
!     Flux surface averaged energy (per particle)
!..................................................................

          sn=sn+tam2(j)*cint2(j)*tcsgm1(j)

!..................................................................
!     Midplane average energy per particle
!..................................................................

          en=en+tam11(j)*cint2(j)*tcsgm1(j)

!..................................................................
!     flux surface averaged parallel and perpendicular energy density
!..................................................................

          wpar_=wpar_+tam4(j)*cint2(j)*tcsgm1(j)
          wperp_=wperp_+tam5(j)*cint2(j)*tcsgm1(j)

 40     continue

!..................................................................
!     If density computed is negative call exit
!     Energies are in Kev or Kev/cm**3 (after multiplication by
!     fions)
!..................................................................

        if (gn.le.0.)  call diagwrng(10)
        energym(k,l_)=en/gn*fions(k)
        enrgypa(k,ls_)=en/gn*fions(k)
        if (l_ .eq. lmdpln_) then
          energy(k,lr_)=sn/hn*fions(k)
!BH080502          wpar(k,lr_)=wpar_*fions(k)/zmaxpsi(lr_)*ergtkev
!BH080502          wperp(k,lr_)=wperp_*fions(k)/zmaxpsi(lr_)*ergtkev
          wpar(k,lr_)=wpar_*fions(k)/hn
          wperp(k,lr_)=wperp_*fions(k)/hn
        endif

!..................................................................
!     At timet=0. scale the
!     distribution function to desired initial density.
!     NOTE: The momentum mesh, x, extends from 0. to 1.  The maximum
!     actual value of momentum/rest mass  or speed is vnorm.
!     A non-normalized mesh would extend from 0. to vnorm.
!     The distribution function in the code, f, will thus differ
!     from the cgs f (f_cgs) by the following f/vnorm**3=f_cgs
!     In what follows we scale f so that its initial MIDPLANE
!     density is reden.
!     Thus if f is integrated out using the code normalized velocity
!     mesh the integral will be density reden at timet 0.
!..................................................................

        if (n.eq.0) then
!%OS
          print *,'diaggnde2: l_= ',l_,'  gn= ',gn
!%OS
          gni=1./gn
          hni=0.0
          if (l_ .eq. lmdpln_) hni=1./hn
          if (l_ .eq. lmdpln_) hnis(k,lr_)=hni

!..................................................................
!     This assures the midplane density of reden for f
!..................................................................

          zfact=hni*zmaxpsi(lr_)*reden(k,lr_)
          if (cqlpmod .eq. "enabled") zfact=gni*denpar(k,ls_)
!BH070414  if (.not. nlrestrt) call dscal(iyjx2,zfact,f(0,0,k,l_),1)
          if (nlrestrt.eq."disabled") &
               call dscal(iyjx2,zfact,f(0:iyjx2-1,0,k,l_),1)

!..................................................................
!     xlndn is the field line density (particles/cm**2)
!     xlndn/zmaxpsi(lr_) will equal the flux surface averaged (volume
!     averaged) density.
!..................................................................

          if (l_ .eq. lmdpln_) then
            xlndn(k,lr_)=reden(k,lr_)*zmaxpsi(lr_)
!BH080502            wperp(k,lr_)=wperp(k,lr_)*reden(k,lr_)*hni*zmaxpsi(lr_)
!BH080502            wpar(k,lr_)=wpar(k,lr_)*reden(k,lr_)*hni*zmaxpsi(lr_)
            wpar(k,lr_)=wpar_*fions(k)/hn
            wperp(k,lr_)=wperp_*fions(k)/hn
          endif

!..................................................................
!     currm is the non explicitly s dependent parallel current.
!     In CQL3D, this is equal to the current at the midplane.
!     note:
!     current(z)=currm(k,l_)*bbpsi(z,lr__)
!     currm(k,lmdpln_)=parallel current at midplane
!
!     curr is the flux surface averaged current.
!     To wit: curr*dA=the total number of particles/sec passing through
!     an area element associated with the cross sectional area of the
!     flux surface in question. Once again curr and currm are
!     PARALLEL current.
!..................................................................


          faccur=reden(k,lr_)*hni*zmaxpsi(lr_)*vnorm*bnumb(k)*charge
          if (cqlpmod .eq. "enabled") faccur=denpar(k,ls_)*vnorm* &
            bnumb(k)*charge*gni
          currm(k,l_)=faccur*cn

!..................................................................
!     psifct will contain the factor needed to modify midplane current
!     to flux surface averaged current.
!..................................................................

!..................................................................
!     old bug next
!     psifct=radmaj/(radmaj**2-rovera(lr_)**2*radmin**2)**.5
!     for eqmod=disabled: psifct=sqrt((1+eps)/(1-eps))
!..................................................................

          psifct=psiovr(lr_)/onovrp(1,lr_)
          if (l_ .eq. lmdpln_) curr(k,lr_)=currm(k,l_)*psifct
          fgni=faccur*psifct/3.e+9
          if (cqlpmod .eq. "enabled") fgni=faccur/3.e+9
          call dscal(jx,fgni,currv(1:jx,k,l_),1)
          call dscal(jx,fgni,currvs(1:jx,k),1)

!..................................................................
!     Standard, non-initialization time step logic follows..
!..................................................................

        else
          faccur=one_*vnorm*bnumb(k)*charge
          psifct=psiovr(lr_)/onovrp(1,lr_)

!..................................................................
!     xlndn(k,lr_) is the field line density
!     reden(k,lr_) is the flux averaged density
!     denpar(k,ls_) is the local in s density
!     currm(k,l_) is the local in s current (s:along the magnetic field)
!     curr(k,lr_) is the flux surface average current at l_=lmdpln_.
!     currt(lr_) is the current summed over all ionic general species.
!     currmt(l_) is the local in s current summed over all ionic g. species
!     currtp(lr_) adds in general electron contribution to currt(l_).
!     currmtp(l_) adds in general local in s electron contribution to currm
!..................................................................

          facpsi=faccur*psifct
          call dscal(jx,facpsi,currv(1:jx,k,l_),1)
          call dscal(jx,facpsi,currvs(1:jx,k),1)
          currm(k,l_)=faccur*cn
          currmtp(l_)=currmtp(l_)+currm(k,l_)
          denpar(k,ls_)=one_*gn
          if (sbdry.eq."periodic" .and. transp.eq."enabled") then
            if (ls_ .eq. 1) then
              denpar(k,lsmax+1)=denpar(k,1)
              enrgypa(k,lsmax+1)=enrgypa(k,1)
            else if (ls_ .eq. lsmax) then
              denpar(k,0)=denpar(k,lsmax)
              enrgypa(k,0)=enrgypa(k,lsmax)
            endif
          endif
          if (l_ .eq. lmdpln_) then
            curr(k,lr_)=currm(k,l_)*psifct
            currtp(lr_)=currtp(lr_)+curr(k,lr_)
            xlndn(k,lr_)=one_*hn
            reden(k,lr_)=one_*hn/zmaxpsi(lr_)
          endif

        endif

 90   continue

!..................................................................
!     Above end of do loop over general (time advanced) species
!     Below subtract electron general species electron current
!     contribution to currt(lr_), currmt(l_).
!..................................................................

      if (kelecg.gt.0) then
        if (l_ .eq. lmdpln_) currt(lr_)=currtp(lr_)-curr(kelecg,lr_)
        currmt(l_)=currmtp(l_)-currm(kelecg,l_)
      else
        if (l_ .eq. lmdpln_) currt(lr_)=currtp(lr_)
        currmt(l_)=currmtp(l_)
      endif

!.......................................................................
!     following variables depend only on midplane values,
!     i.e. are strictly radial quantities
!.......................................................................

      if (l_ .eq. lmdpln_) then

!..................................................................
!     Compute current drive efficiency; Amps/watt
!     bdre(lr_) - ions only;   bdrep(lr_) - ions and electrons (kelecg.ne.
!     sorpwt(lr_) is the total source power (RF+NBI for urf and fr module
!..................................................................

        sorpwt(lr_)=0.
        do 85 kk=1,ngen
          sorpwt(lr_)=sorpwt(lr_)+sorpw_rf(kk,lr_)+sorpw_nbi(kk,lr_) !-YuP
 85     continue

!..................................................................
!     lrzmax.gt.0 flags a multi-flux surface (CQL3D) run.
!     Define dA/DV
!..................................................................

        if (lrzmax.gt.1.and.n.gt.0) then
          areaovol=darea(lr_)/dvol(lr_)
        else
          if (eqmod.eq."enabled") then
            areaovol=1./twopi*onovrp(1,lr_)
          else
            areaovol=1./twopi/radmaj
          endif
        endif

!..................................................................
!     sorpwt(lr_) is the combined nbi and rf power (from urf module)
!..................................................................

        if (n.gt.0) then
          if (sorpwt(lr_).ge.1.e-25) then
            bdre(lr_)=areaovol*currt(lr_)/(sorpwt(lr_)*3.e+9)
            bdrep(lr_)=areaovol*currtp(lr_)/(sorpwt(lr_)*3.e+9)
          else
            bdre(lr_)=0.
            bdrep(lr_)=0.
          endif
        endif

!......................................................................
!     Compute Z-effective
!......................................................................

      zeff(lr_)=0.
      zeff1=0.
      xq=0.
!$$$        if (izeff.eq."ion") then
!$$$          k1=ngen+1
!$$$        else
!$$$          k1=1
!$$$        endif
!$$$        do 80 k=k1,ntotal
!$$$          if (k.eq.kelecg .or. k.eq.kelecm) goto 80
!$$$cBobH990128          if (k.eq.izeff) goto 80
!$$$          xq=xq+1
!$$$          zeff(lr_)=zeff(lr_)+bnumb(k)**2*reden(k,lr_)
!$$$          zeff4(lr_)=bnumb(k)**4*reden(k,lr_)+zeff4(lr_)
!$$$          zeff1=zeff1+bnumb(k)*reden(k,lr_)
!$$$ 80     continue

      if (izeff.eq."backgrnd") then
        do 80 k=1,ntotal
           if (k.eq.kelecg .or. k.eq.kelecm) goto 80
           xq=xq+1
           zeff(lr_)=zeff(lr_)+bnumb(k)**2*reden(k,lr_)
           zeff4(lr_)=bnumb(k)**4*reden(k,lr_)+zeff4(lr_)
           zeff1=zeff1+bnumb(k)*reden(k,lr_)
 80     continue

      elseif (izeff.eq."ion") then
         do kk=1,nionm
            k=kionm(kk)
            xq=xq+1
            zeff(lr_)=zeff(lr_)+bnumb(k)**2*reden(k,lr_)
            zeff4(lr_)=bnumb(k)**4*reden(k,lr_)+zeff4(lr_)
            zeff1=zeff1+bnumb(k)*reden(k,lr_)
         enddo

      else
         write(*,*) 'Problem with izeff specification'
         stop
      endif

      zeff4(lr_)=zeff4(lr_)/xq
      zeff(lr_)=zeff(lr_)/zeff1

!.......................................................................
!     end of strictly radial dependent variables
!.......................................................................

      endif

!.......................................................................
!     Compute the analytic (Cordey-Start) electron current degradation
!     factor for NBI type calculations. This assumes that electrons
!     are not a general species. If they are, the electron contribution
!     was computed above.
!.......................................................................

      if (kelecg.eq.0.and. eleccomp.eq."enabled") then
        factor=diagcfac(eps(lr_),zeff(lr_))
        if (l_ .eq. lmdpln_) currtp(lr_)=currt(lr_)*factor
        currmtp(l_)=currmt(l_)*factor
      endif

      if (l_ .ne. lmdpln_) go to 99

!..................................................................
!     Compute J/P which includes ion and electron currents
!..................................................................

      if (sorpwt(lr_).ge.1.e-25) then
        bdrep(lr_)=areaovol*currtp(lr_)/(sorpwt(lr_)*3.e+9)
      else
        bdrep(lr_)=0.
      endif

!..................................................................
!     Redefine electron density to maintain quasineutrality.
!..................................................................

      s=0.
      if (kelecm.eq.0 .or. kelecg.ne.0) goto 72
      if (qsineut .ne. "disabled") then
        if (qsineut.eq."maxwel") goto 71
        do 74 k=1,ngen
          s=s+xlndn(k,lr_)*zmaxpsii(lr_)*bnumb(k)
 74     continue
 71     do 70 k=ngen+1,ntotal
          if (k.eq.kelecm) goto 70
          s=s+reden(k,lr_)*bnumb(k)
 70     continue
        reden(kelecm,lr_)=s
        locquas="disabled"
      endif
 72   continue

!..................................................................
!     determine density as a function of poloidal angle.
!..................................................................

      call diagdenz

 99   continue
!
      return
      end
end module diaggnde2_mod
