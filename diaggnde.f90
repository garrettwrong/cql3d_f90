module diaggnde_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use cqlcomm_mod
  use diagdenz_mod, only : diagdenz
  use diagcfac_mod, only : diagcfac
  use diaggnde2_mod, only : diaggnde2
  use diagwrng_mod, only : diagwrng
  use iso_c_binding, only : c_double
  use param_mod
  use r8subs_mod, only : dscal
  use soucrit_mod, only : soucrit

  !---END USE

contains

  subroutine diaggnde
    !implicit integer (i-n), real(c_double) (a-h,o-z)
    implicit none
    real(c_double) :: areaovol
    real(c_double) :: cn
    real(c_double) :: en
    real(c_double) :: faccur
    real(c_double) :: facpsi
    real(c_double) :: factor
    real(c_double) :: fgni
    real(c_double) :: gn
    real(c_double) :: gni
    real(c_double) :: hn
    real(c_double) :: hni
    real(c_double) :: psifct
    real(c_double) :: s
    real(c_double) :: sn
    real(c_double) :: wpar_
    real(c_double) :: wperp_
    real(c_double) :: xq
    real(c_double) :: zeff1
    real(c_double) :: zfact

    integer :: i
    integer :: j
    integer :: k
    integer :: kk
    save

!..................................................................
!     Main diagnostic routine; computes densities, currents,
!     energies, z=effective, J/P. Also redefies electron density
!     to maintain quasineutrality if necessary.
!..................................................................

!MPIINSERT_INCLUDE

!.......................................................................

      if (oldiag .ne. "enabled") then
!     use Legendre decomposition for computing the integrals
        call diaggnde2
        return
      endif

!..................................................................
!     Compute integrals for various moments
!     tam1 will contain midplane density moment
!     tam2 will contain field line (flux surface averaged) density
!     tam3 will contain parallel current density - at midplane.
!     tam4 related to parallel energy calc
!     tam5 reltaed to perpendicular energy calc
!..................................................................

      currmt(l_)=0.
      currmtp(l_)=0.

!     For cqlpmod.ne."enabled", default lmdpln_ will take on values
!       lrz:1 and successive calls (for given time step).  l_ also
!       takes on these values, so following if statement is satisfied
!       at each call diaggnde.
!     For cqlpmod.eq."enabled", lmpdpln_ will equal 1 at each call,
!       whereas l_ varies lrors:1, representing distance along B-field.
!       Only a single flux surface is examined, lrz=1, and lr_ is a
!       chosen integer.
      if (l_ .eq. lmdpln_) then
        currt(lr_)=0.
        currtp(lr_)=0.
        curtor(lr_)=0.
        curpol(lr_)=0.
      endif
      do 90 k=1,ngen
        call bcast(tam1,zero,jx)
        call bcast(tam2,zero,jx)
        call bcast(tam3,zero,jx)
        call bcast(tam4,zero,jx)
        call bcast(tam5,zero,jx)
        do 10 j=1,jx
        do 20 i=1,iy
            !if(gone(i,j,k,l_).eq.zero)then !YuP[07-14-2014]
               ! Only count particles that are NOT in the loss cone
               ! (gone array is set to -1 in the loss cone).
               ! Another option would be simply to remove particles
               ! in the loss cone, but tests show that the solution
               ! can become unstable in such case.
               ! So, keep f() non-zero in the loss cone,
               ! but do not add such f() to the integrals/sums over d3v.
              tam1(j)=tam1(j)+f(i,j,k,l_)*cynt2(i,l_)
              tam3(j)=tam3(j)+f(i,j,k,l_)*cynt2(i,l_)*coss(i,l_)
!     include vptb=|cos(th0)| * tau
            tam2(j)=tam2(j)+f(i,j,k,l_)*cynt2(i,l_)* &
              abs(coss(i,lmdpln_))*tau(i,lr_)
            tam4(j)=tam4(j)+f(i,j,k,l_)*cynt2(i,l_)* &
              abs(coss(i,lmdpln_))*tau(i,lr_)*(1.-sinn(i,l_)**2* &
              psiba(i,lr_))
!BH080502            tam5(j)=tam5(j)+f(i,j,k,l_)*cynt2(i,l_)*vptb(i,lr_)*
!BH080502  Error:  only need one vptb= abs(coss)*tau factor
            tam5(j)=tam5(j)+f(i,j,k,l_)*cynt2(i,l_)* &
              abs(coss(i,lmdpln_))*tau(i,lr_)*(sinn(i,l_)**2* &
              psiba(i,lr_))
            !endif
 20     continue ! i
 10     continue ! j

!     write(*,*) 'diaggnde:l_,(tau(i,l_),i=1,iy)', l_,(tau(i,l_),i=1,iy)

        denra(k,l_)=0.
        curra(k,l_)=0.
!-----------------------------
!   denra is the flux surface average runaway density.
!   curra is the flux surface area average parallel runaway current
!   density.
!   fdenra is fraction of density that is runaway.
!   fcurra is fraction of current due to runaway.
!-----------------------------
!
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

          gn=gn+tam1(j)*cint2(j) !== n0 midplane density
          !except, may not be properly scaled yet, at t=0, because
          !distr.func. which is set in finit (for ZOW or HYBRID_FOW)
          !is not yet properly scaled. This is done further in this
          !subroutine.

!..................................................................
!     miplane current density
!     currv is integrand in velocity integral for parallel current,
!       at the midplane.
!     currvs contains partial sums of the currv, up to velocity x.
!     cn is the full integral.
!..................................................................

          currv(j,k,l_)=tam3(j)*x(j)*gammi(j)*cint2(j)*dxi(j)
          if (j.eq.1) then
            currvs(j,k)=dx(1)*currv(j,k,l_)
          else
            currvs(j,k)=currvs(j-1,k)+currv(j,k,l_)*dx(j)
          endif
          cn=cn+currv(j,k,l_)*dx(j)

!..................................................................
!     field line density (particles/cm**2)
!     Re eqsym: tam2 has tau factor, to be divided below by zmaxpsi.
!     tau is v*tau_B, and it's over 1/2 orbit if eqsym.ne.none
!..................................................................
          hn=hn+tam2(j)*cint2(j) !actually 1/2 line-density if eqsym.ne.none

!..................................................................
!     Flux surface averaged energy (per particle)
!..................................................................
          sn=sn+tam2(j)*cint2(j)*tcsgm1(j)

!..................................................................
!     Midplane average energy per particle
!..................................................................

          en=en+tam1(j)*cint2(j)*tcsgm1(j)

!..................................................................
!     flux surface averaged parallel and perpendicular energy density
!     Re eqsym:  wpar_/wperp_ has psiba=int_over_l{dtau/bbpsi} factor
!......................................................................

          wpar_=wpar_+tam4(j)*cint2(j)*tcsgm1(j)
          wperp_=wperp_+tam5(j)*cint2(j)*tcsgm1(j)

 40     continue ! j
!-------------------


!     call routine to determine jxcrit (presently determined as
!       index for amin1(3.*clight,ucrit).
        call soucrit

        if (jxcrit(k,lr_).ge.jx) go to 50
        do j=jxcrit(k,lr_),jx
           denra(k,l_)=denra(k,l_)+tam2(j)*cint2(j)
           curra(k,l_)=curra(k,l_)+currv(j,k,l_)*dx(j)
        enddo
 50     continue
        fdenra(k,l_)=denra(k,l_)/hn
        fcurra(k,l_)=curra(k,l_)/cn

!..................................................................
!     If density computed is negative call exit
!     Energies are in Kev or Kev/cm**3 (after multiplication by
!     fions)
!..................................................................

        if (gn.le.0.d0)  call diagwrng(10)
        energym(k,l_)=en/gn*fions(k)  ! at midplane
        enrgypa(k,ls_)=en/gn*fions(k) ! at midplane
        if (l_ .eq. lmdpln_) then
          energy(k,lr_)=sn/hn*fions(k) ! FSA for ZOW only
          ! at n>0, this definition will be over-written by FOW/ZOW
          ! universal procedure: through reconstruction of local f(R,Z).
!BH080502          wpar(k,lr_)=wpar_*fions(k)/zmaxpsi(lr_)*ergtkev
!BH080502          wperp(k,lr_)=wperp_*fions(k)/zmaxpsi(lr_)*ergtkev
!BH180531:  Following wpar/wperp (diaggnde2 also) needs clarification..
          wpar(k,lr_)=wpar_*fions(k)/hn
          wperp(k,lr_)=wperp_*fions(k)/hn
        endif
!MPIINSERT_IF_RANK_EQ_0
!        WRITE(*,*) &
!       'diaggnde: k,lr_,l_,en,hn,fions(k),energym(k,l_),energy(k,lr_)', &
!                  k,lr_,l_,en,hn,fions(k),energym(k,l_),energy(k,lr_)
!MPIINSERT_ENDIF_RANK

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
!          print *,'diaggnde:  l_= ',l_,'  gn,hn= ',gn,hn
!%OS
          gni=1./gn
          hni=0.0
          !XXXXXX infinity?  YuP: hn and gn are never 0 - they are integrals over distr.func.
          if (l_ .eq. lmdpln_) hni=1./hn
          if (l_ .eq. lmdpln_) hnis(k,lr_)=hni

!..................................................................
!     This assures the midplane density of reden for f
!..................................................................

          zfact=hni*zmaxpsi(lr_)*reden(k,lr_)
          if (cqlpmod .eq. "enabled") zfact=gni*denpar(k,ls_)

!MPIINSERT_IF_RANK_EQ_0
          WRITE(*,*)'----------------------- lr_===', lr_
          WRITE(*,'(a,i4,3e12.4)') &
           'diaggnde_n=0 lr_, reden, gn, sum_ij(gone)', &
                         lr_,reden(k,lr_),gn,sum(gone(1:iy,1:jx,k,lr_))
!MPIINSERT_ENDIF_RANK

          if ( (nlrestrt.ne."disabled") &
               .or. (fpld(1,1).eq.-1.0)) then
             continue
          else
             call dscal(iyjx2,zfact,f(0:iy+1,0:jx+1,k,l_),1)
             !XXXX, why not ?  YuP:Could use f(:,:,k,l_)=f(:,:,k,l_)*zfact Not sure which way is faster
             !f(:,0,k,l_) = f(:,0,k,l_)*zfact
          endif
!          write(*,*)'diaggnde: k,lr_,zfact,reden(k,lr_)=',
!     +                         k,lr_,zfact,reden(k,lr_)

!..................................................................
!     xlndn is the field line density (particles/cm**2)
!     xlndn/zmaxpsi(lr_) will equal the flux surface averaged (volume
!     averaged) density.
!..................................................................

          if (l_ .eq. lmdpln_) then
            xlndn(k,lr_)=reden(k,lr_)*zmaxpsi(lr_) !n=0: an approximation
!BH080502            wperp(k,lr_)=wperp(k,lr_)*reden(k,lr_)*hni*zmaxpsi(lr_)
!BH080502            wpar(k,lr_)=wpar(k,lr_)*reden(k,lr_)*hni*zmaxpsi(lr_)
            wpar(k,lr_)=wpar_*fions(k)/hn
            wperp(k,lr_)=wperp_*fions(k)/hn
          endif

!..................................................................
!     currm(k,l_) is the non-explicitly s-dependent parallel current
!       density, for each general species k.
!     In CQL3D, this is equal to the current density at the outer
!        midplane due to each general species k.
!     note:
!     current(z)=currm(k,l_)*bbpsi(z,lr_)
!     currm(k,lmdpln_)=parallel current at outer midplane
!
!     curr is the flux surface toroidal-area-averaged current density.
!     To wit: curr*dA=the total number of particles/sec passing through
!     an area element associated with the cross sectional area of the
!     flux surface in question. Once again curr and currm are
!     PARALLEL current.  (In non-tokamak situations, it will make
!     more sense to consider Toroidal current times toroidal area.)
!..................................................................

          faccur=reden(k,lr_)*hni*zmaxpsi(lr_)*vnorm*bnumb(k)*charge
          if (cqlpmod .eq. "enabled") faccur=denpar(k,ls_)*vnorm* &
            bnumb(k)*charge*gni
          currm(k,l_)=faccur*cn
          curra(k,l_)=faccur*curra(k,l_)

!..................................................................
!     psifct will contain the factor needed to modify midplane current
!     to flux surface (toroidal-area-)averaged current.
!     Thus,
!       Avg curr = int{d(tor area) * parallel_curr density}/
!                  int{d(tor area)}
!       parallel_curr=parallel_cur_at_min_B_pt * (B/min_B)
!       gives
!      Avg curr = <B/(min_B * R)>/<1/R>, where < >= FSAvg.
!
!..................................................................

!..................................................................
!     old bug next
!     psifct=radmaj/(radmaj**2-rovera(lr_)**2*radmin**2)**.5
!     for eqmod=disabled: psifct=sqrt((1+eps)/(1-eps))
!..................................................................

          psifct=psiovr(lr_)/onovrp(1,lr_)
          if (l_ .eq. lmdpln_) then
             curr(k,lr_)=currm(k,l_)*psifct
             curra(k,l_)=curra(k,l_)*psifct
             denra(k,l_)=one_*denra(k,l_)/zmaxpsi(lr_)
          endif

          fgni=faccur*psifct/3.e+9
          if (cqlpmod .eq. "enabled") fgni=faccur/3.e+9
          call dscal(jx,fgni,currv(1:jx,k,l_),1)
          call dscal(jx,fgni,currvs(1:jx,k),1)

!..................................................................
!     Standard, non-initialization time step (n.ne.0) logic follows..
!..................................................................

        else                                  ! n.ne.0
          faccur=one_*vnorm*bnumb(k)*charge
          psifct=psiovr(lr_)/onovrp(1,lr_)

!..................................................................
!     xlndn(k,lr_) is the field line density
!     reden(k,lr_) is the flux averaged density
!     denpar(k,ls_) is the local in s density
!
!     currm(k,l_) is the local in s current (s:along the magnetic field)
!       (For CQL3D, l_=lmdpln_, this is current density at midplane.
!        For CQLP, l_=lmdpln_=1, is current density at the midplane.)
!     currmtp(l_) sums currm(k,l_) over general species k, including
!                the general local s=0 electron contribution.
!     currmt(l_) is the local s=0 sum of parallel curr. den.
!                over the ionic gen. species,
!                i.e., currmtp(l_) minus any electron general species.
!                contribution .
!     curr(k,lr_) is flux surface (area-)average current
!                density at l_=lmdpln_.
!     currtp(lr_) is sum over general species k of curr(k,lr_)
!     currt(lr_) is sum or flux surface (toroidal-area-)average
!                parallel curr. den. over the ionic gen. species,
!                i.e., currmtp(l_) minus any electron general species.
!                contribution.
!..................................................................

          facpsi=faccur*psifct
!         Scaling units of currv to be statamps/cm**2*vnorm,
!                          currvs to be statamps/cm**2
          call dscal(jx,facpsi,currv(1:jx,k,l_),1)
          call dscal(jx,facpsi,currvs(1:jx,k),1)
          currm(k,l_)=faccur*cn
          curra(k,l_)=faccur*curra(k,l_)
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
            curra(k,l_)=curra(k,l_)*psifct
            currtp(lr_)=currtp(lr_)+curr(k,lr_)
            xlndn(k,lr_)=one_*hn
            reden(k,lr_)=one_*hn/zmaxpsi(lr_)
            denra(k,l_)=one_*denra(k,l_)/zmaxpsi(lr_)
          endif

        endif  ! n.eq.0/else


!..................................................................
!     Breaking parallel current densities down into toroidal
!     and poloidal contributions (Amps/cm**2):
!     curtor(lr_) is parallel current density vector dotted into
!                  toroidal direction, evaluated at outer eq. plane.
!     curpol(lr_) is parallel current density vector dotted into
!                  poloidal direction, evaluated at outer eq. plane.
!..................................................................

        curtor(lr_)=curtor(lr_)+(fpsi(lr_)/rpcon(lr_))/ &
             bmidplne(lr_)*currm(k,l_)/3.e9
        curpol(lr_)=curpol(lr_)+bthr0(lr_)/ &
             bmidplne(lr_)*currm(k,l_)/3.e9

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

!.......................................................................
!     Compute Z-effective
!.......................................................................

      zeff(lr_)=0.
      zeff1=0.
      zeff4(lr_)=0.d0 !Yup[2014-05-27] Initialize to 0.
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

      endif  !l_.eq.lmdpln_

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
!     Save copy of pitch angle integrated f in tem5
!
      call bcast(tem5,zero,jx)
      do j=1,jx
         do i=1,iy
            tem5(j)=tem5(j)+f(i,j,1,l_)*cynt2(i,l_)*cint2(j)* &
                 vptb(i,lr_)/zmaxpsi(lr_)
         enddo
      enddo
      return
      end subroutine diaggnde


end module diaggnde_mod
