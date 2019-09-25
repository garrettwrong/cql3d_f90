module zfreya_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE


contains

      real(c_double) function bkef(vpar, tpar, emzpar)
      use frsubs_mod, only : asimp
      implicit none 
!
!     This routine evaluates the function Ke, which is related to the trans-
!     fer of parallel momentum from fast ions slowing down on electrons.
!     The function is defined by Callen et al., IAEA Tokyo, Vol. I, 645
!     (1974).
!
      real(c_double) :: vpar, tpar, emzpar ! in arg.list
      real(c_double) :: a,b,ep !local, input in called func.

      integer m,m4 !local, input in called func.
      real(c_double) :: zero !local
      !external bkefun ! Now in zfreya_mod module
      real(c_double) :: vcvo, tstcx, emzrat !in /gecom/
      common /gecom/ vcvo, tstcx, emzrat

      vcvo = vpar
      tstcx = tpar
      emzrat = emzpar
      a = 0.
      b = 1.
      ep = 0.1
      m4 = 3
      bkef = asimp(a,b,ep,m,m4,bkefun)
      return
      end function bkef


      real(c_double) function bkefun(y)
      implicit none 
      real(c_double) :: y ! input
      real(c_double) :: zero,v3,arg,pcx,b !local

      real(c_double) :: vcvo, tstcx, emzrat !in /gecom/
      common /gecom/ vcvo, tstcx, emzrat
      
      zero=0.d0
      bkefun = zero
      if(y.eq.zero) return
      v3 = vcvo**3
      arg = (1.+v3)/(y**3+v3)
      pcx = (arg)**(-tstcx/3.)
      b = (y**3*arg)**(emzrat/3.)
      bkefun = y**3*pcx*b/(y**3+v3)
      return
      end function bkefun

      subroutine bproc(atw,ebkev,ibion,ke,kj,mb,nj,nion, &
        ene,en,te,z,sb,nprim,sbcx,sbion)
      implicit none 
      save
!----------------------------------------------------------------
!     inputs
!----------------------------------------------------------------
!     atw(i)             : atomic weight of ion species i
!     ebkev(ib)          : maximum particle energy in beam ib
!     ibion              : ion index of beam species
!     ke                 : maximum number of beam components
!     kj                 : maximum number of radial points
!     mb                 : number of beams
!     nj                 : number of radial points
!     nion               : number of ions
!     ene                : electron density (#/cm**3)
!     en(j,k)            : ion density at r(j) for species k (#/cm**3)
!     te(j)              : electron temperature at r(j) (kev)
!     z(j,k)             : z at r(j) for ion species k
!     sb(j,ic,ib)        : hot ion source at r(j) , beam ib, component ic
!     (#/cm**3-sec)
!----------------------------------------------------------------
!     outputs
!----------------------------------------------------------------
!     sbcx(j,i)          : hot ion source at r(j) due to charge exchange
!     with primary ion species i   (#/cm**3-sec)
!     sbion(j)           : hot ion source at r(j) due to ionization
!     by electrons  (#/cm**3-sec)
!----------------------------------------------------------------
!
!     mean free paths for charge exchange, proton ionization and
!     electron ionization are calculated from the fitted results
!     of freeman and jones, clm-r137, culham (1974).
!
!
      integer ibion,ke,kj,mb,nj,nion,nprim ! input
      real(c_double) :: atw(*),ebkev(*),ene(*),en(kj,*),te(kj)
      real(c_double) :: sb(kj,ke,*),sbcx(kj,2),sbion(*)
      real(c_double) :: z(kj,*)
      real(c_double) :: alogt,teev,expo,sgvxne,ebev,vbeam,veli !local
      real(c_double) :: rpathei,rpath,e,aloge,sigi,sigcx,ebi !local
      real(c_double) :: rpathcx,rpathii,fractcx1,fractcx2 !local
      integer ii,ib,ic !local

      real(c_double) :: cfione(7),cfionp(7),rpathcxn(2)
      integer j,if7,k !local

      data cfione/-3.173850e+01,1.143818e+01, &
        -3.833998,7.046692e-01,-7.431486e-02,4.153749e-03, &
        -9.486967e-05/
      data cfionp/-4.203309e+01,3.557321, &
        -1.045134,3.139238e-01,-7.454475e-02,8.459113e-03, &
        -3.495444e-04/

      do 5000 j=1,nj
        sbcx(j,1)=0.
        sbcx(j,2)=0.
        sbion(j)=0.
        alogt=0.
        teev = 1.e3*te(j)
!990131        if(teev.gt.1.) alogt=alog(teev)
        if(teev.gt.1.) alogt=log(teev)
        if(teev.gt.1.e+05) alogt=11.51
        expo=cfione(7)
        do 2000 ii=1,6
          if7=7-ii
          expo=expo*alogt+cfione(if7)
 2000   continue
        sgvxne=exp(expo)*ene(j)
!cc
        do 4000 ib=1,mb
          do 4001 ic=1,3
            ebev=1.e3* ebkev(ib)/(ic*atw(ibion))
            vbeam=1.384e6*sqrt(ebev)
            veli=1./vbeam
            rpathei=sgvxne*veli
            rpath=rpathei
            rpathcxn(1)=0.
            rpathcxn(2)=0.
            do 3200 k=1,nion
              if(atw(k).gt.3.01) go to 2500
              e = 1.e3*ebkev(ib)/(ic*atw(ibion))
!990131              aloge=alog10(e)
              aloge=log10(e)
              sigcx=.6937e-14*(1.-.155*aloge)**2/ &
                (1.+.1112e-14*e**3.3)
              aloge=aloge*2.302585093-6.907755279
              if(aloge.le.-2.30258) go to 2220
              expo=cfionp(7)
              do 2200 ii=1,6
                if7=7-ii
                expo=expo*aloge+cfionp(if7)
 2200         continue
              sigi=exp(expo)
              go to 2240
 2220         continue
              sigi=0.
 2240         continue
!cc
              rpathcx = sigcx*en(j,k)
              rpathii = sigi*en(j,k)
              go to 3000
!------------------------------------------------------------
!     r.e. olson et al., phys. rev. let. 41, 163 (1978)
!     ionization by impact with impurity
!------------------------------------------------------------
 2500         continue
              ebi=ebkev(ib)/(ic*atw(ibion))
              rpathii=1.0e-17*en(j,k)*46.*z(j,k)* &
                (32.*z(j,k)/ebi)*(1.-exp(-ebi/(32.*z(j,k))))
              rpathcx=0.
!------------------------------------------------------------
 3000         continue
              rpath=rpath+rpathcx+rpathii
              if(k.le.nprim.and.k.le.2) rpathcxn(k)=rpathcx
 3200       continue
            fractcx1=rpathcxn(1)/rpath
            fractcx2=rpathcxn(2)/rpath
            sbcx(j,1)=sbcx(j,1) + fractcx1*sb(j,ic,ib)
            sbcx(j,2)=sbcx(j,2) + fractcx2*sb(j,ic,ib)
            sbion(j)=sbion(j) + (1.-fractcx1-fractcx2)*sb(j,ic,ib)
 4001     continue
 4000   continue
 5000 continue
!cc
      return
      end subroutine bproc
!
!
      real(c_double) function ceef(i,j,te)
      implicit none 
      save
!*****************************************************************
!--   e-impact excitation rate coefficient from vriens & smeets.
!     used only for excitations to final state with n > 2 (i.e. j > 3).
!
      integer i,j !input
      real(c_double) :: te !input
      real(c_double) :: ryd,ceef1 !local
      integer ms,mc,id !local
      parameter (ms=21,mc=35)
      integer ni,nj !local
      
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      
      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug
      
      real(c_double) :: en,dg,ae,be,de1,de2,ge1,ge2 !in /b4/
      common/b4/en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc), &
        de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
        
      data ryd/13.6/
      !YuP was ceef1(ni,nj,te)=1.6e-7*sqrt(te) &
      !YuP[2019-06-17] Why indexes are needed? Removed.
      ceef1=1.6e-7*sqrt(te) &
        *exp(-(ryd/te)*(1./DBLE(ni)**2-1./DBLE(nj)**2)) &
!990131        *(ae(ni,nj)*alog(.3*te/ryd+de2(ni,nj))+be(ni,nj))
        *(ae(ni,nj)*log(.3*te/ryd+de2(ni,nj))+be(ni,nj)) &
!990131        /(te+ge2(ni,nj)*alog(1.+(DBLE(ni)**3)*te/ryd))
        /(te+ge2(ni,nj)*log(1.+(DBLE(ni)**3)*te/ryd))
!
      if((i.ge.j) .or. (j.le.3)) then
        ihxbug = 8
        if(nouthx.gt.0) write(nouthx,3737) i,j
 3737   format(' ??? ceef error: i4= ',i3,'        j= ',i3)
        return
      endif
!
      ni=nfhx(i)
      nj=j-1
      ceef=ceef1 !(ni,nj,te) !YuP[2019-06-17] Why indexes are needed? Removed.
      id=2
      if(ni.eq.1)id=1
      return
      end function ceef
      
      
      subroutine chek(x,xmin,xmax,imin,imax,ick)
      implicit none 
      real(c_double) :: x,xmin,xmax !input
      integer imin,imax !input
      integer ick !output
      ick = 0
      if(imin.eq.0) go to 10
      if(x.ge.xmin) go to 10
      ick = 1
      return
 10   continue
      if(imax.eq.0) return
      if(x.gt.xmax) ick = 1
      return
      end subroutine chek
!
!
      real(c_double) function cief(i,te)
      implicit none 
!*****************************************************************
!--   e-impact ionization rate coefficient from vriens & smeets,
!     phys. rev. a 22, 940 (1980).
!     used only for 2p (i=3) or for n > 2 (i=n+1).
!
      integer i !input
      real(c_double) :: te !input
      integer ms,mc !local
      parameter (ms=21,mc=35)
      
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      
      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug
      
      real(c_double) :: en,dg,ae,be,de1,de2,ge1,ge2 !in /b4/
      common/b4/en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc), &
        de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)

      real(c_double) :: ent !local
      
      if(i.le.2) then
        ihxbug = 7
        if(nouthx.gt.0) write(nouthx,3737) i
 3737   format(' ??? cief error: i4= ',i3)
        return
      endif

      ent=en(i)/te !YuP[2019-06-17] Index i in ent(i) was used as in a statement function
      cief=9.56e-6*exp(-ent) &
        /(te*sqrt(te)*(ent**2.33+4.38*ent**1.72+1.32*ent))
      return
      end function cief

!Kinsey substituted crsecs version below (from onetwo)
!$$$!     ONETWO DIVERGENCE
!$$$      subroutine crsecs(atw, ebkev, ibion, ke, kz, mb, mfm1,
!$$$     *  nion, vbeam, zne, zni, zte, zzi,
!$$$     *  sgvxne, sgxn, sgxnmi,iexcit,
!$$$     *  dtemax,dnemax,dzemax)
!$$$      implicit none 
!$$$      save
!$$$!
!$$$!     ONETWO DIVERGENCE
!$$$
!$$$!     this subroutine calculates charge exchange, proton ionization, and
!$$$!     electron ionization cross sections from the fitted results of
!$$$!     freeman and jones, clm-r137, culham (1974) for iexcit=0.
!$$$!     But for iexcit=-1 W. Stearns simplified formula is used.
!$$$!
!$$$      dimension atw(*), ebkev(*), vbeam(ke,*),
!$$$     *  zne(*), zni(kz,*), zte(*), zzi(kz,*),
!$$$     *  sgvxne(*), sgxn(kz,ke,*), sgxnmi(ke,*)
!$$$ccc
!$$$      dimension cfione(7),cfionp(7)
!$$$ccc
!$$$      data cfione/-3.173850e+01,1.143818e+01,
!$$$     1  -3.833998,7.046692e-01,-7.431486e-02,4.153749e-03,
!$$$     2  -9.486967e-05/
!$$$      data cfionp/-4.203309e+01,3.557321,
!$$$     1  -1.045134,3.139238e-01,-7.454475e-02,8.459113e-03,
!$$$     2  -3.495444e-04/
!$$$ccc
!$$$!      write(*,*) 'cfione(1),cfione(2),cfionp(1),cfionp(1)',
!$$$!     +            cfione(1),cfione(2),cfionp(1),cfionp(1)
!$$$!     ONETWO DIVERGENCE
!$$$      if (iexcit.eq.-1) go to 100
!$$$      do 10 i=1,mfm1
!$$$        alogt=0.
!$$$        teev = 1.e3*zte(i)
!$$$!990131        if(teev.gt.1.) alogt=alog(teev)
!$$$        if(teev.gt.1.) alogt=log(teev)
!$$$        if(teev.gt.1.e+05) alogt=11.51
!$$$        expo=cfione(7)
!$$$        do 11 ii=1,6
!$$$          if7=7-ii
!$$$          expo=expo*alogt+cfione(if7)
!$$$ 11     continue
!$$$        sgvxne(i)=exp(expo)*zne(i)
!$$$ 10   continue
!$$$ccc
!$$$      do 120 ib=1,mb
!$$$        do 121 j=1,3
!$$$          veli=1./vbeam(j,ib)
!$$$          do 122 i=1,mfm1
!$$$            sgxn(i,j,ib)  = sgvxne(i)*veli
!$$$ 122      continue
!$$$ 121    continue
!$$$ 120  continue
!$$$      do 280 ib=1,mb
!$$$        do 281 ion=1,nion
!$$$          do 27 j=1,3
!$$$            if(atw(ion).gt.3.01) go to 25
!$$$            e = 1.e3*ebkev(ib)/(j*atw(ibion))
!$$$!990131            aloge=alog10(e)
!$$$            aloge=log10(e)
!$$$            sigcx=.6937e-14*(1.-.155*aloge)**2/
!$$$     1        (1.+.1112e-14*e**3.3)
!$$$            aloge=aloge*2.302585093-6.907755279
!$$$            if(aloge.le.-2.30258) go to 21
!$$$            expo=cfionp(7)
!$$$            do 22 ii=1,6
!$$$              if7=7-ii
!$$$              expo=expo*aloge+cfionp(if7)
!$$$ 22         continue
!$$$            sigi=exp(expo)
!$$$            go to 23
!$$$ 21         continue
!$$$            sigi=0.
!$$$ 23         continue
!$$$ccc
!$$$            do 24 i=1,mfm1
!$$$              sgxncx = sigcx*zni(i,ion)
!$$$              sgxni  = sigi*zni(i,ion)
!$$$              sgxn(i,j,ib)=sgxn(i,j,ib)+sgxncx+sgxni
!$$$ 24         continue
!$$$            go to 27
!$$$!------------------------------------------------------------
!$$$!     r.e. olson et al., phys. rev. lett. 41, 163 (1978)
!$$$!------------------------------------------------------------
!$$$ 25         ebi=ebkev(ib)/(j*atw(ibion))
!$$$            do 26 i=1,mfm1
!$$$              rpath=1.0e-17*zni(i,ion)*46.*zzi(i,ion)*
!$$$     &          (32.*zzi(i,ion)/ebi)*(1.-exp(-ebi/(32.*zzi(i,ion))))
!$$$              sgxn(i,j,ib)=sgxn(i,j,ib)+rpath
!$$$ 26         continue
!$$$!------------------------------------------------------------
!$$$ 27       continue
!$$$ 281    continue
!$$$ 280  continue
!$$$ccc
!$$$      do 38 ib=1,mb
!$$$        if(ib.gt.1 .and. ebkev(ib).eq.ebkev(1)) go to 36
!$$$        do 35 j=1,3
!$$$          sgxnm=sgxn(1,j,ib)
!$$$          do 31 i=2,mfm1
!$$$            if(sgxnm.lt.sgxn(i,j,ib)) sgxnm=sgxn(i,j,ib)
!$$$ 31       continue
!$$$          sgxnmi(j,ib)=1./sgxnm
!$$$ 35     continue
!$$$        go to 38
!$$$ 36     do 37 j=1,3
!$$$ 37     sgxnmi(j,ib) = sgxnmi(j,1)
!$$$ 38   continue
!$$$      return
!$$$ccc
!$$$!     ONETWO DIVERGENCE
!$$$
!$$$!..................................................................
!$$$!     Stearns formula
!$$$!..................................................................
!$$$
!$$$ 100  continue
!$$$      do 112 ib=1,mb
!$$$        do 111 j=1,3
!$$$          sigeff1=7.e-19*(dzemax+.2)**.65
!$$$          sigeff2=(3.*dnemax/(2.*dtemax+1.))**.118
!$$$          tem=DBLE(j)*atw(ibion)
!$$$          tem1=-(.7+.02*dzemax)
!$$$          sigeff3=(ebkev(ib)/(1000.*tem))**tem1
!$$$          sigeff=sigeff1*sigeff2*sigeff3
!$$$          do 110 i=1,mfm1
!$$$            sgxn(i,j,ib)=sigeff*zne(i)
!$$$ 110      continue
!$$$          sgxnmi(j,ib)=1./sgxn(1,j,ib)
!$$$ 111    continue
!$$$ 112  continue
!$$$      return
!$$$      end
!
!
      subroutine crsecs(iexcit,atw,ebkev,ibion,mb,mfm1,nion,vbeam,zne, &
                        zni,zte,zzi,dtemax,dnemax,dzemax, &
                        sgvxne,sgxn,sgxnmi)
      use param_mod
!
!     This subroutine calculates charge exchange, proton ionization, and
!     electron ionization cross sections from the fitted results of
!     freeman and jones, clm-r137, culham (1974).
!
      implicit none 
!      include 'params.inc'   BH: The file is in Kinsey's cql3d update.
!                                 Not needed except in standalone freya.
!
      integer iexcit, mfm1 !input
      integer mb, nion, ibion !input
      real(c_double) :: dtemax,dnemax,dzemax !input
      integer ib, i,ii,j,k,if7,ion !local
      real(c_double) :: alogt,aloge,teev,expo,ebi,e,veli,sigcx,sigi !local
      real(c_double) :: sgxncx,sgxni,sgxnm,rpath,tem,tem1,tc,dum !local
      real(c_double) :: sigeff,sigeff1,sigeff2,sigeff3
      
      real(c_double) :: atw(kion),ebkev(kb),vbeam(ke,kb),zne(kz),zni(kz,kion), &
          zte(kz), zzi(kz,kion), sgvxne(kz)
      real(c_double) :: sgxn(kcmp1,kz,kbe,ksge), sgxnmi(ke,kb)  ! stand alone code
!      real(c_double) :: sgxn(kcmp1,kz,ke,kb), sgxnmi(ke,kb)
!
      real(c_double) ::  cfione(7), cfionp(7)
!
      data cfione / -3.173850D+01, 1.143818D+01, -3.833998D0, &
          7.046692D-01, -7.431486D-02, 4.153749D-03, -9.486967D-05/
      data cfionp / -4.203309D+01, 3.557321D0, -1.045134D0, &
          3.139238D-01, -7.454475D-02, 8.459113D-03, -3.495444D-04/
!
!---:----1----:----2----:----3----:----4----:----5----:----6----:----7-!
!
!... initializations
!
      do ib=1,kb
       do j=1,ke
         sgxnmi(j,ib)=0.D0
         do i=1,kz
          do k=1,kcmp1
            sgxn(k,i,j,ib)=0.D0
          enddo
         enddo
       enddo
      enddo
      do i=1,kz
       sgvxne(i)=0.D0
      enddo
!
      if (iexcit.eq.-1) go to 200
!
      do 20 i = 1, mfm1
        alogt = 0.
        teev = 1.D3 * zte(i)
        if (teev.gt.1.D0) alogt = dlog(teev)
        if (teev.gt.1.D5) alogt = 11.51D0
        expo = cfione(7)
        do 10 ii = 1, 6
          if7 = 7 - ii
          expo = expo * alogt + cfione(if7)
   10   continue
        sgvxne(i) = dexp(expo) * zne(i)
   20 continue
!
      do 50 ib = 1, mb
        do 40 j = 1, 3
          veli = 1.0 / vbeam(j,ib)
          do 30 i = 1, mfm1
            sgxn(4,i,j,ib) = sgvxne(i) * veli
   30     continue
   40   continue
   50 continue
!
      do 110 ib = 1, mb
        do 100 ion = 1, nion
          do 90 j = 1, 3
            if (atw(ion).le.3.01) then
              e = 1.D3 * ebkev(ib) / (j*atw(ibion))
              aloge = dlog10(e)
              sigcx = .6937D-14 * (1.D0-.155D0*aloge) ** 2 / &
                  (1.D0+.1112D-14*e**3.3D0)
              aloge = aloge * 2.302585093D0 - 6.907755279D0
              if (aloge.gt.-2.30258) then
                expo = cfionp(7)
                do 60 ii = 1, 6
                  if7 = 7 - ii
                  expo = expo * aloge + cfionp(if7)
   60           continue
                sigi = dexp(expo)
              else
                sigi = 0.D0
              endif
!
              do 70 i = 1, mfm1
                sgxncx = sigcx * zni(i,ion)
                sgxni = sigi * zni(i,ion)
                sgxn(4,i,j,ib) = sgxn(4,i,j,ib) + sgxncx + sgxni
   70         continue
!
            else
!------------------------------------------------------------
!  r.e. olson et al., phys. rev. lett. 41, 163 (1978)
!------------------------------------------------------------
              ebi = ebkev(ib) / (j*atw(ibion))
              do 80 i = 1, mfm1
                rpath = 1.0e-17 * zni(i,ion) * 46.D0 * zzi(i,ion) * &
                   (32.D0*zzi(i,ion)/ebi) * &
                   (1.D0-dexp(-ebi/(32.D0*zzi(i,ion))))
                sgxn(4,i,j,ib) = sgxn(4,i,j,ib) + rpath
   80         continue
            endif
!------------------------------------------------------------
   90     continue
  100   continue
  110 continue
!
!         do ib = 1, mb
!           do j = 1, 3
!             do i=1,mfm1
!       write(*,'(3i3,2x,1p1e10.4,a14)') i,j,ib,sgxn(4,i,j,ib), ' sgxn'
!             enddo
!           enddo
!         enddo
!
      do 150 ib = 1, mb
        if (ib.le.1.or.ebkev(ib).ne.ebkev(1)) then
          do 130 j = 1, 3
            sgxnm = sgxn(4,1,j,ib)
            do 120 i = 2, mfm1
              if (sgxnm.lt.sgxn(4,i,j,ib)) sgxnm = sgxn(4,i,j,ib)
  120       continue
            sgxnmi(j,ib) = 1.D0 / sgxnm
  130     continue
        else
          do 140 j = 1, 3
            sgxnmi(j,ib) = sgxnmi(j,1)
  140     continue
        endif
  150 continue
!
      return
!
!..................................................................
!     Stearns formula
!..................................................................
!
 200  continue
      do 212 ib=1,mb
        do 211 j=1,3
          sigeff1=7.e-19*(dzemax+.2)**.65
          sigeff2=(3.*dnemax/(2.*dtemax+1.))**.118
          tem=DBLE(j)*atw(ibion)
          tem1=-(.7+.02*dzemax)
          sigeff3=(ebkev(ib)/(1000.*tem))**tem1
          sigeff=sigeff1*sigeff2*sigeff3
          do 210 i=1,mfm1
            sgxn(4,i,j,ib)=sigeff*zne(i)
 210      continue
          sgxnmi(j,ib)=1./sgxn(4,1,j,ib)
 211    continue
 212  continue
!
!         do ib = 1, mb
!           do j = 1, 3
!             do i=1,mfm1
!       write(*,'(3i3,2x,1p1e10.4,a14)') i,j,ib,sgxn(4,i,j,ib), ' sgxn2'
!             enddo
!           enddo
!         enddo
!
      return
      end subroutine crsecs
!
!
      real(c_double) function cxr (x)
      implicit none 
!
! ----------------------------------------------------------------------
! this function calculates the charge exchange rate for hydrogen atoms
! interacting with protons in units of cm**3/s.
! x is in units of keV for 1.0e-3 .le. x .le. 100, the
! the formula is taken from the paper by r.l. freeman and e.m. jones
! clm-r 137 culham laboratory 1974.
! for x .lt. 1.0e-3, a rate coefficient derived from an analytic average
! over the approximate cross section  sigma = 0.6937e-14*(1.0-0.155*LOG10
! (e/1ev))**2 is used.  this cross section is an approximation to that
! given by riviere, nuclear fusion 11,363(1971).
! ----------------------------------------------------------------------
!
      real(c_double) :: x !input
      real(c_double) :: tc,dum !local
      if (x .lt.   1.0e-3)  go to 10
      if (x .gt. 100     )  go to 20
      tc  = LOG (x)+6.9077553
      dum = 0.2530205e-4-tc*0.8230751e-6
      tc  = -0.1841757e2+tc*(0.528295-tc*(0.2200477-tc*(0.9750192e-1-tc* &
            (0.1749183e-1-tc*(0.4954298e-3+tc*(0.2174910e-3-tc*dum))))))
      tc  = EXP (tc)
      cxr = tc
      return
!
   10 tc  = 0.50654 - 6.7316e-2 * LOG  (x)
      cxr =           3.3340e-7 * SQRT (x) * (tc*tc + 7.454e-3)
      return
!
   20 cxr = 0.0
      return
!
      end function cxr
!
!
!      BLOCK DATA d01bbf_data
!!
!!     Took the following data statements out of subroutine
!!     d01bbf, since the Lahey compiler objected to their
!!     presence in d01bbf.
!!     Their presense in d01bbf is not ANSI F77.
!!     BobH, 010618.
!!
!      implicit none 
!!
!      common / cgh / wgh(70), xgh(70)
!!
!      equivalence (wgh10, wgh(1)), (xgh10, xgh(1))
!      equivalence (wgh16, wgh(11)), (xgh16, xgh(11))
!      equivalence (wgh20, wgh(27)), (xgh20, xgh(27))
!      equivalence (wgh24, wgh(47)), (xgh24, xgh(47))
!!
!      dimension wgh10(10), xgh10(10)
!      dimension wgh16(16), xgh16(16)
!      dimension wgh20(20), xgh20(20)
!      dimension wgh24(24), xgh24(24)
!!
!      common / cgl / wgl(70), xgl(70)
!!
!      equivalence (wgl10, wgl(1)), (xgl10, xgl(1))
!      equivalence (wgl16, wgl(11)), (xgl16, xgl(11))
!      equivalence (wgl20, wgl(27)), (xgl20, xgl(27))
!      equivalence (wgl24, wgl(47)), (xgl24, xgl(47))
!
!      dimension wgl10(10), xgl10(10)
!      dimension wgl16(16), xgl16(16)
!      dimension wgl20(20), xgl20(20)
!      dimension wgl24(24), xgl24(24)
!
!!
!      data wgh10   / &
!        7.6404329e-06,  1.3436457e-03,  3.3874394e-02,  2.4013861e-01, &
!        6.1086263e-01,  6.1086263e-01,  2.4013861e-01,  3.3874394e-02, &
!        1.3436457e-03,  7.6404329e-06/
!
!
!      data xgh10   / &
!        3.4361591e+00,  2.5327317e+00,  1.7566836e+00,  1.0366108e+00, &
!        3.4290133e-01, -3.4290133e-01, -1.0366108e+00, -1.7566836e+00, &
!        -2.5327317e+00, -3.4361591e+00/
!
!
!      data wgh16   / &
!        2.6548075e-10,  2.3209808e-07,  2.7118601e-05,  9.3228401e-04, &
!        1.2880312e-02,  8.3810041e-02,  2.8064746e-01,  5.0792948e-01, &
!        5.0792948e-01,  2.8064746e-01,  8.3810041e-02,  1.2880312e-02, &
!        9.3228401e-04,  2.7118601e-05,  2.3209808e-07,  2.6548075e-10/
!
!
!      data xgh16   / &
!        4.6887389e+00,  3.8694479e+00,  3.1769992e+00,  2.5462022e+00, &
!        1.9517880e+00,  1.3802585e+00,  8.2295145e-01,  2.7348105e-01, &
!        -2.7348105e-01, -8.2295145e-01, -1.3802585e+00, -1.9517880e+00, &
!        -2.5462022e+00, -3.1769992e+00, -3.8694479e+00, -4.6887389e+00/
!
!
!      data wgh20   / &
!        2.2293936e-13,  4.3993410e-10,  1.0860694e-07,  7.8025565e-06, &
!        2.2833864e-04,  3.2437733e-03,  2.4810521e-02,  1.0901721e-01, &
!        2.8667551e-01,  4.6224367e-01,  4.6224367e-01,  2.8667551e-01, &
!        1.0901721e-01,  2.4810521e-02,  3.2437733e-03,  2.2833864e-04, &
!        7.8025565e-06,  1.0860694e-07,  4.3993410e-10,  2.2293936e-13/
!
!
!      data xgh20   / &
!        5.3874809e+00,  4.6036824e+00,  3.9447640e+00,  3.3478546e+00, &
!        2.7888061e+00,  2.2549740e+00,  1.7385377e+00,  1.2340762e+00, &
!        7.3747373e-01,  2.4534071e-01, -2.4534071e-01, -7.3747373e-01, &
!        -1.2340762e+00, -1.7385377e+00, -2.2549740e+00, -2.7888061e+00, &
!        -3.3478546e+00, -3.9447640e+00, -4.6036824e+00, -5.3874809e+00/
!
!
!      data wgh24   / &
!        1.6643685e-16,  6.5846202e-13,  3.0462543e-10,  4.0189712e-08, &
!        2.1582457e-06,  5.6886916e-05,  8.2369248e-04,  7.0483558e-03, &
!        3.7445471e-02,  1.2773962e-01,  2.8617954e-01,  4.2693116e-01, &
!        4.2693116e-01,  2.8617954e-01,  1.2773962e-01,  3.7445471e-02, &
!        7.0483558e-03,  8.2369248e-04,  5.6886916e-05,  2.1582457e-06, &
!        4.0189712e-08,  3.0462543e-10,  6.5846202e-13,  1.6643685e-16/
!
!
!      data xgh24   / &
!        6.0159256e+00,  5.2593829e+00,  4.6256628e+00,  4.0536644e+00, &
!        3.5200068e+00,  3.0125461e+00,  2.5238810e+00,  2.0490036e+00, &
!        1.5842500e+00,  1.1267608e+00,  6.7417111e-01,  2.2441455e-01, &
!        -2.2441455e-01, -6.7417111e-01, -1.1267608e+00, -1.5842500e+00, &
!        -2.0490036e+00, -2.5238810e+00, -3.0125461e+00, -3.5200068e+00, &
!        -4.0536644e+00, -4.6256628e+00, -5.2593829e+00, -6.0159256e+00/
!
!
!      data wgl10   / &
!        3.0844112e-01,  4.0111993e-01,  2.1806829e-01,  6.2087456e-02, &
!        9.5015170e-03,  7.5300839e-04,  2.8259233e-05,  4.2493140e-07, &
!        1.8395648e-09,  9.9118272e-13/
!
!
!      data xgl10   / &
!        1.3779347e-01,  7.2945455e-01,  1.8083429e+00,  3.4014337e+00, &
!        5.5524961e+00,  8.3301527e+00,  1.1843786e+01,  1.6279258e+01, &
!        2.1996586e+01,  2.9920697e+01/
!
!
!
!
!      data wgl16   / &
!        2.0615171e-01,  3.3105785e-01,  2.6579578e-01,  1.3629693e-01, &
!        4.7328929e-02,  1.1299900e-02,  1.8490709e-03,  2.0427192e-04, &
!        1.4844587e-05,  6.8283193e-07,  1.8810248e-08,  2.8623502e-10, &
!        2.1270790e-12,  6.2979670e-15,  5.0504737e-18,  4.1614624e-22/
!
!
!      data xgl16   / &
!        8.7649410e-02,  4.6269633e-01,  1.1410578e+00,  2.1292836e+00, &
!        3.4370866e+00,  5.0780186e+00,  7.0703385e+00,  9.4383143e+00, &
!        1.2214223e+01,  1.5441527e+01,  1.9180157e+01,  2.3515906e+01, &
!        2.8578730e+01,  3.4583399e+01,  4.1940453e+01,  5.1701160e+01/
!
!
!      data wgl20   / &
!        1.6874680e-01,  2.9125436e-01,  2.6668610e-01,  1.6600245e-01, &
!        7.4826065e-02,  2.4964417e-02,  6.2025508e-03,  1.1449624e-03, &
!        1.5574177e-04,  1.5401441e-05,  1.0864864e-06,  5.3301209e-08, &
!        1.7579812e-09,  3.7255024e-11,  4.7675293e-13,  3.3728442e-15, &
!        1.1550143e-17,  1.5395221e-20,  5.2864427e-24,  1.6564566e-28/
!
!
!      data xgl20   / &
!        7.0539890e-02,  3.7212682e-01,  9.1658210e-01,  1.7073065e+00, &
!        2.7491993e+00,  4.0489253e+00,  5.6151750e+00,  7.4590175e+00, &
!        9.5943929e+00,  1.2038803e+01,  1.4814293e+01,  1.7948896e+01, &
!        2.1478788e+01,  2.5451703e+01,  2.9932555e+01,  3.5013434e+01, &
!        4.0833057e+01,  4.7619994e+01,  5.5810796e+01,  6.6524417e+01/
!
!
!      data wgl24   / &
!        1.4281197e-01,  2.5877411e-01,  2.5880671e-01,  1.8332269e-01, &
!        9.8166273e-02,  4.0732478e-02,  1.3226019e-02,  3.3693491e-03, &
!        6.7216256e-04,  1.0446121e-04,  1.2544722e-05,  1.1513158e-06, &
!        7.9608130e-08,  4.0728590e-09,  1.5070082e-10,  3.9177365e-12, &
!        6.8941811e-14,  7.8198004e-16,  5.3501888e-18,  2.0105175e-20, &
!        3.6057659e-23,  2.4518188e-26,  4.0883016e-30,  5.5753458e-35/
!
!
!      data xgl24   / &
!        5.9019852e-02,  3.1123915e-01,  7.6609691e-01,  1.4255976e+00, &
!        2.2925621e+00,  3.3707743e+00,  4.6650837e+00,  6.1815351e+00, &
!        7.9275392e+00,  9.9120980e+00,  1.2146103e+01,  1.4642732e+01, &
!        1.7417993e+01,  2.0491460e+01,  2.3887330e+01,  2.7635937e+01, &
!        3.1776041e+01,  3.6358406e+01,  4.1451720e+01,  4.7153106e+01, &
!        5.3608575e+01,  6.1058531e+01,  6.9962240e+01,  8.1498279e+01/
!
!
!      end BLOCK DATA d01bbf_data
!!
!
!CMG  changed 11/13/2017    subroutine d01bbf(type,inum,weight,abscis,ier)
      subroutine d01bbf_cql3d(type,inum,weight,abscis,ier)
      implicit none 
      save
      character*2 type
      integer inum !input
      integer ier !output
      integer i,index,irel !local
      real(c_double) :: weight(inum), abscis(inum)
      
      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug

      integer num(4)
      integer iorg(4)
!
      real(c_double) :: wgh,xgh !in /cgh/
      common / cgh / wgh(70), xgh(70)
!
      equivalence (wgh10, wgh(1)), (xgh10, xgh(1))
      equivalence (wgh16, wgh(11)), (xgh16, xgh(11))
      equivalence (wgh20, wgh(27)), (xgh20, xgh(27))
      equivalence (wgh24, wgh(47)), (xgh24, xgh(47))
!
      real(c_double) :: wgh10(10), xgh10(10)
      real(c_double) :: wgh16(16), xgh16(16)
      real(c_double) :: wgh20(20), xgh20(20)
      real(c_double) :: wgh24(24), xgh24(24)
!
      real(c_double) :: wgl,xgl !in /cgl/
      common / cgl / wgl(70), xgl(70)
!
      equivalence (wgl10, wgl(1)), (xgl10, xgl(1))
      equivalence (wgl16, wgl(11)), (xgl16, xgl(11))
      equivalence (wgl20, wgl(27)), (xgl20, xgl(27))
      equivalence (wgl24, wgl(47)), (xgl24, xgl(47))

      real(c_double) :: wgl10(10), xgl10(10)
      real(c_double) :: wgl16(16), xgl16(16)
      real(c_double) :: wgl20(20), xgl20(20)
      real(c_double) :: wgl24(24), xgl24(24)
!
      data num /10,16,20,24/
      data iorg /1,11,27,47/

      do 50 i = 1, 13
        index = i
        if(inum.eq.num(i)) go to 60
 50   continue
      ier = 1
      return

 60   continue
      ier = 0
!
      if (type.ne.'gl') go to 180
      do 100 i = 1, inum
        irel = iorg(index) - 1 + i
        weight(i) = wgl(irel)
        abscis(i) = xgl(irel)
 100  continue
      go to 1000
!
 180  continue
      if(type.ne.'gh') go to 300
      do 200 i = 1, inum
        irel = iorg(index) - 1 + i
        weight(i) = wgh(irel)
        abscis(i) = xgh(irel)
 200  continue
      go to 1000
!
 300  continue
      ier = 1
      return
!
 1000 continue
!
      return
      end subroutine d01bbf_cql3d
!
!
      subroutine decomp(ndim,n,a,cond,ipvt,work)
      implicit none 
      integer ndim,n
      real(c_double) :: a(ndim,n),cond,work(n)
      real(c_double) :: zero
      integer ipvt(n)
!
!     decomposes a real matrix by gaussian elimination
!     and estimates the condition of the matrix.
!
!     use solveq to compute solutions to linear systems.
!
!     input..
!
!     ndim = declared row dimension of the array containing  a.
!
!     n = order of the matrix.
!
!     a = matrix to be triangularized.
!
!     output..
!
!     a  contains an upper triangular matrix  u  and a permuted
!     version of a lower triangular matrix  i-l  so that
!     (permutation matrix)*a = l*u .
!
!     cond = an estimate of the condition of  a .
!     for the linear system  a*x = b, changes in  a  and  b
!     may cause changes  cond  times as large in  x .
!     if  cond+1.0 .eq. cond , a is singular to working
!     precision.  cond  is set to  1.0e+32  if exact
!     singularity is detected.
!
!     ipvt = the pivot vector.
!     ipvt(k) = the index of the k-th pivot row
!     ipvt(n) = (-1)**(number of interchanges)
!
!     work space..  the vector  work  must be declared and included
!     in the call.  its input contents are ignored.
!     its output contents are usually unimportant.
!
!     the determinant of a can be obtained on output by
!     det(a) = ipvt(n) * a(1,1) * a(2,2) * ... * a(n,n).
!
      real(c_double) :: ek, t, anorm, ynorm, znorm
      integer nm1, i, j, k, kp1, kb, km1, m
!
      zero=0.d0

      ipvt(n) = 1
      if (n .eq. 1) go to 80
      nm1 = n - 1
!
!     compute 1-norm of a
!
      anorm = 0.0
      do 10 j = 1, n
        t = 0.0
        do 5 i = 1, n
          t = t + abs(a(i,j))
 5      continue
        if (t .gt. anorm) anorm = t
 10   continue
!
!     gaussian elimination with partial pivoting
!
      do 35 k = 1,nm1
        kp1= k+1
!
!     find pivot
!
        m = k
        do 15 i = kp1,n
          if (abs(a(i,k)) .gt. abs(a(m,k))) m = i
 15     continue
        ipvt(k) = m
        if (m .ne. k) ipvt(n) = -ipvt(n)
        t = a(m,k)
        a(m,k) = a(k,k)
        a(k,k) = t
!
!     skip step if pivot is zero
!
        if (t .eq. 0.0) go to 35
!
!     compute multipliers
!
        do 20 i = kp1,n
          a(i,k) = -a(i,k)/t
 20     continue
!
!     interchange and eliminate by columns
!
        do 30 j = kp1,n
          t = a(m,j)
          a(m,j) = a(k,j)
          a(k,j) = t
          if (t .eq. 0.0) go to 30
          do 25 i = kp1,n
            a(i,j) = a(i,j) + a(i,k)*t
 25       continue
 30     continue
 35   continue
!
!     cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)
!     estimate obtained by one step of inverse iteration for the
!     small singular vector.  this involves solving two systems
!     of equations, (a-transpose)*y = e  and  a*z = y  where  e
!     is a vector of +1 or -1 chosen to cause growth in y.
!     estimate = (1-norm of z)/(1-norm of y)
!
!     solve (a-transpose)*y = e
!
      do 50 k = 1, n
        t = 0.0
        if (k .eq. 1) go to 45
        km1 = k-1
        do 40 i = 1, km1
          t = t + a(i,k)*work(i)
 40     continue
 45     ek = 1.0
        if (t .lt. 0.0) ek = -1.0
        if (a(k,k) .eq. 0.0) go to 90
        work(k) = -(ek + t)/a(k,k)
 50   continue
      do 60 kb = 1, nm1
        k = n - kb
        t = 0.0
        kp1 = k+1
        do 55 i = kp1, n
          t = t + a(i,k)*work(k)
 55     continue
        work(k) = t
        m = ipvt(k)
        if (m .eq. k) go to 60
        t = work(m)
        work(m) = work(k)
        work(k) = t
 60   continue
!
      ynorm = 0.0
      do 65 i = 1, n
        ynorm = ynorm + abs(work(i))
 65   continue
!
!     solve a*z = y
!
      call solveq(ndim, n, a, work, ipvt)
!
      znorm = zero
      do 70 i = 1, n
        znorm = znorm + abs(work(i))
 70   continue
!
!     estimate condition
!
      cond = anorm*znorm/ynorm
      if (cond .lt. 1.0) cond = 1.0
      return
!
!     1-by-1
!
 80   cond = 1.0
      if (a(1,1) .ne. 0.0) return
!
!     exact singularity
!
 90   cond = 1.0e+32
      return
      end subroutine decomp
!
!
      real(c_double) function dfhx(beta)
      implicit none 
      save
!*****************************************************************
!--   approximation to the function defined by janev & presnyakov,
!     j. phys. b 13, 4233 (1980).
!
      real(c_double) :: beta !input
      real(c_double) :: dd(38)
      integer ia
      real(c_double) :: beta1,a !local
      
      data dd/ &
        .10450,.121, .138, .157, .175, .200, .229, .260, .300, .339, &
        .367, .389, .402, .410, .398, .376, .346, .317, .285, .255, &
        .227, .205, .185, .168, .152, .138, .124, .110, .099, .089, &
        .079, .070, .062, .054, .047, .041, .035, .02898/
!
      beta1=1./beta
      if(beta1 .le. .2)goto 110
      if(beta1 .ge. 1.e3)goto 120
!990131      a=10.*alog10(beta1)+8.
      a=10.*log10(beta1)+8.
      ia=min0(37,int(a))
      dfhx=dd(ia)+(a-DBLE(ia))*(dd(ia+1)-dd(ia))
      return
 110  dfhx=.5*beta*(1.-1./(8.*beta*sqrt(beta)))*exp(-sqrt(2.*beta))
      return
!990131 120  dfhx=4.*beta*alog(1.4*beta1)
 120  dfhx=4.*beta*log(1.4*beta1)
      return
      end function dfhx
!
!
      subroutine eigen(xeig)
      use frsubs_mod, only : eigrf
      implicit none 
      save
!*****************************************************************
      real(c_double) :: xeig !input
      integer ms,mc,mz,mi !local
      parameter (ms=21,mc=35,mz=1,mi=2 )
      
      integer i,j,ijob,ifail1 !local
      real(c_double) :: eigmin,sum,cond,cpua,cpubb

      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf

      integer ns,nc,numi,numz,iz,izstrp !in /b0/
      real(c_double) :: er0,v0,te,ti,ami,deni,amz,denz,zcor,zsqcor,dene !in /b0/
      common/b0/ns,nc,numi,numz,iz(mz),izstrp(mz), &
        er0,v0,te,ti,ami(mi),deni(mi), &
        amz(mz),denz(mz),zcor(mz),zsqcor(mz),dene

      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug
      
      real(c_double) :: q !in /b7/
      common/b7/q(ms+1,ms+1)
      
      real(c_double) :: eigvr,eigvl,cj !in /b9/
      common/b9/eigvr(ms+1,ms+1),eigvl(ms+1),cj(ms+1)
      
      real(c_double) :: xfrac !in b10
      common/b10/xfrac
      
      real(c_double) :: qq(ms+1,ms+1),eigr(ms+1),eigi(ms+1)
      integer intger(ms+1)
      real(c_double) :: evr(ms+1),evi(ms+1),efr(ms+1,ms+1),efi(ms+1,ms+1)
      real(c_double) :: xn(ms+1),yn(ms+1),work2(3*(ms+1)),coef(2)
      real(c_double) :: w(2*(ms+1)),z(2*(ms+1),ms+1)
      real(c_double) :: cjmat(ms+1,ms+1)
      real(c_double) :: workdc(ms+1)
      integer ipvt(ms+1)
      real(c_double) :: infac(ms+1,ms+1), lexp(ms+1), in(ms+1), dindx(ms+1)
      real(c_double) :: zero,xtest
      data xtest/1./
!
!     on return from eigrf, w contains the complex*16 eigenvalues
!     as w= [wr(1),wi(1), ... ,wr(ns),wi(ns)]
!     on return from eigrf, z contains the complex*16 eigenvectors
!     as z(col1)= [zr(1),zi(1), ... ,zr(ns),zi(ns)], etc.,
!     each column representing an eigenvector
!
!BH080118      call second(cpua)
!

      zero=0.d0

      do 10 i=1,ns
        do 101 j=1,ns
          qq(i,j)=q(i,j)
 101    continue
 10   continue
      if(nouthx.gt.0) then
        write(nouthx,11)
 11     format( 'q _ matrix')
        do 13 i=1,ns
          write(nouthx,12) (q(i,j), j=1,ns)
 12       format(1x,1p7e11.3)
 13     continue
      endif


!@    ifail1=0
!@    call f02aff(qq,ms+1,ns,eigr,eigi,intger,ifail1)
!@    if(ifail1.eq.0) go to 15
      ijob = 1
      call eigrf(qq,ns,ms+1,ijob,w,z,ms+1,work2,ifail1)
      if(ifail1.eq.0)goto 15
      ihxbug = 9
      if(nouthx.gt.0) write(nouthx,3939) ifail1
 3939 format(' ??? eigrf routine error: ifail1= ',i4)
      return
!
 15   continue
      do 18 j=1,ns
        do 17 i=1,ns
          eigvr(i,j)=z(2*i-1,j)
 17     continue
        eigvl(j) = w(j*2-1)
 18   continue
      eigmin=1.e30
      do 20 i=1,2*ns,2
!990131        eigmin=amin1(eigmin,abs(w(i)))
        eigmin=min(eigmin,abs(w(i)))
 20   continue
!@    do 20 i=1,ns
!@    eigmin=amin1(eigmin,abs(eigr(i)))
!@    20 continue
      xeig=0.
      if(eigmin.ne.zero)xeig=v0/eigmin
!
!     determine fraction of 3rd excited state (approximate)
      sum=0
      do 25 i=1,ns
        sum=sum+eigvr(i,1)
 25   continue
      xfrac=eigvr(4,1)/sum
      if(nouthx.eq.0) return

      write(nouthx,1000) xeig
 1000 format(' length (eigenvalue):  ',1pe11.4)
      write(nouthx,1010) (i,eigvl(i),i=1,ns)
 1010 format('    i= ',i2,'         eigenvalue= ',1pe11.4)
      do 1030 j = 1, ns
        write(nouthx,1020) j,(eigvr(i,j),i=1,ns)
 1020   format(' eigenvector[',i2,']= ',/,(10x,1p5e11.4))
 1030 continue

!...  apply initial condition constraint
      do 1200 i = 1, ns
        do 1201 j = 1, ns
          cjmat(i,j) = eigvr(i,j)
 1201   continue
 1200 continue
      cj(1) = 1.0
      do 1300 i = 2, ns
        cj(i) = 0.
 1300 continue
      call decomp (ms+1,ns,cjmat,cond,ipvt,workdc)
      call solveq(ms+1,ns,cjmat,cj,ipvt)

      write(nouthx,1320) (cj(i), i=1,ns)
 1320 format(' bc constants: ',/,(10x,1p5e11.4))

      do 1340 j = 1, ns
        do 1350 i = 1, ns
          lexp(i) = eigvl(i) / v0
          infac(j,i) = cj(i) * eigvr(j,i)
 1350   continue
        write (nouthx,1400) j, (infac(j,i), lexp(i), i = 1, ns)
 1400   format( ' *** i(',i2,')=   ',/, &
          (12x, 1pe15.7,' * exp(',1pe15.7,' * x)') )
 1340 continue
!
!%OS
      print *,' WARNING: xtest not defined in EIGEN'
!%OS
      do 1330 j = 1,ns
        in(j) = 0.
        do 1345 i = 1, ns
          in(j) = in(j) + infac(j,i) * exp(lexp(i)*xtest)
 1345   continue
 1330 continue

      write(nouthx,1450)
 1450 format(///' qnj * ij / v0 ...........',/)
      do 1360 j = 1, ns
        sum = 0.
        do 1370 i = 1, ns
          sum = sum + q(j,i) * in(i)
 1370   continue
        write(nouthx,1500) j, sum/v0
 1500   format(15x,i3,5x,1pe15.7)
 1360 continue

!BH080118      call second(cpubb)
      cpu7 = cpu7 + cpubb - cpua
!
      return
      end subroutine eigen
!
!
      real(c_double) function encapf(vcvo, tstcx)
      implicit none 
      real(c_double) :: vcvo,tstcx !input
      real(c_double) :: v3,tfts !local
!
!     This routine evaluates the function N, which is related to the
!      rate at which fast ions slow down on electrons and thermal ions.
!       This function is defined by Callen et al., IAEA Tokyo,
!       Vol. I, 645(1974).
!
      v3 = vcvo**3
!990131      tfts = alog(1.+1./v3)/3.
      tfts = log(1.+1./v3)/3.
      if(tstcx.gt.0.01) go to 10
      encapf = tfts
      return
 10   encapf = (1.-exp(-tfts*tstcx))/tstcx
      return
      end function encapf
!
!
      subroutine freyorb(atw, codeid, drutpi, drini, droti, &
        ic, iskip, izi, mfm1, norb, pinsid, potsid, pzone, &
        rinsid, rotsid, rzone, ri, zetai, v, zi, &
        ipass, iaxis, ier, izp, wid, wt, zeta,angmtp)
      implicit none 
      real(c_double) :: atw,drutpi,drini,droti !input
      real(c_double) :: pzone,rzone,ri,zetai,v,zi,wid !input
      integer ic,iskip,izi,mfm1,norb,ipass,iaxis,izp !input
      integer ier !output
      save
!
!     This subroutine calculates the following orbit parameters for a
!     fast ion:
!     ipass : 1, particle is passing
!     0, particle is trapped
!     iaxis : 1, particle encircles the axis
!     0, particle does not encircle the axis
!     ier   : 1, routine had difficulty determining orbit
!     0, routine obtained reasonable orbit
!     izp   : outermost zone of orbit
!     wid   : weight giving approximate fraction of time spent
!     by fast ion in each zone
!     zeta  : approximate pitch angle cosine in each zone along
!     orbit
!     The fast ion is assumed to be initially traveling in the
!     co-current direction.
!
      real(c_double) :: pinsid(*), potsid(*), rinsid(*), rotsid(*)
      real(c_double) :: wt(*), zeta(*), angmtp(*)
      real(c_double) :: eoverc,xmassp,zero,one,onem5,onep5,onem3,psim,rm !local
      integer itm,izm,i,itp,nout !local
      real(c_double) :: xmass,rmajor,psiax,psii,pang,rt,psit,pangt,pangx
      real(c_double) :: zetax2,zetax,zetap,zetam,c1,c2,c3,ra,rp,psip
      real(c_double) :: wtx,wtm,wtp,xzp,xzm,zdif
      
      character*8 codeid,flagm,flagp
      data eoverc / 1.6e-20 /
      data xmassp / 1.673e-24 /
!
      zero=0.d0
      one=1.d0
      onem5=1.d-5
      onep5=1.d5
      onem3=1.d-3

!     initialize parameters
!
      flagp = ' '
      flagm = ' '
      ier   = 0
      itm   = 0
      izm   = izi
      psim  = 0.
      rm    = 0.
      wid   = 0.
      do 2 i=1,mfm1
        wt(i) = 0.
        angmtp(i) = 0.
!     angmtp will be used for orbit smearing aof angular momentum
!     this option will be implemented later
    2 zeta(i) = 0.
      wt(izi) = 1.
      zeta(izi) = zetai
!
!     calculate ion mass times c/e; initialize rmajor and psiax
!
      xmass  = atw*xmassp/eoverc
      rmajor = rotsid(1)
      psiax  = potsid(1)
!
!     calculate poloidal magnetic flux at the fast ion birth point
!
      psii = pzone
      if(codeid.ne.'onedee') go to 5
      call inter(droti,0,mfm1,zero,potsid,rzone,psii)
    5 continue
!
!     determine whether fast ion is passing or trapped
!
      pang = xmass*v*zetai*ri - psii
      ipass = 1
      if(abs(zetai).eq.1.) go to 40
      if(zetai.eq.zero) go to 30
      rt = (1.-zetai**2)*ri
      if(rt.gt.rmajor) go to 10
      call inter(drini,0,mfm1,rmajor,pinsid,rt,psit)
      go to 20
 10   call inter(droti,0,mfm1,rmajor,potsid,rt,psit)
 20   pangt = -psit
      if(pang.gt.pangt) go to 40
 30   ipass = 0
 40   continue
!
!     determine whether ion encircles the magnetic axis
!
      iaxis = 1
      zetax2 = 1. - (1.-zetai**2)*ri/rmajor
      if(zetax2.le.0.) go to 50
      zetax = sqrt(zetax2)
      if(ipass.eq.0) zetax = -zetax
      pangx = xmass*v*zetax*rmajor - psiax
      if(ipass.ne.0 .and. pang.lt.pangx) go to 60
      if(ipass.eq.0 .and. pang.gt.pangx) go to 60
 50   iaxis = 0
 60   continue
!
!     calculate poloidal magnetic flux and zone index at outer major
!     radius of orbit
!
      c1 = (1.-zetai**2)*ri
      c2 = 1./(xmass*v)
      c3 = zetai*ri - c2*psii
      if(codeid.ne.'onedee') go to 110
      ra = rmajor + rzone
      go to 115
 110  call inter(drutpi,1,mfm1,psiax,rotsid,psii,ra)
 115  ra = ra + 2.

      call scalit(ra, rp, zero, onep5, 1, 1, itp, 10,onem5,onem3,one,0, &
        c1, c2, c3, droti, mfm1, psip, potsid, rmajor)
      if(rp.lt.ri-0.01 .or. rp.lt.rmajor .or. itp.lt.0) flagp = '*'
      if(flagp.eq.'*') go to 340
      if(codeid.ne.'onedee') go to 150
      xzp = (rp-rmajor)*droti + 1.
      go to 160
 150  xzp = sqrt(psip-psiax)*drutpi + 1.
 160  izp = xzp
      izp = max0(izp,izi)
      izp = min0(izp,mfm1+1)
      zetap = (pang+psip)/(xmass*v*rp)
!
!     determine whether ion is confined or lost; return if it is lost
!
      if(izp.gt.mfm1) go to 350
!
!     calculate poloidal magnetic flux and zone index at inner major
!     radius of orbit
!
      if(iaxis.eq.0) go to 220
      if(codeid.ne.'onedee') go to 210
      ra = rmajor - rzone
      go to 215
 210  call inter(drutpi,1,mfm1,psiax,rinsid,psii,ra)
 215  call scalit(ra, rm, zero, onep5, 1, 1, itm,10,onem5,-onem3,one,0, &
        c1, c2, c3, drini, mfm1, psim, pinsid, rmajor)
      go to 240
 220  ra = rmajor + 0.25*(rp-rmajor)
      call scalit(ra, rm, zero, onep5, 1, 1, itm,10,onem5,onem3,one,0, &
        c1, c2, c3, droti, mfm1, psim, potsid, rmajor)
 240  continue
      if(ipass.eq.0 .and. rm.gt.ri+0.01) flagm = '*'
      if(iaxis.ne.0 .and. rm.gt.rmajor) flagm = '*'
      if(itm.lt.0) flagm = '*'
      if(rm.gt.rp-0.01) flagm = '*'
      if(flagm.eq.'*') go to 340
      if(codeid.ne.'onedee') go to 250
      xzm = abs(rmajor-rm)*droti + 1.
      go to 260
 250  xzm = sqrt(psim-psiax)*drutpi + 1.
 260  izm = xzm
      izm = min0(izm,izi)
      zetam = (pang+psim)/(xmass*v*rm)
!
!     calculate orbit width
!
      wid = (rp-rmajor) - abs(rmajor-rm)
!990131      wid = amax1(wid,0.)
      wid = max(wid,zero)
!
!     calculate weights by zone; assume ion spends equal time in each
!     zone it fully traverses and proportionately less time in the
!     innermost and outermost zones of the orbit
!
      if(izp.eq.izm) go to 350
      wtm = 1. - (xzm - izm)
      wtp = xzp - izp
      wtx = 1./(wtm+wtp+izp-izm-1.)
      wt(izm) = wtm*wtx
      wt(izp) = wtp*wtx
      if(izp.eq.izm+1) go to 320
      do 310 i=izm+1,izp-1
 310  wt(i) = wtx
 320  continue
!
!     calculate pitch angle cosines by zone; assume cosine varies
!     linearly with zone index
!
      zdif = zetap - zetam
      zeta(izm) = zetam + 0.5*wt(izm)*zdif
      zeta(izp) = zetap - 0.5*wt(izp)*zdif
      if(izp.eq.izm+1) go to 350
      do 330 i=izm+1,izp-1
 330  zeta(i) = zeta(i-1) + wt(i)*zdif
      go to 350
!
!     set error flag if error was detected
!
 340  izp = izi
      ier = 1
!
!     print out selected results
!
 350  continue
      if(norb.eq.0) return
      nout = iabs(norb)
      if(ic.eq.1) write(nout,1010)
      if(mod(ic-1,iskip).ne.0) return
      write(nout,1020) ic, ipass, iaxis, zetai, ri, zi, psii, rp, flagp, &
        psip, rm, flagm, psim, itp, itm, izi, izm, izp
      if(norb.lt.0) return
      if(izp.gt.mfm1) return
      write(nout,1030) (wt(i), i=1,mfm1)
      write(nout,1040) (zeta(i), i=1,mfm1)
      return
!
 1010 format(/' orbit parameters'// &
        4x,'ic ip ix',7x,'zetai',10x,'ri',10x,'zi',8x,'psii', &
        10x,'rp',8x,'psip',10x,'rm',8x,'psim', &
        ' itp itm izi izm izp')
 1020 format(i6,2i3,f12.4,2f12.2,e12.3,2(f12.2,a1,e11.3),5i4)
 1030 format(' wt',f9.4,9f12.4/(10f12.4))
 1040 format(' zeta',f7.4,9f12.4/(10f12.4))
      end subroutine freyorb
!
!
      real(c_double) function gef(vpar, tpar)
      use frsubs_mod, only : asimp
      implicit none 
      real(c_double) :: vpar,tpar !input
      real(c_double) :: y,rooty,root13,fact4 !local
      real(c_double) :: term2,term3,term4,a,b,ep ! local
      integer m,m4
      save
!
!     This routine evaluates the function Ge, which is related to the trans-
!     fer of energy from fast ions slowing down on electrons.  The func-
!     tion is defined by Callen et al., IAEA Tokyo, Vol. I, 645 (1974).
!
      !external gefun ! Now in zfreya_mod module
      
      real(c_double) :: vcvo, tstcx, emzrat !in /gecom/
      common /gecom/ vcvo, tstcx, emzrat
      
!     root13 = sqrt(1./3.)
!     fact4  = root13*atan(root13)
      data root13 /0.577350/, fact4 /0.302300/
      if(tpar.gt.0.01) go to 10
      rooty = 1./vpar
      y = rooty**2
!990131      term2 = alog((1.-rooty+y)/(1.+rooty)**2)/(6.*y)
      term2 = log((1.-rooty+y)/(1.+rooty)**2)/(6.*y)
      term3 = root13*atan(root13*(2.*rooty-1.))/y
      term4 = fact4/y
      gef = 2.*(0.5-term2-term3-term4)
      return
 10   vcvo = vpar
      tstcx = tpar
      a = 0.
      b = 1.
      ep = 0.1
      m4 = 3
      gef = asimp(a,b,ep,m,m4,gefun)
      return
      end function gef

      real(c_double) function gefun(y)
      implicit none 
      real(c_double) :: y !input
      
      real(c_double) :: vcvo, tstcx, emzrat !in /gecom/
      common /gecom/ vcvo, tstcx, emzrat

      real(c_double) :: v3,arg,pcx !lovcal
      
      v3 = vcvo**3
      arg = (1.+v3)/(y**3+v3)
      pcx = (arg)**(-tstcx/3.)
      gefun = 2.*y**4*pcx/(y**3+v3)
      return
      end function gefun
!
!
!      BLOCK DATA hexnb_data
!      implicit none 
!      save
!
!!     Data statement of common block variables moved outside
!!     subroutine hexnb, and not ANSI F77 (and Lahey compiler
!!     complains).  BobH, 010618.
!      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
!        cpuii, cpuiz, cpuplf
!
!      data cpu1/0./, cpu2/0./, cpu3/0./, cpu4/0./, cpu5/0./, &
!        cpu6/0./, cpu7/0./, cpu8/0./
!
!      end BLOCK DATA hexnb_data
!!
!
      subroutine hexnb(istarx,iexcix,ilorenx,mstatx,ncontx, &
        er0x,tex,tix,numix,amix,denix, &
        numzx,izx,amzx,izstrx,denzx,bperpx, &
        ldene,ldeni,ldenz,lsvi,lsvz,lsve,lrad, &
        lngh,lngl,louthx,lcor, &
        nsigmav,lambda,hexfrac,ihxerr)
      implicit none 
      integer istarx,iexcix,ilorenx,mstatx,ncontx,numix,numzx !input
      integer ldene,ldeni,ldenz,lsvi,lsvz,lsve,lrad !input
      integer lngh,lngl,louthx,lcor !input
      integer ihxerr !output
      real(c_double) :: er0x,tex,tix,bperpx,hexfrac !input
      save
!*****************************************************************
!     version 21
!--   author:
!     c. d. boley
!     pppl 1984
!--   modified by : (vax / modular)
!     r. m. wieland
!     pppl july,1985
!--   ref.:
!     c. d. boley, r. k. janev, and d. e. post, phys. rev. letters
!     52, 534 (1984).
!--   calling sequence
!     call hexnb(istart,iexcit,ilorent,mstate,ncont,
!     &                   er0,tex,tix,numi,ami,deni,
!     &                   numz,iz,amz,izstrp,denz,bperp,
!     &                   kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,
!     &                   ngh,ngl,nouthx,ncorin,
!     &                   nsigmav,lambda)
!--   input parameters:
!     istart = 1 to initialize rad rates & fine pts       0 otherwise
!     iexcit if =0 then ignore contribution of excited states
!     if =1 then include contributionof excited states.
!     ilorent =0 or 1       whether or not to calculate the maximum
!     principal quantum number ns of the populated
!     states as given by the lorentz ionization
!     limit (1 ==> yes)
!     mstate< parameter ms
!     for ilorent=0, use ns = mstate+1
!     note: the operationalsignificance of ns in the code is to
!     set an upper limit such that any excitations to levels n
!     higher than ns are counted as "ionizations".
!     ncont< parameter mc
!     upper bound to number of continuum states
!     er0    beam energy per amu (ev)
!     te     electron temperature (ev)
!     ti     ion temperature (ev)
!     numi   number of hydrogenic ion species
!     ami    masses of the hydrogenic ions (amu)
!     deni   densities of the hydrogenic ions (cm**-3)
!     numz   number of impurity species
!     if numz gt 0, the coronal data file 'coronb' is read.
!     iz     atomic numbers of the impurities
!     izstrp =0 for coronal equilibrium values of <z> & <zsq>
!     =1 for fully stripped impurities
!     amz    atomic mass number of the impurities [iff izstrp#0]
!     denz   densities of the impurities (cm**-3)
!     bperp  magnetic field perpendicular to beamline (tesla (for ilorent=1))
!     kdene  =0: discard electron reactions
!     =1: include electron reactions (default)
!     kdeni  =0: discard ion (hydrogen) reactions
!     =1: include ion reactions (default)
!     kdenz  =0: discard impurity reactions
!     =1: include impurity reactions (default)
!     ksvi   =0: simple multiplication instead of full sigma-v
!     integral for ions
!     =1: full sigma-v integral for ions
!     ksvz   =0: simple multiplication instead of full sigma-v
!     integral for impurities
!     =1: full sigma-v integral for impurities
!     ksve   =0: sigma-v integrals for electron reactions
!     involving 1s, 2s, and 2p       rates for other
!     electron reactions from vriens & smeets.
!     =1: full sigma-v integrals for electron reactions
!     krad   =0: discard hydrogen line radiation
!     =1: include hydrogen line radiation .
!     ngh        order of gauss-hermite integrations for
!     ion and impurity reactions (default 24)
!     ngl        order of gauss-laguerre integrations for
!     electron reactions
!     permissible values: 10,16,20,24
!     nouthx unit number for messages(to unit nouthx if > 0)
!     (also nouthx>0 gives detailed diagnostic output)
!     ncorin unit number for reading coronb ce data file
!
!--   output:
!     nsigmav n<sigma-v> [sec**-1] : all processes included
!     lambda  mean free path (cm)
!     hexfrac fraction of 3rd excited state
!     ihxerr: error return
!     0     no error
!     1     radiation file coronb not found
!     2     d01bbf error in setting up ngh integration arrays
!     3     d01bbf error in setting up ngl integration arrays
!     4     sezf error: states "i" & "j" are the same!
!     5     sezf error: state "j" is 0
!     6     seef error       states "i" &/or "j" are in error
!     7     cief error
!     8     ceef error
!     9     eigrf error       eigenvalue determination is incorrect
!     10    input parameter error
!--   note:this routine uses imsl libraries.
!
!--   foreign files required:
!     hx2:[wieland.hex]coronb. :  this file contains all the coronal
!     equilibrium data for the various impurities.
!     dsk0:[imsl]imslibs/lib -- the imsl replacements for the nag codes.
!     nag:    f02aff et. al.
!     imsl:   eigrf  et. al.
!     hx2:[wieland.lib]wielib  for diagnostic matrix solvers
!     "decomp" & "solveq" for ax = b .
!     subr second -- a cpu timing routine of your choice
!--   recommended namelist values:
!     kdene = 1
!     kdeni = 1
!     kdenz = 1
!     ksvi  = 0
!     ksvz  = 0
!     ksve  = 0
!     krad  = 1
!     ngh   = 10
!     ngl   = 10
!     iexcit= 1
!     ilorent=0
!     mstate= 4
!     ncont = 30
!.....
!     these result in a cpu time of approx. 7 sec for a 25 point
!     radial profile of n*sigma(rho) with reasonably good accuracy.
!
      integer ms,mc,mz,mi !local
      parameter (ms=21,mc=35)
      parameter (mz=1,mi=2)

      real(c_double) :: amix(mi),denix(mi),denzx(mz),amzx(mz)
      integer izx(mz),izstrx(mz)
        
      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf

      integer ns,nc,numi,numz,iz,izstrp !in /b0/
      real(c_double) :: er0,v0,te,ti,ami,deni,amz,denz,zcor,zsqcor,dene !in /b0/
      common/b0/ns,nc,numi,numz,iz(mz),izstrp(mz), &
        er0,v0,te,ti,ami(mi),deni(mi), &
        amz(mz),denz(mz),zcor(mz),zsqcor(mz),dene
        
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      
      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug
      
      integer iexcit,ilorent,mstate,ncont !in /b8/
      common/b8/iexcit,ilorent,mstate,ncont
      
      real(c_double) :: xfrac !in b10
      common/b10/xfrac

      real(c_double) :: nsigmav, lambda
      
      integer i,kdene,kdeni,kdenz,ksvi,ksve,ksvz,krad,ki,kz !local
      integer ngh,ngl,ncor,nsvigma,iwatch
      real(c_double) :: bperp,xeig
!
      istart = istarx

!--   move vbls in hexnb commons
      iexcit = iexcix
      ilorent = ilorenx
      mstate = mstatx
      ncont = ncontx
      er0 = er0x
      te = tex
      ti = tix
      numi = numix
      if(numi.gt.mi) go to 9999
      do 2 i = 1, mi
        ami(i) = amix(i)
        deni(i) = denix(i)
    2 continue
      if(numz.gt.mz) go to 9999
      numz = numzx
      do 3 i = 1, mz
        iz(i) = izx(i)
        izstrp(i) = izstrx(i)
        amz(i) = amzx(i)
        denz(i) = denzx(i)
    3 continue
      bperp = bperpx
      kdene = ldene
      kdeni = ldeni
      kdenz = ldenz
      ksvi = lsvi
      ksvz = lsvz
      ksve = lsve
      krad = lrad
      ngh = lngh
      ngl = lngl
      nouthx = louthx
      ncor = lcor
!
      ihxbug = 0
      nsvigma = 0.
      lambda = 0.
      if(istart.eq.0)go to 5

      if(mstate+1.gt.ms .or. ncont+1 .gt. mc) go to 9999
!     ncor=31
      call hradin(ncor,numz,iz,izstrp,amz)
      if(ihxbug.gt.0) go to 20
      call hxinit
!
    5 continue
!
      v0=(1.3841e6)*sqrt(er0)
      call hxradi(te,zcor,zsqcor,numz,iz,izstrp,iwatch)
      dene=0.
      do 17 ki=1,numi
        dene=dene+deni(ki)
 17   continue
      do 18 kz=1,numz
        dene=dene+zcor(kz)*denz(kz)
 18   continue
!
!--   determine calculational mode: include excitations or not?
!
      goto (21,22),iexcit+1
!
!     no excitations included
 21   continue
      nc=1
!     if(nouthx.gt.0) write(nouthx,1000)
 1000 format(/,' no excitations')
      goto 23
!
!     include excitations
 22   continue
      nc=ncont+1
!     if(nouthx.gt.0) write(nouthx,1004) ncont
 1004 format(/,' excitations, with continuum at n=',i3)
!
 23   continue
      call lorent(v0,bperp,nc,ns)
!     if(nouthx.gt.0) write(nouthx,1005)ns
 1005 format(/,' max princ qn= ',i3)

      call hxsvi
      if(ihxbug.gt.0) go to 20

      call hxsve
      if(ihxbug.gt.0) go to 20

      call matri

      call eigen(xeig)
      if(ihxbug.gt.0) go to 20

      lambda = xeig
      nsigmav = v0 / xeig
      hexfrac=xfrac

!
 20   continue
      ihxerr = ihxbug
      istart=0
      return

 9999 continue
      ihxbug = 10
      if(nouthx.gt.0) write(nouthx,3939) numi,numz,mstate,ncont
 3939 format(' ??? inconsistency between input vbls and upper ', &
        'limits as defined by parameters :',/, &
        3x,'numi= ',i4,/ &
        3x,'numz= ',i4,/ &
        3x,'mstate= ',i4,/ &
        3x,'ncont= ',i4)
      go to 20

      end subroutine hexnb
!
!
      subroutine hradin(ncor,numz,nz,izstrp,amz)
      implicit none 
!*****************************************************************
!--   ifnumz ne 0, read the post radiation tables (file 'coronb').
      integer ncor,numz !input
      integer mz,mi !local
      parameter (mz=1,mi=2)

      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf

      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug

      real(c_double) :: arad !in locrad
      common/locrad/arad(8,15,mz)

      character*8 spec(mz)
      real(c_double) :: amz(mz)
      integer izstrp(mz),nz(mz)
      character*8 dum,dum2
      
      integer izsum,i,imp,j,k,iz,ia
      real(c_double) :: cpua,cpubb
!
!--   arrangement of arad(j,k,i) for given i:
!
!     1    2    3    4    5    6    7    8    9    10   11   12   13   14  15
!     1  t1   t2   t3   t4   0.   t1   t2   t3   t4   0.   t1   t2   t3   t4  0.
!     2  t2   t3   t4   t5   0.   t2   t3   t4   t5   0.   t2   t3   t4   t5  0.
!     3 a(0) a(0) a(0) a(0)  0.  b(0) b(0) b(0) b(0)  0.  c(0) c(0) c(0) c(0) 0.
!     4  .    .    .    .    .    .    .    .    .    .    .    .    .    .   .
!     5  .    .    .    .    .    .    .    .    .    .    .    .    .    .   .
!     6  .    .    .    .    .    .    .    .    .    .    .    .    .    .   .
!     7  .    .    .    .    .    .    .    .    .    .    .    .    .    .   .
!     8 a(5) a(6) a(5) a(5)  0.  b(5) b(5) b(5) b(5)  0.  c(5) c(5) c(5) c(5) 0.
!
!BH080118      call second(cpua)
      if(numz.eq.0) return
      izsum = 0
      do 5 i = 1,numz
        izsum = izsum + izstrp(i)
    5 continue
      if(izsum.eq.numz) return

      do 60 imp=1,numz
        do 20 k=1,15
          do 10 j=1,8
            arad(j,k,imp)=0.0
 10       continue
 20     continue
!**   find desired element in data table.
!
        open(unit=ncor,file='coronb',status='old')
!     if(length.eq.-1) go to 70
        do 30 j=1,10000
          read(ncor,1010,err=70,end=70)dum,dum2,iz,ia
 69       continue
          if(iz.eq.nz(imp)) go to 40
          read(ncor,1020) dum
          read(ncor,1020) dum
 30     continue
!**   found desired element.
 40     spec(imp)=dum
        amz(imp)=ia
!**   read data for this specie.
        k=1
 50     read(ncor,1020)(arad(j,k,imp),j=1,8),dum
        k=k+1
!990131        if(dum.eq.2hzb.and.k.lt.6) k=6
        if(dum.eq.'zb'.and.k.lt.6) k=6
!990131        if(dum.eq.2hzs.and.k.lt.11) k=11
!990131        if(k.ge.15.and.dum.ne.2hzs) go to 55
        if(dum.eq.'zs'.and.k.lt.11) k=11
        if(k.ge.15.and.dum.ne.'zs') go to 55
        go to 50
!
 55     continue
        close(unit=ncor)
!
 60   continue
!BH080118      call second(cpubb)
      cpu8 = cpu8 + cpubb - cpua
      return
!**   could not find desired element.
 70   if(nouthx.gt.0) write(nouthx,1030) nz(imp),ncor
      ihxbug = 1
 1010 format(1a2,15x,a4,2i4)
 1020 format(2e10.3,3e15.6,/,3e15.6,/,12x,1a2)
 1030 format(1x,9helement #,i10,19h not found on unit ,i5)
      return
      end subroutine hradin
!
!
      real(c_double) function hxgf(n,y)
      implicit none 
      integer n !input
      real(c_double) :: y !input
      real(c_double) :: rn, g0_rn, g1_rn, g2_rn !local
!*****************************************************************
      !YuP[2019-06-17] Replaced statement functions with scalars, below
      !g0(r)=.9935+.2328/r-.1296/(r*r)
      !g1(r)=-(.6282-.5598/r+.5299/(r*r))/r
      !g2(r)=(.3887-1.181/r+1.470/(r*r))/(r*r)
!
      rn=DBLE(n)
      if(n.eq.1)goto 11
      if(n.eq.2)goto 12
      g0_rn=.9935+.2328/rn-.1296/(rn*rn)
      g1_rn=-(.6282-.5598/rn+.5299/(rn*rn))/rn
      g2_rn=(.3887-1.181/rn+1.470/(rn*rn))/(rn*rn)
      hxgf= g0_rn +g1_rn/y +g2_rn/(y*y)
      return
!
 11   hxgf=1.1330-.4059/y+.07014/(y*y)
      return
!
 12   hxgf=1.0785-.2319/y+.02947/(y*y)
      return
      end function hxgf
!
!
      subroutine hxinit
      implicit none 
      save
!*****************************************************************
      integer ms,mc !local
      parameter (ms=21,mc=35)
      
      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf
        
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      
      real(c_double) :: f,ar !in /b3/
      common/b3/f(ms,mc),ar(ms+1,ms+1)
      
      real(c_double) :: en,dg,ae,be,de1,de2,ge1,ge2 !in /b4/
      common/b4/en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc), &
        de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
        
      integer iexcit,ilorent,mstate,ncont !in /b8/
      common/b8/iexcit,ilorent,mstate,ncont
      
      real(c_double) :: be1(mc)
      real(c_double) :: ryd,rrate,an1,en1,an2,en12,s,tot,ajsq,frac2s !local
      real(c_double) :: cpua,cpubb
      integer i,j,ni,nj,n1,n2
      data ryd/13.6/
!
!--   total radiation rate from level n2 to level n1:
      !YuP[2019-06-17] Replaced the statement function by a scalar, below
      !rrate(n2,n1)=(8.0323e9)*(((1./DBLE(n1))**2-(1./DBLE(n2))**2) &
      !  *(DBLE(n1)/DBLE(n2)))**2*f(n1,n2)
!
!BH080118      call second(cpua)
!
!--   tabulate oscillator strengths:
      call hxosc
!
      dg(1)=1.
      dg(2)=1.
      dg(3)=3.
      en(1)=ryd
      en(2)=ryd/4.
      en(3)=ryd/4.
      do 7 i=4,mc+1
        dg(i)=(DBLE(i-1))**2
        en(i)=ryd/dg(i)
    7 continue
!
      do 8 n1=1,mc
        an1=DBLE(n1)
!990131        be1(n1)=1.4*alog(an1)/an1-.7/an1-.51/(an1*an1)+1.16/(an1**3)
        be1(n1)=1.4*log(an1)/an1-.7/an1-.51/(an1*an1)+1.16/(an1**3) &
          -.55/(an1**4)
    8 continue
!
      do 9 n1=1,mstate
        an1=DBLE(n1)
        en1=ryd/(an1*an1)
        do 10 n2=n1+1,mc
          an2=DBLE(n2)
          en12=en1-ryd/(an2*an2)
          ae(n1,n2)=2.*ryd*f(n1,n2)/en12
          be(n1,n2)=4.*ryd*ryd*(1./(en12*en12)+4.*en1/(3.*en12**3) &
            +be1(n1)*en1*en1/(en12**4))/(an2**3)
          de1(n1,n2)=exp(-be(n1,n2)/ae(n1,n2))-.4*en12/ryd
          s=an2-an1
          de2(n1,n2)=exp(-be(n1,n2)/ae(n1,n2))+.06*s*s/(an1*an1*an2)
          ge1(n1,n2)=ryd*(8.+23.*(s/an1)**2) &
            /(8.+1.1*an2*s+.8/(s*s)+.4*(s-1.)*sqrt(an2*an2*an2/s))
          ge2(n1,n2)=ryd*(3.+11.*(s/an1)**2) &
            /(6.+1.6*an2*s+.3/(s*s)+.8*(s-.6)*sqrt(an2*an2*an2/s))
 10     continue
    9 continue
!
!--   radiation rates:
      do 205 i=1,mstate+1
        do 206 j=1,mstate+1
          ar(i,j)=0.
 206    continue
 205  continue
      if(krad.eq.0)go to 1000
!
      do 210 j=4,mstate+1
        nj=j-1
          !YuP[2019-06-17] Replacing the statement function:
          n1=1
          n2=nj
          !rrate(n2,n1)=
          rrate= (8.0323e9)*(((1./DBLE(n1))**2-(1./DBLE(n2))**2) &
                *(DBLE(n1)/DBLE(n2)))**2*f(n1,n2)
        ar(j,1)=rrate  !(nj,1)
        do 211 i=4,j-1
          ni=i-1
          !YuP[2019-06-17] Replacing the statement function:
          n1=ni
          n2=nj
          !rrate(n2,n1)=
          rrate= (8.0323e9)*(((1./DBLE(n1))**2-(1./DBLE(n2))**2) &
                *(DBLE(n1)/DBLE(n2)))**2*f(n1,n2)
          ar(j,i)=rrate !(nj,ni)
 211    continue
 210  continue
!
!--   2s to 1s:
      ar(2,1)=0.
!
!--   2p to 1s:
          !YuP[2019-06-17] Replacing the statement function:
          n1=1
          n2=2
          !rrate(n2,n1)=
          rrate= (8.0323e9)*(((1./DBLE(n1))**2-(1./DBLE(n2))**2) &
                *(DBLE(n1)/DBLE(n2)))**2*f(n1,n2)
      ar(3,1)=(4./3.)*rrate !(2,1)
!
!--   2p to 2s:
      ar(3,2)=0.
!
      do 220 j=4,mstate+1
        nj=j-1
        ajsq=(DBLE(nj))**2
!
!--   total radiation to 2s+2p:
          !YuP[2019-06-17] Replacing the statement function:
          n1=2
          n2=nj
          !rrate(n2,n1)=
          rrate= (8.0323e9)*(((1./DBLE(n1))**2-(1./DBLE(n2))**2) &
                *(DBLE(n1)/DBLE(n2)))**2*f(n1,n2)
        tot=rrate !(nj,2)
!
!--   fraction to 2s:
        frac2s=12.*(ajsq-1.)*(ajsq-4.)
        frac2s=frac2s/(frac2s+ajsq*(ajsq-4.)+32.*ajsq*(ajsq-1.))
!
        ar(j,2)=frac2s*tot
        ar(j,3)=(1.-frac2s)*tot
 220  continue
!
 1000 continue
!BH080118      call second(cpubb)
      cpu1 = cpubb - cpua
!
      return
      end subroutine hxinit
!
!
      subroutine hxosc
      implicit none 
!*****************************************************************
!--   calculate oscillator strengths.
!
      integer ms,mc !local
      parameter (ms=21,mc=35)
      
      real(c_double) :: f,ar !in /b3/
      common/b3/f(ms,mc),ar(ms+1,ms+1)

      integer iexcit,ilorent,mstate,ncont !in /b8/
      common/b8/iexcit,ilorent,mstate,ncont

      real(c_double) :: zero,one,pi,const,ai,aj,y
      integer i,j
!
      zero=0.d0
      one=1.d0
!      pi=acos(-one)
      pi=atan2(zero,-one)
      const=32./(sqrt(27.)*pi)
!
      do 10 i=1,mstate
        ai=DBLE(i)
        do 11 j=i+1,mc
          aj=DBLE(j)
          y=1.-(ai/aj)**2
          f(i,j)=(const*ai/((aj*y)**3))*hxgf(i,y)
 11     continue
 10   continue
      return
      end subroutine hxosc
!
!
      subroutine hxradi(teev,z,zsq,numimp,iz,izstrp,iwatch)
      implicit none 
      save
!*****************************************************************
!--   evaluate coronal z, zsq, and radiation.
      integer numimp,iwatch !input
      real(c_double) :: teev !input
      
      integer mz,mi !local
      parameter (mz=1,mi=2)

      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf

      real(c_double) :: arad !in locrad
      common/locrad/arad(8,15,mz)

      real(c_double) :: rad(mz),z(mz),zsq(mz),b(3)
      integer iz(mz),izstrp(mz)
      
      integer j,jj,jp,ihigh,imp,k,kk
      
      real(c_double) :: erad,zero,tekev,factor,t1,tl,bb,cc,dd
      real(c_double) :: cpua,cpubb
      
      data erad/1./
      zero=0.d0
!BH080118      call second(cpua)
      iwatch=0
      ihigh=0
      tekev=teev*1.0e-3
      do 60 imp=1,numimp
        if(izstrp(imp).eq.1) go to 51
!**   find temperature region. test for less than 5 intervals.
        do 10 j=1,5
          if(arad(2,j,imp).eq.zero) go to 10
          jj=j
          if(tekev.ge.arad(1,j,imp).and.tekev.lt.arad(2,j,imp)) go to 30
 10     continue
!**   temperature out of range.
        iwatch=1
        if(tekev.ge.arad(1,1,imp)) go to 20
!**   temperature too low. return zero power.
!     rad(imp)=0.0
!     z(imp)=0.0
!     zsq(imp)=0.0
!**   interpolation for low te (cdb):
!990131        t1=alog10(arad(1,1,imp))
        t1=log10(arad(1,1,imp))
        do 15 j=1,3
          jp=5*j-4
          b(j)=arad(3,jp,imp)+t1*(arad(4,jp,imp)+t1*(arad(5,jp,imp) &
            +t1*(arad(6,jp,imp)+t1*(arad(7,jp,imp) &
            +t1*arad(8,jp,imp)))))
 15     continue
        factor=(tekev/arad(1,1,imp))**erad
        rad(imp)=factor*(10.**b(1))
        z(imp)=factor*b(2)
        zsq(imp)=factor*b(3)
        go to 50
!**   temperature too high. compute value for tmax, then extrapolate.
!990131 20     tl=alog10(arad(2,jj,imp))
 20     tl=log10(arad(2,jj,imp))
        ihigh=1
!990131 30     if(ihigh.eq.0) tl=alog10(tekev)
 30     if(ihigh.eq.0) tl=log10(tekev)
!**   polonomial fits bb=sum(a(ik)*log(te)**k) etc.
        bb=0.0
        cc=0.0
        dd=0.0
        do 40 kk=1,6
          k=7-kk
          bb=bb*tl+arad((k+2),jj,imp)
          cc=cc*tl+arad((k+2),(jj+5),imp)
          dd=dd*tl+arad((k+2),(jj+10),imp)
 40     continue
        rad(imp)=0.0
        z(imp)=cc
        zsq(imp)=dd
        if(bb.lt.-38.0) go to 50
        rad(imp)=10.0**bb
!**   for te>tmax extrapolate as for pure bremsstrahlung.
        if(ihigh.eq.1) rad(imp)=rad(imp)*sqrt(tekev/arad(2,4,imp))
        ihigh=0
 50     continue
        go to 60

!...  fully stripped option [ izstrp(imp)=1 ]
 51     z(imp)   = iz(imp)
        zsq(imp) = iz(imp)**2

 60   continue
!
!BH080118      call second(cpubb)
      cpu2 = cpu2 + cpubb - cpua
!
      return
      end subroutine hxradi
!
!
      subroutine hxsve
      implicit none 
!*****************************************************************
!--   calculate the rate coefficients for collisions with electrons.
!
      real(c_double) :: zero,one,sinh1,ve,sqpi,facte,pi,arg !local
      real(c_double) :: cpua,cpubb !local
      integer i,j,kk,nsp1,ifail !local
      integer ms,mc,mz,mi !local
      parameter (ms=21,mc=35,mz=1,mi=2 )
      
      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf

      integer ns,nc,numi,numz,iz,izstrp !in /b0/
      real(c_double) :: er0,v0,te,ti,ami,deni,amz,denz,zcor,zsqcor,dene !in /b0/
      common/b0/ns,nc,numi,numz,iz(mz),izstrp(mz), &
        er0,v0,te,ti,ami(mi),deni(mi), &
        amz(mz),denz(mz),zcor(mz),zsqcor(mz),dene
      
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      
      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug
      
      real(c_double) :: en,dg,ae,be,de1,de2,ge1,ge2 !in /b4/
      common/b4/en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc), &
        de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
        
      real(c_double) :: cii,cei,ciz,cez,cie,cee,ccxi,ccxz !in /b6/
      common/b6/cii(ms+1,mi),cei(ms+1,ms+2,mi),ciz(ms+1,mz), &
        cez(ms+1,ms+2,mz),cie(ms+1),cee(ms+1,ms+2), &
        ccxi(ms+1,mi),ccxz(ms+1,mz)
        
      real(c_double) :: et(mc+1)
      real(c_double) :: xgl(64),wgl(64)
      character*2 typent
!
      sinh1(arg)=sinh(arg)/arg
!
!BH080118      call second(cpua)
!
!--   initialize:
      if(istart.eq.0)goto 5
      zero=0.d0
      one=1.d0
!      pi=acos(-one)
      pi=atan2(zero,-one)
      sqpi=sqrt(pi)
!
!     e-impact thresholds (temp.)
      et(1)=(13.723/13.6)*en(1)
      do 1 i=2,mc+1
        et(i)=en(i)
    1 continue
!
      ifail=0
      typent = 'gl'
!CMG  changed 11/13/2017         call d01bbf(typent,ngl,wgl,xgl,ifail)
      call d01bbf_cql3d(typent,ngl,wgl,xgl,ifail)
!@    call d01bbf(d01bax,0.,1.,0,ngl,wgl,xgl,ifail)
      if(ifail.eq.0)goto 5
      if(nouthx.gt.0) write(nouthx,3939) ifail
 3939 format(' ??? d01bbf [ngl]  call error: ifail= ',i4)
      ihxbug = 3
      return
!
    5 continue
      nsp1=ns+1
!
      do 205 i=1,ns
        cie(i)=0.
        do 206 j=i+1,nsp1
          cee(i,j)=0.
 206    continue
 205  continue
      ve=5.931e7*sqrt(te)
!
      if(kdene.eq.0)go to 1000
      goto(220,210),ksve+1
!
!--   gauss-laguerre integrations:
 210  continue
      facte=(2./sqpi)*ve*exp(-v0*v0/(ve*ve))
      do 211 kk=1,ngl
        do 212 i=1,ns
          cie(i)=cie(i) &
            +wgl(kk)*(xgl(kk)+et(i)/te) &
            *sinh1(2.*v0*sqrt(xgl(kk)+et(i)/te)/ve) &
            *sief(i,et(i)+te*xgl(kk)) &
            *facte*exp(-et(i)/te)
          do 213 j=i+1,ns
            cee(i,j)=cee(i,j) &
              +wgl(kk)*(xgl(kk)+(et(i)-et(j))/te) &
              *sinh1(2.*v0*sqrt(xgl(kk)+(et(i)-et(j))/te)/ve) &
              *seef(i,j,et(i)-et(j)+te*xgl(kk)) &
              *facte*exp(-(et(i)-et(j))/te)
 213      continue
          do 214 j=nsp1,nc
            cee(i,nsp1)=cee(i,nsp1) &
              +wgl(kk)*(xgl(kk)+(et(i)-et(j))/te) &
              *sinh1(2.*v0*sqrt(xgl(kk)+(et(i)-et(j))/te)/ve) &
              *seef(i,j,et(i)-et(j)+te*xgl(kk)) &
              *facte*exp(-(et(i)-et(j))/te)
 214      continue
 212    continue
 211  continue
      go to 1000
!
 220  continue
!--   gauss-laguerre averages for ionization of 1s, ionization of 2s,
!     and excitations among 1s, 2s, and 2p:
      facte=(2./sqpi)*ve*exp(-v0*v0/(ve*ve))
      do 221 kk=1,ngl
        do 222 i=1,min0(2,ns)
          cie(i)=cie(i) &
            +wgl(kk)*(xgl(kk)+et(i)/te) &
            *sinh1(2.*v0*sqrt(xgl(kk)+et(i)/te)/ve) &
            *sief(i,et(i)+te*xgl(kk)) &
            *facte*exp(-et(i)/te)
          do 2221 j=i+1,min0(3,ns)
            cee(i,j)=cee(i,j) &
              +wgl(kk)*(xgl(kk)+(et(i)-et(j))/te) &
              *sinh1(2.*v0*sqrt(xgl(kk)+(et(i)-et(j))/te)/ve) &
              *seef(i,j,et(i)-et(j)+te*xgl(kk)) &
              *facte*exp(-(et(i)-et(j))/te)
 2221     continue
 222    continue
 221  continue
!
!--   rates from vriens & smeets for other reactions:
!     ionizations--
      do 223 i=3,ns
        cie(i)=cief(i,te)
 223  continue
!
!     excitations--
      do 224 i=1,ns
        do 225 j=i+1,ns
          if(j.le.3)goto 225
          cee(i,j)=ceef(i,j,te)
 225    continue
        do 226 j=nsp1,nc
          cee(i,nsp1)=cee(i,nsp1)+ceef(i,j,te)
 226    continue
 224  continue
!
 1000 continue
!BH080118      call second(cpubb)
      cpu5 = cpu5 + cpubb - cpua
!
      return
      end subroutine hxsve
!
!
      subroutine hxsvi
      implicit none 
      save
!*****************************************************************
!--   calculate the rate coefficients for collisions with ions
!     and impurities.
!
      real(c_double) :: zero,one,pi,sqpi,xp,xm,s,ei,ez,zz !local
      real(c_double) :: facti,factz,er0old !local
      real(c_double) :: cpua1,cpua2,cpua,cpubb !local
      integer i,j,kk,ki,kz,nsp1,ngh1,ifail,ister0 !local
      integer ms,mc,mz,mi !local
      parameter (ms=21,mc=35,mz=1,mi=2)
      
      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf

      integer ns,nc,numi,numz,iz,izstrp !in /b0/
      real(c_double) :: er0,v0,te,ti,ami,deni,amz,denz,zcor,zsqcor,dene !in /b0/
      common/b0/ns,nc,numi,numz,iz(mz),izstrp(mz), &
        er0,v0,te,ti,ami(mi),deni(mi), &
        amz(mz),denz(mz),zcor(mz),zsqcor(mz),dene

      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl

      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug
      
      real(c_double) :: cii,cei,ciz,cez,cie,cee,ccxi,ccxz !in /b6/
      common/b6/cii(ms+1,mi),cei(ms+1,ms+2,mi),ciz(ms+1,mz), &
        cez(ms+1,ms+2,mz),cie(ms+1),cee(ms+1,ms+2), &
        ccxi(ms+1,mi),ccxz(ms+1,mz)

      real(c_double) :: vi(mi),vz(mz)
      real(c_double) :: xgh(64),wgh(64),w1(32)
      character*2 typent
      data ister0/1/
!
!--   initialize:
!BH080118      call second(cpua)
!
      if(istart.eq.0)goto 5
      zero=0.d0
      one=1.d0
!      pi=acos(-one)
      pi=atan2(zero,-one)
      sqpi=sqrt(pi)
!
      ngh1=(ngh+1)/2
      do 3 i=1,ngh1
        w1(i)=2.
    3 continue
      if(2*ngh1-ngh .eq. 1)w1(ngh1)=1.
      ifail = 0
      typent = 'gh'
!CMG  changed 11/13/2017      call d01bbf(typent,ngh,wgh,xgh,ifail)
      call d01bbf_cql3d(typent,ngh,wgh,xgh,ifail)
!@    call d01bbf(d01baw,0.,1.,0,ngh,wgh,xgh,ifail)
      if(ifail.eq.0)goto 5
      if(nouthx.gt.0) write(nouthx,3939) ifail
 3939 format(' ??? d01bbf [ngh]  call error: ifail= ',i4)
      ihxbug = 2
      return
!
    5 continue
      nsp1=ns+1
      print *,'****************************************'
      print *,' warning er0old not defined in zfreya'
      print *,'****************************************'
      if(er0.ne.er0old) ister0=1
!
!--   ion reactions-----------------------------------------------
!
!BH080118      call second(cpua1)
!
!--   first zero out arrays:
      if( (ister0.eq.1 .or. istart.eq.1) .or. ksvi.eq.1) then
        do 10 ki=1,numi
          do 11 i=1,ns
            cii(i,ki)=0.
            ccxi(i,ki)=0.
            do 12 j=i+1,nsp1
              cei(i,j,ki)=0.
 12         continue
 11       continue
 10     continue
      endif
!
      do 13 ki=1,numi
        vi(ki)=1.3841e6*sqrt(ti/ami(ki))
 13   continue

      if(kdeni.eq.0)goto 100
      if(ksvi.eq.0)goto 40
!
!--   gauss-hermite integrations:
      do 20 kk=1,ngh1
        do 21 ki=1,numi
          xp=(xgh(kk)+v0/vi(ki))**2
          xm=(xgh(kk)-v0/vi(ki))**2
          s=1.
          if(xgh(kk)-v0/vi(ki) .lt. 0.)s=-1.
          ei=ti/ami(ki)
          do 22 i=1,ns
            cii(i,ki)=cii(i,ki)+wgh(kk)*w1(kk) &
              *(xp*siif(i,ei*xp)-s*xm*siif(i,ei*xm))
            ccxi(i,ki)=ccxi(i,ki)+wgh(kk)*w1(kk) &
              *(xp*scxif(i,ei*xp)-s*xm*scxif(i,ei*xm))
            do 23 j=i+1,ns
              cei(i,j,ki)=cei(i,j,ki)+wgh(kk)*w1(kk) &
                *(xp*sezf(i,j,one,ei*xp)-s*xm*sezf(i,j,one,ei*xm))
 23         continue
            do 24 j=nsp1,nc
              cei(i,nsp1,ki)=cei(i,nsp1,ki)+wgh(kk)*w1(kk) &
                *(xp*sezf(i,j,one,ei*xp)-s*xm*sezf(i,j,one,ei*xm))
 24         continue
 22       continue
 21     continue
 20   continue
!
      do 30 ki=1,numi
        facti=vi(ki)*vi(ki)/(2.*sqpi*v0)
        do 31 i=1,ns
          cii(i,ki)=facti*cii(i,ki)
          ccxi(i,ki)=facti*ccxi(i,ki)
          do 32 j=i+1,nsp1
            cei(i,j,ki)=facti*cei(i,j,ki)
 32       continue
 31     continue
 30   continue
      goto 100
!
!--   simple multiplication instead of maxwellian averaging:
!--   cii & cei are independent of plasma parameters in this
!--   approximation       so as long as er0 doesn't change, use
!--   the old values.
!--   ister0=1 ==> restart everything
!--   =0 ==> cruise in no-update mode
 40   continue
!--   redo everything if er0 changes
      if( (ister0.eq.1 .or. istart.eq.1) .or. ksvi.eq.1) then
        do 41 ki=1,numi
          do 411 i=1,ns
            cii(i,ki)=siif(i,er0)*v0
            ccxi(i,ki)=scxif(i,er0)*v0
            do 42 j=i+1,ns
              cei(i,j,ki)=sezf(i,j,one,er0)*v0
 42         continue
            do 43 j=nsp1,nc
              cei(i,nsp1,ki)=cei(i,nsp1,ki)+sezf(i,j,one,er0)*v0
 43         continue
 411      continue
 41     continue
      endif
 100  continue
!BH080118      call second(cpua2)
      cpuii = cpuii + cpua2 - cpua1
!
!--   impurity reactions------------------------------------------
!
!--   zero out arrays:
      do 110 kz=1,numz
        vz(kz)=1.3841e6*sqrt(ti/amz(kz))
        do 111 i=1,ns
          ciz(i,kz)=0.
          ccxz(i,kz)=0.
          do 112 j=i+1,nsp1
            cez(i,j,kz)=0.
 112      continue
 111    continue
 110  continue
!
      if(kdenz.eq.0)go to 1000
      if(ksvz.eq.0)goto 140
!
!--   gauss-hermite integrations:
      do 120 kk=1,ngh1
        do 121 kz=1,numz
          xp=(xgh(kk)+v0/vz(kz))**2
          xm=(xgh(kk)-v0/vz(kz))**2
          s=1.
          if(xgh(kk)-v0/vz(kz) .lt. 0.)s=-1.
          ez=ti/amz(kz)
          zz=zcor(kz)
          do 122 i=1,ns
            ciz(i,kz)=ciz(i,kz)+wgh(kk)*w1(kk) &
              *(xp*sizf(i,zz,ez*xp)-s*xm*sizf(i,zz,ez*xm))
            ccxz(i,kz)=ccxz(i,kz)+wgh(kk)*w1(kk) &
              *(xp*scxzf(i,zz,ez*xp)-s*xm*scxzf(i,zz,ez*xm))
            do 123 j=i+1,ns
              cez(i,j,kz)=cez(i,j,kz)+wgh(kk)*w1(kk) &
                *(xp*sezf(i,j,zz,ez*xp)-s*xm*sezf(i,j,zz,ez*xm))
 123        continue
            do 124 j=nsp1,nc
              cez(i,nsp1,kz)=cez(i,nsp1,kz)+wgh(kk)*w1(kk) &
                *(xp*sezf(i,j,zz,ez*xp)-s*xm*sezf(i,j,zz,ez*xm))
 124        continue
 122      continue
 121    continue
 120  continue
!
      do 130 kz=1,numz
        factz=vz(kz)**2/(2.*sqpi*v0)
        do 131 i=1,ns
          ciz(i,kz)=factz*ciz(i,kz)
          ccxz(i,kz)=factz*ccxz(i,kz)
          do 132 j=i+1,nsp1
            cez(i,j,kz)=factz*cez(i,j,kz)
 132      continue
 131    continue
 130  continue
      go to 1000
!
!--   simple multiplications:
 140  continue
      do 141 kz=1,numz
        do 142 i=1,ns
          ciz(i,kz)=sizf(i,zcor(kz),er0)*v0
          ccxz(i,kz)=scxzf(i,zcor(kz),er0)*v0
          do 143 j=i+1,ns
            cez(i,j,kz)=sezf(i,j,zcor(kz),er0)*v0
 143      continue
          do 144 j=nsp1,nc
            cez(i,nsp1,kz)=cez(i,nsp1,kz)+sezf(i,j,zcor(kz),er0)*v0
 144      continue
 142    continue
 141  continue
!
 1000 continue
!BH080118      call second(cpubb)
      cpu4 = cpu4 + cpubb - cpua
      cpuiz = cpuiz + cpubb - cpua2
!
      ister0=0
      er0old = er0
      return
      end subroutine hxsvi
!
!
!$$$ Replaced by Kinsey, when adding ADAS cross-sections to old (1990's)
!$$$ CQL3D/NFREYA  (June-July 2013).
!$$$      subroutine inject(atw, codeid, drutpi, droti, dri, dzi,
!$$$     *  elongi, ib, ie, ke, kz, kbe, ki, mfm1, mim1, mjm1,
!$$$     *  newpar, psiax, psi,
!$$$     *  r, rmajor, rin, rmax, sgxn, sgxnmi,
!$$$     *  x0, y0, z0, vx0, vy0, vz0, vbeam, z, zax, zmin, zmax,
!$$$     *  izone, pzone, rzone, rpos, xpos, ypos, zpos)
!$$$      implicit integer (i-n), real*8 (a-h,o-z)
!$$$      save
!$$$!
!$$$!     this subroutine follows the particle from the pivot point
!$$$!     into, through, or around the plasma.
!$$$!
!$$$      dimension psi(ki,*), sgxn(4,kz,kbe,*), sgxnmi(ke,*),
!$$$     *  r(*), z(*), vbeam(ke,*)
!$$$      character codeid*8
!$$$
!$$$      integer iflag
!$$$      data       seed0    /0.0d0/
!$$$
!$$$!      iflag=0
!$$$!
!$$$!%OS
!$$$cbh   (930920) print *,' WARNING: isol not defined in inject'
!$$$!%OS
!$$$cbh   (930920) if(newpar.eq.0 .and. isol.eq.0) go to 140
!$$$      if(newpar.eq.0) go to 100
!$$$!
!$$$!     calculate times for particle to enter and exit toroidal box
!$$$!     surrounding plasma
!$$$!
!$$$      call timtor(rin, rmax, x0, y0, z0, vx0, vy0, vz0, zmin, zmax,
!$$$     *  tenter, texit)
!$$$!      write(*,*) 'inject:rin,rmax,x0,y0,z0,vx0,vy0,vz0,zmin,zmax',
!$$$!     +     ',tenter,texit',
!$$$!     +     rin,rmax,x0,y0,z0,vx0,vy0,vz0,zmin,zmax,tenter,texit
!$$$
!$$$!      write(*,*)'inject: tenter',tenter
!$$$      if(tenter.le.-1.e10) go to 140
!$$$!      write(*,*)'inject: tenter',tenter
!$$$!
!$$$!     advance particle to edge of box
!$$$!
!$$$      x0=x0+vx0*tenter
!$$$      y0=y0+vy0*tenter
!$$$      z0=z0+vz0*tenter
!$$$!
!$$$!     set coordinates and time for entering box
!$$$!
!$$$ 100  continue
!$$$      xpos=x0
!$$$      ypos=y0
!$$$      zpos=z0
!$$$      tt=tenter
!$$$!      write(*,*)'inject:xpos,ypos,zpos,tt',xpos,ypos,zpos,tt
!$$$!
!$$$!     follow particle into plasma
!$$$!
!$$$ 110  continue
!$$$!990131      dfac=-alog(ranf())
!$$$!000114*Not portable*      dfac=-log(drand(iflag))
!$$$      dfac=-log(RANDOM_my(seed0))
!$$$      tstep=dfac*sgxnmi(ie,ib)/vbeam(ie,ib)
!$$$      tt=tt+tstep
!$$$!      write(*,*) 'dfac,ie,ib,sgxnmi(ie,ib),vbeam(ie,ib),tstep',
!$$$!     +             dfac,ie,ib,sgxnmi(ie,ib),vbeam(ie,ib),tstep
!$$$      if(tt.ge.texit) go to 140
!$$$      xpos=xpos+vx0*tstep
!$$$      ypos=ypos+vy0*tstep
!$$$      zpos=zpos+vz0*tstep
!$$$      rpos=sqrt(xpos**2+ypos**2)
!$$$!      write(*,*)'inject:xpos,ypos,zpos,rpos',xpos,ypos,zpos,rpos
!$$$
!$$$!
!$$$!     determine zone in which particle collides for 'onedee' geometry
!$$$!
!$$$      if(codeid.ne.'onedee') go to 120
!$$$      rzone2 = (rpos-rmajor)**2 + (elongi*(zpos-zax))**2
!$$$      rzone = sqrt(rzone2)
!$$$      izone = rzone*droti + 1.
!$$$      go to 130
!$$$!
!$$$!     determine zone in which particle collides for general geometry;
!$$$!     use bilinear interpolation away from magnetic axis and
!$$$!     biquadratic interpolation near the axis
!$$$!
!$$$ 120  i=(rpos-r(1))*dri+1.
!$$$      j=(zpos-z(1))*dzi+1.
!$$$      if(i.gt.mim1) i=mim1
!$$$      if(j.gt.mjm1) j=mjm1
!$$$!990131      psix  = amin1(psi(i,j), psi(i+1,j), psi(i,j+1), psi(i+1,j+1))
!$$$      psix  = min(psi(i,j), psi(i+1,j), psi(i,j+1), psi(i+1,j+1))
!$$$      ptest = (psix-psiax)*(drutpi/mfm1)**2
!$$$      if(ptest.lt.0.02) go to 124
!$$$      area1=(rpos-r(i))*(zpos-z(j))
!$$$      area2=(r(i+1)-rpos)*(zpos-z(j))
!$$$      area3=(r(i+1)-rpos)*(z(j+1)-zpos)
!$$$      area4=(rpos-r(i))*(z(j+1)-zpos)
!$$$      pzone=(area3*psi(i,j)+area4*psi(i+1,j)+area1*psi(i+1,j+1)
!$$$     1  +area2*psi(i,j+1))*dri*dzi
!$$$      go to 126
!$$$ 124  call pfit(psi(i-1,j-1), r(i-1), z(j-1), rpos, zpos, ki, pzone,
!$$$     *  dum, dum)
!$$$!990131 126  pzone = amax1(pzone,psiax)
!$$$ 126  pzone = max(pzone,psiax)
!$$$      izone = sqrt(pzone-psiax)*drutpi + 1.
!$$$!
!$$$!     if particle has psuedo-collision, continue following particle;
!$$$!     if particle has real collision, return
!$$$!
!$$$ 130  if(izone.gt.mfm1) go to 110
!$$$!990131      if(ranf().gt.sgxn(izone,ie,ib)*sgxnmi(ie,ib))
!$$$!000114*Not portable*      if(drand(iflag).gt.sgxn(izone,ie,ib)*sgxnmi(ie,ib))
!$$$      if(RANDOM_my(seed0).gt.sgxn(4,izone,ie,ib)*sgxnmi(ie,ib))
!$$$     1  go to 110
!$$$      return
!$$$!
!$$$!     set flag to indicate that particle hit wall
!$$$!
!$$$ 140  izone = mfm1+1
!$$$      return
!$$$      end
!
!
      subroutine inter(dxi,iroot,nxm1,x0,ytab,x,y)
      implicit none 
!
!     this subroutine performs a fast linear interpolation (or extrapo-
!     lation) for y vs. x for either of two cases:
!     1.  if iroot.eq.0, then x is uniformly spaced with
!     spacing 1/dxi,
!     2.  if iroot.ne.0, then sqrt(x) is uniformly spaced with
!     spacing 1/dxi.
!
      real(c_double) :: ytab(*)
      real(c_double) :: dxi,x0,x,y !input
      integer iroot,nxm1 !input
      integer i !local
      real(c_double) :: xx,wt,xtab !local

      xtab = abs(x-x0)
      if(iroot.ne.0) xtab = sqrt(xtab)
      xx = xtab*dxi + 1.
      i  = xx
      i  = min0(i,nxm1)
      wt = xx - i
      y  = (1.-wt)*ytab(i) + wt*ytab(i+1)
      return
      end subroutine inter
!
!
      subroutine itorb(c1,c2,c3,drtabi,mfm1,psival,ptab,rmajor,r,fcap)
      implicit none 
      real(c_double) :: ptab(*)
      real(c_double) :: c1,c2,c3,drtabi,psival,rmajor,r !input
      real(c_double) :: fcap !output
      integer mfm1 !input
!
      call inter(drtabi,0,mfm1,rmajor,ptab,r,psival)
      fcap = r**2 - c1*r - c2**2*psival**2 - 2.*c2*c3*psival - c3**2
      fcap = fcap/r
      return
      end subroutine itorb
!
!
      subroutine logint(x,y)
      implicit none 
      real(c_double) :: x !input
      real(c_double) :: y !output
      integer mdat,mdatm,im,i0,ip !local
      real(c_double) :: ylog,ylogm,ylog0,ylogp,dm,d0,dp,d0m,d0p !local
      real(c_double) :: dp0,dpm,dm0,dmp,facm,fac0,facp !local
      save
!cc
!cc   interpolates y(x) quadratically and logarithmically
!cc
      real(c_double) :: xdat(15),ydat(15)
!cc
      data xdat/4.,6.,8.,10.,20.,30.,40.,60.,80.,100., &
        200.,300.,400.,600.,800./
      data ydat/8.95e-01,8.75e-01,8.70e-01,8.65e-01,8.20e-01, &
        7.25e-01,6.25e-01,4.40e-01,2.90e-01,1.90e-01,2.40e-02, &
        5.25e-03,1.20e-03,1.60e-04,5.40e-05/
!cc
      mdat=15
      mdatm=mdat-1
      do 10 i0=2,mdatm
        if(xdat(i0).ge.x) go to 11
 10   continue
 11   if(i0.gt.mdatm) i0=mdatm
      im=i0-1
      ip=i0+1
!990131      ylogp=alog(ydat(ip))
!990131      ylog0=alog(ydat(i0))
!990131      ylogm=alog(ydat(im))
      ylogp=log(ydat(ip))
      ylog0=log(ydat(i0))
      ylogm=log(ydat(im))
      dm=x-xdat(im)
      d0=x-xdat(i0)
      dp=x-xdat(ip)
      d0m=xdat(i0)-xdat(im)
      dp0=xdat(ip)-xdat(i0)
      dpm=xdat(ip)-xdat(im)
      dm0=-d0m
      d0p=-dp0
      dmp=-dpm
      facm=d0*dp/(dm0*dmp)
      fac0=dm*dp/(d0m*d0p)
      facp=dm*d0/(dpm*dp0)
      ylog=facm*ylogm+fac0*ylog0+facp*ylogp
      y=exp(ylog)
!cc
      return
      end subroutine logint
!
!
      subroutine lorent(v0,bperp,nc,ns)
      implicit none 
      real(c_double) :: v0,bperp !input
      integer nc,ns !input
      integer i,il,i1,n,ni,ncrit !local
      real(c_double) :: fl,w0,almin,almax,expl,zero,one,ep,eight !local
      real(c_double) :: an3,ani3,el !local
      real(c_double) :: cpua,cpubb !local
      save
!*****************************************************************
      integer ms !local
      parameter (ms=21)

      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf
        
      real(c_double) :: al !in /b5/
      common/b5/al(ms+1)
      
      integer iexcit,ilorent,mstate,ncont !in /b8/
      common/b8/iexcit,ilorent,mstate,ncont
      
      real(c_double) :: sl(30)
      data w0/4.1341e16/,almin/1.e-10/,almax/1.e15/,expl/4./
!     data sl/
!     *1.00000000e+00,2.50000000e-01,2.77777778e-02,1.73611111e-03,
!     *6.94444444e-05,1.92901235e-06,3.93675989e-08,6.15118733e-10,
!     *7.59405843e-12,7.59405843e-14,6.27608135e-16,4.35838982e-18,
!     *2.57892889e-20,1.31578005e-22,5.84791131e-25,2.28434036e-27,
!     *7.90429189e-30,2.43959626e-32,6.75788439e-35,1.68947110e-37,
!     *3.83100022e-40,7.91528971e-43,1.49627405e-45,2.59769800e-48,
!     *4.15631680e-51,6.14839763e-54,8.43401595e-57,1.07576734e-59,
!     *1.27915260e-62,1.42128067e-65/
      data sl/ &
        1.00000000e+00,2.50000000e-01,2.77777778e-02,1.73611111e-03, &
        6.94444444e-05,1.92901235e-06,3.93675989e-08,6.15118733e-10, &
        7.59405843e-12,7.59405843e-14,6.27608135e-16,4.35838982e-18, &
        2.57892889e-20,1.31578005e-22,5.84791131e-25,2.28434036e-27, &
        7.90429189e-30,2.43959626e-32,6.75788439e-35,1.68947110e-37, &
        10*0.0/
!
      fl(n,an3,el)=((4./(an3*el))**(2*n-1)) &
        *exp(-2./(3.*an3*el))/an3
!     input:
!     v0
!     bperp
!     nc
!     output:
!     ns
!     al(i),i=1,ns
!
!BH080118      call second(cpua)
!
      zero=0.d0
      if(nc.ge.2)goto 1
      ns=1
      al(1)=0.
      go to 1000
!
    1 continue
      if(ilorent.eq.0) then
        ns=mstate+1
        do 2 il=1,ns
          al(il)=0.
 2      continue
        go to 1000
      endif
!
      if(bperp.eq.zero)goto 3
!
!--   ep= abs(v x b)/(electric field at 1st bohr radius)
      ep=v0*bperp/5.1417e13
      ncrit=.5/(ep**(1./expl))
      if(ncrit.le.1)ns=1
      if(ncrit.ge.2)ns=ncrit+1
      if(ns.le.(mstate+1))goto 4
!
    3 continue
      ns=mstate+1
!
    4 continue
      do 10 i=1,ns
        if(i.ge.4)goto 13
        goto(11,12,12),i
!990131 11     al(1)=fl(1,1.,ep)
 11     one=1.d0
        al(1)=fl(1,one,ep)
        goto 10
!990131 12     al(i)=fl(2,8.,ep)*.5
 12     eight=8.
        al(i)=fl(2,eight,ep)*.5
        goto 10
 13     ni=i-1
        ani3=(DBLE(ni))**3
        al(i)=fl(ni,ani3,ep)*sl(ni)
 10   continue
!
!--   normalize rates       zero out rates which are too small:
      do 15 i=1,ns
        al(i)=al(i)*w0
        if(al(i) .le. almin)al(i)=0.
 15   continue
!
!--   if some rates are too high, truncate system:
      do 20 i=1,ns
        i1=i
        if(al(i) .gt. almax)goto 21
 20   continue
      go to 1000
 21   continue
      ns=i1-1
 1000 continue
!BH080118      call second(cpubb)
      cpu3 = cpu3 + cpubb - cpua
!
      return
      end subroutine lorent
!
!
      real(c_double) function polyf(c,n,x)
      implicit none
      integer n !input
!@    real function polyf(c,n,x)
!*****************************************************************
!--   evaluate the polynomial c(1)+c(2)*x+...+c(n)*x**(n-1).

      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf
        
      real(c_double) :: x,c(n) !input
      
      integer j !local
      real(c_double) :: polysv !local
      real(c_double) :: cpua,cpubb !local
      
!BH080118      call second(cpua)
      polyf=c(n)
      do 10 j=1,n-1
        polyf=polyf*x+c(n-j)
        polysv = polyf
 10   continue
!BH080118      call second(cpubb)
      cpuplf = cpuplf + cpubb - cpua
      return
      end function polyf
!
!
      subroutine matri
      implicit none 
      real(c_double) :: cpua,cpubb !local
      integer i,j,jp,nsp1,ki,kz !local
!*****************************************************************
!--   form the matrix of rate coefficients which occurs
!     in the system of rate equations.  evaluate deexcitation via
!     detailed balance.
!--   vectorized.
!
      integer ms,mc,mz,mi !local
      parameter (ms=21,mc=35,mz=1,mi=2 )

      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf
        
      integer ns,nc,numi,numz,iz,izstrp !in /b0/
      real(c_double) :: er0,v0,te,ti,ami,deni,amz,denz,zcor,zsqcor,dene !in /b0/
      common/b0/ns,nc,numi,numz,iz(mz),izstrp(mz), &
        er0,v0,te,ti,ami(mi),deni(mi), &
        amz(mz),denz(mz),zcor(mz),zsqcor(mz),dene
        
      real(c_double) :: f,ar !in /b3/
      common/b3/f(ms,mc),ar(ms+1,ms+1)
      
      real(c_double) :: en,dg,ae,be,de1,de2,ge1,ge2 !in /b4/
      common/b4/en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc), &
        de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
        
      real(c_double) :: al !in /b5/
      common/b5/al(ms+1)
      
      real(c_double) :: cii,cei,ciz,cez,cie,cee,ccxi,ccxz !in /b6/
      common/b6/cii(ms+1,mi),cei(ms+1,ms+2,mi),ciz(ms+1,mz), &
        cez(ms+1,ms+2,mz),cie(ms+1),cee(ms+1,ms+2), &
        ccxi(ms+1,mi),ccxz(ms+1,mz)
        
      real(c_double) :: q !in /b7/
      common/b7/q(ms+1,ms+1)
      
      real(c_double) :: dge(ms+1),dgi(ms+1)
!
!cc*   dump c arrays
!     write(60,3310) cii
!c3310 format(/,' cii= ',/,(1x,1p5e11.3))
!     write(60,3315) ccxii
!c3315 format(/,' ccxii= ',/,(1x,1p5e11.3))
!     write(60,3320) cei
!c3320 format(/,' cei= ',/,(1x,1p5e11.3))
!     write(60,3330) ciz
!c3330 format(/,' ciz= ',/,(1x,1p5e11.3))
!     write(60,3335) ccxz
!c3335 format(/,' ccxz= ',/,(1x,1p5e11.3))
!     write(60,3340) cez
!c3340 format(/,' cez= ',/,(1x,1p5e11.3))
!     write(60,3350) cie
!c3350 format(/,' cie= ',/,(1x,1p5e11.3))
!     write(60,3360) cee
!c3360 format(/,' cee= ',/,(1x,1p5e11.3))

!BH080118      call second(cpua)
!
      nsp1=ns+1
      do 10 i=1,ns
        do 11 j=1,ns
          q(i,j)=0.
 11     continue
 10   continue
!
!--   detailed-balance factors (note that en's are positive):
      do 15 i=1,ns
        dge(i)=dg(i)*exp(en(i)/te)
        dgi(i)=dg(i)*exp(en(i)/ti)
 15   continue
!
!--   rates due to collisions with electrons
!     radiation rates       lorentz rates:
!
      do 20 j=1,ns
        q(j,j)= &
          -cie(j)*dene &
          -al(j)
        do 21 jp=j+1,nsp1
          q(j,j)=q(j,j) &
            -cee(j,jp)*dene
 21     continue
        do 22 jp=1,j-1
          q(j,j)=q(j,j) &
            -(dge(jp)/dge(j))*cee(jp,j)*dene &
            -ar(j,jp)
 22     continue
        do 23 i=j+1,ns
          q(i,j)=q(i,j) &
            +cee(j,i)*dene
 23     continue
        do 24 i=1,j-1
          q(i,j)=q(i,j) &
            +(dge(i)/dge(j))*cee(i,j)*dene &
            +ar(j,i)
 24     continue
 20   continue
!
!--   add rates due to collisions with ions:
      do 30 ki=1,numi
        do 301 j=1,ns
          q(j,j)=q(j,j) &
            -cii(j,ki)*deni(ki)
          do 31 jp=j+1,nsp1
            q(j,j)=q(j,j) &
              -cei(j,jp,ki)*deni(ki)
 31       continue
          do 32 jp=1,j-1
            q(j,j)=q(j,j) &
              -(dgi(jp)/dgi(j))*cei(jp,j,ki)*deni(ki)
 32       continue
          do 33 i=j+1,ns
            q(i,j)=q(i,j) &
              +cei(j,i,ki)*deni(ki)
 33       continue
          do 34 i=1,j-1
            q(i,j)=q(i,j) &
              +(dgi(i)/dgi(j))*cei(i,j,ki)*deni(ki)
 34       continue
 301    continue
 30   continue
!
!--   add rates due to collisions with impurities:
      do 40 kz=1,numz
        do 401 j=1,ns
          q(j,j)=q(j,j) &
            -ciz(j,kz)*denz(kz)
          do 41 jp=j+1,nsp1
            q(j,j)=q(j,j) &
              -cez(j,jp,kz)*denz(kz)
 41       continue
          do 42 jp=1,j-1
            q(j,j)=q(j,j) &
              -(dgi(jp)/dgi(j))*cez(jp,j,kz)*denz(kz)
 42       continue
          do 43 i=j+1,ns
            q(i,j)=q(i,j) &
              +cez(j,i,kz)*denz(kz)
 43       continue
          do 44 i=1,j-1
            q(i,j)=q(i,j) &
              +(dgi(i)/dgi(j))*cez(i,j,kz)*denz(kz)
 44       continue
 401    continue
 40   continue
!
!BH080118      call second(cpubb)
      cpu6 = cpu6 + cpubb - cpua
!
      return
      end subroutine matri
!
!
      subroutine newgrid(xold,yold,nold,xnew,ynew,nnew)
      implicit none 
      integer nold,nnew !input
      integer iold,inew,ipm1 !local
      real(c_double) :: sex,s,tol,del !local
      save
!     ---------------------------------------------------------------
!     convert to new grid
!     ---------------------------------------------------------------
      real(c_double) :: xold(nold),yold(nold),xnew(nnew),ynew(nnew)
      data tol/1.0e-20/
!
      iold=1
      inew=1
      sex=(yold(nold)-yold(nold-1))/(xold(nold)-xold(nold-1))
!     ---------------------------------------------------------------
!     interpolate to new grid that is compute ynew(xnew(inew))
!     ---------------------------------------------------------------
 2000 if(inew.gt.nnew) return
 2100 if(iold.gt.nold) go to 2200
      if(xnew(inew).gt.xold(1)) go to 2105
      ynew(inew)=yold(1)
      inew=inew+1
      go to 2000
 2105 continue
      del=reldif(xold(iold),xnew(inew))
      if(del.gt.tol) go to 2110
      ynew(inew)=yold(iold)
      inew=inew+1
      go to 2000
 2110 if(xold(iold).gt.xnew(inew)) go to 2120
      iold=iold+1
      go to 2100
 2120 ipm1=iold-1
      s=(yold(iold)-yold(ipm1))/(xold(iold)-xold(ipm1))
      ynew(inew)=yold(ipm1)+(xnew(inew)-xold(ipm1))*s
      inew=inew+1
      go to 2000
 2200 ynew(inew)=yold(nold)+sex*(xnew(inew)-xold(nold))
      if(ynew(inew).lt. 0.0) ynew(inew)=0.0
      inew=inew+1
      go to 2000
      end subroutine newgrid
      
      
! jakub urban 110708: changed to integer function
!      real*8 function nfhx(i)
      integer function nfhx(i)
      implicit none 
      integer i !input
!*****************************************************************
      nfhx=i-1
      if(i.le.2)nfhx=i
      return
      end function nfhx
!
!
      subroutine nubplt
      implicit none 
!
!     record debug information
!
!     ONETWO DIVERGENCE
      return
 8100 format(i6)
 8120 format(6e12.5)
 8000 format('****continue****')
 8005 format(5e12.3,i10)
 8010 format(6e12.3)
 8020 format(i6)
      end subroutine nubplt
      

      subroutine pfit(p,x,y,xv,yv,nx,f,dfdx,dfdy)
      implicit none 
      integer nx !input
      real(c_double) :: xv,yv,f,dfdx,dfdy !input
      integer i,j !local
      real(c_double) :: a1,a2,a3,a4,b1,b2,b3,b4 !local
      real(c_double) :: d2fdx2,d2fdy2 !local
      real(c_double) :: p(nx,*),x(*),y(*) !input
      real(c_double) :: cx(4),cy(4),cxp(4),cyp(4)
      real(c_double) :: c2xp(4),c2yp(4)
!
!-----------------------------------------------------------------------
!
!     bi-quadratic interpolator
!
!     input
!     1.  p     - effectivly a 4x4 matrix of function values
!     but for generality the first dimension of
!     p is nx.
!     2.  x     - associates a grid to the first dimension of p
!     3.  y     - associates a grid to the second dimension of p
!     4.  xv    - location at which bi-quadratic function is evaluated
!     5.  yv    - location at which bi-quadratic function is evaluated
!     6.  nx    - first dimension of p in calling program
!
!     output
!     1.  f     - interpolated value of p at (xv,yv)
!     2.  dfdx  - interpolated value of x-partial of p at (xv,yv)
!     3.  dfdy  - interpolated value of y-partial of p at (xv,yv)
!
!-----------------------------------------------------------------------
!
!
      a1=(x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))
      a2=(x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))
      a3=(x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))
      a4=(x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))
!
      cx(1)=(xv-x(2))*(xv-x(3))*(xv-x(4))/a1
      cx(2)=(xv-x(1))*(xv-x(3))*(xv-x(4))/a2
      cx(3)=(xv-x(1))*(xv-x(2))*(xv-x(4))/a3
      cx(4)=(xv-x(1))*(xv-x(2))*(xv-x(3))/a4
!
      cxp(1)=((xv-x(3))*(xv-x(4)) &
        +      (xv-x(2))*(xv-x(3)) &
        +      (xv-x(2))*(xv-x(4)))/a1
      cxp(2)=((xv-x(3))*(xv-x(4)) &
        +      (xv-x(1))*(xv-x(3)) &
        +      (xv-x(1))*(xv-x(4)))/a2
      cxp(3)=((xv-x(1))*(xv-x(2)) &
        +      (xv-x(1))*(xv-x(4)) &
        +      (xv-x(2))*(xv-x(4)))/a3
      cxp(4)=((xv-x(1))*(xv-x(2)) &
        +      (xv-x(1))*(xv-x(3)) &
        +      (xv-x(2))*(xv-x(3)))/a4
!
      c2xp(1)=2.0*(3.0*xv-x(2)-x(3)-x(4))/a1
      c2xp(2)=2.0*(3.0*xv-x(1)-x(3)-x(4))/a2
      c2xp(3)=2.0*(3.0*xv-x(1)-x(2)-x(4))/a3
      c2xp(4)=2.0*(3.0*xv-x(1)-x(2)-x(3))/a4
!
!
      b1=(y(1)-y(2))*(y(1)-y(3))*(y(1)-y(4))
      b2=(y(2)-y(1))*(y(2)-y(3))*(y(2)-y(4))
      b3=(y(3)-y(1))*(y(3)-y(2))*(y(3)-y(4))
      b4=(y(4)-y(1))*(y(4)-y(2))*(y(4)-y(3))
!
      cy(1)=(yv-y(2))*(yv-y(3))*(yv-y(4))/b1
      cy(2)=(yv-y(1))*(yv-y(3))*(yv-y(4))/b2
      cy(3)=(yv-y(1))*(yv-y(2))*(yv-y(4))/b3
      cy(4)=(yv-y(1))*(yv-y(2))*(yv-y(3))/b4
!
      cyp(1)=((yv-y(3))*(yv-y(4)) &
        +      (yv-y(2))*(yv-y(3)) &
        +      (yv-y(2))*(yv-y(4)))/b1
      cyp(2)=((yv-y(3))*(yv-y(4)) &
        +      (yv-y(1))*(yv-y(3)) &
        +      (yv-y(1))*(yv-y(4)))/b2
      cyp(3)=((yv-y(1))*(yv-y(2)) &
        +      (yv-y(1))*(yv-y(4)) &
        +      (yv-y(2))*(yv-y(4)))/b3
      cyp(4)=((yv-y(1))*(yv-y(2)) &
        +      (yv-y(1))*(yv-y(3)) &
        +      (yv-y(2))*(yv-y(3)))/b4
      c2yp(1)=2.0*(3.0*yv-y(2)-y(3)-y(4))/b1
      c2yp(2)=2.0*(3.0*yv-y(1)-y(3)-y(4))/b2
      c2yp(3)=2.0*(3.0*yv-y(1)-y(2)-y(4))/b3
      c2yp(4)=2.0*(3.0*yv-y(1)-y(2)-y(3))/b4
!
!
!
      f=0.0
      dfdx=0.0
      dfdy=0.0
      d2fdx2=0.0
      d2fdy2=0.0
      do 10 i=1,4
        do 11 j=1,4
          f=f+cx(i)*cy(j)*p(i,j)
          dfdx=dfdx+cxp(i)*cy(j)*p(i,j)
          dfdy=dfdy+cx(i)*cyp(j)*p(i,j)
          d2fdx2=d2fdx2+c2xp(i)*cy(j)*p(i,j)
          d2fdy2=d2fdy2+c2yp(j)*cx(i)*p(i,j)
 11     continue
 10   continue
      if (nx.gt.0) return
      dfdx=dfdx/d2fdx2
      dfdy=dfdy/d2fdy2
!
!
      return
      end subroutine pfit
!
!
      subroutine polfit(ideg,npnts,x,y,coeff,ier)
      use frsubs_mod, only : leqt1fl
      implicit none 
      integer ideg,npnts !input
      integer ier !output
      integer i,j,k,jstart,jend,idegp1,nsums !local
      real(c_double) :: s !local
      real(c_double) :: x(101),y(101),coeff(7),sums(203)
      real(c_double) :: a(49),b(7),wkarea(7)
      ier=0
      idegp1=ideg+1
      if(ideg.ge.(npnts+1)) ier=1
      if(npnts.gt.101) ier=1
      if(ier.eq.1) return
!---------------------------------------------------------------------
      nsums=2*ideg
      do 2400 i=0,nsums
        s=0.0
        do 2200 j=1,npnts
          s=s+x(j)**i
 2200   continue
        sums(i+1)=s
 2400 continue
      k=0
      do 2800 i=1,idegp1
        jstart=i
        jend=i+idegp1-1
        do 2600 j=jstart,jend
          k=k+1
          a(k)=sums(j)
 2600   continue
 2800 continue
      do 3200 i=0,ideg
        s=0.0
        do 3000 j=1,npnts
          s=s+y(j)*(x(j)**i)
 3000   continue
        b(i+1)=s
 3200 continue
!---------------------------------------------------------------
!     solve set of simultaneous equations
!---------------------------------------------------------------
      call leqt1fl(a,1,idegp1,idegp1,b,0,wkarea,ier)
      do 3400 i=1,idegp1
        coeff(i)=b(i)
 3400 continue
      return
      end subroutine polfit

      real(c_double) function polval(xval,coeff,ideg)
      implicit none 
      integer ideg !input
      real(c_double) :: xval !input
      real(c_double) :: yval !local
      integer i !local
      real(c_double) :: coeff(*)
      yval=coeff(ideg+1)
      do 2000 i=ideg+1,2,-1
        yval=yval*xval+coeff(i-1)
 2000 continue
      polval=yval
      return
      end function polval
!
!
      subroutine postnub
      implicit none 
!     ONETWO DIVERGENCE
      return
      end subroutine postnub
!
!
      subroutine prenub
      implicit none 
!     this subroutine calculates certain flux surface information needed
!     by freya
!     ONETWO DIVERGENCE
      return
      end subroutine prenub
!
!
      real(c_double) function ranorm(iseed)
      implicit none 
      integer iseed !input
      save
!     generate one normal (0,1) pseudo random number using a method
!     based on central limit theorem and power residue method
!     output becomes truly normal as k (below) goes to infinity.
!     general formula for the deviate is
!     y=(sum of x(i),i=1,k) -k/2. )/ sqrt(k/12.)
!     where x(i) are are uniformally distributed on 0,1.
!     method is borrowed thru ibm. they use k=12.
!
!     drand(iflag) is [0.,1.] random number; iflag is real*4.
!     iflag=0 gives next number in sequence.  See Absoft f77.
!     For Cray f90: Use ranf
!     Portable: (from GA portlib.f, 000114, bh) Use RANDOM_my(seed0)
      integer i,k !local
      real(c_double) :: rtkd12,seed0,a,y !local
      parameter (k=12,rtkd12=1.0000000)

!      integer iflag
!      iflag=0
      data       seed0 /0.0d0/  !I.E. Continue with intialized random_my

!     rtkd12=sqrt(k/12.)
!
      a=0.
      do 10 i=1,k
!     y=ranf(iseed)
!990131        y=ranf()
!000114*Not portable*         y=drand(iflag)
      y=RANDOM_my(seed0)
 10   a=a+y
      ranorm=(a-0.5*k)/rtkd12
      return
      end function ranorm
      
      
      real(c_double) function reldif(x1,x2)
      implicit none 
      real(c_double) :: x1,x2 !input
      real(c_double) del,xabs !local
      del=abs(x1-x2)
      xabs=abs(x1)
      if(xabs.gt.1.0e-20) del=del/xabs
      reldif=del
      return
      end function reldif
!
!
      subroutine revers(x,work,n)
      implicit none 
      integer n !input
      integer i,k !local
      real(c_double) :: work(n),x(n)
!
      do 2100 i=1,n
        work(i)=x(i)
 2100 continue
      k=0
      do 2200 i=n,1,-1
        k=k+1
        x(k)=work(i)
 2200 continue
      return
      end subroutine revers
!
!
      subroutine rotate(naptr, iatype, aheigh, awidth, alen, bhofset, &
        bvofset, cangv, cangh, ib, isourc, &
        costp, sintp, costpp, sintpp, blenp, &
        nsourc, sangv, sangh, rpivot, zpivot, &
        mlost, x0, y0, z0, vx0, vy0, vz0)
      implicit none 
      integer naptr,ib,isourc,nsourc,mlost !input
      real(c_double) :: x0,y0,z0,vx0,vy0,vz0 !input
      integer iap,i !local
      real(c_double) :: zero,alen1,alen2,tpvt,temp,zcos,zsin !local
!cc
!     this subroutine advances a particle from source to the pivot point,
!     and transforms coordinates.
!
!     translation(s):  Particle to each aperature aligned along the souce
!     axis.  Return with mlost=1 if particle hits aperature.
!     rotations:       about y-axis so z-axis is aligned vertically;
!     about new z-axis so x-axis parallel to beamline axis.
!     translation:     coordinates to intersection of beamline axis with
!     source exit plane.
!     translation:     particles to each of beamline axis centered aperatures,
!     checking for mlost=1 condition.
!     translation:     particle and then coordinate axis to pivot point.
!     rotation:        x-axis through pivot point and torus center.
!     translation:     origin to torus center.
!
!     ONETWO DIVERGENCE
      integer nap !local
      parameter(nap=10)
      real(c_double) :: aheigh(nap,*), awidth(nap,*), &
        alen(nap,*), bhofset(*), bvofset(*), blenp(*), cangv(*), &
        cangh(*), sangv(*), sangh(*), rpivot(*), zpivot(*), &
        costp(*), sintp(*), costpp(*),  sintpp(*)
      integer iatype(nap,*)
!cc
      zero=0.d0

!BH110314:  Define iap for nsourc=1 case
      if(nsourc.eq.1)  then
         iap=0
         go to 19
      endif

!     Move particle to each of source-axis centered aperatures, and
!     test for particle passage.
      alen1=0.0
      alen2=0.0
!
      do 10  i=1,naptr
!     iatype.le.4 are source centered aperatures.
        if(iatype(i,ib).gt.4)  go to 11
        alen1=alen(i,ib)-alen2
        alen2=alen(i,ib)
        tpvt=-alen1/vx0
        x0=x0+vx0*tpvt
        y0=y0+vy0*tpvt
        z0=z0+vz0*tpvt
!
        mlost=0
        go to (12,13,14,15),  iatype(i,ib)
 12     if((y0**2+z0**2).le.(0.5*awidth(i,ib))**2) go to 10
        go to 16
 13     if(abs(y0).le.0.5*awidth(i,ib).and. &
          abs(z0).le.0.5*aheigh(i,ib))  go to 10
        go to 16
 14     if(abs(z0).le.0.5*aheigh(i,ib))  go to 10
        go to 16
 15     if(abs(y0).le.0.5*awidth(i,ib))  go to 10
 16     mlost=1
        return
 10   iap=i
 11   continue
!
!     Source center is taken to lie at vertical distance
!     bvofset and horiz. distance bhofset from beam line axis.
!     Source centerline intesects beam line axis at distance
!     bleni from source.  Coords. are changed from source
!     coordinates to those with x-axis along beamline axis
!     directed in +R direction, z in the vertical direction, and
!     origin in source exit plane.
!
!     Rotate about y-axis of source coordinate system so x-axis
!     is vertical and in exit plane defined by beam axis
!     and source centers.
!      write(*,*) 'isourc = ',isourc
!      write(*,*) 'costp,sintp = ',costp(ib),sintp(ib)
!      write(*,*) 'bhofset,bvofset = ',bhofset(ib),bvofset(ib)
!      write(*,*) 'cangh,cangv = ', cangh(ib),cangv(ib)
!      write(*,*) 'sangh,sangv = ', sangh(ib),sangv(ib)
      temp=x0*costp(ib)-isourc*z0*sintp(ib)
      z0=isourc*x0*sintp(ib)+z0*costp(ib)
      x0=temp
      temp=vx0*costp(ib)-isourc*vz0*sintp(ib)
      vz0=isourc*vx0*sintp(ib)+vz0*costp(ib)
      vx0=temp
!
!     Rotate about z-axis of source system so y-axis lies
!     in source exit plane.
      temp=x0*costpp(ib)-isourc*y0*sintpp(ib)
      y0=isourc*x0*sintpp(ib)+y0*costpp(ib)
      x0=temp
      temp=vx0*costpp(ib)-isourc*vy0*sintpp(ib)
      vy0=isourc*vx0*sintpp(ib)+vy0*costpp(ib)
      vx0=temp
!
!     Translate coordinate axes to beamline axis.
      x0=x0
      y0=y0+isourc*bhofset(ib)
      z0=z0+isourc*bvofset(ib)
!
!
!     Translate particle to beamline centered aperatures and set
!     mlost=1 if particle lost.
 19   continue
      if(iap.eq.naptr)  go to 29
      alen1=0.0
      alen2=abs(x0)
      do 20  i=iap+1,naptr
        alen1=alen(i,ib)-alen2
        if(alen1.lt.0.)  stop 'zfreya: rotate'
        alen2=alen(i,ib)
        tpvt=-alen1/vx0
        x0=x0+vx0*tpvt
        y0=y0+vy0*tpvt
        z0=z0+vz0*tpvt
        mlost=0
        go to (22,23,24,25,26),  iatype(i,ib)-4
 22     if((y0**2+z0**2).le.(0.5*awidth(i,ib))**2) go to 20
        go to 27
 23     if(abs(y0).le.0.5*awidth(i,ib).and. &
          abs(z0).le.0.5*aheigh(i,ib))  go to 20
        go to 27
 24     if(abs(z0).le.0.5*aheigh(i,ib))  go to 20
        go to 27
 25     if(abs(y0).le.0.5*awidth(i,ib))  go to 20
        go to 27
!     DIII-D  special case:
 26     if(abs(z0).gt.24.75)  go to 27
        if(abs(y0).le.13.45)  go to 20
        if(abs(y0).gt.21.7)  go to 27
        if(abs(z0).gt.24.75*((abs(y0)-13.45)/(-8.25)+1))  go to 27
        go to 20
 27     mlost=1
        return
 20   continue
 29   continue
!
!     translate particle and then coordinate axes to pivot point
!
      tpvt=-(blenp(ib)+x0)/vx0
      x0=0.0
      y0=y0+vy0*tpvt
      z0=z0+vz0*tpvt
!cc
      if(sangv(ib).eq.zero) go to 30
      zcos=cangv(ib)
      zsin=sangv(ib)
      temp=x0*zcos+z0*zsin
      z0=-x0*zsin+z0*zcos
      x0=temp
      temp=vx0*zcos+vz0*zsin
      vz0=-vx0*zsin+vz0*zcos
      vx0=temp
 30   continue
!cc
      if(sangh(ib).eq.0.) go to 40
      zcos=cangh(ib)
      zsin=sangh(ib)
      temp=x0*zcos+y0*zsin
      y0=-x0*zsin+y0*zcos
      x0=temp
      temp=vx0*zcos+vy0*zsin
      vy0=-vx0*zsin+vy0*zcos
      vx0=temp
 40   continue
!cc
      x0 = x0 + rpivot(ib)
      z0 = z0 + zpivot(ib)
!cc
      return
      end subroutine rotate
!
!
      subroutine scalit(xa,xb,xmin,xmax,imin,imax,ifail,maxit,toler &
        ,pertf,trans,iprint &
        ,c1,c2,c3,drtabi,mfm1,psival,ptab,rmajor)
      implicit none 
      integer imin,imax,ifail,maxit,iprint,mfm1 !input/output
      real(c_double) :: xa,xb,xmin,xmax,c1,c2,c3,toler,pertf,trans !input
      real(c_double) :: drtabi,psival,rmajor !input
      integer icall,ick,iter !local
      real(c_double) :: test,dxa,xbmag,dxb,slope,fc,dxmag,deltax
      save
!
!
!*****************************************************************
!
!**   subroutine for scalar iteration.
!     scalit is a non-subscripted simplification of vectit, a more
!     general code for iteration with vectors.
!
!**   method - an extension of the secant method.  the procedure
!     resembles newton's method.  after each iteration a linear
!     extrapolation to the root is made.  the slope used for
!     extrapolation is computed by two methods.  the first
!     method is to perturb x a small amount, giving a good
!     approximation to the tangent.  this is accurate, but it
!     does not lead to productive iterations.  the second
!     method is to do one more ordinary iteration and calculate
!     a secant.  this is more productive and sufficiently
!     accurate provided one is already quite close to the
!     solution.  a transition from method 1 to method 2 is made
!     when the test quantity (dx/x) is less than a prescribed
!     value = 'trans'.  the method is quadratically convergent
!     and works even when the uncorrected iteration is
!     divergent.  as with newton's method it is not guaranteed
!     to work unless the first guess is 'sufficiently close' to
!     the desired root.  a false guess correction is provided
!     to force the root (if one exists) to lie within a given
!     interval.
!
!**   summary of arguments
!
!     name      definition
!
!     xa      - first guess - destroyed during calculation.
!     xb      - final guess (converged value if ifail > 0).
!     xmin    - allowable lower limit to x if imin.ne.0
!     xmax    - allowable upper limit to x if imax.ne.0
!     imin    - if = 0, ignore xmin criterion.
!     imax    - if = 0, ignore xmax criterion.
!     ifail   - results indicator -
!     > 0  - converged within ifail cycles,
!     = 0  - abnormal exit,
!     = -1 - failed to converge in maxit cycles.
!     maxit   - upper limit to number of iteration cycles.
!     note - 2 local iterations per cycle.
!     maxit = 10 should be sufficient.
!     toler   - convergence criterion for ratio of  abs(dx/x).
!     pertf   - perturbation factor.  deltax = pertf*abs(x).
!     see above.  typically  0.1
!     trans   - transition value of (dx/x) for solution method.
!     typically .01
!     iprint  - if ne.0, debug print of all x iterations
!     is provided in the printout.
!**   subroutines required by scalit
!
!     fortran name        description
!
!     chek      -  checks to see if a given x is between xmax and xmin.
!
!     itorb     -  a user-supplied external subroutine that defines
!     the basic uncorrected iteration. given x, calculates
!
!**   use of the program scalit
!
!     user supplies calling program.  calling program initializes xa.
!     calling program then calls scalit.  entire iteration up to point
!     of convergence is done in scalit.  scalit returns converged value
!     of x in xb.
!     basic uncorrected iteration must also be defined by user in an
!     external subroutine named 'itorb(x,dx)'.  itorb has arguments
!     x  - the input value of a guess x (must be saved)
!     dx - the output value of a change in x
!
!     p. d. smith   5/10/74
!
!***********************************************************************
!
      real(c_double) :: ptab(*)
      data icall/0/
      icall = icall + 1
      iter = 0
      test = 1.0e6
      if(iprint.ne.0) write(6,800) icall
      if(iprint.ne.0) write(6,801) iter,xa
 800  format(/' debug print for call',i4,' to scalit'/ &
        ,' x value at end of each iteration printed below',/)
 801  format(' cycle',i3,3x,e15.8)
!
!*    ** now begin iteration/extrapolation cycles.
!
      ifail = 0
 200  continue
!
!*    ** start of an iteration cycle
      iter = iter + 1
!
!*    ** check to see if maximum iteration criterion is exceeded
      if(iter.le.maxit) go to 210
      ifail = -1
      return
 210  continue
!
!*    ** o k to do another cycle.
!
      if(test.lt.trans) go to 150
!
!*    ** use perturbation method to estimate slope.
!*    ** whenever the iteration is still not close to the solution.
!
!*    ** first iteration.
      call itorb(c1,c2,c3,drtabi,mfm1,psival,ptab,rmajor,xa,dxa)
      xb = xa + dxa
      xbmag = abs(xb)
      deltax = pertf*xbmag
!
!*    ** perturb xa
      xa = xa + deltax
!
!*    ** iterate using perturbed xa.
      call itorb(c1,c2,c3,drtabi,mfm1,psival,ptab,rmajor,xa,dxb)
!
!*    ** change xa back to initial value.
      xa = xa - deltax
!
!*    ** calculate slope.
      slope = (dxb - dxa)/deltax
!
!*    ** negative inverse of slope is correction factor.
      fc = -1./slope
 150  continue
!
!*    ** now begin extrapolation/iteration cycle, assuming that the
!*    ** difference method will work.  if it doesn't work, the solution
!*    ** will be cycled back to the perturbation method.
!
!
!*    ** first extrapolate previous xa and dxa to get a value for xb.
      dxb = fc*dxa
      xb = xa + dxb
      if(iprint.ne.0) write(6,801) iter,xb
!
!*    ** check limits of xb
      call chek(xb,xmin,xmax,imin,imax,ick)
      if(ick.eq.0) go to 220
!
!*    ** overshot boundary. increase effectiveness of damping.
      fc = .5*fc
!
!*    ** start iteration over again.
      go to 200
 220  continue
!
!*    ** extrapolation was successful.  continue
!
!*    ** check convergence and return if test is passed.
      xbmag = abs(xb)
      dxmag = abs(dxb)
      if(dxmag.le.0.)    go to 240
      if(xbmag.le.0.) go to 250
      test = dxmag/xbmag
      if(test.lt.trans) go to 230
!
!*    ** solution is not sufficiently close.  must use perturbation.
!*    ** update extrapolation prior to next perturbation.
      xa = xb
 225  continue
      go to 200
!
 230  continue
      if(test.gt.toler)  go to 250
 240  continue
!
!*    ** converged
      ifail = iter
      return
 250  continue
!
!*    ** no convergence - continue with calculations.
!
!*    ** perform an iteration with xb to get dxb.
      call itorb(c1,c2,c3,drtabi,mfm1,psival,ptab,rmajor,xb,dxb)
!
!*    ** calculate new secant slope
      slope = (dxb - dxa)/(xb - xa)
!
!*    ** update xa and dxa prior to next iteration/extrapolation.
      xa = xb
      dxa = dxb
!
!*    ** calculate another correction factor.
      fc = -1./slope
!
!*    ** go back and do another iteration cycle
      go to 200
      end subroutine scalit
!
!
      real(c_double) function scxif(i,er)
      implicit none 
      integer i !input
      real(c_double) :: er !input
      real(c_double) :: aer,siif1,ansq,xtilda !local
!*****************************************************************
!--   estimated cross sections for chg. exchange .
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
!
      !real(c_double) :: polyf ! Now in this module
!
!
      if(i.ge.2)goto 200
!
!--   1s:
!990131      aer=alog(er)
      aer=log(er)
!
!     charge exchange.
!     below 1.e6 ev, use riviere's fit       above 1.e6 ev, use
!     power-law extrapolation.
!
      if(er .gt. 1.e6)goto 105
      siif1=(.6937e-14)*(1.-.06732*aer)**2/(1.+.1112e-14*er**3.3)
      goto 110
 105  siif1=(4.8363e-22)/((1.e-6)*er)**3.3
!
 110  scxif=siif1
      return
!
 200  continue
      ansq=(nfhx(i))**2
      xtilda = ansq*er/9260.
      scxif = 0.
      if(xtilda.lt.1.) scxif=1.465e-11*xtilda*(ansq/er)
      return
      end function scxif
!
!
      real(c_double) function scxzf(i,z,er)
      implicit none 
      integer i !input
      real(c_double) :: z,er !input
      integer ni,id !local
      real(c_double) :: ansq,xtilda
!*****************************************************************
!--   cross section for electron loss via collision with impurity.
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
!
      ni=nfhx(i)
      id=2
      if(ni.eq.1)id=1
!
      ansq=(DBLE(ni))**2
      xtilda = ansq*er/(31250.*z)
      scxzf = 0.
      if(xtilda.lt.1.) scxzf=1.465e-11*xtilda*(ansq*z*z/er)
      return
      end function scxzf
!
!
      real(c_double) function sdaccf(ni,nj,z,er)
      implicit none 
      integer ni,nj !input
      real(c_double) :: z,er !input
      real(c_double) :: omega,alam,beta !local
!*****************************************************************
!--   dacc cross section.
!     janev & presnyakov, j. phys. b 13, 4233 (1980).
!
      integer ms,mc !local
      parameter (ms=21,mc=35)
      real(c_double) :: f,ar !in /b3/
      common/b3/f(ms,mc),ar(ms+1,ms+1)
!
      omega=.5*(1./((DBLE(ni))**2)-1./((DBLE(nj))**2))
      alam=sqrt(f(ni,nj)/(2.*omega))
      beta=z*alam*omega/(er/24982.)
      sdaccf=(1.7595e-16)*(z*alam/omega)*dfhx(beta)
      return
      end function sdaccf
!
!
      real(c_double) function seef(i,j,er)
      implicit none 
      integer i,j !input
      real(c_double) :: er !input
      integer ni,nj,id !local
      real(c_double) :: zero,ryd,emin12,emax12,emin13,emax13,aer !local
      save
!*****************************************************************
!--   cross sections for excitation from i to j due to electron impact.
      integer ms,mc !local
      parameter (ms=21,mc=35)
      
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      
      integer nouthx,istart,ihxbug !in /b2/
      common/b2/nouthx,istart,ihxbug
      
      real(c_double) :: en,dg,ae,be,de1,de2,ge1,ge2 !in /b4/
      common/b4/en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc), &
        de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
!
      !real(c_double) :: polyf ! Now in this module
      real(c_double) :: ae12(8)
      real(c_double) :: ae13(9)

      data ryd/13.6/
!
!--   data for 1s to 2s (e impact):
      data emin12/10.2/, emax12/1140.5/
      data ae12/ &
        -1.5138513657d-17, 1.8231019028d-16,-1.9274654811d-16, &
        8.9291530932d-17,-2.2108041618d-17, 3.0448025264d-18, &
        -2.2039298214d-19, 6.5491238735d-21/
!@    real ae12
!@    dimension ae12(8)
!@    data ae12/
!@    *-1.5138513657e-17, 1.8231019028e-16,-1.9274654811e-16,
!@    * 8.9291530932e-17,-2.2108041618e-17, 3.0448025264e-18,
!@    *-2.2039298214e-19, 6.5491238735e-21/
!
!--   data for 1s to 2p (e impact)
      data emin13/10.2/, emax13/747.4/
      data ae13/ &
        -2.1078372920d-14, 3.7548968459d-14,-2.8788463637d-14, &
        1.2394689683d-14,-3.2696005067d-15, 5.4068250737d-16, &
        -5.4759059881d-17, 3.1084887079d-18,-7.5815695055d-20/
!@    real ae13
!@    dimension ae13(9)
!@    data ae13/
!@    *-2.1078372920e-14, 3.7548968459e-14,-2.8788463637e-14,
!@    * 1.2394689683e-14,-3.2696005067e-15, 5.4068250737e-16,
!@    *-5.4759059881e-17, 3.1084887079e-18,-7.5815695055e-20/
!
      if(i.ge.j)goto 300
      seef=0.
!990131      aer=alog(er)
      aer=log(er)
      if((i.ge.3) .or. (j.ge.4))goto 250
      goto(210,240),i
!
 210  continue
      goto(300,220,230),j
!
!--   1s to 2s:
 220  continue
      if(er .le. emin12)return
      if(er .le. 11.1)  goto 221
      if(er .ge. emax12)goto 222
      seef=polyf(ae12,8,aer)
      goto 223
!--   linear interpolation for region just above threshold:
 221  continue
      seef=((er-10.2)/.9)*1.6297e-17
      goto 223
!--   asymptotic:
 222  continue
      seef=(5.312e-16)/er
 223  continue
      return
!
!--   1s to 2p:
 230  continue
      if(er .le. emin13)return
      if(er .ge. emax13)goto 232
      seef=polyf(ae13,9,aer)
      goto 233
!--   asymptotic:
 232  continue
      seef=(2.65571e-15/er)*(aer-2.120681)
 233  continue
      return
!
!--   2s to 2p:
 240  continue
!990131      seef=(8.617e-15/er)*alog(1.14e4*er)
      seef=(8.617e-15/er)*log(1.14e4*er)
      return
!
 250  continue
      ni=nfhx(i)
      nj=j-1
!--   from vriens & smeets (eq.(14)):
!990131   seef=(1.75947e-16*ryd)*(ae(ni,nj)*alog(.5*er/ryd+de1(ni,nj))
      seef=(1.75947e-16*ryd)*(ae(ni,nj)*log(.5*er/ryd+de1(ni,nj)) &
        +be(ni,nj))/(er+ge1(ni,nj))
!990131      seef=amax1(seef,0.)
      zero=0.d0
      seef=max(seef,zero)
      id=2
      if(ni.eq.1)id=1
      return
!
 300  continue
      ihxbug = 6
      if(nouthx.gt.0) write(nouthx,3939) i,j
 3939 format(' ??? seef error: i4,j are in error       i4=',i3, &
        /,'   j = ',i3)
      return
      end function seef
!
!
      subroutine setrz(ndim,rmin,rmax,zmin,zmax,dr,dz,nr,nz)
      implicit none 
      integer nr,nz !output
      real(c_double) :: rmin,rmax,zmin,zmax,dr,dz !input/output
!-----------------------------------------------------------------------
!     this subroutine sets up (r,z) grid parameters for freya when the
!     optional two-dimensional deposition calculation is used.  default
!     parameters correspond to original eqdisk grid
!-----------------------------------------------------------------------
      integer ndim(2)
      nr=ndim(1)
      nz=ndim(2)
      dr=(rmax-rmin)/(nr-1)
      dz=(zmax-zmin)/(nz-1)
      return
      end subroutine setrz
!
!
      real(c_double) function sezf(i,j,z,er)
      implicit none 
      integer i,j !input
      real(c_double) :: z,er !input
      real(c_double) :: emin12,emax12,emin13,emax13,zero,one,aer,sezfsv !local
      integer ni,n12,nj,id
      save
!*****************************************************************
!--   cross sections for excitation from i to j due to ion impact.
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl

      common/b2/nouthx,istart,ihxbug
      integer nouthx,istart,ihxbug !in /b2/
!
      !real(c_double) :: polyf ! Now in this module
      real(c_double) :: ai12(11)
      real(c_double) :: ai13(12)
!
!--   data for 1s to 2s (p impact):
      data emin12/1.e3/, emax12/1.e6/
      data ai12/ &
        -3.3812171941d+05, 3.5976726159d+05,-1.7037301759d+05, &
        4.7282644608d+04,-8.5167211591d+03, 1.0405267504d+03, &
        -8.7343818297d+01, 4.9754307346d+00,-1.8412279745d-01, &
        3.9984716057d-03,-3.8707719188d-05/
!
!--   data for 1s to 2p (p impact):
      data emin13/1.e3/, emax13/1.e6/
      data ai13/ &
        -1.3490069287d+06, 1.4573274888d+06,-7.0815407013d+05, &
        2.0432442357d+05,-3.8900004212d+04, 5.1319758650d+03, &
        -4.7884757131d+02, 3.1608287540d+01,-1.4469255104d+00, &
        4.3759824250d-02,-7.8716881911d-04, 6.3824434435d-06/
!
      if(i.ge.j) then
        ihxbug = 4
        if(nouthx.gt.0) write(nouthx,3939) i,j
 3939   format(' ??? sezf input error: i4=',i3,'        j=',i3)
        return
      endif

      if(z .gt. 1.01)goto 300
!990131      aer=alog(er)
      aer=log(er)
      if((i.ge.3) .or. (j.ge.4))goto 250
      goto(210,240),i
!
 210  continue
      goto(400,220,230),j
!
!--   1s to 2s:
 220  continue
      if(er .ge. emin12)goto 221
      sezf=2.8875e-18*(er/emin12)**.7093
      goto 223
 221  continue
      if(er .ge. emax12)goto 222
      sezf=exp(polyf(ai12,11,aer))
!@    sezf=exp(polyf(ai12,11,aer))
      goto 223
 222  continue
      sezf=1.9564e-12/er
 223  continue
      sezfsv = sezf
      return
!
!--   1s to 2p:
 230  continue
      if(er .ge. emin13)goto 231
      sezf=1.4053e-17*(er/emin13)**.95695
      goto 233
 231  continue
      if(er .ge. emax13)goto 232
!990131      sezf=exp(polyf(ai13,12,aer))
!@    sezf=exp(polyf(ai13,12,aer))
!BH120202      twelve=12
!BH120202      sezf=exp(polyf(ai13,twelve,aer))
      n12=12
      sezf=exp(polyf(ai13,n12,aer))
      goto 233
 232  continue
!990131      sezf=(6.7085e-13/er)*alog(5701.79*er)
      sezf=(6.7085e-13/er)*log(5701.79*er)
 233  continue
      sezfsv = sezf
      return
!
!--   2s to 2p:
!990131 240  sezf=(1.584e-10/er)*alog(.62*er)
 240  sezf=(1.584e-10/er)*log(.62*er)
!990131      sezf=amax1(sezf,0.)
      zero=0.d0
      sezf=max(sezf,zero)
      sezfsv = sezf
      return
!
 250  continue
      ni=nfhx(i)
      id=2
      if(ni.eq.1)id=1
      nj=j-1
      sezf=sdaccf(ni,nj,one,er)
      sezfsv = sezf
      return
!
!--   impurity scattering:
 300  continue
      if((i.eq.1).and.(j.eq.2))goto 312
      if((i.eq.1).and.(j.eq.3))goto 313
      sezf=0.
      sezfsv = sezf
      if((i.eq.2).and.(j.eq.3))return
      ni=nfhx(i)
      nj=j-1
      id=2
      if(ni.eq.1)id=1
      sezf=sdaccf(ni,nj,z,er)
      sezfsv = sezf
      return
 312  continue
      sezf=.25*sdaccf(1,2,z,er)
      sezfsv = sezf
      return
 313  continue
      sezf=.75*sdaccf(1,2,z,er)
      sezfsv = sezf
      return
!
 400  ihxbug = 5
      if(nouthx.gt.0) write(nouthx,3941)
!     ONETWO DIVERGENCE
 3941 format(" ??? sezf error: state j is 0")
      return
      end function sezf
!
!
      real(c_double) function sief(i,er)
      implicit none 
      integer i !input
      real(c_double) :: er !input
      real(c_double) :: ryd,emin1,emax1,emin2,emax2,aer,ansq,ery !local
      real(c_double) :: sbea !YuP: not sure if every compiler 
                             !will understand that this is a statement function
      save
!*****************************************************************
!--   cross sections for ionization due to e impact.
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      
      !real(c_double) :: polyf ! Now in this module
      real(c_double) :: ae1(7)
      real(c_double) :: ae2(9)
!
      data ryd/13.6/
!
!--   data for e-impact ionization of 1s:
      data emin1/13.723/, emax1/620./
      data ae1/ &
        2.3036652148d-15,-3.4298197580d-15, 1.9465132340d-15, &
        -5.4457519508d-16, 8.0929651995d-17,-6.1527279210d-18, &
        1.8889193736d-19/
!@    data ae1/
!@    * 2.3036652148e-15,-3.4298197580e-15, 1.9465132340e-15,
!@    *-5.4457519508e-16, 8.0929651995e-17,-6.1527279210e-18,
!@    * 1.8889193736e-19/
!
!--   data for e-impact ionization of 2s:
      data emin2/3.4/, emax2/386./
      data ae2/ &
        8.4162370468d-15,-2.3545910382d-14, 2.4728342689d-14, &
        -1.2718639429d-14, 3.6382722533d-15,-6.0044378263d-16, &
        5.5114285405d-17,-2.4253587346d-18, 3.0573876445d-20/
!@    data ae2/
!@    * 8.4162370468e-15,-2.3545910382e-14, 2.4728342689e-14,
!@    *-1.2718639429e-14, 3.6382722533e-15,-6.0044378263e-16,
!@    * 5.5114285405e-17,-2.4253587346e-18, 3.0573876445e-20/
!
!--   be across section (ansq=n**2, ery=e/ryd):
      ! statement function:
      sbea(ansq,ery)= &
        (3.519e-16)*(5.*ansq/3.-1./ery-2./(3.*ansq*ery*ery)) &
        /(ery+3.25/ansq)
!
      sief=0.
      if(i.ge.3)goto 130
!990131      aer=alog(er)
      aer=log(er)
      goto(110,120),i
!
!--   1s:
 110  continue
      if(er .le. emin1)return
      if(er .ge. emax1)goto 112
      sief=polyf(ae1,7,aer)
      goto 113
 112  continue
      sief=(1.3563e-15/er)*(aer+1.823647)
 113  continue
      return
!
!--   2s:
 120  continue
      if(er .le. emin2)return
      if(er .ge. emax2)goto 122
      sief=polyf(ae2,9,aer)
      goto 123
 122  continue
      sief=(8.195137e-15/er)*(aer-.9445204)
 123  continue
      return
!
!--   2p and higher:
 130  continue
      ansq=(DBLE(i-1))**2
      if(er .le. ryd/ansq)return
      sief=sbea(ansq,er/ryd)
      return
      end function sief
!
!
      real(c_double) function siif(i,er)
      implicit none 
      integer i !input
      real(c_double) :: er !input
      real(c_double) :: aer,siif1,siif2,fnum1,fnum2,ansq !local
      save
!*****************************************************************
!--   cross sections for chg. exchange and ion-impact ionization.
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
!
      !real(c_double) :: polyf ! Now in this module
!
!--   data for ion-impact ionization of h(1s):
      real(c_double) :: ai1(6)
      real(c_double) :: afj(7)
      data ai1/ &
        -4.410635d+02, 1.846170d+02,-3.429509d+01, 3.217263d+00, &
        -1.512004d-01, 2.825854d-03/
!@    real ai1
!@    real(c_double) :: ai1(6)
!@    data ai1/
!@    *-4.410635e+02, 1.846170e+02,-3.429509e+01, 3.217263e+00,
!@    *-1.512004e-01, 2.825854e-03/
!
!--   data for ion-impact ionization of h(1s) (freeman & jones):
      data afj/ &
        -.4203309d+2, .3557321d+1,-.1045134d+1, .3139238d+0, &
        -.7454475d-1,.8459113d-2,-.3495444d-3/
!@    data afj/
!@    *-.4203309e+2, .3557321e+1,-.1045134e+1, .3139238e+0,
!@    * -.7454475e-1,.8459113e-2,-.3495444e-3/
!
      if(i.ge.2)goto 200
!
!--   1s:
!990131      aer=alog(er)
      aer=log(er)
!
!     charge exchange.
!     below 1.e6 ev, use riviere's fit       above 1.e6 ev, use
!     power-law extrapolation.
!
      if(er .gt. 1.e6)goto 105
      siif1=(.6937e-14)*(1.-.06732*aer)**2/(1.+.1112e-14*er**3.3)
      goto 110
 105  siif1=(4.8363e-22)/((1.e-6)*er)**3.3
!
!     ion-impact ionization.
!     at low energies or high energies, use riviere's fits
!     at intermediate energies, use fit to rkj's curve.
!
 110  continue
!
!     freeman and jones option:
!     siif2=exp(polyf(afj,7,alog(.001*er)))
!     goto 120
!     111 continue
!
!990131
      fnum1=807.4
      fnum2=154400.
      if(er .lt. fnum1)goto 113
      if(er .gt. fnum2)goto 114
      siif2=exp(polyf(ai1,6,aer))
!@    siif2=exp(polyf(ai1,6,aer))
      goto 120
 113  continue
      siif2=exp(-80.206+8.156*aer-.3784*aer*aer)
      goto 120
 114  siif2=(1.56346e-12/er)*(aer-1.792160)
 120  continue
      siif=siif1+siif2
      return
!
 200  continue
      ansq=(nfhx(i))**2
      siif=(1.465e-11)*(ansq/er)*(1.-exp(-ansq*er/9260.))
      return
      end function siif
!
!
      real(c_double) function sizf(i,z,er)
      implicit none 
      integer i !input
      real(c_double) :: z,er !input
      integer ni,id !local
      real(c_double) :: ansq !local
!*****************************************************************
!--   cross section for electron loss via collision with impurity.
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl !in /b1/
      common/b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
!
      ni=nfhx(i)
      id=2
      if(ni.eq.1)id=1
!
      ansq=(DBLE(ni))**2
      sizf=(1.465e-11)*(ansq*z*z/er)*(1.-exp(-ansq*er/(31250.*z)))
      return
      end function sizf
!
!
      !YuP[2019-06-17] This subroutine has a conflict for definition of emzrat
      !as an array with common block /gecom/ where emzrat is a scalar.
      !But this subroutine is not used anyway, so - commenting it.
!      subroutine slow1(atw,atwf,ene,en,enn,ezero,ibcur,ibcx, &
!        ifirst,kj,nj,nion,te,zsq,zsqf, &
!        bke,bki,ecrit,emzrat,encap,ge,gi,taus,fionx)
!      implicit integer (i-n), real*8 (a-h,o-z)
!      save
!!     ------------------------------------------------------------------
!!     This subroutine evaluates the transfer functions N,Ge,Gi,Ke, and
!!     Ki defined by Callen et al., IAEA Tokyo, Vol. I, 645 (1974) to
!!     describe fast ion slowing down.  Ge and Gi are the fractions of
!!     the fast ion energy transferred to the electrons and thermal ions,
!!     respectively.  The remaining energy fraction is assumed lost due
!!     to secondary charge exchange, i.e., charge exchange between fast
!!     ions and thermal neutrals.  Similarly, Ke and Ki are the fractions
!!     of the fast ion parallel momentum transferred to electrons and
!!     thermal ions.  Moreover, N, Ge, and Ke are proportional to the
!!     slowing-down times for particles, energy, and parallel momentum
!!     and are used in subroutine slow2.  The fast ions may be beam ions
!!     from neutral beam heating or alpha particles from fusion.
!!     --------------------------------------------------------------------
!      real(c_double) :: atw(*), ene(*), en(kj,*), enn(kj,*), te(*), zsq(kj,*)
!      real(c_double) :: bke(*), bki(*), ecrit(*), emzrat(*), encap(*), &
!        ge(*), gi(*), taus(*)
!      data pi, rot2pi / 3.14159, 2.50663 /
!      data xmasse / 9.110d-28 /
!      data xmassp / 1.673d-24 /
!      data charg4 / 5.308d-38 /
!!
!!     calculate constants
!!
!      econst = 0.5*(4.5*pi*xmassp/xmasse)**0.333333
!      tconst = atwf*xmassp/(zsqf*xmasse)
!!
!!     begin loop over mesh points
!!
!      do 100 j=1,nj
!        if(ifirst.eq.0) go to 20
!!
!!     calculate ecrit, the "critical energy" at which fast ions transfer
!!     equal energy to electrons and thermal ions; also calculate
!!     emzrat = mi*<Z>/(mf*[Z])
!!
!        sum1 = 0.
!        sum2 = 0.
!        do 10 k=1,nion
!          sum1 = sum1 + en(j,k)*zsq(j,k)
! 10     sum2 = sum2 + en(j,k)*zsq(j,k)/atw(k)
!        ecrit(j) = econst*te(j)*atwf*(sum2/ene(j))**0.666667
!        emzrat(j) = sum1/(atwf*sum2)
!!
!!     calculate taus, the Spitzer momentum exchange time for electron-ion
!!     collisions
!!
!        teerg = 1.6e-9*te(j)
!!990131        xlam  = 24. - alog(sqrt(ene(j))/(1.e3*te(j)))
!        xlam  = 24. - log(sqrt(ene(j))/(1.e3*te(j)))
!        ztaue = sqrt(xmasse*teerg**3) &
!          /(1.33333*rot2pi*ene(j)*charg4*xlam)
!        taus(j) = tconst*ztaue
!!
!!     calculate vrat = sqrt(ecrit/ezero) and taurat = taus/taucx,
!!     where taucx is the charge exchange time
!!
! 20     vrat = sqrt(ecrit(j)/ezero)
!        taurat = 0.
!        if(ibcx.eq.0) go to 50
!        if(atwf.gt.3.) go to 50
!        taurat = taus(j) * (enn(j,1)+enn(j,2)) * cxr(ezero/atwf)
!!
!!     calculate encap, ge, and gi
!!
! 50     encap(j) = encapf(vrat,taurat)
!        ge(j) = gef(vrat,taurat)
!        if(ge(j).lt.0.0) ge(j) = 0.0
!        gi(j) = 1. - (1.+0.5*taurat)*ge(j)
!        if(gi(j).lt.0.0) gi(j) = 0.0
!!     fionx allows testing of non-classical slowing down
!        if(fionx.eq.0.0) go to 70
!        ft=fionx*gi(j)
!        if(ft.gt.ge(j)) ft=ge(j)
!        gi(j)=gi(j) + ft
!        ge(j)=ge(j) - ft
! 70     continue
!!
!!     calculate bke and bki if ibcur.ne.0
!!
!        if(ibcur.eq.0) go to 100
!        bke(j) = bkef(vrat,taurat,emzrat)
!        bki(j) = 1. - (1.+taurat)*bke(j)
!!
!!     end loop over mesh points
!!
! 100  continue
!      return
!      end subroutine slow1
!
!
      subroutine slow2(bke,dtt,encap,enfsav,ge,ibcur,nj, &
        ppfsav,qfsav,sfsav,spfsav,taus,wfsav, &
        enf,ppf,qf,sf,spf,taupf,tauppf,tauef,wf)
      implicit none 
      real(c_double) :: dtt !input
      integer ibcur,nj !input
      integer j !local
!
!     This subroutine models fast ion slowing down.  The fast ions may
!     be beam ions from neutral beam heating or alpha particles from
!     fusion.
!
      real(c_double) :: bke(*), encap(*), enfsav(*), ge(*), ppfsav(*), qfsav(*), &
        sfsav(*), spfsav(*), taus(*), wfsav(*)
      real(c_double) :: enf(*), ppf(*), qf(*), sf(*), spf(*), taupf(*), &
        tauppf(*), tauef(*), wf(*)
!
!     begin loop over mesh points
!
      do 100 j=1,nj
!
!     calculate slowing down times
!
        taupf(j) = taus(j)*encap(j)
        tauppf(j) = taus(j)*bke(j)
        tauef(j) = 0.5*taus(j)*ge(j)
!
!     calculate fast ion particle density and delayed particle source
!
        enf(j) = enfsav(j)*exp(-dtt/taupf(j)) &
          + sfsav(j)*taupf(j)*(1.-exp(-dtt/taupf(j)))
        sf(j) = enf(j)/taupf(j)
!
!     calculate fast ion energy density and delayed energy source
!
        wf(j) = wfsav(j)*exp(-dtt/tauef(j)) &
          + qfsav(j)*tauef(j)*(1.-exp(-dtt/tauef(j)))
        qf(j) = wf(j)/tauef(j)
!
!     calculate density and delayed source of fast ion parallel momentum
!
        if(ibcur.eq.0) go to 100
        ppf(j) = ppfsav(j)*exp(-dtt/tauppf(j)) &
          + spfsav(j)*tauppf(j)*(1.-exp(-dtt/tauppf(j)))
        spf(j) = ppf(j)/tauppf(j)
!
!     end loop over mesh points
!
 100  continue
      return
      end subroutine slow2
!
!
!BH990903      subroutine smooth(x,y,n,ideg,nuse,nupdate)
      subroutine smooth1(x,y,n,ideg,nuse,nupdate)
      implicit none 
      integer n,ideg,nuse,nupdate !input
      integer i,i1,i2,ier !local
      real(c_double) :: big,small !local
      real(c_double) :: x(n),y(n),coeff(8)
!
      do 1800 i=1,8
        coeff(i)=0.0
 1800 continue
      big=-1.0
      do 2000 i=1,n
        if(y(i).gt.big) big=y(i)
 2000 continue
      small=big*.001
      i1=0
 2100 i1=i1+1
      i2=i1+nuse-1
      if(i2.le.n) go to 2120
!-----------------------------------------------------
!     use last set of coefficients to finish with
!-----------------------------------------------------
      do 2110 i=i1,n
        y(i)=polval(x(i),coeff,ideg)
 2110 continue
      go to 2150
!-----------------------------------------------------
 2120 if(y(i1).le.small .and. y(i1+1).le.small) go to 2100
      call polfit(ideg,nuse,x(i1),y(i1),coeff,ier)
      do 2140 i=i1,i1+nupdate-1
        y(i)=polval(x(i),coeff,ideg)
 2140 continue
      go to 2100
 2150 do 2160 i=1,n
        if(y(i).lt. 0.0) y(i)=0.0
 2160 continue
      return
      end subroutine smooth1
!
!
      subroutine solveq(ndim, n, a, b, ipvt)
      implicit none 
!
      integer ndim, n, ipvt(n)
      real(c_double) :: a(ndim,n),b(n)
!
!     solution of linear system, a*x = b .
!     do not use if decomp has detected singularity.
!
!     input..
!
!     ndim = declared row dimension of array containing a .
!
!     n = order of matrix.
!
!     a = triangularized matrix obtained from decomp .
!
!     b = right hand side vector.
!
!     ipvt = pivot vector obtained from decomp .
!
!     output..
!
!     b = solution vector, x .
!
      integer kb, km1, nm1, kp1, i, k, m
      real(c_double) :: t
!
!     forward elimination
!
      if (n .eq. 1) go to 50
      nm1 = n-1
      do 20 k = 1, nm1
        kp1 = k+1
        m = ipvt(k)
        t = b(m)
        b(m) = b(k)
        b(k) = t
        do 10 i = kp1, n
          b(i) = b(i) + a(i,k)*t
 10     continue
 20   continue
!
!     back substitution
!
      do 40 kb = 1,nm1
        km1 = n-kb
        k = km1+1
        b(k) = b(k)/a(k,k)
        t = -b(k)
        do 30 i = 1, km1
          b(i) = b(i) + a(i,k)*t
 30     continue
 40   continue
 50   b(1) = b(1)/a(1,1)
      return
      end subroutine solveq
!
!
!BH130915:  Kinsey replaces sorspt with sorspt1
!
!$$$      subroutine sorspt(bshape, bheigh, bwidth, bhfoc, bvfoc,
!$$$     *  bhdiv, bvdiv, ib, ie, isourc, ke,
!$$$     *  nsourc, sfrac1, vbeam,
!$$$     *  x0, y0, z0, vx0, vy0, vz0)
!$$$      implicit integer (i-n), real*8 (a-h,o-z)
!$$$      save
!$$$ccc
!$$$ccc   generates a particle at the injector surface with
!$$$ccc   coordinates and velocities x0,y0,z0,vx0,vy0,vz0
!$$$ccc   These coordinates are attached to the source center, with the
!$$$ccc   x-direction perp to the source along the source centerline.
!$$$ccc
!$$$      character*8 bshape
!$$$      dimension bshape(*), bheigh(*), bwidth(*), bhfoc(*), bvfoc(*),
!$$$     *  bhdiv(*), bvdiv(*), sfrac1(*), vbeam(ke,*)
!$$$      integer iflag
!$$$!     drand(iflag) is real*8 [0.,1.] random number.
!$$$!     iflag=0 gives next number in the sequence.  See Absoft f77.
!$$$!     For Cray F90: use ranf.
!$$$
!$$$      data pio180/0.017453293/
!$$$      data rt2/1.414213562/
!$$$      data       seed0 /0.0d0/
!$$$ccc
!$$$
!$$$
!$$$      x0=0.
!$$$      isourc=1
!$$$      if(nsourc.eq.1)  go to 10
!$$$!     two sources
!$$$!     sfrac1 is fraction of source current coming from source 1 (upper
!$$$!     or rightmost source, assuming positive offsets).
!$$$!990131      if(ranf().gt.sfrac1(ib))  isourc=-1
!$$$
!$$$!000114*Not portable*      if(RANDOM_my(seed0).gt.sfrac1(ib))  isourc=-1
!$$$      rv=RANDOM_my(seed0)
!$$$      if(rv.gt.sfrac1(ib))  isourc=-1
!$$$ 10   if(bshape(ib).ne."circ") go to 20
!$$$ccc
!$$$!990131 12   y0=ranf()-0.5
!$$$!990131      z0=ranf()-0.5
!$$$
!$$$!000114*Not portable* 12   y0=drand(iflag)-0.5
!$$$!000114*Not portable*      z0=drand(iflag)-0.5
!$$$ 12   rv=RANDOM_my(seed0)
!$$$      y0=rv-0.5
!$$$      rv=RANDOM_my(seed0)
!$$$      z0=rv-0.5
!$$$
!$$$
!$$$      zsqu=y0**2+z0**2
!$$$      if(zsqu.gt.0.25) go to 12
!$$$      y0=y0*bwidth(ib)
!$$$      z0=z0*bwidth(ib)
!$$$!      write(*,*)'sorspt: rv,ib,bwidth(ib)',rv,ib,bwidth(ib)
!$$$      go to 30
!$$$ccc
!$$$ 20   continue
!$$$!990131      y0=bwidth(ib)*(ranf()-0.5)
!$$$!990131      z0=bheigh(ib)*(ranf()-0.5)
!$$$!000114*Not portable*      y0=bwidth(ib)*(drand(iflag)-0.5)
!$$$!000114*Not portable*      z0=bheigh(ib)*(drand(iflag)-0.5)
!$$$      rv=RANDOM_my(seed0)
!$$$      y0=bwidth(ib)*(rv-0.5)
!$$$      rv=RANDOM_my(seed0)
!$$$      z0=bheigh(ib)*(rv-0.5)
!$$$
!$$$ccc
!$$$ccc
!$$$ccc
!$$$ccc
!$$$ 30   continue
!$$$!
!$$$!     Special coding for DIII-D Long Pulse Sources.
!$$$      if(bshape(ib).ne."rect-lps")  go to 32
!$$$      if(abs(z0).gt.12.)  go to 31
!$$$!
!$$$!     particle on central 2 modules
!$$$      vdx=-1.
!$$$      vdy=0.0
!$$$      vdz=0.0
!$$$      go to 33
!$$$!
!$$$!     particle on upper or lower modules
!$$$ 31   vdx=-1.
!$$$      vdy=0.0
!$$$      vdz=-z0/bvfoc(ib)
!$$$      go to 33
!$$$!
!$$$ 32   vdx=-1.
!$$$      vdy=-y0/bhfoc(ib)
!$$$      vdz=-z0/bvfoc(ib)
!$$$ 33   vsqrt=1./sqrt(vdx**2+vdy**2+vdz**2)
!$$$      vx0=vbeam(ie,ib)*vdx*vsqrt
!$$$      vy0=vbeam(ie,ib)*vdy*vsqrt
!$$$      vz0=vbeam(ie,ib)*vdz*vsqrt
!$$$ccc
!$$$      thz=ranorm(1)*bvdiv(ib)/rt2
!$$$      thy=ranorm(1)*bhdiv(ib)/rt2
!$$$      vz0=vz0+thz*pio180*vx0
!$$$      vy0=vy0+thy*pio180*vx0
!$$$!      write(*,*)'sorspt: ie,ib,bvfoc(ib),y0,z0,vbeam(ie,ib)',
!$$$!     1                   ie,ib,bvfoc(ib),y0,z0,vbeam(ie,ib)
!$$$      return
!$$$      end
!
!
      subroutine sorspt1(nbshape,bheigh,bwidth,bhfoc,bvfoc, &
                        bhdiv,bvdiv,ib,ie,isourc,nsourc,sfrac1,vbeam, &
                        x0,y0,z0,vx0,vy0,vz0)
      use param_mod
!... generates a particle at the injector surface with
!    coordinates and velocities x0,y0,z0,vx0,vy0,vz0
!    These coordinates are attached to the source center, with the
!    x-direction perp to the source along the source centerline.
!
      implicit none 
      integer ib,ie,isourc,nsourc !input
      real(c_double) :: x0,y0,z0,vx0,vy0,vz0 !output
      real(c_double) :: seed0,rt2,pio180 !local
      real(c_double) :: vdx,vdy,vdz,thy,thz,zsqu,vsqrt !local
!      include 'params.inc'
!
      real(c_double) ::  bheigh(kb), bwidth(kb), bhfoc(kb), bvfoc(kb), &
              bhdiv(kb), bvdiv(kb), sfrac1(kb), vbeam(ke,kb)
      character*8 nbshape(kb)
!
      data pio180 /0.017453293D0/
      data rt2 /1.414213562D0/
      data seed0 /0.0/
!
      x0 = 0.D0
      isourc = 1
      if (nsourc.eq.1) goto 10
!
!     two sources (DIII-D case)
!     sfrac1 is fraction of source current coming from source 1 (upper
!     or rightmost source, assuming positive offsets).

!        if (RANF().gt.sfrac1(ib)) isourc = -1
      if (RANDOM_my(seed0) .gt. sfrac1(ib)) isourc = -1
   10 if (nbshape(ib).ne.'circ') go to 20
!
!   12 y0 = RANF() - 0.5D0
   12 y0 = RANDOM_my(seed0) - 0.5D0
!        z0 = RANF() - 0.5D0
      z0 = RANDOM_my(seed0) - 0.5D0
      zsqu = y0 ** 2 + z0 ** 2
      if (zsqu.gt.0.25) go to 12
      y0 = y0 * bwidth(ib)
      z0 = z0 * bwidth(ib)
      go to 30
!
  20  continue
      y0 = bwidth(ib) * (RANDOM_my(seed0)-0.5D0)
      z0 = bheigh(ib) * (RANDOM_my(seed0)-0.5D0)
!
!  Special coding for DIII-D Long Pulse Sources
!
   30 if (nbshape(ib) .ne. 'rect-lps')  go to 32
      if (   ABS (z0) .gt.  12.0     )  go to 31
!
!  particle on central 2 modules
!
      vdx = -1.0
      vdy =  0.0
      vdz =  0.0
      go to 33
!
!  particle on upper or lower modules
!
   31 vdx = -1.D0
      vdy = 0.D0
      vdz = -z0 / bvfoc(ib)
      go to 33
!
   32 vdx = -1.0
      vdy = -y0 / bhfoc(ib)
      vdz = -z0 / bvfoc(ib)
!
  33  vsqrt = 1.D0 / dsqrt(vdx**2+vdy**2+vdz**2)
      vx0 = vbeam(ie,ib) * vdx * vsqrt
      vy0 = vbeam(ie,ib) * vdy * vsqrt
      vz0 = vbeam(ie,ib) * vdz * vsqrt
!
      thz = ranorm(1) * bvdiv(ib) / rt2
      thy = ranorm(1) * bhdiv(ib) / rt2
      vz0 = vz0 + thz * pio180 * vx0
      vy0 = vy0 + thy * pio180 * vx0
!
      return
      end subroutine sorspt1
!
!
      subroutine timtor(rin, rmax, x0, y0, z0, vx0, vy0, vz0, &
        zmin, zmax, tenter, texit)
      implicit none 
      real(c_double) :: rin,rmax,x0,y0,z0,vx0,vy0,vz0 !input/output
      real(c_double) :: zmin,zmax,tenter,texit !input/output
      integer i,iin,isol !local
      real(c_double) :: zero,aa,bb,cc,tt,xx,yy,zz,rr,arg,sqr !local
!
!     This subroutine calculates the times for a particle to enter and
!     exit a toroidal box surrounding the plasma starting from the
!     point (x0,y0,z0) and moving with velocity (vx0,vy0,vz0).
!
      real(c_double) :: edge(4),timsol(6)
!
      zero=0.d0

!     specify edges of toroidal box surrounding plasma
!

!      write(*,*) 'x0,y0,z0,vx0,vy0,vz0,tenter,texit',
!     +     x0,y0,z0,vx0,vy0,vz0,tenter,texit

      edge(1)=rin
      edge(2)=rmax
      edge(3)=zmin
      edge(4)=zmax
!
!     find times for particle to intersect inside and outside of box
!
      isol=0
      aa=vx0**2+vy0**2
      if(aa.eq.zero) go to 50
      bb=2.*(vx0*x0+vy0*y0)
      do 40 i=1,2
        cc=x0**2+y0**2-edge(i)**2
        arg=bb**2-4.*aa*cc
        if(arg) 40,30,10
 10     sqr=sqrt(arg)
        tt=(-bb-sqr)/(2.*aa)
        zz=z0+vz0*tt
        if(zz.lt.zmin) go to 20
        if(zz.gt.zmax) go to 20
        isol=isol+1
        timsol(isol)=tt
 20     tt=(-bb+sqr)/(2.*aa)
        zz=z0+vz0*tt
        if(zz.lt.zmin) go to 40
        if(zz.gt.zmax) go to 40
        isol=isol+1
        timsol(isol)=tt
        go to 40
 30     tt=-bb/(2.*aa)
        zz=z0+vz0*tt
        if(zz.lt.zmin) go to 40
        if(zz.gt.zmax) go to 40
        isol=isol+1
        timsol(isol)=tt
 40   continue
!
!     find times for particle to intersect top and bottom of box
!
 50   if(vz0.eq.zero) go to 70
      do 60 i=3,4
        tt=(edge(i)-z0)/vz0
        xx=x0+vx0*tt
        yy=y0+vy0*tt
        rr=sqrt(xx**2+yy**2)
        if(rr.lt.rin) go to 60
        if(rr.gt.rmax) go to 60
        isol=isol+1
        timsol(isol)=tt
 60   continue
 70   continue
!
!     return if particle misses box
!
      tenter = -1.e10
      if(isol.eq.0) return
!
!     calculate times to enter and exit box
!
      tenter=timsol(1)
      iin=1
      do 80 i=2,isol
        if(timsol(i).ge.tenter) go to 80
        iin=i
        tenter=timsol(i)
 80   continue
      texit=1.e10
      do 90 i=1,isol
        if(i.eq.iin) go to 90
        if(timsol(i).lt.texit) texit=timsol(i)
 90   continue
!      write(*,*) 'x0,y0,z0,vx0,vy0,vz0,tenter,texit',
!     +     x0,y0,z0,vx0,vy0,vz0,tenter,texit

      return
      end subroutine timtor
!
!
      subroutine tozone(x,y,n,xz,yz,nz,nout,ncrt)
      implicit none 
      integer n,nz,nout,ncrt !input/output
      real(c_double) :: tol,yzlast,del,s !local
      integer ip,iz,ipm1 !local
      save
!     ---------------------------------------------------------------
!     convert from point values to zone values
!     ---------------------------------------------------------------
      real(c_double) :: x(*),y(*),xz(*),yz(*)
      data tol/1.0e-20/
!
      ip=1
      iz=1
      yzlast=y(n)
!     ---------------------------------------------------------------
!     interpolate to new grid that is compute yz(xz(iz))
!     ---------------------------------------------------------------
 2000 if(iz.gt.nz) go to 2210
 2100 if(ip.gt.n) go to 9200
      del=reldif(x(ip),xz(iz))
      if(del.gt.tol) go to 2110
      yz(iz)=y(ip)
      iz=iz+1
      go to 2000
 2110 if(x(ip).gt.xz(iz)) go to 2120
      ip=ip+1
      go to 2100
 2120 ipm1=ip-1
      s=(y(ip)-y(ipm1))/(x(ip)-x(ipm1))
      yz(iz)=y(ipm1)+(xz(iz)-x(ipm1))*s
      iz=iz+1
      go to 2000
!     ---------------------------------------------------------------
!     convert to zone
!     ---------------------------------------------------------------
 2210 continue
      do 2300 iz=1,nz-1
        yz(iz)=.5*(yz(iz)+yz(iz+1))
 2300 continue
      yz(nz)=.5*(yz(nz)+yzlast)
      return
!     ---------------------------------------------------------------
!     fatal errors
!     ---------------------------------------------------------------
 9200 continue
      write(nout,8000)
      write(ncrt,8000)
 8000 format(' fatal error in subroutine tozone',//, &
        ' xz(iz) or xz(izp1) is not contained in interval x(1),x(n)')
      stop
      end subroutine tozone

      real(c_double) function RANDOM_my(seed)
!
! --- basic uniform (0,1) pseudo-random number generator (portable) ----
!     From GA portlib.f (Jan. 14, 2000, bh).
!
!     BH131015 comment:
!     If called with seed.ne.0.d0, it will start a random number
!     sequence with the given seed.  Subsequent calls with seed=0.0d0
!     will get the ensuing RNs in a sequence.
!
      implicit none
!
      save         saved_seed
      real(c_double) :: seed, saved_seed, &
             big,                  big_minus_1
      data   big /2147483648.0d0/, big_minus_1 /2147483647.0d0/
!
      if (seed .ne. 0.0d0)  saved_seed = seed
      saved_seed = MOD (16807.0d0*saved_seed, big_minus_1)
      RANDOM_my  = saved_seed / big
      return
!
      end function RANDOM_my

!.......................................................................
!BH130915:  Rest of this file down through spline_12 subroutine are
!BH130915:  additions per Kinsey upgrade of freya stopping.
!.......................................................................

      subroutine inject_old (atw,codeid,drutpi,droti,dri,dzi,elongi, &
          ib,ie,mfm1,mim1,mjm1,newpar,psiax,psi,r,rmajor,rin,rmax,sgxn, &
          sgxnmi,x0,y0,z0,vx0,vy0,vz0,vbeam,z,zax,zmin,zmax,izone,pzone, &
          rzone,rpos,xpos,ypos,zpos)
      use param_mod !contains kb,kbe,ksge,ke,kz,ki,kion,kj,kcmp1
      implicit none 
      integer ib,ie,mfm1,mim1,mjm1,nebin,newpar,nout,izone,myid !input
      real(c_double) :: debin,drutpi,droti,dri,ds1,dzi !input
      real(c_double) :: elongi,psiax,rmajor,x0,y0,z0 !input
      real(c_double) :: vx0,vy0,vz0,smax,tenter,texit !input
      real(c_double) :: pzone,rzone,rpos,xpos,ypos,zpos !input
      real(c_double) :: dt1,delt,x1,y1,z1,r1,p1,vdotu,vrel1,vrel2,eova !local
      real(c_double) :: seed0,dum,psix,rzone2,tstep,dfac
      real(c_double) :: usq,area1,area2,area3,area4,tt,xnorm,ptest
      real(c_double) :: smin_time,smin_step
      integer i,j,ne1,n1,ir,iz,izone1,index,isol !local

!     This subroutine follows the particle from the pivot point
!     into, through, or around the plasma. It is the older version
!     of sub. inject to be used when crsecs (iexcit=0) is used.
!
!      include 'params.inc'
!      save
!
      real(c_double) :: psi(ki,kj), sgxn(kcmp1,kz,kbe,ksge), sgxnmi(ke,kb), &
           r(ki), z(kj), vbeam(ke,kb)
      real(c_double) :: atw(kion), rin, rmax, zax, zmin, zmax
      character*8 codeid
      data  seed0  /0.0d0/
!
      if (newpar.ne.0.or.isol.ne.0) then !YuP[2019-06-18]: isol not defined?
        if (newpar.ne.0) then
!
!... Calculate times for particle to enter and exit toroidal box
!    surrounding plasma:
!
!      write(*,*)
!      write(*,*) 'inside subroutine inject...',newpar
!      write(*,*) 'atw = ',atw
!      write(*,*) 'codeid = ',codeid
!      write(*,*) 'debin = ',debin
!      write(*,*) 'drutpi = ',drutpi
!      write(*,*) 'droti = ',droti
!      write(*,*) 'dri = ',dri
!      write(*,*) 'dzi = ',dzi
!      write(*,*) 'elongi = ',elongi
!      write(*,*) 'rin,rmax = ',rin,rmax
!      write(*,*) 'rmajor = ',rmajor
!      write(*,*) 'psiax = ',psiax
!      write(*,*) 'x0 = ',x0
!      write(*,*) 'y0 = ',y0
!      write(*,*) 'z0 = ',z0
!      write(*,*) 'ib,ie,kbe = ',ib,ie,kbe
!      write(*,*) 'ke,kz,ki = ',ke,kz,ki
!      write(*,*) 'mfm1,mim1,mjm1 = ',mfm1,mim1,mjm1
!      write(*,*) 'nebin = ',nebin
!      write(*,*) 'newpar = ',newpar
!
          call timtor(rin,rmax,x0,y0,z0,vx0,vy0,vz0,zmin,zmax,tenter, &
              texit)
!      write(*,*)
!      write(*,*) 'after timtor ...'
!      write(*,*) 'rin = ',rin
!      write(*,*) 'rmax = ',rmax
!      write(*,*) 'x0 = ',x0
!      write(*,*) 'y0 = ',y0
!      write(*,*) 'z0 = ',z0
!      write(*,*) 'vx0 = ',vx0
!      write(*,*) 'vy0 = ',vy0
!      write(*,*) 'vz0 = ',vz0
!      write(*,*) 'zmin = ',zmin
!      write(*,*) 'zmax = ',zmax
!      write(*,*) 'tenter = ',tenter
          if (tenter.le.-1.e10) go to 20
!
!... advance particle to edge of box
!
          x0 = x0 + vx0 * tenter
          y0 = y0 + vy0 * tenter
          z0 = z0 + vz0 * tenter
!
        endif ! newpar.ne.0
!
!... set coordinates and time for entering box
!
        xpos = x0
        ypos = y0
        zpos = z0
        tt = tenter
!
!... follow particle into plasma
!
   10   continue
!        dfac = -dlog(ranf())
        dfac = -LOG(RANDOM_my (seed0))
        tstep = dfac * sgxnmi(ie,ib) / vbeam(ie,ib)
        tt = tt + tstep
!        write(*,*) 'dfac = ',dfac
!        write(*,*) 'vbeam = ',vbeam(ie,ib)
!        write(*,*) 'sgxnmi = ',sgxnmi(ie,ib)
!        write(*,*) 'tstep = ',tstep
!        write(*,*) 'tt,texit = ',tt,texit
        if (tt.lt.texit) then
          xpos = xpos + vx0 * tstep
          ypos = ypos + vy0 * tstep
          zpos = zpos + vz0 * tstep
          rpos = dsqrt(xpos**2+ypos**2)
!
!... determine zone in which particle collides for 'onedee' geometry
!
          if (codeid.eq.'onedee') then
            rzone2 = (rpos-rmajor) ** 2 + (elongi*(zpos-zax)) ** 2
            rzone = SQRT (rzone2)
            izone = rzone * droti + 1.
          else
!
!... Determine zone in which particle collides for general geometry;
!    use bilinear interpolation away from magnetic axis and
!    biquadratic interpolation near the axis
!
            i = (rpos-r(1)) * dri + 1.
            j = (zpos-z(1)) * dzi + 1.
            if (i.gt.mim1) i = mim1
            if (j.gt.mjm1) j = mjm1
            psix = dmin1(psi(i,j),psi(i+1,j),psi(i,j+1),psi(i+1,j+1))
            ptest = (psix-psiax) * (drutpi/mfm1)
            if (ptest.ge.0.02) then
              area1 = (rpos-r(i)) * (zpos-z(j))
              area2 = (r(i+1)-rpos) * (zpos-z(j))
              area3 = (r(i+1)-rpos) * (z(j+1)-zpos)
              area4 = (rpos-r(i)) * (z(j+1)-zpos)
              pzone = (area3*psi(i,j)+area4*psi(i+1,j)+area1* &
                  psi(i+1,j+1)+area2*psi(i,j+1)) * dri * dzi
            else
              call pfit(psi(i-1,j-1),r(i-1),z(j-1),rpos,zpos,ki,pzone, &
                  dum,dum)
            endif
            pzone = MAX(pzone,psiax)
            izone = SQRT (pzone-psiax) * drutpi + 1.
!            write(*,*) 'izone-inject = ',izone
          endif ! codeid
!
!... If particle has psuedo-collision, continue following particle;
!    If particle has real collision, return
!
          if (izone.gt.mfm1) go to 10
!          if (ranf().gt.sgxn(izone,ie,ib)*sgxnmi(ie,ib)) go to 10
          if (RANDOM_my(seed0) .gt. &
             sgxn(4,izone,ie,ib)*sgxnmi(ie,ib))  go to 10
          return
!
        endif ! tt.lt.texit
!
      endif  ! newpar.ne.0.or.isol.ne.0
!
!... set flag to indicate that particle hit wall
!
   20 izone = mfm1 + 1
!
      return
      end subroutine inject_old


      subroutine inject1 (atw,codeid,debin,drutpi,droti,dri,ds1,dzi, &
                         elongi,ib,ie,kb,kbe,ksge,ke,kz,ki,mfm1,mim1, &
                         mjm1,nebin,newpar,nout,psiax,psi,r,rmajor,rin, &
                         rmax,sgxn,sgxnloc,sgxnmi,x0,y0,z0,vx0,vy0,vz0, &
                         vbeam,z,zangrot,zax,zmin,zmax,izone, &
                         pzone,rzone,rpos,xpos,ypos,zpos,myid,tenter, &
                         smax,texit)
!
      implicit none 
      integer ib,ie,kb,kbe,ksge,ke,kz,ki,mfm1,mim1,mjm1,nebin,newpar,nout,izone,myid !input
      real(c_double) :: atw,debin,drutpi,droti,dri,ds1,dzi !input
      real(c_double) :: elongi,psiax,rmajor,rin,rmax,x0,y0,z0 !input
      real(c_double) :: vx0,vy0,vz0,zax,zmin,zmax,smax,tenter,texit !input
      real(c_double) :: pzone,rzone,rpos,xpos,ypos,zpos !input
      real(c_double) :: dt1,delt,x1,y1,z1,r1,p1,vdotu,vrel1,vrel2,eova !local
      real(c_double) :: seed0,dum,psix,rzone2,tstep,dfac
      real(c_double) :: usq,area1,area2,area3,area4,tt,xnorm,ptest
      real(c_double) :: smin_time,smin_step
      integer i,j,ne1,n1,ir,iz,izone1,index !local
      
!
! ----------------------------------------------------------------------
!  This subroutine follows the particle from the pivot point into,
!  through, or around the plasma.  note:  bilinear interpolation is
!  used to obtain psi as a function of track length along a
!  collisionless neutral trajectory (rotating discharges only).
!  bicubic spline interpolation was tested, but found to provide no
!  appreciable increase in accuracy.
! ----------------------------------------------------------------------
!
      real(c_double) :: cvec(200)
      real(c_double) :: eova_wk(1) !YuP: for proper input to getsgxn
      integer izone_wk(1)  !YuP: for proper input to getsgxn
      
      data      cvec &
           /  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, &
             11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, &
             21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, &
             31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, &
             41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, &
             51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, &
             61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, &
             71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, &
             81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, &
             91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0,100.0, &
            101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,110.0, &
            111.0,112.0,113.0,114.0,115.0,116.0,117.0,118.0,119.0,120.0, &
            121.0,122.0,123.0,124.0,125.0,126.0,127.0,128.0,129.0,130.0, &
            131.0,132.0,133.0,134.0,135.0,136.0,137.0,138.0,139.0,140.0, &
            141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0, &
            151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0, &
            161.0,162.0,163.0,164.0,165.0,166.0,167.0,168.0,169.0,170.0, &
            171.0,172.0,173.0,174.0,175.0,176.0,177.0,178.0,179.0,180.0, &
            181.0,182.0,183.0,184.0,185.0,186.0,187.0,188.0,189.0,190.0, &
            191.0,192.0,193.0,194.0,195.0,196.0,197.0,198.0,199.0,200.0/
!
      integer ktk !local
      parameter (ktk = 100)
      character*8  codeid
      integer i1(ktk)
      real(c_double) :: sgxn(4,kz,kbe,ksge), sgxnloc(kbe), &
             sgxnmi(ke,kb), vbeam(ke,kb)
      real(c_double) :: r(*), z(*), psi(ki,*)
      real(c_double) :: zangrot(kz), e1(ktk), sgxntab(ktk)
      data   seed0    /0.0/
      real(c_double) :: kevperg, xmassp
      data xmassp / 1.67262310e-24 /    ! proton   mass   (g)
      data kevperg /6.2415064e+08/
!
!     the following assumes that data from previous
!     particle is saved. this ok only if passed through
!     argument list. HSJ
!
      do i=1,ktk
        i1(i)=0
        e1(i)=0.D0
        sgxntab(i)=0.D0
      enddo
      do i=1,kbe
        sgxnloc(i)=0.D0
      enddo
!
!      write(*,*)
!      write(*,*) 'inside subroutine inject...',newpar
!      write(*,*) 'atw = ',atw
!      write(*,*) 'codeid = ',codeid
!      write(*,*) 'debin = ',debin
!      write(*,*) 'drutpi = ',drutpi
!      write(*,*) 'droti = ',droti
!      write(*,*) 'dri = ',dri
!      write(*,*) 'ds1 = ',ds1
!      write(*,*) 'dzi = ',dzi
!      write(*,*) 'elongi = ',elongi
!      write(*,*) 'x0 = ',x0
!      write(*,*) 'y0 = ',y0
!      write(*,*) 'z0 = ',z0
!      write(*,*) 'ib,ie,kbe = ',ib,ie,kbe
!      write(*,*) 'ke,kz,ki = ',ke,kz,ki
!      write(*,*) 'mfm1,mim1,mjm1 = ',mfm1,mim1,mjm1
!      write(*,*) 'nebin = ',nebin
!      write(*,*) 'newpar = ',newpar
      if (newpar .eq. 0)  go to 100
!
! calculate times for particle to enter and exit toroidal box surrounding plasma
!
      call timtor (rin,rmax,x0,y0,z0,vx0,vy0,vz0,zmin,zmax,tenter,texit)
!      write(*,*)
!      write(*,*) 'after timtor ...'
!      write(*,*) 'rin = ',rin
!      write(*,*) 'rmax = ',rmax
!      write(*,*) 'x0 = ',x0
!      write(*,*) 'y0 = ',y0
!      write(*,*) 'z0 = ',z0
!      write(*,*) 'vx0 = ',vx0
!      write(*,*) 'vy0 = ',vy0
!      write(*,*) 'vz0 = ',vz0
!      write(*,*) 'zmin = ',zmin
!      write(*,*) 'zmax = ',zmax
!      write(*,*) 'tenter = ',tenter
      if (tenter .le. -1.0e10)  go to 140
!
! advance particle to edge of box
!
      x0 = x0 + vx0*tenter
      y0 = y0 + vy0*tenter
      z0 = z0 + vz0*tenter
!
      if (nebin .ne. 0) then
!
! ----------------------------------------------------------------------
! follow collisionless neutral trajectory to obtain minimum mean free path.
! this is required to account for toroidal rotation.
!
! ALSO:
!    Only include those energy groups that fall in a valid flux zone,
!    (izone1 .le. mfm1) for the search in subroutine GETSGXN. The array
!    of energy values e1 is now ne1 elements long rather than n1.
!                                            Daniel Finkenthal    9-6-95
! ----------------------------------------------------------------------
!
      smax = vbeam(ie,ib) * (texit-tenter)
      n1   = 2.0 + smax/ds1
      if (n1 .gt. ktk) then
        n1  = ktk
        ds1 = smax / FLOAT (n1-1)
        write (nout, 200) ds1
      endif
      dt1 = smax / (FLOAT (n1-1) * vbeam(ie,ib))
      ne1 = 0
!
      if (codeid .eq. 'onedee') then
        do i=1,n1
          delt   = (cvec(i) - 1.0) * dt1
          x1     = x0 + delt*vx0
          y1     = y0 + delt*vy0
          z1     = z0 + delt*vz0
          r1     = SQRT (x1**2 + y1**2)
          ir     = (r1-r(1))*dri + 1.0
          iz     = (z1-z(1))*dzi + 1.0
          ir     = MIN0 (ir,mim1)
          iz     = MIN0 (iz,mjm1)
          if (ir .le. 0 .or. iz .le. 0)go to 10  !HSJ 1/28/2000
          p1     = SQRT ((r1-rmajor)**2+(elongi*(z1-zax))**2)
          izone1 = p1*droti + 1.0
          if (izone1 .le. mfm1) then
            ne1     = ne1 + 1
            i1(ne1) = izone1
            usq     = ( x1**2 + y1**2 ) * zangrot(i1(ne1))**2
            vdotu   = (x1*vy0 - y1*vx0) * zangrot(i1(ne1))
            vrel1   = vbeam(ie,ib)**2 + usq - 2.0*vdotu
            e1(ne1) = 0.5 * xmassp * vrel1 * kevperg
            e1(ne1) = ABS (e1(ne1))
          endif
        enddo
      else
        do i=1,n1
          delt  = (cvec(i) - 1.0) * dt1
          x1    = x0 + delt*vx0
          y1    = y0 + delt*vy0
          z1    = z0 + delt*vz0
          r1    = SQRT (x1**2+y1**2)
          ir    = (r1-r(1))*dri + 1.0
          iz    = (z1-z(1))*dzi + 1.0
          ir    = MIN0 (ir,mim1)
          iz    = MIN0 (iz,mjm1)
          if (ir .le. 0 .or. iz .le. 0)go to 10  !HSJ 1/28/2000
          area1 = (r1-r(ir))*(z1-z(iz))
          area2 = (r(ir+1)-r1)*(z1-z(iz))
          area3 = (r(ir+1)-r1)*(z(iz+1)-z1)
          area4 = (r1-r(ir))*(z(iz+1)-z1)
          p1    = (area3*psi(ir,iz)   + area4*psi(ir+1,iz) &
                +  area2*psi(ir,iz+1) + area1*psi(ir+1,iz+1))*dri*dzi
          p1     =  MAX (p1,psiax)
          izone1 = SQRT (p1-psiax)*drutpi + 1.0
          if (izone1 .le. mfm1) then
            ne1     = ne1 + 1
            i1(ne1) = izone1
            usq     = (x1**2+y1**2)*zangrot(i1(ne1))**2
            vdotu   = (x1*vy0-y1*vx0)*zangrot(i1(ne1))
            vrel1   = vbeam(ie,ib)**2 + usq - 2.0*vdotu
            e1(ne1) = 0.5 * xmassp * vrel1 * kevperg
            e1(ne1) = ABS (e1(ne1))
          endif
        enddo
      endif
!
!      write(*,*)
!       write(*,*) 'calling getsgxn ...'
!       write(*,*) 'ne1 = ',ne1
!      write(*,*) 'e1(ne1) = ',e1(ne1)
!      write(*,*) 'i1(ne1) = ',i1(ne1)
!      write(*,*) 'vrel1 = ',vrel1
!      write(*,*) 'usq = ',usq
!      write(*,*) 'vdotu = ',vdotu
!      write(*,*) 'zangrot  = ',zangrot(i1(ne1))
!      write(*,*) 'ib = ',ib
!      write(*,*) 'ie = ',ie
!      write(*,*) 'nebin = ',nebin
!      write(*,*) 'debin = ',debin,' in inject'
!      write(*,*) 'sgxn(4,22,3,2)-pre = ',sgxn(4,22,3,2)
!
 10   call getsgxn (e1, i1, ne1, ktk, ib, ie, sgxn, nebin, debin, kbe, &
                    ksge, kz, nout, sgxntab)
!      write(*,*) 'after call to getsgxn ...'
!      write(*,*) 'ne1 = ',ne1
!      do i=1,ne1
!       write(*,*) i, sgxntab(i), ' sgxntab'
!      enddo
      xnorm         = amaxaf (sgxntab, 1, ne1)
      sgxnmi(ie,ib) = 1.0 / xnorm
!      write(*,*) 'sgxnmi(ie,ib) = ',sgxnmi(ie,ib)
      end if
!      stop
!
! ----------------------------------------------------------------------
! inject neutral into plasma
! ----------------------------------------------------------------------
!
!... set coordinates and time for entering box
!
  100 xpos      = x0
      ypos      = y0
      zpos      = z0
      tt        = tenter
      izone     = mfm1 + 1     ! initially neutral is outside the plasma
      smin_step = 0.1                           ! 0.1 cm min step or
      smin_step = MIN (smin_step, smax/1000.0)  ! make scale-independent
      smin_time = smin_step / SQRT (vx0**2 + vy0**2 + vz0**2)
!
!... follow particle into plasma
!
  110 dfac  = -LOG (RANDOM_my (seed0))
!
!     if neutral is not yet in the plasma (izone ge mf) then take steps
!     of minimum size 1 mm until we enter the plasma. we could find the
!     exact plasma boundary but that would be overkill. with 1mm step
!     size we certainly are within any physics scales we could resolve
!     near the plasma edge. this avoids wasting a lot of steps outside
!     the plasma and will be a significant savings if the cross sections
!     are large. dfac is 1/exponentially distributed, so for large steps
!     don't modify tstep. --------------------------- 27 Oct 95 ---- HSJ
!
      if (izone .le. mfm1) then                     !  inside the plasma
        tstep =      dfac * sgxnmi(ie,ib) / vbeam(ie,ib)
      else
        tstep = MAX (dfac * sgxnmi(ie,ib) / vbeam(ie,ib), &
                     smin_time)                     ! outside the plasma
      endif
!
      tt    = tt  + tstep
      if (tt .ge. texit)  go to 140
      xpos = xpos + vx0*tstep
      ypos = ypos + vy0*tstep
      zpos = zpos + vz0*tstep
      rpos = SQRT (xpos**2 + ypos**2)
!
!  determine zone in which particle collides for 'onedee' geometry
!
      if (codeid .eq. 'onedee') then
        rzone2 = (rpos-rmajor)**2 + (elongi*(zpos-zax))**2
        rzone  = SQRT (rzone2)
        izone  = rzone*droti + 1.0
      else
!
!  determine zone in which particle collides for general geometry;
!     use bilinear interpolation away from magnetic axis,
!     and biquadratic interpolation near the axis.
!
        i     = (rpos-r(1))*dri+1.0
        j     = (zpos-z(1))*dzi+1.0
        i     = MIN0 (i,mim1)
        j     = MIN0 (j,mjm1)
        psix  = MIN  (psi(i,j),psi(i+1,j),psi(i,j+1),psi(i+1,j+1))
        ptest = (psix-psiax)*(drutpi/mfm1)**2
        if (ptest .ge. 0.02) then
          area1 = (rpos-r(i))*(zpos-z(j))
          area2 = (r(i+1)-rpos)*(zpos-z(j))
          area3 = (r(i+1)-rpos)*(z(j+1)-zpos)
          area4 = (rpos-r(i))*(z(j+1)-zpos)
          pzone = (area3*psi(i,j) + area4*psi(i+1,j) &
                + area1*psi(i+1,j+1) + area2*psi(i,j+1))*dri*dzi
        else
         !YuP160330: Avoiding out of bounds occuring near R=0 for mirror geometry
         if (i.ne.1) then
         call pfit(psi(i-1,j-1), r(i-1), z(j-1), rpos, zpos, ki, pzone, &
         dum, dum)
         else
         call pfit(psi(1,j-1), r(1), z(j-1), rpos, zpos, ki, pzone, &
         dum, dum)
         endif

        endif
        pzone =  MAX (pzone,psiax)
        izone = SQRT (pzone-psiax)*drutpi + 1.0
!        write(*,*) '** izone = ',izone,psiax
      endif
!
!... if particle has psuedo-collision, continue following particle.
!    if particle has real collision, return.
!
      if (izone .gt. mfm1)  go to 110 ! the particle is inside the box..
!                                     ! but still outside the plasma
      if (nebin .ne. 0) then
        usq   = (rpos*zangrot(izone))**2
        vdotu = (xpos*vy0-ypos*vx0)*zangrot(izone)
        vrel2 = vbeam(ie,ib)**2 + usq - 2.0*vdotu
        eova  = 0.5 * xmassp*vrel2*kevperg
        eova  = ABS (eova)
!        write(*,*) 'eova = ',eova
!        write(*,*) 'vrel2 = ',vrel2
        eova_wk(1)=eova !YuP: as an array of size 1
        izone_wk(1)=izone !YuP: as an array of size 1
        call getsgxn (eova_wk, izone_wk, 1, ktk, ib, ie, sgxn, nebin, &
                      debin, kbe, ksge, kz, nout, sgxnloc)
      else
        index      = ke*(ib-1) + ie
        sgxnloc(1) = sgxn(1,izone,index,1)
        sgxnloc(2) = sgxn(2,izone,index,1)
        sgxnloc(3) = sgxn(3,izone,index,1)
        sgxnloc(4) = sgxn(4,izone,index,1)
      end if
!     if (RANF   (     ) .gt. sgxnloc(4) * sgxnmi(ie,ib))  go to 110
      if (RANDOM_my (seed0) .gt. sgxnloc(4) * sgxnmi(ie,ib))  go to 110
      return
!
!... set flag to indicate that particle hit wall
!
  140 izone = mfm1 + 1
  200 format (' WARNING from subroutine INJECT1:'                  / &
              '         maximum number of grid elements exceeded' / &
              '         increasing ds1 to ', e10.3, ' cm')
      return
!
      end subroutine inject1

!
!
      subroutine nbsgxn (iexcit,namep,namei,mb,mfm1,nebin,nprim,nimp, &
                         nion,atw_beam,atw,ebkev,ebfac,ibion,vbeam, &
                         zne,zni,zte,zti,zzi,debin, &
                         dtemax,dnemax,dzemax,hxfrac,sgxn,sgxnmi)
!
! ----------------------------------------------------------------------
!
!  this subroutine calculates the neutral beam attenuation array, sgxn.
!
!     input:
!          iexcit         - atomic excitation model
!                           0=Freeman-Jones, 5=ADAS, 6=Boley parameterization
!          namei          - name of ith impurity ion species
!          nprim, nimp
!          atw(k)         - atomic mass of ion species k
!          ebkev(mb)      - maximum energy of mb-th neutral beamline
!                           (keV)
!          ebfac          - factor defining upper bound on energy bin
!                           range,  > 1.0
!          atw_beam         mass no. of beam
!
!          mb             - number of beamlines modeled
!          mfm1           - number of flux zones minus one
!          nebin          - number of energy bins (rotation case only)
!          vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
!                           beamline (cm/sec)
!          zne(mfm1)      - local electron density (cm**-3)
!          zni(mfm1,nion) - local density of nion-th ion species
!                           (cm**-3)
!          zzi(mfm1,nion) - local average charge state of nion-th ion
!                           species
!          dtemax,dnemax,dzemax - used in crsecs Stearn formula
!     ONETWO stuff:
!     input from common /io/:
!          ncrt,nout,nouthx,ncorin
!     input from common /ions/:
!          namei, atw
!     input from common /nub3/:
!          iexcit,ilorent,mstate,ncont,kdene,kdeni,kdenz,ksvi,ksvz,ksve,
!          krad,ngh,ngl,znipm,atwpm,iz,zti,izstrp
!     input from common /numbrs/:
!          nprim,nimp,nion
!     END ONETWO stuff
!
!     output:
!          sgxn(i,j,k,l)
!               i - mode index
!                   = 1, fraction of reactions producing electrons;
!                   = 2, fraction of reactions producing species 1 ion;
!                   = 3, fraction of reactions producing species 2 ion;
!                   = 4, total inverse mean free path;
!               j - FREYA zone index
!               k - beam index
!                   k = 3*(ib-1) + ie, where ib is the beam index and ie
!                                      is the beam energy group index
!               l - index for relative energy.  bins are equispaced in
!                   delta(energy) between 0 and max(ebkev(ib))*ebfac,
!                   with bin width given by
!                             delta(energy) = ebmax*ebfac/nebin.
!                   ebfac .ge. 1.0 and nebin are user supplied namelist
!                   variables.
!          sgxnmi(ie,mb)
!                 - minimum inverse mean free path for ie-th energy
!                   group of mb-th beamline.  calculated for stationary
!                   target case only.  otherwise, calculated for each
!                   neutral trajectory in subroutine INJECT.
!          hxfrac(ie,mb)
!                 - fractional density of n = 3 excited state neutral.
!                   calculated in stationary target case only.
!          debin
!                 - width of energy bins (keV/amu)
! ----------------------------------------------------------------------
!
      use param_mod
      implicit none 
      integer iexcit,nebin,nimp !input
!
!      include 'params.inc'  ! both of these .inc's used in stand alone
!      include 'fileiou.inc'
!
      real(c_double) :: debin,dtemax,dnemax,dzemax !input
      real(c_double) :: tiev,teev,rmfp,ebmax,ebin,bperp,eova,vbin,hexfrac
      integer ksrc, i,j,k,ie,ib,ind,istart,nout,nouthx,ncrt,ihxerr !local
      parameter (ksrc = 3)
      character*8 namei(kimp), namep(kprim)   ! added for cql3d
      integer mb, mfm1, nion, nprim
      real(c_double) ::  ebkev(kb), vbeam(ke,kb)
      real(c_double) ::  hxfrac(ke,kb),sgxnmi(ke,kb), sgxne(kz)
      real(c_double) ::  sgvxne(kz), sgxn(kcmp1,kz,kbe,ksge)  ! stand alone code
!      real(c_double) ::  sgvxne(kz), sgxn(kcmp1,kz,ke,kb)
! here, atw(kion)=atw(5), atwpm(kprim)=atwpm(3), atwim(kimp)=atwim(2)
! and iz(kimp)=iz(2) in ONETWO
! nimpeqdsk, nprimedsk, and nioneqdsk are read in in prenub.f
! kcm = 3, kcmp1 = kcm + 1 in param.i in Onetwo
!
      integer ibion, iz(kion)
      real(c_double) :: ebfac, atw_beam
      real(c_double) :: atw(kion), atwpm(kion), atwim(kion)
      real(c_double) :: znipm(kion), zniim(kion)
      real(c_double) :: zne(kz), zni(kz,kion), zte(kz), zti(kz), zzi(kz,kion)
!
! ----------------------------------------------------------------------
! Original Freeman-Jones routine (iexcit.le.0)
!   -1 for Stearn formula
!    0 for Freeman-Jones
! ----------------------------------------------------------------------
      if (iexcit.le.0) then
        write(*,*) 'Calling Freeman-Jones routine...'
        call crsecs(iexcit,atw,ebkev,ibion,mb,mfm1,nion,vbeam, &
                    zne,zni,zte,zzi,dtemax,dnemax,dzemax, &
                    sgvxne,sgxn,sgxnmi)
!
!         do ib = 1, mb
!           do j = 1, 3
!             do i=1,mfm1
!        write(*,'(3i3,2x,1p1e10.4,a12)') i,j,ib,sgxn(4,i,j,ib), ' sgxn'
!             enddo
!           enddo
!         enddo
!        stop
        return
      endif
! ----------------------------------------------------------------------
! Use ADAS routines (iexcit = 5)
! ----------------------------------------------------------------------
      if (iexcit .eq. 5) then
        write(*,*) 'Calling adassgxn routine...'
        call adassgxn (namep,namei,mb,mfm1,nprim,nimp,nion,atw, &
                       ebkev,ebfac,ibion,nebin,vbeam,zne,zni,zte, &
                       zti,zzi,debin,sgxn,sgxnmi,atw_beam)
!
!        do ib = 1, mb
!          do j = 1, 3
!            do i=1,mfm1
!             write(*,'(3i4,2x,1p1e10.4,a16)')
!     &         i,j,ib,sgxn(4,i,j,ib), ' sgxn-nbsgxn'
!            enddo
!          enddo
!        enddo
!        stop
!
        return
      endif
!
! ---------------------------------------------------- HSJ-2/5/98 ------
! Boley's parametrization including MSI effects (iexcit=6)
! ----------------------------------------------------------------------
       if (iexcit .eq. 6) then
!         call wrap_xboley (atw,zzi,ebkev,ibion,mb,mfm1,nebin,ebfac,
!     &                     debin,zne,zte,zni,nion,sgxn,sgxnmi,atw_beam)
         return
       endif
!
! ----------------------------------------------------------------------
! HEXNB (iexcit.ne.0)
! ----------------------------------------------------------------------
!
      if (iexcit .ne. 0) then
        if (nouthx .gt. 0) &
        open (unit = nouthx, file = 'hexnbout', status = 'UNKNOWN')
        istart = 1
!
!       separate primary and impurity ions
!
        do j=1,nion
          if (j .le. nprim) then
            atwpm(j) = atw(j)
          else
            k = j - nprim
            atwim(k) = atw(j)
          endif
        enddo
!
!       get atomic number of impurities
!
        do i=1,nimp
          if(trim(namei(i)).eq.'he' .or. &
             trim(namei(i)).eq.'HE' .or. &
             trim(namei(i)).eq.'He')  iz(i) =  2

          if(trim(namei(i)).eq.'b' .or. &
             trim(namei(i)).eq.'B' )  iz(i) =  5  ! YuP added [2015]

          if(trim(namei(i)).eq.'c' .or. &
             trim(namei(i)).eq.'C' )  iz(i) =  6

          if(trim(namei(i)).eq.'o' .or. &
             trim(namei(i)).eq.'O' )  iz(i) =  8

          if(trim(namei(i)).eq.'si' .or. &
             trim(namei(i)).eq.'SI' .or. &
             trim(namei(i)).eq.'Si')  iz(i) = 14

          if(trim(namei(i)).eq.'ar' .or. &
             trim(namei(i)).eq.'AR' .or. &
             trim(namei(i)).eq.'Ar')  iz(i) = 18

          if(trim(namei(i)).eq.'cr' .or. &
             trim(namei(i)).eq.'CR' .or. &
             trim(namei(i)).eq.'Cr')  iz(i) = 24

          if(trim(namei(i)).eq.'fe' .or. &
             trim(namei(i)).eq.'FE' .or. &
             trim(namei(i)).eq.'Fe')  iz(i) = 26

          if(trim(namei(i)).eq.'ni' .or. &
             trim(namei(i)).eq.'NI' .or. &
             trim(namei(i)).eq.'Ni')  iz(i) = 28

          if(trim(namei(i)).eq.'kr' .or. &
             trim(namei(i)).eq.'KR' .or. &
             trim(namei(i)).eq.'Kr')  iz(i) = 36

          if(trim(namei(i)).eq.'mo' .or. &
             trim(namei(i)).eq.'MO' .or. &
             trim(namei(i)).eq.'Mo')  iz(i) = 42

          if(trim(namei(i)).eq.'w' .or. &
             trim(namei(i)).eq.'W' )  iz(i) = 74
        enddo
!
!... no Lorentz ionization limit
!
        bperp = 0.0
      endif
!
! ----------------------------------------------------------------------
! stationary plasma case
! ----------------------------------------------------------------------
!
      if (nebin .eq. 0) then
        do ib=1,mb
          do ie=1,ke
            sgxnmi(ie,ib) = 0.D0
            ind  = ke*(ib-1) + ie
            eova = 1.D3*ebkev(ib)/(FLOAT (ie)*atw_beam)
            vbin = vbeam(ie,ib)
            do i=1,mfm1
              teev = 1.D3*zte(i)
              if (iexcit .ne. 0) then
                tiev = 1.D3*zti(i)
                do j=1,nion
                  if (j .le. nprim) then
                    znipm(j) = zni(i,j)
                  else
                    k = j - nprim
                    zniim(k) = zni(i,j)
                  endif
                enddo
              endif
              sgxne(ind) = fsgxne(vbin,teev,zne(i))
              sgxn(1,i,ind,1) = 0.D0
              do k=1,nion
                sgxn(1,i,ind,1) = sgxn(1,i,ind,1) &
                                + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
                if ((k .le. nprim) .and. (k .le. 2.0)) &
                  sgxn(k+1,i,ind,1) = fsgxncx(atw(k),eova,zni(i,k))
              enddo
              sgxn(1,i,ind,1) = sgxn(1,i,ind,1) + sgxne(ind)
              sgxn(4,i,ind,1) = sgxn(1,i,ind,1) + sgxn(2,i,ind,1) &
                              + sgxn(3,i,ind,1)
              sgxn(1,i,ind,1) = sgxn(1,i,ind,1)/sgxn(4,i,ind,1)
              sgxn(2,i,ind,1) = sgxn(2,i,ind,1)/sgxn(4,i,ind,1)
              sgxn(3,i,ind,1) = sgxn(3,i,ind,1)/sgxn(4,i,ind,1)
              if (iexcit .ne. 0) then
!                call hexnb (istart, 1, ilorent, mstate, ncont, eova,
!     .                      teev, tiev, nprim, atwpm, znipm, nimp, iz,
!     .                      atwim, izstrp, zniim, bperp, kdene, kdeni,
!     .                      kdenz, ksvi, ksvz, ksve, krad, ngh, ngl,
!     .                      nouthx, ncorin, rerate, rmfp, hexfrac,
!     .                      ihxerr)
                if (ihxerr .ne. 0) write (ncrt, 1000) ihxerr
                if (ihxerr .eq. 1) then
                  write (ncrt, 1010)
                  write (nout, 1010)
                  call STOP ('subroutine NBSGXN: problem #1', 45)
                end if
                sgxn(4,i,ind,1) = 1.D0/rmfp
                if (i .eq. 1) hxfrac(ie,ib) = hexfrac
              endif
              sgxnmi(ie,ib) = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
            enddo
          enddo
        enddo
        do ib=1,mb
          do ie=1,ke
            sgxnmi(ie,ib) = 1.D0 / sgxnmi(ie,ib)
          enddo
        enddo
      else
!
! ----------------------------------------------------------------------
! rotating plasma case
! ----------------------------------------------------------------------
!
        ebmax = 0.D0
        do ib=1,mb
          ebmax = MAX (ebmax,ebkev(ib))
        enddo
        ebmax = ebmax/atw_beam
        debin = ebmax*ebfac/FLOAT (nebin)
!
        do 340 i=1,mfm1
        teev = 1.D3 * zte(i)
        if (iexcit .ne. 0) then
          tiev = 1.D3 * zti(i)
          do j=1,nion
            if (j .le. nprim) then
              znipm(j) = zni(i,j)
            else
              zniim(j-nprim) = zni(i,j)
            endif
          enddo
        endif
!
!       electron impact ionization independent of rotation speed
!
        do ib=1,mb
          do ie=1,ke
            ind = ke*(ib-1) + ie
            sgxne(ind) = fsgxne(vbeam(ie,ib),teev,zne(i))
          enddo
        enddo
!
        do 340 j=1,nebin
        ebin = FLOAT (j) * debin
        eova = 1.D3 * ebin
        sgxn(1,i,1,j) = 0.D0
        do k=1,nion
          sgxn(1,i,1,j) = sgxn(1,i,1,j) &
                        + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
          if ((k .le. nprim) .and. (k .le. 2.0)) &
            sgxn(k+1,i,1,j) = fsgxncx(atw(k),eova,zni(i,k))
        enddo
        if (iexcit .ne. 0) then
!          call hexnb (istart, 1, ilorent, mstate, ncont, eova,
!     .                teev, tiev, nprim, atwpm, znipm, nimp, iz,
!     .                atwim, izstrp, zniim, bperp, kdene, kdeni,
!     .                kdenz, ksvi, ksvz, ksve, krad, ngh, ngl, nouthx,
!     .                ncorin, rerate, rmfp, hexfrac, ihxerr)
          if (ihxerr .ne. 0) write (ncrt, 1000) ihxerr
          if (ihxerr .eq. 1) then
            write (ncrt, 1010)
            write (nout, 1010)
            call STOP ('subroutine NBSGXN: problem #2', 46)
          endif
        endif
        do 340 k=ind,1,-1
          sgxn(1,i,k,j) = sgxn(1,i,1,j) + sgxne(k)
          sgxn(4,i,k,j) = sgxn(1,i,k,j) + sgxn(2,i,1,j) + sgxn(3,i,1,j)
          sgxn(1,i,k,j) = sgxn(1,i,k,j)/sgxn(4,i,k,j)
          sgxn(2,i,k,j) = sgxn(2,i,1,j)/sgxn(4,i,k,j)
          sgxn(3,i,k,j) = sgxn(3,i,1,j)/sgxn(4,i,k,j)
          if (iexcit .ne. 0)  sgxn(4,i,k,j) = 1.D0/rmfp
  340   continue
      endif
!
 1000 format (' subroutine NBSGXN reports a HEXNB return code of ', i5)
 1010 format (' ERROR: execution terminated - file "coronb" not found')
!
      return
!
      end subroutine nbsgxn
      
!
      real(c_double) function fsgxncx (atw, e, zni)
      implicit none 
      real(c_double) :: atw, e, zni !input
      real(c_double) :: aloge,sigcx
!
! ----------------------------------------------------------------------
! this subprogram calculates inverse mean free path due to charge exchange
! from the fitted results of freeman and jones, clm-r137, culham (1974).
!
!     input:
!             atw - atomic weight of target ion
!             e   - relative energy of impinging neutral (ev/amu)
!             zni - density of target ion (cm**-3)
! ----------------------------------------------------------------------
!
      if (atw .gt. 3.01) then
        sigcx = 0.0
      else
        aloge = LOG10 (e)
        sigcx = 0.6937e-14 * (1.0 - 0.155*aloge)**2 / &
                             (1.0 + 0.1112e-14*e**3.3)
      endif
      fsgxncx = sigcx*zni
!
      return
      end function fsgxncx
      

      real(c_double) function fsgxne (vb, te, zne)
      implicit none 
      real(c_double) :: vb, te, zne !input
      real(c_double) :: alogt,expo !local
!
! ----------------------------------------------------------------------
! this subprogram evaluates local inverse mean free path for electron impact
! ionization from fitted results of freeman and jones, clm-r137,culham (1974).
!
!     input:
!             vb  - speed of impinging neutral (cm/sec)
!             te  - target electron temperature (ev)
!             zne - target electron density (cm**-3)
! ----------------------------------------------------------------------
!
      real(c_double) :: cfione(7)
      data      cfione /-3.173850e+01,  1.143818e+01, -3.833998    , &
                         7.046692e-01, -7.431486e-02,  4.153749e-03, &
                        -9.486967e-05/
!
      alogt = 0.0
      if (te .gt. 1.0    )  alogt = LOG (te)
      if (te .gt. 1.0e+05)  alogt = 11.51
      expo = (((((cfione(7) *alogt + cfione(6))*alogt + cfione(5))*alogt &
                + cfione(4))*alogt + cfione(3))*alogt + cfione(2))*alogt &
                + cfione(1)
      fsgxne = EXP (expo) * zne / vb
      return
!
      end function fsgxne
      
!
      real(c_double) function fsgxni (atw, eova, zni, zzi)
      implicit none 
      real(c_double) :: atw, eova, zni, zzi !input
      real(c_double) :: aloge,sigi,expo,ekev !local
!
! ----------------------------------------------------------------------
!  this subprogram calculates inverse mean free path due to proton and
!     impurity impact ionization.  proton impact ionization cross
!     sections from fitted results of freeman and jones,clm-r137, culham
!     (1974).  impurity impact ionization cross sections from r.e. olson
!     et al., Phys. Rev. Lett. 41, 163 (1978).
!
!     input:
!             atw  - atomic weight of target ion
!             eova - relative energy of impinging neutral (ev/amu)
!             zni  - density of target ion (cm**-3)
!             zzi  - average charge state of target ion
! ----------------------------------------------------------------------
!
      real(c_double) :: cfionp(7)
      data      cfionp /-4.203309e+01,  3.557321    , -1.045134, &
                         3.139238e-01, -7.454475e-02,  8.459113e-03, &
                        -3.495444e-04/
!
      if (atw .le. 3.01) then
        aloge = LOG10 (eova)
        aloge = aloge * 2.302585093 - 6.907755279
        if (aloge .le. -2.30258) then
          sigi = 0.0
        else
          expo = (((((cfionp(7) *aloge + cfionp(6))*aloge &
                    + cfionp(5))*aloge + cfionp(4))*aloge &
                    + cfionp(3))*aloge + cfionp(2))*aloge &
                    + cfionp(1)
          sigi = EXP (expo)
        endif
        fsgxni = sigi*zni
      else
        ekev   = 1.0e-3*eova
        fsgxni = 1.0e-17*zni*46.0*zzi*(32.0*zzi/ekev)* &
                    (1.0 - EXP (-ekev/(32.0*zzi)))
      endif
      return
!
      end function fsgxni
!
!
      subroutine getsgxn (e, iz, ns, ktk, ib, ie, sgxn, nbins, &
                          debin, kbe, ksge, kz, nout, sgxnloc)
!
      implicit none 
!
! ----------------------------------------------------------------------
! this subroutine locates n*sigma in the lookup table, sgxn
! ----------------------------------------------------------------------
!
! --- input through argument list:
!          debin          - energy bin width     (keV/amu)
!          e(ns)          - energy to be located (keV/amu)
!          ib             - beam line index
!          ie             - beam energy group index
!          iz(ns)         - flux zone index
!          nbins          - total number of energy bins
!          ns             - total number of points to be evaluated
!          sgxn(i,iz,j,k) - array of neutral stopping data
!                           i  - data type index
!                                =1, fraction of interactions producing
!                                    electrons;
!                                =2, fraction of interactions producing
!                                    neutral species 1;
!                                =3, fraction of interactions producing
!                                    neutral species 2;
!                                =4, inverse mean free path (cm**-1)
!                           iz - flux zone index
!                           j  - beam/energy index, j = 3*(ib-1)+ie
!                           k  - energy bin index
!
! --- output through argument list:
!          sgxnloc(i)     - local neutral stopping data (cm**-1)
!                           i - data type index (see above)
!                           note:  iftns>1, only the inverse mean free
!                           path is evaluated (i = 4) and ns data points
!                           are passed.  this facilitates faster
!                           execution when evaluating the max imfp along
!                           the collisionless neutral trajectory (used
!                           in rotating discharges only).  otherwise
!                           more detailed information is passed in the
!                           first four elements only (ns = 1).
! ----------------------------------------------------------------------
!
!   Created:    Glenn Sager?               Date:  ??-???-????
!
!   changes from original version:
!      1) The code will now terminate if the maximum rotational energy
!         bin is exceeded. Detailed Error message provide to both the
!         OUTONE file (nout) and the standard output device (ncrt = 6).
!      2) Fixed some math in array index calculations.
!      3) Today is the 50th Aniversary of Hiroshima. NEVER AGAIN!
!
!                                  D.F. Finkenthal 6-AUG-95
!
! ----------------------------------------------------------------------
!
      real(c_double) :: e(ns) !YuP[2019-06-17]: was e(ktk)   BUG?
      real(c_double) :: sgxn(4,kz,kbe,ksge), sgxnloc(kbe)  ! stand alone code
      integer iz(ns) !YuP[2019-06-17]: was iz(ktk)   BUG?
      integer imaxa(200)
      integer ns,ktk,ib,ie,nbins,kbe,ksge,nout !input/output
      integer i,ind,ncrt,ibin,kz !local
      real(c_double) :: debin !input
      real(c_double) :: emax !local
      
      integer,save :: imax
!
      ncrt = 6
!
      ind = 3 * (ib - 1) + ie
      if (ns .gt. 1) then
        do i=1,ns
          ibin       = e(i) / debin + 1.0
          imaxa(i)   = ibin
          ibin       = MIN0 (ibin, nbins)
          sgxnloc(i) = sgxn (4, iz(i), ind, ibin)
        end do
        imax = maxaf (imaxa, 1, ns)
      else
        ibin       = e(1) / debin + 1.0
        imax       = MAX0 (ibin, imax )
        ibin       = MIN0 (ibin, nbins)
        sgxnloc(1) = sgxn (1, iz(1), ind, ibin)
        sgxnloc(2) = sgxn (2, iz(1), ind, ibin)
        sgxnloc(3) = sgxn (3, iz(1), ind, ibin)
        sgxnloc(4) = sgxn (4, iz(1), ind, ibin)
      end if
!      write(*,*) 'sgxn(4,22,3,2) = ',sgxn(4,22,3,2)
!      if(ie.eq.2) stop
!
      if (imax .gt. nbins) then
        emax = amaxaf (e, 1, ns)
        write (nout, 110)  nbins * debin, emax
        write (ncrt, 110)  nbins * debin, emax
        call STOP ('subroutine GETSGXN: max bin energy exceeded', 51)
      end if
!
 110  format (/ &
       ' ERROR in GETSGXN: maximum rotational energy bin exceeded'     / &
       ' The highest calculated energy bin is now ', f10.2, ' keV/amu' / &
       ' but GETSGXN tried to look up a value of  ', f10.2, ' keV/amu' / &
       ' The FE_TK input parameter must be increased to overcome', &
       ' this problem!'                                                / &
       ' ONETWO will be terminated to avoid further problems.')
!
      return
!
      end subroutine getsgxn
!
!
      subroutine adassgxn (namep,namei,mb,mfm1,nprim,nimp,nion,atw, &
                           ebkev,ebfac,ibion,nebin,vbeam,zne,zni,zte, &
                           zti0,zzi,debin,sgxn,sgxnmi,atw_beam)
!
      use param_mod
      implicit none 
!
!      character what_id*45
!      save      what_id
!      data      what_id/"@(#)cray321.f      12 Feb 98  General Atomics"/
!
! ----------------------------------------------------------------------
!
!  ADASSGXN calculates the effective neutral beam stopping cross sections
!  using the JET Atomic Data and Analysis Structure (ADAS). The cross
!  sections are returned in array sgxn for each beam, beam energy
!  component, and FREYA flux zone.
!
!  The plasma ions are assumed to be fully stripped. Only ion species
!  H, He, B, Be, C, O, N, and Ne are available from ADAS. Any other ion
!  will terminate the code. This can be improved later.
!
!  The routine was made to be as compatible as possible with the existing
!  code. A call is made to the original cross section package in order to
!  determine relative deposition fractions only. The HEXNB package is
!  totally avoided.
!
!  Reference:   Finkenthal, Ph.D. Thesis, UC-Berkeley, 1994
!
!  Created  :   06-jul-1995    by  D. Finkenthal
! ----------------------------------------------------------------------
!
!  Input:
!
!     nprim,nimp,nion
!     atw
!     ebkev(mb)      - full energy of mb-th neutral beamline
!                      (keV)
!     ebfac          - factor defining upper bound on energy bin
!                      range,  .gt. 1.0
!     ibion          - index of beam ion species
!                      (e.g. atwb = atw(ibion))
!                      if ibion = -1 beam is dt mixture
!     atw_beam         (use atw_beam for mass in this case)
!     mb             - number of beamlines modeled
!     mfm1           - number of flux zones minus one
!     nebin          - number of energy bins (rotation case only)
!     vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
!                      beamline (cm/sec)
!     zne(mfm1)      - local electron density (cm**-3)
!     zni(mfm1,nion) - local density of nion-th ion species
!                      (cm**-3)
!     zte(mfm1)      - local electron temperature (KeV)
!     zti0(mfm1)     - local ion temperature (KeV)
!     zzi(mfm1,nion) - local average charge state of nion-th ion
!                           species
!     ONETWO stuff:
!     input from common /io/:
!          ncrt,nout
!     input from common /ions/:
!          namei, atw
!     input from common /numbrs/:
!          nprim,nimp,nion
!     END ONETWO stuff
!
!  Output:
!
!      sgxn(i,j,k,l)
!        i - mode index
!            = 1, fraction of reactions producing electrons;
!            = 2, fraction of reactions producing species 1 ion;
!            = 3, fraction of reactions producing species 2 ion;
!            = 4, total inverse mean free path;
!        j - FREYA zone index
!        k - beam index
!            k = 3*(ib-1) + ie, where ib is the beam index and ie
!                is the beam energy group index
!        l - index for relative energy.  bins are equispaced in
!            delta(energy) between 0 and max(ebkev(ib))*ebfac,
!            with bin width given by
!                   delta(energy) = ebmax*ebfac/nebin.
!            ebfac.ge.1.0 and nebin are user given namelist variables.
!      sgxnmi(ie,mb)
!         - minimum inverse mean free path for ie-th energy
!           group of mb-th beamline.  calculated for stationary
!           target case only.  otherwise, calculated for each
!           neutral trajectory in subroutine INJECT.
!      debin
!         - width of energy bins (keV/amu)
!           Used only for rotatation case.
! ----------------------------------------------------------------------
!
!       include 'params.inc'

!
      character*8 namei(kimp), namep(kprim)   ! added for cql3d
      integer mb, mfm1, nion, nprim
      real(c_double) :: ebkev(kb),vbeam(ke,kb)
      real(c_double) :: zne(kz),zni(kz,kion),zte(kz),zti0(kz),zzi(kz,kion)
      real(c_double) :: sgxn(kcmp1,kz,kbe,ksge),sgxnmi(ke,kb)  ! stand alone code
!      real(c_double) :: sgxn(kcmp1,kz,ke,kb),sgxnmi(ke,kb)
!
      real(c_double) :: ebfac, debin, atw_beam, atw(kion)
      real(c_double) :: sgxeff(3), cnz(kz,20), zeffx(kz), izatom(kion)
      real(c_double) :: atwb,tiev,ebin,ebmax,ecol,sgold,sgxnscale
      integer nimp,ibion,nebin !input
      integer i,j,k,ib,ie,init,ierr,izbm,nout,jreff,ncrt,ind
      save init, izatom, atwb, izbm
      data init /0/
!
      do i=1, kz
        do k = 1, 20
          cnz(i,k) = 0.D0
        enddo
      enddo
!
      do i=1,kz
       do j=1,ke
         do k=1,kb
          sgxn(1,i,j,k)=0.D0
          sgxn(2,i,j,k)=0.D0
          sgxn(3,i,j,k)=0.D0
          sgxn(4,i,j,k)=0.D0
        enddo
       enddo
      enddo
!
      do i=1,3
       sgxeff(i)=0.D0
      enddo
!
      if (init .eq. 0) then
!
! ----------------------------------------------------------------------
! determine atomic number of primary ions
! ----------------------------------------------------------------------
!
      do i=1,nprim
        izatom(i) = 1 ! for H,D,T
        if (trim(namep(i)).eq.'he'   .or. &
            trim(namep(i)).eq.'HE'   .or. &
            trim(namep(i)).eq.'He') izatom(i) = 2
      enddo
!
! ----------------------------------------------------------------------
! determine atomic number of impurity ions
! ----------------------------------------------------------------------
!
      if (nimp .eq. 0)  go to 3430
      do i=1,nimp
        k = nprim + i
        izatom(k) = 0
        if(trim(namei(i)).eq.'he' .or. &
           trim(namei(i)).eq.'HE' .or. &
           trim(namei(i)).eq.'He')  izatom(k) =  2

        if(trim(namei(i)).eq.'b' .or. &
           trim(namei(i)).eq.'B' )  izatom(k) =  5  ! YuP added [2015]

        if(trim(namei(i)).eq.'c' .or. &
           trim(namei(i)).eq.'C' )  izatom(k) =  6

        if(trim(namei(i)).eq.'o' .or. &
           trim(namei(i)).eq.'O' )  izatom(k) =  8

        if(trim(namei(i)).eq.'si' .or. &
           trim(namei(i)).eq.'SI' .or. &
           trim(namei(i)).eq.'Si')  izatom(k) = 14

        if(trim(namei(i)).eq.'ar' .or. &
           trim(namei(i)).eq.'AR' .or. &
           trim(namei(i)).eq.'Ar')  izatom(k) = 18

        if(trim(namei(i)).eq.'cr' .or. &
           trim(namei(i)).eq.'CR' .or. &
           trim(namei(i)).eq.'Cr')  izatom(k) = 24

        if(trim(namei(i)).eq.'fe' .or. &
           trim(namei(i)).eq.'FE' .or. &
           trim(namei(i)).eq.'Fe')  izatom(k) = 26

        if(trim(namei(i)).eq.'ni' .or. &
           trim(namei(i)).eq.'NI' .or. &
           trim(namei(i)).eq.'Ni')  izatom(k) = 28

        if(trim(namei(i)).eq.'kr' .or. &
           trim(namei(i)).eq.'KR' .or. &
           trim(namei(i)).eq.'Kr')  izatom(k) = 36

        if(trim(namei(i)).eq.'mo' .or. &
           trim(namei(i)).eq.'MO' .or. &
           trim(namei(i)).eq.'Mo')  izatom(k) = 42

        if(trim(namei(i)).eq.'w' .or. &
           trim(namei(i)).eq.'W' )  izatom(k) = 74
      enddo
!
! Get beam species
!
 3430 if (ibion .gt. 0) then
        atwb = atw(ibion)
        izbm = izatom(ibion)
      else
        atwb = atw_beam
        izbm = 1
      endif
!
! ----------------------------------------------------------------------
!... Check to make sure that the impurities requested are compatable
!    with ADAS (i.e., H, He, B, Be, C, O, N, or Ne). Terminate with
!    error if an invalid impurity is listed.
! ----------------------------------------------------------------------
!
      do k=1,nion
        if (izatom(k) .gt. 10) then
          write (ncrt, 1010)
          write (nout, 1010)
          call STOP ('subroutine ADASSGXN: unallowed impurity', 181)
        end if
 1010 format &
         (' *** Execution Terminated:'                                 / &
          '     The ADAS database only contains atomic cross-sections' / &
          '     for fully-stripped H, He, B, Be, C, O, N, and Ne ions.'/ &
          '     You must restrict your choice of impurity species to'  / &
          '     these ions.')
      end do
!
      init = 1
      endif    ! init
!
! ----------------------------------------------------------------------
!... Get original cross sections-Used to determine relative deposition
!    fractions only. Total cross section is determined using ADAS.
! ----------------------------------------------------------------------
!
      call nbsgold (mb,mfm1,nebin,nprim,nimp,nion,atw,ebkev, &
                    ebfac,ibion,vbeam,zne,zni,zte,zzi, &
                    debin,sgxn,sgxnmi,atw_beam)
!
! ----------------------------------------------------------------------
!... Set up the cnz arrays (concentrations of plasma ions) and Zeffx
!    zni and zne are the (FREYA zone) densities of electron and ions:
! ----------------------------------------------------------------------
!
      do i=1, mfm1
        do k = 1, nion
          cnz(i,k) = 0.D0
          if (izatom(k) .eq.  1)  cnz(i, 1) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  2)  cnz(i, 2) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  4)  cnz(i, 4) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  5)  cnz(i, 5) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  6)  cnz(i, 6) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  7)  cnz(i, 7) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  8)  cnz(i, 8) = zni(i,k)/zne(i)
          if (izatom(k) .eq. 10)  cnz(i,10) = zni(i,k)/zne(i)
!          write(*,*) 'cnz = ',i,k,cnz(i,k)
        enddo
!
! The Zeff(i=1,...mfm1) array has been stored in zzi(i,nion+1)
        zeffx(i) = zzi(i,nion+1)
!
      enddo
!
! ----------------------------------------------------------------------
!     stationary plasma case
! ----------------------------------------------------------------------
!
      if (nebin .eq. 0) then
        ierr = 0
        do i=1,mfm1
          do ib=1,mb
            tiev = 1.D3 * zti0(i)
            ecol = 1.D3 * ebkev(ib) / (atwb)
            call adasqh6 (tiev,ecol,izbm,0,zne(i),zeffx(i), &
                          cnz(i,2),cnz(i,4),cnz(i,5),cnz(i,6), &
                          cnz(i,7),cnz(i,8),cnz(i,10), &
                          sgxeff,ierr)
            if (ierr .eq. 1) then
              call STOP ('subroutine ADASSGXN: problem #1', 183)
            end if
            do ie = 1,ke
              ind  = ke*(ib-1) + ie
              sgxn(4,i,ind,1) = zne(i)*sgxeff(ie)
              sgxnmi(ie,ib)   = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
            enddo
          enddo
        enddo
!
      else
!
! ----------------------------------------------------------------------
!     rotating plasma case
! ----------------------------------------------------------------------
!
        ebmax = 0.D0
        do ib=1,mb
          ebmax = MAX (ebmax,ebkev(ib))
        end do
        ebmax = ebmax/atwb
        debin = ebmax*ebfac/FLOAT (nebin)
!
        do i=1,mfm1
          tiev  = 1.D3*zti0(i)
          jreff = 0
          do j=1,nebin
            ebin = FLOAT (j) * debin
            ecol =     1.D3 *  ebin
!
!... Only call ADAS if beam energy (ecol) is greater than 5 keV/amu.
!    Skip over the energy bins that are less than 5 keV/amu. These
!    low energy bins will be calculated using the old cross sections
!    scaled to the adas data at some reference energy.
!
            if (ecol .ge. 5000.0) then
              call adasqh6 (tiev,ecol,izbm,1,zne(i),zeffx(i), &
                            cnz(i,2),cnz(i,4),cnz(i,5),cnz(i,6), &
                            cnz(i,7),cnz(i,8),cnz(i,10), &
                            sgxeff,ierr)
              if (ierr .eq. 1) then
                call STOP ('subroutine ADASSGXN: problem #2', 184)
              endif
!
!... Calculate a scaling factor from the first available (ge. 5.0)
!    ADAS bin energy This will be used to scale the old cross
!    sections to fit in with ADAS for the low energy bins.
!
              if (jreff .eq. 0) then
                jreff     = j
                sgxnscale = zne(i)*sgxeff(1)/sgxn(4,i,1,j)
              endif
              do k=1,ke*mb
                sgxn(4,i,k,j) = zne(i)*sgxeff(1)
              enddo
!
            endif  ! if ecol
!
          enddo ! j-loop
!
!... Now scale the old cross sections to match the first available ADAS
!    datapoint for the energy bins below 5.0 keV/amu
!
          if (jreff .gt. 1) then
            do  j=1,jreff-1
             do k=1,ke*mb
               sgold=sgxn(4,i,k,j)
               sgxn(4,i,k,j) = sgxn(4,i,k,j)*sgxnscale
             enddo
            enddo
          endif
!
        enddo  ! i-loop
      endif  ! if nebin
!
!      do ib=1,10
!       do j=1,kbe
!         write(*,'(2i3,2x,1p4e12.4,a12)') ib,j,sgxn(4,ib,j,1),
!     &        sgxn(4,ib,j,2),sgxn(4,ib,j,3),sgxn(4,ib,j,4),' sgxn'
!       enddo
!      enddo
!
      return
      end subroutine adassgxn
!
      !

            real(c_double) function seval (n, u, x, y, b, c, d)
!
      implicit none 
!.......................................................................
!      include 'fileiou.inc'
!      BH: explicitly including:
!-----Include file, FILEIOU.INC, 03/01/92
!-----This file includes the unit numbers of the different I/O files
!-----for the NFREYA project, as well as the title:

      integer ntty, nin, nout, nmhd, ngrf
      common/iou/ntty, nin, nout, nmhd, ngrf ! only here 
      character*80 title
      common/titl/title ! only here
!.......................................................................
!
      integer n
      real(c_double) :: u, x(n), y(n), b(n), c(n), d(n)
!
! this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) < u < x(i+1), using Horner's rule
!
!  if  u .lt. x(1)  then  i = 1  is used.
!  if  u .ge. x(n)  then  i = n  is used.

!  ---------------------------HSJ Modification 1/30/98------------------
!  The above two rules are not acceptible in the context in which
!  this routine is used in the ONETWO code. (actually they are not
!  acceptable in any context in my opinion.
!  Interpolating outside the cubic spline fit is totally bogus.)
!  I have changed the rules to read
!   if u .gt. x(n) then use Y(n) as the most reasonable
!   approximation if you dont actually want to stop the code.
!   similarly if u .lt. x(1) then use y(1) as the returned value.
!  using the first or last point may not be much better in some circumstances
!  but at least we will get physically relevant values !!!!!!!!!!!HSJ
!  Note that the switch extend_seval must be set to 1 to get this
!  behavior. Otherwise an error stop is taken when trying to interpoalte
!  outside the tables.
! -----------------------------------------------------------------------
!
!  input..
!
!    n     = the number of data points
!    u     = the abscissa at which the spline is to be evaluated
!    x,y   = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
      integer i, j, k
      real(c_double) ::  dx
      integer extend_seval
      data    i/1/, extend_seval/1/ ! YuP[07-2016] changed extend_seval from 0 to 1
      !YuP: Had problem with extend_seval=0 in mirror machine runs:
      ! Code quits with a message:
      !Subroutine Seval,a spline fit evaluater,
      !has detected that out of bounds interpolation
      !is occuring. The value at which the cubic spline
      !is to be evaluated is     4.16666667E+03
      !the interpolating table extends from     5.00000000E+03 to     1.00000000E+06
      !You can set extend_seval=1 in the first namelist
      !in inone if you want to ignore this problem
      !subroutine SEVAL: Bad Interpolation
!
!
!  added  1/29/98  HSJ
!
      if (u .lt. x(1)) then
         if(extend_seval .eq. 1)then
           seval=y(1)
           return
         else
           write(*,100)u,x(1),x(n)
            write(nout,100)u,x(1),x(n)
!           call STOP('subroutine SEVAL: Bad Interpolation', 1226)
           write(nout,*)'subroutine SEVAL: Bad Interpolation'
           STOP
         endif
      else if(u .gt. x(n)) then
         if(extend_seval .eq. 1)then
           seval=y(n)
           return
         else
           write(*,100)u,x(1),x(n)
            write(nout,100)u,x(1),x(n)
!           call STOP('subroutine SEVAL: Bad Interpolation', 1226)
           write(nout,*)'subroutine SEVAL: Bad Interpolation'
           STOP
         endif
      endif
!
      if (i .ge. n)  i = 1
      if (u .lt. x(i))  go to 10
      if (u .le. x(i+1))  go to 30
!
!  binary search
!

  10  i = 1
      j = n+1
  20  k = (i+j)/2
      if (u .lt. x(k))  j = k
      if (u .ge. x(k))  i = k
      if (j .gt. i+1 )  go to 20
!
!  evaluate spline
!
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
!
 100  format('Subroutine Seval,a spline fit evaluater,',/, &
             ' has detected that out of bounds interpolation',/, &
             ' is occuring. The value at which the cubic spline',/, &
             ' is to be evaluated is ',1pe18.8,/, &
             ' the interpolating table extends from ',1pe18.8, &
             ' to ',1pe18.8,/, &
             ' You can set extend_seval=1 in the first namelist',/, &
             ' in inone if you want to ignore this problem')
!
    end function seval
    
      subroutine adasqh6 (ti,ecol,bmz,iflag,ne,zeff,conche, &
                          concbe,concb,concc,concn,conco,concne, &
                          qrat,ierr)
      use param_mod
!
!  Calculate and return the rate coefficient for beam stopping
!  This routine is a modified version of L.D. Horton's (JET)
!  QHIOCH6 routine used to access and evaluate the ADAS beam
!  stopping cross sections from ion-specific effective data files.
!
!   changes from QHIOCH6:
!      1) The need for the NAG library has been eliminated using
!         the spline routines SPLINE and SEVAL.
!      2) Only the full-energy component can be calculated by
!         setting iflag=1. This allows lower energies to be
!         read from the ion-specific files.
!
!                                  D.F. Finkenthal 6-JUL-95
!
!   changes from QHIOCH5:
!      1) modified to force rereading of input files when beam
!         species has changed from the last call (using ipass)
!
!                                  L.D. Horton    2 June 92
!
!   changes from QHIOCH4:
!      1) modified to accept beam species as input and to then
!         calculate the beam stopping for either hydrogen or
!         helium beams
!
!                                  L.D. Horton   16 August 91
!
!   changes from QHIOCH3:
!      1) incorporated Fritsch's new calculation for excitation by
!         Be4+ and He2+
!      2) improved calculation of Maxwellian averages in bundled-n
!         code
!      3) proper inclusion of beam energy dependence
!          - stopping rate is now read from a matrix on beam energy
!            and target density; only the temperature dependence is
!            done separately
!      4) included boron, nitrogen and neon as input concentrations
!          - the code skips all zero concentrations
!
!                                  L.D. Horton   31 July 91
!
!   also:  back to Wilhelm's 3 energies-at-a-time so that all spline
! fits can be done at once. Only if the beam energy has changed by
! more than 1% will the density splines be redone.  Since only one
! energy is used per bank, this means that the spline will be done
! only twice per call of ATTS4.
!
!                                  L.D. Horton   14 August 91
!
! ti        : REAL   : ion temperature in eV
! ecol      : REAL   : collision (=beam) energy in eV/amu
! bmz       : INTEGER: beam nuclear charge
! iflag     : INTEGER: flag for full energy calculation only
! ne        : REAL   : electron density in m**-3
! conche    : REAL   : relative Helium concentration   ( 0 < CHE < 1)
! concbe    : REAL   : relative Beryllium concentration
! concb     : REAL   : relative Boron concentration
! concc     : REAL   : relative Carbon concentration
! concn     : REAL   : relative Nitrogen concentration
! conco     : REAL   : relative Oxygen concentration
! concne    : REAL   : relative Neon concentration
! qrat(3)   : REAL   : requested cross section m**2 for full, half,
!                      and third energy beam components
!
!
!      USE param
!      USE io
!      USE ext_prog_info, only : nchars_12,onetwo_xsct
!
      implicit none 
!      include 'params.inc'
!      include 'io.i'

!
! ipass : file read switch.  Reread files if beam species has changed
!
      integer    ipass, nchars_12
      character*256 code_xsct
      character  dsn(8,2)*35
      data       ipass/0/
      data       dsn  / 'data/adas/h_h1.dat'  , &
                        'data/adas/h_he2.dat' , &
                        'data/adas/h_be4.dat' , &
                        'data/adas/h_b5.dat'  , &
                        'data/adas/h_c6.dat'  , &
                        'data/adas/h_n7.dat'  , &
                        'data/adas/h_o8.dat'  , &
                        'data/adas/h_ne10.dat', &
                        'data/adas/he_h1.dat' , &
                        'data/adas/he_he2.dat', &
                        'data/adas/he_be4.dat', &
                        'data/adas/he_b5.dat' , &
                        'data/adas/he_c6.dat' , &
                        'data/adas/he_n7.dat' , &
                        'data/adas/he_o8.dat' , &
                        'data/adas/he_ne10.dat'/
!
! Physics Constants
!
      real(c_double) :: amu, eV
      parameter (amu = 1.6605e-24)
      parameter (eV  = 1.6022e-12)
!
! Local variables
!
!      integer    bmz, iflag, nebeam
      integer    bmz, iflag
      integer,save :: nebeam
      real(c_double) ::     ti, ecol, ne, qrat(3)
      real(c_double) ::     conche, concbe, concb, concc, concn, conco, concne
!****  real(c_double) ::     conch
      integer    nsp, maxe, maxn, maxt, ierr
      parameter (nsp = 8)              ! 8 different ion species
      parameter (maxe = 15, maxn = 10, maxt = 10)
      integer    isp, z(nsp), neb(nsp), ndens(nsp), ntemp(nsp), i, j, k
      real(c_double) ::     eb(maxe,nsp), dens(maxn,nsp), temp(maxt,nsp)
      real(c_double) ::     tref(nsp), ebref(nsp), denref(nsp), svref(nsp)
      real(c_double) ::     sven(maxe,maxn,nsp), svt(maxt,nsp)
      character  line*80
      data       z /1, 2, 4, 5, 6, 7, 8, 10/
!
!      real(c_double) ::     ti8, ecol8, ne8, qrat8(3)
      real(c_double) ::     ti8, ne8, qrat8(3)
      real(c_double),    save  :: ecol8
      real(c_double) ::     conc(nsp)
!
      integer ifail
      real(c_double) :: be(maxe,maxn,nsp),ce(maxe,maxn,nsp),de(maxe,maxn,nsp)
      real(c_double) :: bt(maxt,nsp),ct(maxt,nsp),dt(maxt,nsp)
      real(c_double) :: bn(maxn,nsp,3),cn(maxn,nsp,3),dn(maxn,nsp,3)
      real(c_double) :: svintn(maxn,nsp,3),svintt,svtot(nsp),svtcor(nsp)
!
      real(c_double) :: zeffm1, vbeam, zeff8
      real(c_double) :: zeff !input
      integer ncrt,nout,nunadas,ierror !local
      save z, svref, neb, ndens, tref, eb, dens, sven, ntemp
      save ebref, denref, temp, svt
!

!BH171118      code_xsct is the name of directory where ADAS directory
!BH171118      data/ is to be found. We first choose to either put this
!BH171118      directory in the local execution directory, or put a
!BH171118      link to it therein. A copy of adas_dir data comes with
!BH171118      the code distribution.
!BH171118      Alternative, cql3d will search /usr/local/ for adas_dir.
      code_xsct="./adas_dir/"

      ierr    = 0
      qrat(1) = 0.0
      qrat(2) = 0.0
      qrat(3) = 0.0
      ti8     = ti
      zeff8   = zeff
      ne8     = ne * 1.0e-13
!****  conc(1) = conch
      conc(2) = conche
      conc(3) = concbe
      conc(4) = concb
      conc(5) = concc
      conc(6) = concn
      conc(7) = conco
      conc(8) = concne
      nchars_12 = LENGTH(code_xsct)
!
! open and read input file only once
! ** need to save data
!
      if (ipass .ne. bmz) then
        write (ncrt, 1100)
        write (nout, 1100)
!
        ecol8 = ecol
!
        do isp=1,nsp
          call getioun(nunadas,nunadas)
          print *,'file =',code_xsct(1:nchars_12)//dsn(isp,bmz) !jmp.den
          open (unit = nunadas, &
             file = code_xsct(1:nchars_12)//dsn(isp,bmz),status='OLD', &
             action='read',iostat=ierror)
          if (ierror.ne.0) then
             call getioun(nunadas,nunadas)
             print *,'file = ','/usr/local/adas_dir/'//dsn(isp,bmz) !jmp.den
             open (unit = nunadas, &
             file='/usr/local/adas_dir/'//dsn(isp,bmz),status='OLD', &
             action='read',iostat=ierror)
             if (ierror.ne.0) then
                WRITE(*,*)
                WRITE(*,*) 'Error finding adas_dir: It is neither in'
                WRITE(*,*) 'the execution directory, or in /usr/local .'
                WRITE(*,*) 'Refer to zfreya.f for more information'
                stop
             endif
          endif
          read (nunadas,  1000) z(isp),svref(isp)
          read (nunadas, '(a)') line
          read (nunadas,  1001) neb(isp),ndens(isp),tref(isp)
          read (nunadas, '(a)') line
          read (nunadas,  1002) (eb(j,isp),j=1,neb(isp))
          read (nunadas,  1002) (dens(k,isp),k=1,ndens(isp))
          read (nunadas, '(a)') line
          do k=1,ndens(isp)
            read (nunadas, 1002) (sven(j,k,isp),j=1,neb(isp))
          end do
          read  (nunadas, '(a)') line
          read  (nunadas,  1003) ntemp(isp),ebref(isp),denref(isp)
          read  (nunadas, '(a)') line
          read  (nunadas,  1002) (temp(j,isp),j=1,ntemp(isp))
          read  (nunadas, '(a)') line
          read  (nunadas,  1002) (svt(j,isp),j=1,ntemp(isp))
          call giveupus_cql3d(nunadas)
          close (nunadas)
!
          do k = 1, ndens(isp)
            dens(k,isp) = dens(k,isp) * 1.0e-13
          enddo
!
! spline the data on energy and temperature
!
          do k=1,ndens(isp)
            ifail = 0
            call spline_12 (neb(isp),eb(1,isp),sven(1,k,isp), &
                         be(1,k,isp),ce(1,k,isp),de(1,k,isp))
            if (ifail .ne. 0) then
              write (6, '(a)')  ' spline error #1 in ADASQH6'
              ierr = 1
              return
            endif
          enddo
!
          ifail = 0
          call spline_12 (ntemp(isp),temp(1,isp),svt(1,isp), &
                       bt(1,isp),ct(1,isp),dt(1,isp))
          if (ifail .ne. 0) then
            write (6, '(a)')  ' spline error #2 in ADASQH6'
            ierr = 1
            return
          endif
!
! Determine if only the full energy component is requested
! Default is that all three beam energy components are to be determined
! Only the full energy component is interesting for He beams
!
          nebeam = 3
          if (iflag .gt. 0)  nebeam = 1
          if (  bmz .eq. 2)  nebeam = 1
!
! spline on density for each requested energy component (1 or all 3)
!
          do i=1,nebeam
            do k=1,ndens(isp)
              svintn(k,isp,i) = seval(neb(isp),ecol8/i,eb(1,isp), &
                   sven(1,k,isp),be(1,k,isp),ce(1,k,isp),de(1,k,isp))
            enddo
!
            ifail = 0
            call spline_12(ndens(isp),dens(1,isp),svintn(1,isp,i), &
                        bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
!
            if (ifail .ne. 0) then
              write (6, '(a)')  ' spline error #3 in ADASQH6'
              ierr = 1
              return
            endif
          enddo
        enddo
!
 1000   format (i5, 8x, d9.3)
 1001   format (2i5, 7x, d9.3)
 1002   format (8(1x, d9.3) / 8(1x, d9.3))
 1003   format (i5, 7x, d9.3, 7x, d9.3)
 1100   format (/ &
       ' Using ADAS, the effective beam stopping cross sections:' / &
       ' ***** Opening and Reading the Atomic Data Tables *****'  /)
!
        ipass = bmz
!
      else
!
!  redo density splines only if beam energy has changed by more than 1%
!
        if (ABS (ecol8-ecol) / ecol8 .gt. 0.01) then
          do isp=1,nsp
            ecol8 = ecol
            do i=1,nebeam
              do k=1,ndens(isp)
                svintn(k,isp,i) = seval(neb(isp),ecol8/i,eb(1,isp), &
                      sven(1,k,isp),be(1,k,isp),ce(1,k,isp),de(1,k,isp))
              enddo
              ifail = 0
              call spline_12 (ndens(isp),dens(1,isp),svintn(1,isp,i), &
                           bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
              if (ifail .ne. 0) then
                write (6, '(a)')  ' spline error #4 in ADASQH6'
                ierr = 1
                return
              endif
            enddo
          enddo
        endif
      endif
!
! calculate correction to requested temperature
!
      do isp=1,nsp
        if (ti8 .le. temp(1,isp)) then
          svtcor(isp) = svt(1,isp)/svref(isp)
        else
          svintt      = seval(ntemp(isp),ti8,temp(1,isp), &
                        svt(1,isp),bt(1,isp),ct(1,isp),dt(1,isp))
          svtcor(isp) = svintt/svref(isp)
        endif
      enddo
!
! scale the input concentrations to match the required zeff
!
           zeffm1 = 0.0
           do isp=2,nsp
             zeffm1 = zeffm1 + z(isp)*(z(isp)-1)*conc(isp)
           enddo
           if (zeffm1 .gt. 1.0e-5) then
             do isp=2,nsp
               conc(isp) = (zeff8-1.0) / zeffm1 * conc(isp)
             enddo
           endif
           conc(1) = 1.0
           do isp=2,nsp
             conc(1) = conc(1) - z(isp)*conc(isp)
           enddo
!
! evaluate at three energy components
!
      do i=1,nebeam
!
! interpolate to requested density
!
        do isp=1,nsp
          if (ne8 .le. dens(1,isp)) then
            svtot(isp) = seval(ndens(isp),dens(1,isp),dens(1,isp), &
                  svintn(1,isp,i),bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
          else
            svtot(isp) = seval(ndens(isp),ne8,dens(1,isp), &
                  svintn(1,isp,i),bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
          endif
        enddo
!
!  construct the total stopping cross section
!
        qrat8(i) = 0.0
        do isp=1,nsp
          qrat8(i) = qrat8(i) + svtot(isp)*svtcor(isp)*z(isp)*conc(isp)
        enddo
!
!  divide by the beam speed to get a cross section and convert to m**2
!
        vbeam   = SQRT (2.0 * ecol8 / i * ev / amu)
        qrat(i) = qrat8(i) / vbeam
      enddo
!
      return
      end subroutine adasqh6
!
!
      subroutine nbsgold (mb,mfm1,nebin,nprim,nimp,nion,atw,ebkev, &
                          ebfac,ibion,vbeam,zne,zni,zte,zzi,debin, &
                          sgxn,sgxnmi,atw_beam)
      use param_mod
!
!
! ----------------------------------------------------------------------
!
!     This subroutine calculates the neutral beam attenuation array sgxn
!     using the old cross sections based on Freeman & Jones (1972).
!     No excited state effects are included. These cross sections have
!     been found to be outdated and are not to used except for purposes
!     of comparison to the newer ADAS data.
!
!     Created:     12-JUL-1995     Daniel Finkenthal
!
!     input:
!
!          ebkev(mb)      - maximum energy of mb-th neutral beamline
!                           (keV)
!          ebfac          - factor defining upper bound on energy bin
!                           range,  > 1.0
!          ibion          - index of beam ion species
!                           (e.g., atwb = atw(ibion))
!          mb             - number of beamlines modeled
!          mfm1           - number of flux zones minus one
!          nebin          - number of energy bins (rotation case only)
!          vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
!                           beamline (cm/sec)
!          zne(mfm1)      - local electron density (cm**-3)
!          zni(mfm1,nion) - local density of nion-th ion species
!                           (cm**-3)
!          zzi(mfm1,nion) - local average charge state of nion-th ion
!                           species
!     input from common /io/:
!          ncrt,nout,nouthx,ncorin
!     input from common /ions/:
!          namei, atw
!     input from common /numbrs/:
!          nprim,nimp,nion
!
!     output:
!          sgxn(i,j,k,l)
!               i - mode index
!                   = 1, fraction of reactions producing electrons;
!                   = 2, fraction of reactions producing species 1 ion;
!                   = 3, fraction of reactions producing species 2 ion;
!                   = 4, total inverse mean free path;
!               j - FREYA zone index
!               k - beam index
!                   k = 3*(ib-1) + ie, where ib is the beam index and ie
!                                      is the beam energy group index
!               l - index for relative energy.  bins are equispaced in
!                   delta(energy) between 0 and max(ebkev(ib))*ebfac,
!                   with bin width given by
!                             delta(energy) = ebmax*ebfac/nebin.
!                   ebfac .ge. 1.0 and nebin are user supplied namelist
!                   variables.
!          sgxnmi(ie,mb)
!                 - minimum inverse mean free path for ie-th energy
!                   group of mb-th beamline.  calculated for stationary
!                   target case only.  otherwise, calculated for each
!                   neutral trajectory in subroutine INJECT.
!          debin
!                 - width of energy bins (keV/amu)
! ----------------------------------------------------------------------
!
!
      implicit none 
!      include 'params.inc'   ! for stand alone code
!
      integer mb, mfm1, ibion, nebin, nprim, nimp, nion
      real(c_double) :: ebfac,debin,atw(kion),ebkev(kb),vbeam(ke,kb)
      real(c_double) :: zne(kz),zni(kz,kion),zte(kz),zzi(kz,kion)
      real(c_double) :: sgxn(kcmp1,kz,kbe,ksge),sgxnmi(ke,kb)  ! stand alone code
      real(c_double) :: sgxne(kbe)
!      real(c_double) :: sgxn(kcmp1,kz,ke,kb),sgxnmi(ke,kb)
!      real(c_double) :: sgxne(ke)
      real(c_double) :: atw_beam !input
      integer i,j,k,ib,ie,ind !local
      real(c_double) :: eova,vbin,teev,ebmax,ebin !local
!
      do ib=1,kb   !YuP[07-2016] Bug: was 1,ke
        do ie=1,ke !YuP[07-2016] Bug: was 1,kb
          sgxnmi(ie,ib) = 0.D0
        enddo
      enddo
!
      do i=1,kz
       do j=1,ke
         do k=1,kb
          sgxn(1,i,j,k)=0.D0
          sgxn(2,i,j,k)=0.D0
          sgxn(3,i,j,k)=0.D0
          sgxn(4,i,j,k)=0.D0
        enddo
       enddo
      enddo
!
! ----------------------------------------------------------------------
! stationary plasma case
! ----------------------------------------------------------------------
!
      if (nebin .eq. 0) then
        do 220 ib=1,mb
        do 220 ie=1,ke
        sgxnmi(ie,ib) = 0.0
        ind  = ke*(ib-1) + ie
        eova = 1.0e3*ebkev(ib)/(FLOAT (ie)*atw_beam)
        vbin = vbeam(ie,ib)
        do 220 i=1,mfm1
        teev = 1.0e3 * zte(i)
        sgxne(ind) = fsgxne(vbin,teev,zne(i))
        sgxn (1,i,ind,1) = 0.0
        do k=1,nion
          sgxn(1,i,ind,1) = sgxn(1,i,ind,1) &
                          + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
          if ((k .le. nprim) .and. (k .le. 2.0)) &
            sgxn(k+1,i,ind,1) = fsgxncx(atw(k),eova,zni(i,k))
        end do
        sgxn(1,i,ind,1) = sgxn(1,i,ind,1) + sgxne(ind)
        sgxn(4,i,ind,1) = sgxn(1,i,ind,1) + sgxn(2,i,ind,1) &
                        + sgxn(3,i,ind,1)
        sgxn(1,i,ind,1) = sgxn(1,i,ind,1)/sgxn(4,i,ind,1)
        sgxn(2,i,ind,1) = sgxn(2,i,ind,1)/sgxn(4,i,ind,1)
        sgxn(3,i,ind,1) = sgxn(3,i,ind,1)/sgxn(4,i,ind,1)
        sgxnmi(ie,ib) = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
  220   continue
        do   ib=1,mb
          do ie=1,ke
            sgxnmi(ie,ib) = 1.0 / sgxnmi(ie,ib)
          end do
        end do
      else
!
! ----------------------------------------------------------------------
! rotating plasma case
! ----------------------------------------------------------------------
!
        ebmax = 0.0
        do ib=1,mb
           ebmax = MAX (ebmax,ebkev(ib))
        end do
        ebmax = ebmax/atw_beam
        debin = ebmax*ebfac/FLOAT (nebin)
!
        do 340 i=1,mfm1
        teev = 1.0e3*zte(i)
!       electron impact ionization independent of rotation speed
        do 320 ib=1,mb
        do 320 ie=1,ke
        ind = ke*(ib-1) + ie
        sgxne(ind) = fsgxne(vbeam(ie,ib),teev,zne(i))
  320   continue
        do 340 j=1,nebin
        ebin = FLOAT (j)*debin
        eova = 1.0e3*ebin
        sgxn(1,i,1,j) = 0.0
        do 330 k=1,nion
        sgxn(1,i,1,j) = sgxn(1,i,1,j) &
                      + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
        if ((k .le. nprim) .and. (k .le. 2.0)) &
          sgxn(k+1,i,1,j) = fsgxncx(atw(k),eova,zni(i,k))
  330   continue
        do 340 k=ind,1,-1       ! (ind=ke*(ib-1) +ie,beam energy index)
          sgxn(1,i,k,j) = sgxn(1,i,1,j) + sgxne(k)
          sgxn(4,i,k,j) = sgxn(1,i,k,j) + sgxn(2,i,1,j) + sgxn(3,i,1,j)
          sgxn(1,i,k,j) = sgxn(1,i,k,j)/sgxn(4,i,k,j)
          sgxn(2,i,k,j) = sgxn(2,i,1,j)/sgxn(4,i,k,j)
          sgxn(3,i,k,j) = sgxn(3,i,1,j)/sgxn(4,i,k,j)
  340   continue
!
      endif
!
      return
      end subroutine nbsgold
!
!
      real(c_double) function AMAXAF(array, first, last)
!
! --- replacement for this LIBMATH function in UNICOS ------------------
!
      implicit none

      character what_id*45
      save      what_id
      data      what_id/"@(#)portlib.f      17 Nov 98  General Atomics"/
!
      integer   first, last, index
      real(c_double) ::    array(last), maximum
!
      maximum = array(first)
!
      do index=first,last
        if (array(index) .gt. maximum)  maximum = array(index)
      end do
!
      AMAXAF = maximum
      return
!
      end function AMAXAF
!
!
      integer function MAXAF (array, first, last)
!
! --- replacement for this LIBMATH function in UNICOS ------------------
!
      implicit none
!
      integer   first, last, index
      integer   array(last), maximum
!
      maximum = array(first)
!
      do index=first,last
        if (array(index) .gt. maximum)  maximum = array(index)
      end do
!
      MAXAF = maximum
      return
!
      end function MAXAF
!
!
      integer function LENGTH (string)
!
! --- return trimmed (no trailing spaces or tabs) length of a string ---
!
      implicit none
!
      character  last*1, string*(*)
!
      LENGTH = LEN (string)
      if (LENGTH .ne. 0) then
        last = string(LENGTH:LENGTH)
        do while (last .eq. ' ' .or. last .eq. '        ')
          LENGTH = LENGTH - 1
          if (LENGTH .ne. 0) then
            last = string(LENGTH:LENGTH)
          else
            last = '#'
          end if
        end do
      end if
      return
!
      end function LENGTH
!
!
      subroutine getioun (unitv, unitreq)
!  getioun sets io-unit variable unitv to the value unitret.
!  Depending on freeus_cql3d, unitret either = input variable unitreq or is output.
      integer unitv, unitreq, unitret
      unitret = unitreq
      call freeus_cql3d(unitret)
      unitv = unitret
      return
      end subroutine getioun
!
      subroutine freeus_cql3d(unitno)
      !YuP[2017] renamed to avoid conflicts with TRANSP
!
!     $Id: dummies.f,v 1.4 2000/05/31 21:10:41 freemanj Exp $
!
! dummy replacement for Basis's freeus, which returns a free unit-
! number in its argument.
      integer unitno
      return
      end subroutine freeus_cql3d
!
      subroutine giveupus_cql3d(unitno)
      !YuP[2017] renamed to avoid conflicts with TRANSP
! dummy replacement for Basis's giveupus, which gives up io unit-no. unitno.
      integer unitno
      return
      end subroutine giveupus_cql3d
!
!
      subroutine STOP (message, status)
!
! --- print specified message and then exit with the specified status --
!
      implicit none
!
      integer   status
      character message*(*)
!
      if (status .eq. 0) then
        write (*, '(/ '' STOP in '', a /)') &
        message
      else
        write (*, '(/ '' STOP in '', a / '' EXIT status is'', i4 /)') &
        message, status
      end if
      call EXIT (status)
!
      end subroutine STOP
!
!
      subroutine spline (n, x, y, b, c, d)
!
      implicit none 
!
      integer n
      real(c_double) ::  x(n), y(n), b(n), c(n), d(n)
!
!  the coefficients b(i), c(i), and d(i), i = 1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n .ge. 2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
      integer nm1, ib, i
      real(c_double) ::  t
!
      nm1 = n - 1
      if (n .lt. 2)  return
      if (n .lt. 3)  go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i=2, nm1
         d(i)   = x(i+1) - x(i)
         b(i)   = 2.D0 * (d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i)   = c(i+1) - c(i)
      end do
!
!  end conditions.
!  third derivatives at x(1) and x(n) obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.D0
      c(n) = 0.D0
      if (n .eq. 3)  go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
   15 do i=2, n
        t    = d(i-1)/b(i-1)
        b(i) = b(i) - t*d(i-1)
        c(i) = c(i) - t*c(i-1)
      end do
!
!  back substitution
!
      c(n) = c(n)/b(n)
      do ib=1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.0*c(n))
      do i=1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.D0*c(i)
      end do
      c(n) = 3.D0*c(n)
      d(n) = d(n-1)
      return
!
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.D0
      d(1) = 0.D0
      b(2) = b(1)
      c(2) = 0.D0
      d(2) = 0.D0
      return
!
      end subroutine spline
      
!

!
!
      subroutine spline_12 (n, x, y, b, c, d)
!
      implicit none 
!
      integer n
      real(c_double) ::  x(n), y(n), b(n), c(n), d(n)
!
!  the coefficients b(i), c(i), and d(i), i = 1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n .ge. 2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
      integer nm1, ib, i
      real(c_double) ::  t
!
      nm1 = n - 1
      if (n .lt. 2)  return
      if (n .lt. 3)  go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i=2, nm1
         d(i)   = x(i+1) - x(i)
         b(i)   = 2.0 * (d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i)   = c(i+1) - c(i)
      end do
!
!  end conditions.
!  third derivatives at x(1) and x(n) obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      if (n .eq. 3)  go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
   15 do i=2, n
        t    = d(i-1)/b(i-1)
        b(i) = b(i) - t*d(i-1)
        c(i) = c(i) - t*c(i-1)
      end do
!
!  back substitution
!
      c(n) = c(n)/b(n)
      do ib=1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.0*c(n))
      do i=1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.0*c(i)
      end do
      c(n) = 3.0*c(n)
      d(n) = d(n-1)
      return

   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.0
      d(1) = 0.0
      b(2) = b(1)
      c(2) = 0.0
      d(2) = 0.0
      return

      end subroutine spline_12


end module zfreya_mod

!=======================================================================
!=======================================================================
      BLOCK DATA d01bbf_data
      use iso_c_binding, only : c_double
!
!     Took the following data statements out of subroutine
!     d01bbf, since the Lahey compiler objected to their
!     presence in d01bbf.
!     Their presense in d01bbf is not ANSI F77.
!     BobH, 010618.
!
      implicit none 
      
      real(c_double) :: wgh,xgh !in /cgh/
      common / cgh / wgh(70), xgh(70)
!
      equivalence (wgh10, wgh(1)), (xgh10, xgh(1))
      equivalence (wgh16, wgh(11)), (xgh16, xgh(11))
      equivalence (wgh20, wgh(27)), (xgh20, xgh(27))
      equivalence (wgh24, wgh(47)), (xgh24, xgh(47))
!
      real(c_double) :: wgh10(10), xgh10(10)
      real(c_double) :: wgh16(16), xgh16(16)
      real(c_double) :: wgh20(20), xgh20(20)
      real(c_double) :: wgh24(24), xgh24(24)
!
      real(c_double) :: wgl,xgl !in /cgl/
      common / cgl / wgl(70), xgl(70)
!
      equivalence (wgl10, wgl(1)), (xgl10, xgl(1))
      equivalence (wgl16, wgl(11)), (xgl16, xgl(11))
      equivalence (wgl20, wgl(27)), (xgl20, xgl(27))
      equivalence (wgl24, wgl(47)), (xgl24, xgl(47))

      real(c_double) :: wgl10(10), xgl10(10)
      real(c_double) :: wgl16(16), xgl16(16)
      real(c_double) :: wgl20(20), xgl20(20)
      real(c_double) :: wgl24(24), xgl24(24)

!
      data wgh10   / &
        7.6404329e-06,  1.3436457e-03,  3.3874394e-02,  2.4013861e-01, &
        6.1086263e-01,  6.1086263e-01,  2.4013861e-01,  3.3874394e-02, &
        1.3436457e-03,  7.6404329e-06/


      data xgh10   / &
        3.4361591e+00,  2.5327317e+00,  1.7566836e+00,  1.0366108e+00, &
        3.4290133e-01, -3.4290133e-01, -1.0366108e+00, -1.7566836e+00, &
        -2.5327317e+00, -3.4361591e+00/


      data wgh16   / &
        2.6548075e-10,  2.3209808e-07,  2.7118601e-05,  9.3228401e-04, &
        1.2880312e-02,  8.3810041e-02,  2.8064746e-01,  5.0792948e-01, &
        5.0792948e-01,  2.8064746e-01,  8.3810041e-02,  1.2880312e-02, &
        9.3228401e-04,  2.7118601e-05,  2.3209808e-07,  2.6548075e-10/


      data xgh16   / &
        4.6887389e+00,  3.8694479e+00,  3.1769992e+00,  2.5462022e+00, &
        1.9517880e+00,  1.3802585e+00,  8.2295145e-01,  2.7348105e-01, &
        -2.7348105e-01, -8.2295145e-01, -1.3802585e+00, -1.9517880e+00, &
        -2.5462022e+00, -3.1769992e+00, -3.8694479e+00, -4.6887389e+00/


      data wgh20   / &
        2.2293936e-13,  4.3993410e-10,  1.0860694e-07,  7.8025565e-06, &
        2.2833864e-04,  3.2437733e-03,  2.4810521e-02,  1.0901721e-01, &
        2.8667551e-01,  4.6224367e-01,  4.6224367e-01,  2.8667551e-01, &
        1.0901721e-01,  2.4810521e-02,  3.2437733e-03,  2.2833864e-04, &
        7.8025565e-06,  1.0860694e-07,  4.3993410e-10,  2.2293936e-13/


      data xgh20   / &
        5.3874809e+00,  4.6036824e+00,  3.9447640e+00,  3.3478546e+00, &
        2.7888061e+00,  2.2549740e+00,  1.7385377e+00,  1.2340762e+00, &
        7.3747373e-01,  2.4534071e-01, -2.4534071e-01, -7.3747373e-01, &
        -1.2340762e+00, -1.7385377e+00, -2.2549740e+00, -2.7888061e+00, &
        -3.3478546e+00, -3.9447640e+00, -4.6036824e+00, -5.3874809e+00/


      data wgh24   / &
        1.6643685e-16,  6.5846202e-13,  3.0462543e-10,  4.0189712e-08, &
        2.1582457e-06,  5.6886916e-05,  8.2369248e-04,  7.0483558e-03, &
        3.7445471e-02,  1.2773962e-01,  2.8617954e-01,  4.2693116e-01, &
        4.2693116e-01,  2.8617954e-01,  1.2773962e-01,  3.7445471e-02, &
        7.0483558e-03,  8.2369248e-04,  5.6886916e-05,  2.1582457e-06, &
        4.0189712e-08,  3.0462543e-10,  6.5846202e-13,  1.6643685e-16/


      data xgh24   / &
        6.0159256e+00,  5.2593829e+00,  4.6256628e+00,  4.0536644e+00, &
        3.5200068e+00,  3.0125461e+00,  2.5238810e+00,  2.0490036e+00, &
        1.5842500e+00,  1.1267608e+00,  6.7417111e-01,  2.2441455e-01, &
        -2.2441455e-01, -6.7417111e-01, -1.1267608e+00, -1.5842500e+00, &
        -2.0490036e+00, -2.5238810e+00, -3.0125461e+00, -3.5200068e+00, &
        -4.0536644e+00, -4.6256628e+00, -5.2593829e+00, -6.0159256e+00/


      data wgl10   / &
        3.0844112e-01,  4.0111993e-01,  2.1806829e-01,  6.2087456e-02, &
        9.5015170e-03,  7.5300839e-04,  2.8259233e-05,  4.2493140e-07, &
        1.8395648e-09,  9.9118272e-13/


      data xgl10   / &
        1.3779347e-01,  7.2945455e-01,  1.8083429e+00,  3.4014337e+00, &
        5.5524961e+00,  8.3301527e+00,  1.1843786e+01,  1.6279258e+01, &
        2.1996586e+01,  2.9920697e+01/




      data wgl16   / &
        2.0615171e-01,  3.3105785e-01,  2.6579578e-01,  1.3629693e-01, &
        4.7328929e-02,  1.1299900e-02,  1.8490709e-03,  2.0427192e-04, &
        1.4844587e-05,  6.8283193e-07,  1.8810248e-08,  2.8623502e-10, &
        2.1270790e-12,  6.2979670e-15,  5.0504737e-18,  4.1614624e-22/


      data xgl16   / &
        8.7649410e-02,  4.6269633e-01,  1.1410578e+00,  2.1292836e+00, &
        3.4370866e+00,  5.0780186e+00,  7.0703385e+00,  9.4383143e+00, &
        1.2214223e+01,  1.5441527e+01,  1.9180157e+01,  2.3515906e+01, &
        2.8578730e+01,  3.4583399e+01,  4.1940453e+01,  5.1701160e+01/


      data wgl20   / &
        1.6874680e-01,  2.9125436e-01,  2.6668610e-01,  1.6600245e-01, &
        7.4826065e-02,  2.4964417e-02,  6.2025508e-03,  1.1449624e-03, &
        1.5574177e-04,  1.5401441e-05,  1.0864864e-06,  5.3301209e-08, &
        1.7579812e-09,  3.7255024e-11,  4.7675293e-13,  3.3728442e-15, &
        1.1550143e-17,  1.5395221e-20,  5.2864427e-24,  1.6564566e-28/


      data xgl20   / &
        7.0539890e-02,  3.7212682e-01,  9.1658210e-01,  1.7073065e+00, &
        2.7491993e+00,  4.0489253e+00,  5.6151750e+00,  7.4590175e+00, &
        9.5943929e+00,  1.2038803e+01,  1.4814293e+01,  1.7948896e+01, &
        2.1478788e+01,  2.5451703e+01,  2.9932555e+01,  3.5013434e+01, &
        4.0833057e+01,  4.7619994e+01,  5.5810796e+01,  6.6524417e+01/


      data wgl24   / &
        1.4281197e-01,  2.5877411e-01,  2.5880671e-01,  1.8332269e-01, &
        9.8166273e-02,  4.0732478e-02,  1.3226019e-02,  3.3693491e-03, &
        6.7216256e-04,  1.0446121e-04,  1.2544722e-05,  1.1513158e-06, &
        7.9608130e-08,  4.0728590e-09,  1.5070082e-10,  3.9177365e-12, &
        6.8941811e-14,  7.8198004e-16,  5.3501888e-18,  2.0105175e-20, &
        3.6057659e-23,  2.4518188e-26,  4.0883016e-30,  5.5753458e-35/


      data xgl24   / &
        5.9019852e-02,  3.1123915e-01,  7.6609691e-01,  1.4255976e+00, &
        2.2925621e+00,  3.3707743e+00,  4.6650837e+00,  6.1815351e+00, &
        7.9275392e+00,  9.9120980e+00,  1.2146103e+01,  1.4642732e+01, &
        1.7417993e+01,  2.0491460e+01,  2.3887330e+01,  2.7635937e+01, &
        3.1776041e+01,  3.6358406e+01,  4.1451720e+01,  4.7153106e+01, &
        5.3608575e+01,  6.1058531e+01,  6.9962240e+01,  8.1498279e+01/


      end BLOCK DATA d01bbf_data
!
!=======================================================================
!=======================================================================
      BLOCK DATA hexnb_data
      use iso_c_binding, only : c_double
      implicit none 
      save

!     Data statement of common block variables moved outside
!     subroutine hexnb, and not ANSI F77 (and Lahey compiler
!     complains).  BobH, 010618.
      real(c_double) :: cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf !in /cpub/ 
      common/cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8, &
        cpuii, cpuiz, cpuplf

      data cpu1/0./, cpu2/0./, cpu3/0./, cpu4/0./, cpu5/0./, &
        cpu6/0./, cpu7/0./, cpu8/0./

      end BLOCK DATA hexnb_data
!
