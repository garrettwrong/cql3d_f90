module ainalloc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast

  !---END USE

!
!

contains

  subroutine ainalloc
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!MPIINSERT_INCLUDE

!dir$ nobounds

!..................................................................
!     One measure of whether the following allocations are
!     successful is whether code is entered, then exited without
!     a code stopping fault (e.g. Segmentation fault). See istat_tot
!     comments below.
!..................................................................
      write(*,*)'ainalloc:  Entering ainalloc'

      ! YuP-101220: Moved allocation of cqlb()-cqlf() to vlh and vlf

!..................................................................
!     Allocate storage of main arrays:
!     The following counters of storage are vestigial, to previously
!     used cray pointer system, and linking pointers to one large
!     array, dum(1:lndum).   Maintain it here, just to keep track
!     of amount of storage used here, and for possible cross-checking.
!..................................................................

      lnyxgrs=(iyp1+1)*(jxp1+1)*ngen*lrors
      lyxgsp2=(iyp1+1)*(jxp1+1)*ngen*(lrors+2)
      lnyxgrz=(iyp1+1)*(jxp1+1)*ngen*setup0%lrz
      lnyxnrs=iy*jx*ngen*lrors
      lnyxnrs2=iy*jx*ngen*lrors*2
      lnyxrs=iy*jx*lrors
      lnyxrz=iy*jx*setup0%lrz
      lnefnorz=nefitera*nonch*setup0%lrz
      lnypx2=(iy+1)*(jxp1+1)*lrors
      lny2xp=(iyp1+1)*(jx+1)*lrors
      lniyp=(iy+1)*lrors
      lniy=iy*lrors
      lniyrmx=iy*setup0%lrzmax
      lniyrsp=iy*(lrors+1)
      lnlzmx=lza*setup0%lrzmax
      lnoncrs=nonch*lrors
      lnoncsx=nonch*lrors
      lnoncmx=nonch*setup0%lrzmax
      lnlfield=lfielda*setup0%lrzmax
      lnincz=incza*setup0%lrzmax
      lninczp=inczpa*setup0%lrzmax
      lnjx=jx*setup0%lrzmax
      lnlzmxp=lz*(mx+1)*setup0%lrzmax
      lniylzlra=iy*lz*setup0%lrzmax
      lniylzlr=iy*lz*setup0%lrz
      lnxngrs=jx*ngen*lrors
      lnxngrs5=jx*ngen*lrors*5
      lnyngmx=iy*ngen*setup0%lrzmax
      lnlzgne=lz*(ngen+1)*negyrga*setup0%lrzmax
      lnylzm2=iy*lz*(mx+1+1)*setup0%lrzmax
      lnylzm=iy*lz*(mx+1)*setup0%lrz
      lnylzm1=iy*lz*(mx+1)*setup0%lrzmax
      lnnonng=nonch*ngen*lrors
      lnnonz=nonch*setup0%lrz
      lnj=jx
      lnj1=jxp1
      lni=iy
      lnjp=jpxy
      lnip=ipxy
      lnxm=jx*(10+2*mx)
      lnijp=iy*(jx+1)
      lnipj=(iy+1)*jx
      lnij=iy*jx
      lnipjp=(iyp1+1)*(jxp1+1)
      lnjf=jfl
      lnjfp=jfl+1
      lnjxjf=jx*jfl
      lniyjxlz=iy*jx*lz
      lni0lzlr=(i0param+1)*lz*setup0%lrz

!     Temporary add (bobh 960724)
      lnipjpxy=ipxy*jpxy

      lnjxng=jx*ngen
      lnjxpng=(jx+1)*ngen
      lnmxpi=(mx+1)*iy
      lnjmxp=jx*(mx+1)
      lnjxngi=jx*ngen*12
      lnrhs=1
!BH070525      if (implct.eq."enabled") lnrhs=iyjx
!BH080425       if (implct.eq."enabled") lnrhs=iyjx*lrza
      if (implct.eq."enabled") lnrhs=iyjx
      if (soln_method.eq.'it3dv' .or. soln_method.eq.'it3drv') &
           lnrhs=iyjx*lrza
      lnjnsl=jx*ngen*nsoa*setup0%lrzmax
      lnjmpn=jx*(msxr+1)*nena*2

!     New pointer use in this code: fully pointered variables, currv_*,
!       pwrrf,pwrrfs
      lnjxnglr=jx*ngen*lrors

      lndum= 2*lyxgsp2+4*lnyxgrs+3*lnyxgrz+ &
        7*lnyxnrs+2*lnyxnrs2+2*lnyxrs+3*lnyxrz+lnefnorz+lnypx2+lny2xp+ &
        4*lniyp+5*lniy+12*lniyrmx+3*lniy+3*lniyrmx+lniyrsp+8*lnlzmx+ &
        4*lnoncrs+lnoncsx+3*lnoncmx+8*lnlfield+2*lnincz+2*lninczp+1*lnjx &
        +2*lnlzmxp+6*lniylzlra+lnxngrs+lnxngrs5+lnyngmx+lnlzgne+lnylzm2+ &
        lnylzm+5*lnylzm1+3*lnnonng+8*lnnonz+59*lnj+4*lnj1+10*lni+ &
        5*lnjp+lnip+lnxm+4*lnijp+4*lnipj+2*lnij+8*lnipjp+ &
        4*lnjxng+lnjxpng+lnmxpi+2*lnjmxp+lnjxngi+lnrhs+lnjnsl+lnjmpn+ &
        2*lnjf+5*lnjfp+4*lnjxjf+3*lniyjxlz+lni0lzlr+2*lnjxnglr


!****************SCChiu   2/8/95
!    1  3*lnoncrs+lnoncsx+2*lnoncmx+8*lnlfield+2*lnincz+2*lninczp+2*lnjx
!..................................................................
!     Just for record keeping (vestigial to cray pointer system)
      if(relativ .eq. "fully") then
        lnjxmx5=jx*(mx+5)
        lndum=lndum+2*lnjxmx5+lnj
      endif
!..................................................................
!..................................................................
!     Just for record keeping (vestigial to cray pointer system)
      if(ndeltarho.ne."disabled" .or. lossmode(1).eq.'simplban' &
         .or. lossmode(1).eq.'simplbn1' &
         .or. lossmode(1).eq.'simplbn2') then
        lnrzt=nr_delta*nz_delta*nt_delta
        lndum=lndum+lniylzlr+2*lnrzt+nr_delta+nz_delta+nt_delta &
              +nr_delta*nz_delta
      endif
!..................................................................
!..................................................................
!     Just for record keeping (vestigial to cray pointer system)
      if(tavg.ne."disabled") then
        lndum=lndum+lnyxgrs
      endif
!..................................................................

!_cray   call hpalloc(dumptr,lndum,ierr,0)
!_pc     Using malloc for unix library libU77.a;  Allocates bytes.
!_pc     dumptr=malloc(8*lndum)
!_pc     if (dumptr.eq.0) stop 'Problem with malloc(lndum)'


!..................................................................
!     Check on the allocate by adding up istat.
!     If not zero, problem has occurred.
!BH100830:  But, this didn't work with gfortran when nraya parameter
!BH100830:  is  made excessively large.  Code gets Segmentation fault
!BH100830:  [presumably in urfalloc.f].  But this indicates that
!BH100830:  checking on istat is not sufficient to catch a problem
!BH100830:  of using too much memory.

!..................................................................
      istat_tot=0
      allocate(f(0:iy+1,0:jx+1,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      f = zero
      allocate(fxsp(0:iy+1,0:jx+1,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      fxsp = zero
      allocate(f_(0:iy+1,0:jx+1,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      f_ = zero

      allocate(spasou(0:iy+1,0:jx+1,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      spasou = zero
      allocate(velsou(0:iy+1,0:jx+1,ngen,0:lrors+1),STAT=istat)
      istat_tot=istat_tot+istat
      velsou = zero
      allocate(velsou2(0:iy+1,0:jx+1,ngen,0:lrors+1),STAT=istat)
      istat_tot=istat_tot+istat
      velsou2 = zero

      allocate(source(0:iy+1,0:jx+1,ngen,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      source = zero

      allocate(gone(0:iy+1,0:jx+1,ngen,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      gone = zero

      allocate(egylosa(0:iy+1,0:jx+1,ngen,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(egylosa,zero,SIZE(egylosa))

      allocate(i0tran(i0param+1,lz,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(i0tran,0,SIZE(i0tran))

      allocate(cal(iy,jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cal,zero,SIZE(cal))
      allocate(cbl(iy,jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cbl,zero,SIZE(cbl))
      allocate(ccl(iy,jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ccl,zero,SIZE(ccl))
      allocate(cdl(iy,jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cdl,zero,SIZE(cdl))
      allocate(cel(iy,jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cel,zero,SIZE(cel))
      allocate(cfl(iy,jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cfl,zero,SIZE(cfl))
      allocate(eal(iy,jx,ngen,2,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(eal,zero,SIZE(eal))
      allocate(ebl(iy,jx,ngen,2,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ebl,zero,SIZE(ebl))

      allocate(scal(iyjx*ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(scal,zero,SIZE(scal))
      allocate(cet(iy,jx,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cet,zero,SIZE(cet))
      allocate(cex(iy,jx,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cex,zero,SIZE(cex))
      allocate(synca(iy,jx,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(synca,zero,SIZE(synca))
      allocate(syncd(iy,jx,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(syncd,zero,SIZE(syncd))
      allocate(taulos(iy,jx,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(taulos,zero,SIZE(taulos))
      allocate(psi0bar(setup0%lrz,0:nstop),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psi0bar,one,SIZE(psi0bar))
      allocate(delecfld0(1:setup0%lrz,1:nstop),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(delecfld0,zero,SIZE(delecfld0))
      allocate(elecfldn(0:setup0%lrz+1,0:nstop,0:nampfmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(elecfldn,zero,SIZE(elecfldn))
      allocate(delecfld0n(1:setup0%lrz,1:nstop,1:nampfmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(delecfld0n,zero,SIZE(delecfld0n))
!     nefitera=10, presently.
!BH171231      allocate(elecn(1:setup0%lrz,0:nstop,nefitera),STAT=istat)
!BH171231      Fixing write of elecn in tdoutput elecn(,nch(1),)
      allocate(elecn(1:setup0%lrz,0:nstop+1,nefitera),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(elecn,zero,SIZE(elecn))

      allocate(di(0:iy,0:jxp1,1:ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(di,zero,SIZE(di))  ! Maybe set to 0.5 ?
      allocate(dj(0:iyp1,0:jx,1:ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dj,zero,SIZE(dj))  ! Maybe set to 0.5 ?
!BH090810:  Think dyp5,dym5,eyp5, and eym5 are over-dimensioned
!BH090810:  from i=0, when only need i=1.
!HB090826:  NO, at least for dyp5, through hfi(i-1,j,k,l_) in impchk, i=1.
!HB090826:    Although multiplied by df(i,j) which is 0. for i=0,
!BH090826:    still the dyp5(0,j) needs to exist.
      allocate(dym5(1:iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dym5,zero,SIZE(dym5))
      allocate(dyp5(0:iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dyp5,zero,SIZE(dyp5))
      allocate(eym5(1:iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(eym5,zero,SIZE(eym5))
      allocate(eyp5(0:iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(eyp5,zero,SIZE(eyp5))
      allocate(y(iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(y,zero,SIZE(y))
      allocate(dy(iy,lrors),dyi(iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dy,zero,SIZE(dy))
      call bcast(dyi,zero,SIZE(dyi))
      allocate(yptb(iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(yptb,zero,SIZE(yptb))
      allocate(coss(iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(coss,zero,SIZE(coss))
      allocate(cynt2(iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cynt2,zero,SIZE(cynt2))
      allocate(batot(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(batot,zero,SIZE(batot))
      allocate(lmax(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(lmax,0,SIZE(lmax))
      allocate(vpint(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(vpint,zero,SIZE(vpint))
      allocate(psiiv(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psiiv,zero,SIZE(psiiv))
      allocate(psiba(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psiba,zero,SIZE(psiba))
      allocate(psisq(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psisq,zero,SIZE(psisq))
      allocate(psicu(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psicu,zero,SIZE(psicu))
      allocate(psiqu(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psiqu,zero,SIZE(psiqu))
      allocate(bavpd(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(bavpd,zero,SIZE(bavpd))
      allocate(bavdn(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(bavdn,zero,SIZE(bavdn))
      allocate(psiir(iy,setup0%lrzmax),STAT=istat)
      call bcast(psiir,zero,SIZE(psiir))
      istat_tot=istat_tot+istat
      allocate(vderb(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(vderb,zero,SIZE(vderb))
      allocate(sinn(iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sinn,zero,SIZE(sinn))
      allocate(tann(iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tann,zero,SIZE(tann))
      allocate(ymid(iy,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ymid,zero,SIZE(ymid))
      allocate(tau(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tau,zero,SIZE(tau))
      allocate(vptb(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(vptb,zero,SIZE(vptb))
      allocate(zboun(iy,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(zboun,zero,SIZE(zboun))
      allocate(idx(iy,0:lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(idx,0,SIZE(idx))
      !---- (lza,setup0%lrzmax) arrays (at pol. angle grid) ---!
      allocate(imax(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(imax,0,SIZE(imax))
      allocate(dz(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dz,zero,SIZE(dz))
      allocate(pol(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pol,zero,SIZE(pol))
      allocate(solrz(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(solrz,zero,SIZE(solrz))
      allocate(solzz(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(solzz,zero,SIZE(solzz))
      allocate(thtab(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(thtab,zero,SIZE(thtab))
      allocate(z(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(z,zero,SIZE(z))
      allocate(zmid(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(zmid,zero,SIZE(z))
      allocate(bbpsi(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(bbpsi,zero,SIZE(bbpsi))
      allocate(bpolz(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(bpolz,zero,SIZE(bbpsi))
      allocate(btorz(lza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(btorz,zero,SIZE(bbpsi))

      allocate(consnp(nonch,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(consnp,zero,SIZE(consnp))
      allocate(ptime(nonch,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ptime,zero,SIZE(ptime))
!BH120315: Prevent out of bounds reference:
!BH120315:      allocate(sptzrp(nonch,lrors),STAT=istat)
      allocate(sptzrp(nonch,max(lrors,setup0%lrzmax)),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sptzrp,zero,SIZE(sptzrp))
      allocate(pefld(nonch,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pefld,zero,SIZE(pefld))
      allocate(rovsp(nonch,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(rovsp,zero,SIZE(rovsp))
      allocate(restp(nonch,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(restp,zero,SIZE(restp))
      allocate(restnp(nonch,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(restnp,zero,SIZE(restnp))
      allocate(vpov(nonch,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(vpov,zero,SIZE(vpov))
      allocate(es(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(es,zero,SIZE(es))
      allocate(bpsi(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(bpsi,zero,SIZE(bpsi))
      allocate(d2bpsi(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(d2bpsi,zero,SIZE(d2bpsi))
      allocate(d2solrz(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(d2solrz,zero,SIZE(d2solrz))
      allocate(d2solzz(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(d2solzz,zero,SIZE(d2solzz))
      allocate(d2bpolz(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(d2bpolz,zero,SIZE(d2bpolz))
      allocate(d2btorz(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(d2btorz,zero,SIZE(d2btorz))
      allocate(d2thtpol(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(d2thtpol,zero,SIZE(d2thtpol))
      allocate(d2es(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(d2es,zero,SIZE(d2es))
      allocate(thtpol(lfielda,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(thtpol,zero,SIZE(thtpol))
      allocate(esfi(incza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(esfi,zero,SIZE(esfi))
      allocate(psiesfi(incza,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psiesfi,zero,SIZE(psiesfi))
      allocate(psifi(inczpa,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psifi,zero,SIZE(psifi))
      allocate(espsifi(inczpa,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(espsifi,zero,SIZE(espsifi))
      allocate(soupp(jx,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(soupp,zero,SIZE(soupp))
      allocate(waa(lz,0:mx,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(waa,zero,SIZE(waa))
      allocate(wbb(lz,0:mx,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(wbb,zero,SIZE(wbb))
      allocate(cosz(iy,lz,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cosz,zero,SIZE(cosz))
      allocate(dtau(iy,lz,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dtau,zero,SIZE(dtau))
      allocate(sinz(iy,lz,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sinz,zero,SIZE(sinz))
      allocate(tanz(iy,lz,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tanz,zero,SIZE(tanz))
      allocate(yz(iy,lz,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(yz,zero,SIZE(yz))
      allocate(tot(iy,lz,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tot,zero,SIZE(tot))
      allocate(vflux(jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(vflux,zero,SIZE(vflux))
      allocate(f_aveth(jx,ngen,lrors,5),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(f_aveth,zero,SIZE(f_aveth))
      allocate(sincosba(iy,ngen,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sincosba,zero,SIZE(sincosba))
      allocate(densz(lz,ngen+1,negyrga,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(densz,zero,SIZE(densz))

      mmx=mx
      if (softxry.ne.'disabled') mmx= max(msxr, mmx)
      if (sigmamod.ne.'disabled') mmx= max(mmsv, mmx)
!BH Problem if msxr.lt.mx, or mmsv.lt.mx
!BH      if (softxry.ne.'disabled') then
!BH         mmx= max(msxr, mmx)
!BH      elseif (sigmamod.ne.'disabled') then
!BH         mmx= max(mmsv, mmx)
!BH      else
!BH         mmx=mx
!BH      endif
      mmxp1=mmx+1

      allocate(ss(iy,lz,0:mmxp1,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ss,zero,SIZE(ss))
      allocate(dcofleg(iy,lz,0:mmx,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dcofleg,zero,SIZE(dcofleg))
      allocate(dpcosz(iy,lz,0:mmx,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dpcosz,zero,SIZE(dpcosz))
      allocate(ssy(iy,lz,0:mx,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ssy,zero,SIZE(ssy))
      allocate(ssyy(iy,lz,0:mx,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ssyy,zero,SIZE(ssyy))
      allocate(ssyi(iy,lz,0:mx,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ssyi,zero,SIZE(ssyi))
      allocate(ssyyy(iy,lz,0:mx,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ssyyy,zero,SIZE(ssyyy))

      allocate(pcurr(nonch,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pcurr,zero,SIZE(pcurr))
      allocate(pcurrm(nonch,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pcurrm,zero,SIZE(pcurrm))

      allocate(pdens(nonch,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pdens,zero,SIZE(pdens))
      allocate(pdenm(nonch,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pdenm,zero,SIZE(pdenm))

      allocate(pengy(nonch,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pengy,zero,SIZE(pengy))
      allocate(pengym(nonch,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pengym,zero,SIZE(pengym))

      allocate(pdenra(nonch,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pdenra,zero,SIZE(pdenra))
      allocate(pcurra(nonch,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pcurra,zero,SIZE(pcurra))
      allocate(pfdenra(nonch,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pfdenra,zero,SIZE(pfdenra))
      allocate(pfcurra(nonch,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pfcurra,zero,SIZE(pfcurra))
      allocate(pucrit(nonch,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pucrit,zero,SIZE(pucrit))
      allocate(peoe0(nonch,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(peoe0,zero,SIZE(peoe0))
      allocate(psrc(nonch,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(psrc,zero,SIZE(psrc))
      allocate(peoed(nonch,setup0%lrz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(peoed,zero,SIZE(peoed))

         jfljx=max(jfl,jx) !
      allocate(cint2(jfljx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cint2,zero,SIZE(cint2))
      allocate(dx(jfljx),dxi(jfljx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dx,zero,SIZE(dx))
      call bcast(dxi,zero,SIZE(dxi))

      allocate(ifp(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(ifp,0,SIZE(ifp))
      allocate(sg(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sg,zero,SIZE(sg))
      allocate(sgx(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sgx,zero,SIZE(sgx))
      allocate(sgxx(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sgxx,zero,SIZE(sgxx))
      allocate(sh(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sh,zero,SIZE(sh))
      allocate(shx(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(shx,zero,SIZE(shx))
      allocate(shxx(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(shxx,zero,SIZE(shxx))
      allocate(shxxx(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(shxxx,zero,SIZE(shxxx))
      allocate(tam1(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam1,zero,SIZE(tam1))
      allocate(tam2(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam2,zero,SIZE(tam2))
      allocate(tam3(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam3,zero,SIZE(tam3))
      allocate(tam4(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam4,zero,SIZE(tam4))
      allocate(tam5(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam5,zero,SIZE(tam5))
      allocate(tam6(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam6,zero,SIZE(tam6))
      allocate(tam7(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam7,zero,SIZE(tam7))
      allocate(tam8(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam8,zero,SIZE(tam8))
      allocate(tam9(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam9,zero,SIZE(tam9))
      allocate(tam10(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam10,zero,SIZE(tam10))
      allocate(tam11(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam11,zero,SIZE(tam11))
      allocate(tam12(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam12,zero,SIZE(tam12))
      allocate(tam13(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam13,zero,SIZE(tam13))
      allocate(tam14(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam14,zero,SIZE(tam14))
      allocate(tam15(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam15,zero,SIZE(tam15))
      allocate(tam16(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam16,zero,SIZE(tam16))
      allocate(tam17(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam17,zero,SIZE(tam17))
      allocate(tam18(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam18,zero,SIZE(tam18))
      allocate(tam19(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam19,zero,SIZE(tam19))
      allocate(tam20(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam20,zero,SIZE(tam20))
      allocate(tam21(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam21,zero,SIZE(tam21))
      allocate(tam22(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam22,zero,SIZE(tam22))
      allocate(tam23(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam23,zero,SIZE(tam23))
      allocate(tam24(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam24,zero,SIZE(tam24))
      allocate(tam25(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam25,zero,SIZE(tam25))
      allocate(tam26(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam26,zero,SIZE(tam26))
      allocate(tam27(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam27,zero,SIZE(tam27))
      allocate(tam28(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam28,zero,SIZE(tam28))
      allocate(tam29(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam29,zero,SIZE(tam29))
      allocate(tam30(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tam30,zero,SIZE(tam30))
      allocate(x(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(x,zero,SIZE(x))
      allocate(xmidpt(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xmidpt,zero,SIZE(xmidpt))
      allocate(xi(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xi,zero,SIZE(xi))
      allocate(xsq(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xsq,zero,SIZE(xsq))
      allocate(x3i(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(x3i,zero,SIZE(x3i))
      allocate(x2i(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(x2i,zero,SIZE(x2i))
      allocate(xcu(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xcu,zero,SIZE(xcu))
      allocate(xcenter(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xcenter,zero,SIZE(xcenter))
      allocate(xcensq(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xcensq,zero,SIZE(xcensq))
      allocate(xcent3(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xcent3,zero,jx)

      allocate(uoc(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(uoc,zero,SIZE(uoc))
      allocate(enerkev(jx,ngena),STAT=istat) !YuP[2018-01-08] added 2nd index (k)
      istat_tot=istat_tot+istat
      call bcast(enerkev,zero,SIZE(enerkev))
      allocate(gamma(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(gamma,zero,SIZE(gamma))
      allocate(gamsqr(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(gamsqr,zero,SIZE(gamsqr))
      allocate(gamcub(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(gamcub,zero,SIZE(gamcub))
      allocate(gammi(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(gammi,zero,SIZE(gammi))
      allocate(gamm2i(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(gamm2i,zero,SIZE(gamm2i))
      allocate(gamm1(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(gamm1,zero,SIZE(gamm1))
      allocate(tcsgm1(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tcsgm1,zero,SIZE(tcsgm1))
      allocate(gamefac(jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(gamefac,zero,SIZE(gamefac))
      allocate(ident(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(ident,0,SIZE(ident))
      allocate(temc1(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(temc1,zero,SIZE(temc1))
      allocate(temc2(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(temc2,zero,SIZE(temc2))
      allocate(temc3(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(temc3,zero,SIZE(temc3))
      allocate(temc4(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(temc4,zero,SIZE(temc4))
      allocate(itemc1(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(itemc1,0,SIZE(itemc1))
      allocate(itemc2(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(itemc2,0,SIZE(itemc2))
      allocate(l_lower(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(l_lower,0,SIZE(l_lower))
      allocate(lpt(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(lpt,0,SIZE(lpt))
      allocate(mun(iy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(mun,zero,SIZE(mun)) !real(c_double)
      allocate(fll(jpxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(fll,zero,SIZE(fll))
      allocate(xpar(jpxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xpar,zero,SIZE(xpar))
      allocate(rheads(jpxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(rheads,zero,SIZE(rheads))
      allocate(dfvlle(jpxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dfvlle,zero,SIZE(dfvlle))
      allocate(dfvlli(jpxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dfvlli,zero,SIZE(dfvlli))
      allocate(xperp(ipxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xperp,zero,SIZE(xperp))
      allocate(xl(1:jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xl,zero,SIZE(xl))
      allocate(jmaxxl(1:jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(jmaxxl,0,SIZE(jmaxxl))
      allocate(xlm(0:jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xlm,zero,SIZE(xlm))
      allocate(dxl(0:jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dxl,zero,SIZE(dxl))
      allocate(fl(0:jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(fl,zero,SIZE(fl))
      allocate(fl1(0:jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(fl1,zero,SIZE(fl1))
      allocate(fl2(0:jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(fl2,zero,SIZE(fl2))
      allocate(ppars(jx,jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ppars,zero,SIZE(ppars))
      allocate(pprps(jx,jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pprps,zero,SIZE(pprps))
      allocate(faci(jx,jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(faci,zero,SIZE(faci))
      allocate(pparea(jx,jfl),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pparea,zero,SIZE(pparea))
      allocate(wtfl0(iy,jx,lz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(wtfl0,zero,SIZE(wtfl0))
      allocate(wtflm(iy,jx,lz),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(wtflm,zero,SIZE(wtflm))
      allocate(jflbin(iy,jx,lz),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(jflbin,0,SIZE(jflbin))
      allocate(xm(jx,-5-mx:4+mx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xm,zero,SIZE(xm))
      allocate(dbb(iy,0:jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dbb,zero,SIZE(dbb))
      allocate(dff(0:iy,jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dff,zero,SIZE(dff))
      allocate(cah(iy,jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cah,zero,SIZE(cah))
      allocate(cthta(iy,jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(cthta,zero,SIZE(cthta))

      allocate(gon(0:iy+1,0:jx+1),STAT=istat)
      istat_tot=istat_tot+istat
      gon = zero

      allocate(so(0:iy+1,0:jx+1),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(so,zero,SIZE(so))

      allocate(currv(jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(currv,zero,SIZE(currv))
      allocate(currvs(jx,ngen),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(currvs,zero,SIZE(currvs))
      allocate(pwrrf(jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pwrrf,zero,SIZE(pwrrf))
      allocate(tal(jx,ngen),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tal,zero,SIZE(tal))
      allocate(tbl(jx,ngen),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tbl,zero,SIZE(tbl))
      allocate(tfl(jx,ngen),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tfl,zero,SIZE(tfl))
      allocate(pwrrfs(jx,ngen,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pwrrfs,zero,SIZE(pwrrfs))
      allocate(pleg(0:mmx,iy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(pleg,zero,SIZE(pleg))
      ! feta(:,:) is used both in mx-loops and msxr-loops
      allocate(feta(jx,0:mmx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(feta,zero,SIZE(feta))
      allocate(fetb(jx,0:mmx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(fetb,zero,SIZE(fetb))

      allocate(wflux(jx,ngen,0:11),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(wflux,zero,SIZE(wflux))
      allocate(rhs(lnrhs),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(rhs,zero,SIZE(rhs))
      allocate(sovt(jx,ngen,nsoa,setup0%lrzmax),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sovt,zero,SIZE(sovt))
      allocate(sigsxr(jx,0:msxr,nena,2),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sigsxr,zero,SIZE(sigsxr))
      allocate(item1(iyjx2),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(item1,0,SIZE(item1))
      allocate(item2(iyjx2),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(item2,0,SIZE(item2))
      allocate(item3(iyjx2),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(item3,0,SIZE(item3))
      allocate(item4(iyjx2),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(item4,0,SIZE(item4))
      allocate(item5(iyjx2),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(item5,0,SIZE(item5))
      allocate(item6(iyjx2),STAT=istat)
      istat_tot=istat_tot+istat
      call ibcast(item6,0,SIZE(item6))
      allocate(dxm5(jxp1),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dxm5,zero,SIZE(dxm5))
      allocate(exm5(jxp1),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(exm5,zero,SIZE(exm5))
!BH080910      dxp5=>dxm5(2:) !due different start extent
!BH080910      exp5=>exm5(2:) !due different start extent
      allocate(dxp5(0:jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dxp5,zero,SIZE(dxp5))
      allocate(exp5(0:jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(exp5,zero,SIZE(exp5))


      allocate(ix1(0:mx),ix2(0:mx),ix3(0:mx),ix4(0:mx),STAT=istat)
      allocate(ix5(0:mx),ix6(0:mx),ix7(0:mx),ix8(0:mx),STAT=istat)
      allocate(tom1(0:mmxp1),tom2(0:mmxp1),STAT=istat)
      allocate(tom3(0:mmxp1),tom4(0:mmxp1),STAT=istat)
      allocate(fctrl(0:2*mx+2),STAT=istat)
      allocate(choose(0:2*mx+2,0:mx+2),STAT=istat)
      allocate(cog(0:mx,15),STAT=istat)
      allocate(pm(0:mmxp1,lrors),STAT=istat)

      call ibcast(ix1,0,SIZE(ix1))
      call ibcast(ix2,0,SIZE(ix2))
      call ibcast(ix3,0,SIZE(ix3))
      call ibcast(ix4,0,SIZE(ix4))
      call ibcast(ix5,0,SIZE(ix5))
      call ibcast(ix6,0,SIZE(ix6))
      call ibcast(ix7,0,SIZE(ix7))
      call ibcast(ix8,0,SIZE(ix8))
      call bcast(tom1,zero,SIZE(tom1))
      call bcast(tom2,zero,SIZE(tom2))
      call bcast(tom3,zero,SIZE(tom3))
      call bcast(tom4,zero,SIZE(tom4))
      call bcast(fctrl,zero,SIZE(fctrl))
      call bcast(choose,zero,SIZE(choose))
      call bcast(cog,zero,SIZE(cog))
      call bcast(pm,zero,SIZE(pm))

      allocate(tamt1(2,jx,-2:mx+2,-2:mx+2),STAT=istat) !!!YuP: extended range
      istat_tot=istat_tot+istat
      call bcast(tamt1,zero,SIZE(tamt1))
      allocate(tamt2(2,jx,-2:mx+2,-2:mx+2),STAT=istat) !!!YuP: extended range
      istat_tot=istat_tot+istat
      call bcast(tamt2,zero,SIZE(tamt2))
      allocate(da(iy,0:jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(da,zero,SIZE(da))
      allocate(db(iy,0:jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(db,zero,SIZE(db))
      allocate(dc(iy,0:jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dc,zero,SIZE(dc))
      allocate(dd(0:iy,jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(dd,zero,SIZE(dd))
      allocate(de(0:iy,jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(de,zero,SIZE(de))
      allocate(df(0:iy,jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(df,zero,SIZE(df))

      allocate(ca(iy,jx),STAT=istat) !-YuP: changed from ca=>da
      istat_tot=istat_tot+istat
      call bcast(ca,zero,iyjx)
      allocate(cb(iy,jx),STAT=istat) !-YuP: changed from cb=>db
      istat_tot=istat_tot+istat
      call bcast(cb,zero,iyjx)
      allocate(cc(iy,jx),STAT=istat) !-YuP: changed from cc=>dc
      istat_tot=istat_tot+istat
      call bcast(cc,zero,iyjx)
      allocate(cd(iy,jx),STAT=istat) !-YuP: changed from cd=>dd
      istat_tot=istat_tot+istat
      call bcast(cd,zero,iyjx)
      allocate(ce(iy,jx),STAT=istat) !-YuP: changed from ce=>de
      istat_tot=istat_tot+istat
      call bcast(ce,zero,iyjx)
      allocate(cf(iy,jx),STAT=istat) !-YuP: changed from cf=>df
      istat_tot=istat_tot+istat
      call bcast(cf,zero,iyjx)

      allocate(bqlm(iy,jx),STAT=istat)
      call bcast(bqlm,zero,iyjx)

      iyjx2l=max(iy+2,setup0%lrz)*(jx+2)
      allocate(tem1(iyjx2l),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tem1,zero,SIZE(tem1))
      allocate(tem2(iyjx2l),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tem2,zero,SIZE(tem2))
      allocate(tem3(iyjx2l),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tem3,zero,SIZE(tem3))
      allocate(tem4(iyjx2l),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tem4,zero,SIZE(tem4))
      allocate(tem5(iyjx2l),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tem5,zero,SIZE(tem5))
      allocate(tem6(iyjx2l),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(tem6,zero,SIZE(tem6))

      allocate(egg(iy,jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(egg,zero,SIZE(egg))
      allocate(fgg(iy,jx),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(fgg,zero,SIZE(fgg))

      allocate(xhead(jpxy,ipxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xhead,zero,SIZE(xhead))
      allocate(xtail(jpxy,ipxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xtail,zero,SIZE(xtail))
      allocate(ytail(jpxy,ipxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(ytail,zero,SIZE(ytail))
      allocate(yhead(jpxy,ipxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(yhead,zero,SIZE(yhead))
      allocate(fpn(jpxy,ipxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(fpn,zero,SIZE(fpn))

      allocate(temp1(0:iy+1,0:jx+1),STAT=istat)
      istat_tot=istat_tot+istat
      allocate(temp2(0:iy+1,0:jx+1),STAT=istat)
      istat_tot=istat_tot+istat
      allocate(temp3(0:iy+1,0:jx+1),STAT=istat)
      istat_tot=istat_tot+istat
      allocate(temp4(0:iy+1,0:jx+1),STAT=istat)
      istat_tot=istat_tot+istat
      allocate(temp5(0:iy+1,0:jx+1),STAT=istat)
      istat_tot=istat_tot+istat
      allocate(temp6(0:iy+1,0:jx+1),STAT=istat)
      istat_tot=istat_tot+istat
      temp1 = zero
      temp2 = zero
      temp3 = zero
      temp4 = zero
      temp5 = zero
      temp6 = zero

      allocate(xllji(jpxy,ipxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xllji,zero,SIZE(xllji))
      allocate(xppji(jpxy,ipxy),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(xppji,zero,SIZE(xppji))

      if (relativ .eq. "fully") then
         allocate(gamman(jx,-2:mx+2),STAT=istat)
         istat_tot=istat_tot+istat
         call bcast(gamman,zero,SIZE(gamman))
         allocate(alphan(jx,-2:mx+2),STAT=istat)
         istat_tot=istat_tot+istat
         call bcast(alphan,zero,SIZE(alphan))
         allocate(asnha(jx),STAT=istat)
         istat_tot=istat_tot+istat
         call bcast(asnha,zero,SIZE(asnha))
      endif

      if (ndeltarho.ne."disabled".or.lossmode(1).eq.'simplban' &
           .or. lossmode(1).eq.'simplbn1' &
           .or. lossmode(1).eq.'simplbn2') then
         !YuP[2017-11-21] Need lossmode(k), checking all k?
         allocate(deltarho(iy,lz,setup0%lrzmax),STAT=istat)
         call bcast(deltarho,zero,SIZE(deltarho))
         allocate(deltarhop(nt_delta,lz,setup0%lrzmax),STAT=istat)
         call bcast(deltarhop,zero,SIZE(deltarhop))
         allocate(deltarz(nr_delta,nz_delta,nt_delta),STAT=istat)
         call bcast(deltarz,zero,SIZE(deltarz))
         allocate(r_delta(nr_delta),STAT=istat)
         call bcast(r_delta,zero,SIZE(r_delta))
         allocate(z_delta(nz_delta),STAT=istat)
         call bcast(z_delta,zero,SIZE(z_delta))
         allocate(t_delta(nt_delta),STAT=istat)
         call bcast(t_delta,zero,SIZE(t_delta))
         allocate(delta_bdb0(nr_delta,nz_delta),STAT=istat)
         call bcast(delta_bdb0,zero,SIZE(delta_bdb0))
      endif

      if (tavg.ne."disabled") then
         allocate(favg(0:iy+1,0:jx+1,ngen,lrors),STAT=istat)
         istat_tot=istat_tot+istat
            favg = zero
      endif

      allocate(pentr(nonch,ngen,-1:15,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      pentr = zero

      allocate(constp(nonch,lrors),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(constp,zero,SIZE(constp))

      allocate(sigmtt(nonch,4),sigftt(nonch,4),STAT=istat)
      istat_tot=istat_tot+istat
      call bcast(sigmtt,zero,SIZE(sigmtt))
      call bcast(sigftt,zero,SIZE(sigftt))

      allocate(sgaint(8,ngen,lrors),STAT=istat)
      call bcast(sgaint,zero,8*ngen*lrors)
      allocate(entr(ngen,-1:15,lrors),STAT=istat)
      entr = zero

      allocate(xlndnz(ngen+1,negyrga),STAT=istat)
      call bcast(xlndnz,zero,(ngen+1)*negyrga)
      allocate(sounor(ngen,nsoa,lz,setup0%lrz),STAT=istat)
      call bcast(sounor,zero,ngen*nsoa*lz*setup0%lrz)

      call bcast(sgain,zero,8*ngena)

      allocate(truncd(jx),STAT=istat)
      call bcast(truncd,one,jx)

!      allocate(eflux_r_npa(nen_npa,nv_npa,4*setup0%lrz),STAT=istat)
!      allocate(rho_npa(nv_npa,4*setup0%lrz),STAT=istat)
!      allocate(s_npa(nv_npa,4*setup0%lrz),STAT=istat)
!      allocate(ds_npa(nv_npa,4*setup0%lrz),STAT=istat)

      !YuP[07-2016] added: for neutron flux diagnostics (along view lines)
!      allocate(flux_fus_f(4,nv_fus),STAT=istat)
!      allocate(flux_fus_m(4,nv_fus),STAT=istat)
!      allocate(flux_neutron_f(nv_fus), STAT=istat)
!      allocate(flux_neutron_m(nv_fus), STAT=istat)
!      allocate(flux_rad_fus_f(4,nv_fus,4*setup0%lrz),STAT=istat)
!      allocate(flux_rad_fus_m(4,nv_fus,4*setup0%lrz),STAT=istat)
!      allocate(rho_fus(nv_fus,4*setup0%lrz),STAT=istat)
!      allocate(s_fus(nv_fus,4*setup0%lrz),STAT=istat)
!      allocate(ds_fus(nv_fus,4*setup0%lrz),STAT=istat)

!     Check that allocations were OK
      if (istat_tot.ne.0) then
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'ainalloc.f:  Problem with allocation'
         WRITE(*,*)'ainalloc.f:  Reduce param.h paramaters?'
         WRITE(*,*)'ainalloc.f:  Stopping'
!MPIINSERT_ENDIF_RANK
         STOP
      endif


!dir$ bounds

!..................................................................
!     Sucessful ainalloc
!..................................................................
      write(*,*)'ainalloc:  Leaving ainalloc'

      return
    end subroutine ainalloc

end module ainalloc_mod
