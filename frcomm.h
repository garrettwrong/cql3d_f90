!     frcomm.h
!***********************************************************************
!     BEGIN NFREYA (fr..) common blocks
!***********************************************************************

      include 'frname_decl.h'

!................................................................
      integer nion,nneu,nk,nkt,nj,nbctim,njs,ialignf1
      real(c_double) :: zero,one,two,half,pi
      
      common /numbrs/ nion,nneu,nk,nkt,nj,nbctim,njs,ialignf1
      common /numbrs/ zero,one,two,half,pi
!     ONETWO DIVERGENCE

!................................................................

      integer ibeam,ibion,inubplt,mfm1
      real(c_double) :: bion(ke,kb),bneut(ke,kb)
      real(c_double) :: beamof,ebeam(ke,kb)
      real(c_double) :: ennub(k_),fap(ke,kb),fwall(ke,kb)
      real(c_double) :: hibr(k_,ke,kb),hibrz(kz,ke,kb)
      real(c_double) :: hicmz(kz,ke,kb,3)
      real(c_double) :: ftrapfi(kz,ke,kb),ftrapfit(k_,ke,kb)
      real(c_double) :: angmpf(k_,ke,kb),angmpz(kz,ke,kb)
      real(c_double) :: rowpsi(k_)
      real(c_double) :: pbeam(ke,kb),tenub(k_)
      real(c_double) :: sb(k_,ke,kb),sbcx(k_,2),sbion(k_)
      real(c_double) :: spb(k_,ke,kb),qb(k_,ke,kb),qbf(k_,ke,kb)
      real(c_double) :: zzi(kz,kion),zne(kz),zni(kz,kion)
      real(c_double) :: zte(kz),psivol(kz),freyr(kf)
      real(c_double) :: zeffctv(kz)
      
      common /nub/ bion,bneut,beamof,ebeam
      common /nub/ ennub,fap,fwall
      common /nub/ hibr,hibrz
      common /nub/ hicmz
      common /nub/ ftrapfi,ftrapfit
      common /nub/ angmpf,angmpz
      common /nub/ ibeam,ibion,inubplt,mfm1
      common /nub/ rowpsi
      common /nub/ pbeam,tenub
      common /nub/ sb,sbcx,sbion
      common /nub/ spb,qb,qbf
      common /nub/ zzi,zne,zni
      common /nub/ zte,psivol,freyr,zeffctv
!................................................................

      real(c_double) :: bencap(k_,ke,kb),fbe(k_,ke,kb),fbi(k_,ke,kb)
      real(c_double) :: bke(k_,ke,kb),bki(k_,ke,kb)
      real(c_double) :: ecrit(k_),emzrat(k_),enbeam(k_),enbs(k_)
      real(c_double) :: enb(k_,ke,kb),enbsav(k_,ke,kb)
      real(c_double) :: enbav(k_,ke,kb),enbav0(k_),enbav1(k_)
      real(c_double) :: forb(ke,kb),fb11(ke,kb),fb10(ke,kb)
      real(c_double) :: fb01(ke,kb),fb00(ke,kb),fber(ke,kb)
      real(c_double) :: hdep(k_,ke,kb),hdepz(kz,ke,kb)
      real(c_double) :: ppb(k_,ke,kb),ppbsav(k_,ke,kb),ppbav(k_,ke,kb)
      real(c_double) :: pinsid(kf),potsid(kf),rinsid(kf),rotsid(kf)
      real(c_double) :: qbsav(k_,ke,kb)
      real(c_double) :: sbsav(k_,ke,kb),spbsav(k_,ke,kb)
      real(c_double) :: taupb(k_,ke,kb),tauppb(k_,ke,kb)
      real(c_double) :: taueb(k_,ke,kb),taus(k_),wbeam(k_)
      real(c_double) :: wb(k_,ke,kb),wbsav(k_,ke,kb),wbav(k_,ke,kb)
      real(c_double) :: wb11(ke,kb),wb10(ke,kb)
      real(c_double) :: wb01(ke,kb),wb00(ke,kb)
      real(c_double) :: zetaz(kz,ke,kb),zeta(k_,ke,kb)

      common /nub2/ bencap,fbe,fbi
      common /nub2/ bke,bki
      common /nub2/ ecrit,emzrat,enbeam,enbs
      common /nub2/ enb,enbsav
      common /nub2/ enbav,enbav0,enbav1
      common /nub2/ forb,fb11,fb10
      common /nub2/ fb01,fb00,fber
      common /nub2/ hdep,hdepz
      common /nub2/ ppb,ppbsav,ppbav
      common /nub2/ pinsid,potsid,rinsid,rotsid
      common /nub2/ qbsav
      common /nub2/ sbsav,spbsav
      common /nub2/ taupb,tauppb
      common /nub2/ taueb,taus,wbeam
      common /nub2/ wb,wbsav,wbav
      common /nub2/ wb11,wb10
      common /nub2/ wb01,wb00
      common /nub2/ zetaz,zeta

      pointer xpts,ypts,zpts,rpts  ! (1:npart)
!     ONETWO DIVERGENCE
      pointer vx,vy,vz  ! (1:npart)
      real(c_double) :: xpts(:),ypts(:),zpts(:),rpts(:)
      real(c_double) :: vx(:),vy(:),vz(:)
      common /xyzpts/ xpts,ypts,zpts,rpts,vx,vy,vz
!................................................................

      integer ncont,iz(kimp)
      common /nub3/ ncont,iz
      real(c_double) :: znipm(kprim),atwpm(kprim)
      real(c_double) :: atwim(kimp),zniim(kimp),zti(kz)
      common /nub3/znipm,atwpm,atwim,zniim,zti
!     nub3 is is used for variables related to hexnb routine
!     note that nouthx and ncorin,also required for hexnb,
!     have been added to block io.

!................................................................

      integer ncrt,nin,nout,nqik,neqplt,ntrplt,nitre,ngreen
      integer nbplt,neq,nsvsol,nscr,nrguess,nwguess
      integer nmix,ntweak,nyok,nupel,nitrex
      integer ialign25
      integer iprt,ialign26
      integer mprt,ialign27
      integer jprt,ialign28
      integer jflux,jcoef
      integer jsourc,jbal,jtfus,ihead,ineu,inub,irfcalc
      integer ialign30,ifred,ialign31,nterow
      integer jterow(10),ilastp,itimav,ialign32
      integer ncorin,iotoray
      real(c_double) :: eqdskin_fr,guessin,guessout
      real(c_double) :: timprt,prtlst(10)
      real(c_double) :: timplt,pltlst(30)
      real(c_double) :: banktime,extime
      real(c_double) :: vid,ddebug(50)
      common /io/ ncrt,nin,nout,nqik,neqplt,ntrplt,nitre,ngreen
      common /io/ nbplt,neq,nsvsol,nscr,nrguess,nwguess
      common /io/ nmix,ntweak,nyok,nupel,nitrex
      common /io/ ialign25,eqdskin_fr,guessin,guessout
      common /io/ iprt,ialign26,timprt,mprt,ialign27,prtlst
      common /io/ jprt,ialign28,timplt,pltlst,jflux,jcoef
      common /io/ jsourc,jbal,jtfus,ihead,ineu,inub,irfcalc
      common /io/ ialign30,ifred,ialign31,banktime,extime,nterow
      common /io/ jterow,ilastp,itimav,ialign32,vid,ddebug
      common /io/ ncorin,iotoray
!BH070410  mplot removed since not used and possible conflict
!BH070410  with other use of this variable name.

!................................................................

      character*8 namep,namei
      common /ions/ namep(kprim),namei(kimp)
      integer namen(2),nameu(kkq)
      common /ions/ namen,nameu
      real(c_double) :: atw(kion),dzdtim(k_,kimp)
      real(c_double) :: z(k_,kion),zsq(k_,kion),dzdte(k_,kion)
      real(c_double) :: zeff_(k_),rfatmwt
      common /ions/ atw,dzdtim
      common /ions/ z,zsq,dzdte
      common /ions/ zeff_,rfatmwt

!................................................................
!     ONETWO DIVERGENCE

      real(c_double) :: p(ki,kj),xxx(ki),yyy(kj)
      common/mhd1/ p,xxx,yyy
!................................................................




