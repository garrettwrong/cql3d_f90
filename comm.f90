!     comm
!
!.......................................................................
!     There are a number of arrays intended for temporary storage
!     in a subroutine, or for simple passes of data to
!     associated subroutines:
!
!BH180527:  The equivalences listed for tem[1-6]/temp[1-6] are
!BH180527:  have been removed (at some past time).
!     iyjx2=(iy+2)*(jx+2)
!     tem1(iyjx2) --> tem6(iyjx2)
!     temp1(0:iyp1,0:jxp1) --> temp6(0:iyp1,0:jxp1), 
!      "equivalenced" to tem1(iyjx2) --> tem6(iyjx2).
!     item1(iyjx2) --> item6(iyjx2)
!      "equivalenced" to tem1(iyjx2) --> tem6(iyjx2).
!     iyjx2=(iy+2)*(jx+2)
!     MOREOVER: We assume temp[1-6] are in contiguous storage,
!               so we can reference the whole six arrays through tem1.
!BH180527: Dimension of tem[1-6] modified to iyjx2l=max(iy+2,lrz)*(jx+2),
!BH180527: to handle possible situation in netcdfrw2.
!
!     tam1(jx) --> tam30(jx)
!     temc1(iy) --> temc4(iy)
!     tz1(lza) --> tz2(lza)
!     tr(0:lrza)
!     tr1(0:lrza) --> tr5(0:lrza)
!     itemc1(iy) --> itemc2(iy)
!     urftmp(nrayelts*nrayn*5)
!.......................................................................
!
!
!.......................................................................
!     Add in type,size,common declarations for namelist variables.
!     This needs to precede the rest of comm.h declarations so
!     that name_decl.h can be used by itself in aindflt.f.
!     aindflt.f will also be used for setting defaults in the
!     SWIM project Integrated Plasma Simulation (IPS) modules.
!.......................................................................

module comm_mod
  use iso_c_binding, only : c_float             !REAL*4
  use iso_c_binding, only : c_double            !REAL*8
  use iso_c_binding, only : c_double_complex    !COMPLEX*16
  use param_mod
  implicit none
  public

  !.......................................................................
  !     nml variables that take on the values assigned to parameters.
  !.......................................................................
  integer :: iy,jx, &
       lfield,lz,lrz, &
       lrzmax,ls,lsmax, &
       mx, &
       nbctime,negyrg,ngen,nmax

  !.......................................................................
  !     SCALAR INPUT FOR NAMELIST SETUP...
  !.......................................................................

  character(len=8) :: chang, &
       eqmod,eleccomp,f4d_out,tavg, &
       iactst,ineg,idrop,idskf,idskrf,ichkpnt,implct, &
       lbdry0,locquas,lrzdiff,lsdiff,taunew, &
       machine,meshy,manymat, &
       netcdfvecal,netcdfvecc,netcdfvece,netcdfvecrf,netcdfvecs, &
       noplots, &
       psimodel,pltpowe,pltend,pltinput,pltlim,pltrdc, &
       pltrst,plturfb,pltvflu,pltra,pltfvs,pltd,pltprpp,pltfofv,pltlos, &
       pltdn,pltvecal,pltvecc,pltvecrf,pltvece,pltstrm,pltflux, &
       pltsig,pltdnpos, &
       profpsi, &
       qsineut,trapmod,scatmod, &
       relativ,sigmamod,soln_method, &
       symtrap,syncrad,bremsrad, &
       tandem,gamafac, &
       yreset, &
       nmlstout

  character(len=256) :: netcdfnm
  character(len=256) :: mnemonic
  character(len=8) :: iuser,izeff,netcdfshort

  integer :: colmodl


  !common /readscal/ &
  real(c_double) :: btor,bth, &
       contrmin,constr,&
       deltabdb,droptol,dtr,dtr0,xsink, &
       esink,ephicc,esfac,eoved,enorm,enorme,enormi,&
       gsla,gslb,gamaset
  integer :: isigtst,isigsgv1,isigsgv2, &
       kenorm,lfil, &
       nnspec

  complex(c_double_complex) :: vlfemin(nmodsa),vlfeplus(nmodsa)

  !common /readscal/ &
  integer :: nstop,ncoef,nchec,ncont,nrstrt,nstps,nfpld, &
       noncntrl,nonel,noffel,nonvphi,noffvphi,nonloss,noffloss, &
       numby,lnwidth
  real(c_double) :: pltmag, &
       xprpmax, &
       trapredc,scatfrac, &
       radmaj,radmin,rmirror, &
       sigvi,sigvcx, &
       brfac,brfac1,brfacgm3, &
       tfac,tfacz,temp_den, &
       veclnth, &
       vnorm

  real(c_double) :: npwr(0:ntotala),negy(ngena),ntorloss(ngena), &
       mpwr(0:ntotala),megy(ngena),mtorloss(ngena), &
       mpwrzeff,npwrzeff,mpwrvphi,npwrvphi, &
       mpwrelec,npwrelec,mpwrxj,npwrxj

  !common /readscal/ &
   real(c_double) ::   xfac,xpctlwr,xpctmdl,xlwr,xmdl,xsinkm, &
       ylower,yupper, &
       elecscal,enescal,zeffscal,vphiscal,tescal,tiscal,bctimescal
  !..................................................................
  !     VECTOR DIMENSIONED NAMELIST SETUP COMMON BLOCK.....
  !..................................................................

  character(len=8) :: ibox(3), &
       kpress(ntotala),kfield(ntotala), &
       lbdry(ngena),lossmode(ngena), &
       regy(ngena), &
       torloss(ngena), &
       difus_type(ngena), difus_io(ngena)

  character(len=256) :: lossfile(ngena)

  !common /readvec/ &
  real(c_double) :: bnumb(ntotala), &
       fmass(ntotala), &
       ioutput(2),isigmas(6),pltflux1(7), &
       enloss(ngena), &
       gamegy(ngena), &
       paregy(ngena),peregy(ngena),pegy(ngena), &
       zeffin(0:njenea),vphiplin(0:njenea), &
       ennl(npaproca),ennb(npaproca),ennscal(npaproca), &
       bctime(nbctimea),dtr1(ndtr1a), &
       zeffc(nbctimea),zeffb(nbctimea), &
       elecc(nbctimea),elecb(nbctimea), &
       vphic(nbctimea),vphib(nbctimea), &
       xjc(nbctimea),xjb(nbctimea), &
       totcrt(nbctimea), &
       nplot(nplota),nsave(nsavea), &
       tavg1(ntavga),tavg2(ntavga)
  ! XXX bug
  integer :: nondtr1(ndtr1a)

  !..................................................................
  !     TWO DIMENSIONAL NAMELIST SETUP COMMON BLOCK.....
  !..................................................................

  !common /readarr/ &
  real(c_double) :: tauloss(3,ngena)
  character(len=8) :: kspeci(2,ntotala)
  real(c_double) ::  fpld(10,ngena), &
       redenc(nbctimea,ntotala),redenb(nbctimea,ntotala), &
       tempc(nbctimea,ntotala),tempb(nbctimea,ntotala)



  !**********************************************************************
  !     Variables in common block diskx have as their last dimension lrza.
  !     Thus they are dimensioned with respect to the radial coordinate.
  !     If lrza=1  ==> single flux surface CQL run.
  !     Note that input variable lrz can be less than lrza. For some of
  !     the larger arrays we allocate space using lrz rather than
  !     dimensioning with lrza to save some space.
  !
  !     VECTORS
  !
  !**********************************************************************
  !
  !common /diskx/ &
  real(c_double) :: elecfld(0:lrza), &
       tbnd(lrorsa), &
       rovera(lrza), &
       zmax(0:lrza)
  integer :: lrindx(0:lrorsa),lsindx(0:lrorsa)

  !..................................................................
  !     2-D ARRAYS
  !..................................................................

  !common /diskx/ &
  real(c_double) :: reden(ntotala,0:lrza), &
       temp(ntotala,0:lrza)

  !common /diskx/ &
  real(c_double) :: tauegy(ngena,0:lrza),eparc(ngena,0:lrza), &
       eperc(ngena,0:lrza),simpbfac, &
       isoucof,faccof

  !..................................................................
  !     3-D ARRAYS
  !..................................................................

  !..................................................................
  !     4-D ARRAYS
  !..................................................................

  !common /diskx/ &
  real(c_double) :: eegy(negyrga,2,ngena,lrza)
  integer :: jegy(negyrga,2,ngena,lrza)

  !*****************************************************************
  !     BEGIN arrays for analytic ion source (sou..) routines
  !*****************************************************************

  !common /diskx/ &
  integer ::  mpwrsou(0:ngena),npwrsou(0:ngena)
  real(c_double) ::asor(ngena,nsoa,lrza)

  !common /params/ &
  integer :: nso

  character(len=8) :: pltso, soucoord, knockon, komodel, flemodel

  !common /readscal/ &
  integer ::  nsou
  integer :: nkorfn,nonko,noffko
  real(c_double) :: soffvte, soffpr, &
       xlfac,xlpctlwr,xlpctmdl,xllwr,xlmdl
  integer :: jfl

  !common /readarr/ &
  integer :: nonso(ngena,nsoa),noffso(ngena,nsoa)
  real(c_double) :: sellm1(ngena,nsoa),sellm2(ngena,nsoa),seppm1(ngena,nsoa), &
       seppm2(ngena,nsoa),sem1(ngena,nsoa),sem2(ngena,nsoa), &
       sthm1(ngena,nsoa),scm2(ngena,nsoa),szm1(ngena,nsoa), &
       szm2(ngena,nsoa)


  !*****************************************************************
  !     BEGIN arrays for rf package..(rf...,vlh[B,...,vlf...) routines
  !*****************************************************************
  character(len=8) :: urfmod, &
       vlfmod,vlfbes,vlfnpvar,vlhmod,vprprop,vlhplse,vlhprprp(nmodsa), &
       rfread,rdcmod,rdc_clipping,rdc_netcdf
  character(len=256) :: rffile(nmodsa)
  character(len=256) :: rdcfile(nmodsa)

  !common /params/ &
  integer :: nrf
  !common/readscal/
  real(c_double) :: vlhmodes,vdalp,vlh_karney, &
       vlhpon,vlhpoff, &
       vlfmodes, &
       rdc_upar_sign

  !common /readvec/ &
  integer :: nonrf(ngena),noffrf(ngena)
  real(c_double) :: dlndau(nmodsa),vparmin(nmodsa),vparmax(nmodsa), &
       vprpmin(nmodsa),vprpmax(nmodsa), &
       vlhpolmn(nmodsa),vlhpolmx(nmodsa), &
       vlffreq(nmodsa),vlfharms(nmodsa),vlfharm1(nmodsa), &
       vlfnp(nmodsa),vlfdnp(nmodsa),vlfddnp(nmodsa), &
       vlfpol(nmodsa),vlfdpol(nmodsa),vlfddpol(nmodsa), &
       vlfnperp(nmodsa),vlfdnorm(nmodsa), &
       vlfparmn(nmodsa),vlfparmx(nmodsa), &
       vlfprpmn(nmodsa),vlfprpmx(nmodsa)


  !..................................................................
  !     arrays in input
  !..................................................................

  !common/arr3d/ &
  real(c_double) :: acoefne(4),acoefte(4),ryain(njenea),elecin(njenea), &
       enein(njenea,ntotala),tein(njenea),tiin(njenea), &
       ennin(njenea,npaproca)
  integer :: irzplt(lrorsa)
  real(c_double) :: rya(0:lrza+1), &
       thet1(nva),thet2(nva),thet1_npa(nva),thet2_npa(nva)
  integer :: nplt3d(nplota)

  !common/arr3d/ &
  real(c_double) :: enein_t(njenea,ntotala,nbctimea),tein_t(njenea,nbctimea), &
       tiin_t(njenea,nbctimea),zeffin_t(njenea,nbctimea), &
       elecin_t(njenea,nbctimea),xjin_t(njenea,nbctimea), &
       vphiplin_t(njenea,nbctimea),  &
       ennin_t(njenea,nbctimea,npaproca) !neutrals,impurities,etc.

  !common/arr3d/ &
  real(c_double) :: sellm1z(ngena,nsoa,0:lrza),sellm2z(ngena,nsoa,0:lrza), &
       seppm1z(ngena,nsoa,0:lrza),sem1z(ngena,nsoa,0:lrza), &
       sem2z(ngena,nsoa,0:lrza),sthm1z(ngena,nsoa,0:lrza), &
       scm2z(ngena,nsoa,0:lrza),szm1z(ngena,nsoa,0:lrza), &
       seppm2z(ngena,nsoa,0:lrza), &
       szm2z(ngena,nsoa,0:lrza),asorz(ngena,nsoa,0:lrza)


  !****************************************************************
  !     BEGIN arrays for 3-d (td..) driver.
  !****************************************************************

  !common/params/ ndifus_io_t  !Max to be nbctimea
  integer :: ndifus_io_t

  !..................................................................
  !     scalars in input
  !..................................................................

  character(len=8) :: bootst,bootcalc,bootupdt, &
       iprone,iprote,iproti,iprozeff,iprovphi,iproelec,ipronn,iprocur, &
       tmdmeth,partner,pinch,plt3d,pltvs,radcoord, &
       relaxtsp,rzset, &
       ndeltarho,softxry,npa_diag,atten_npa, &
       transp,adimeth, &
       efswtch,efswtchn,efiter,efflag
  ! bug
  character(len=8) :: npa_process(npaproca)

  character(len=256) :: difus_io_file

  !common /s3d/ &
  real(c_double) :: difusr,advectr,&
       enmin,enmax,bootsign,fds
  integer :: kfrsou, &
       mmsv,msxr,njene,njte,njti,nonboot,jhirsh, &
       nrskip,nen,nv,nen_npa,nv_npa,npaproc, &
       nr_delta,nz_delta,nt_delta, &
       nr_f4d,nz_f4d,nv_f4d,nt_f4d
  real(c_double) :: rfacz,roveram,relaxden
  real(c_double) :: enmin_npa,enmax_npa,fds_npa, &
       nonadi,curr_edge,efrelax,efrelax1,currerr

  !common/ar3d/ &
  real(c_double) :: difus_rshape(8),difus_vshape(4),difin(njenea), &
       rd(nva),thetd(nva),x_sxr(nva),z_sxr(nva), &
       rd_npa(nva),thetd_npa(nva),x_npa(nva),z_npa(nva)
  real(c_double) :: difus_io_drrscale(nbctimea,ngena),  &
       difus_io_drscale(nbctimea,ngena), &
       difus_io_t(nbctimea)

  !******************************************************************
  !     BEGIN arrays for EQUILIBRIUM MODEL (eq..) (NON-CIRCULAR CROSS
  !     SECTIONS).
  !******************************************************************

  character(len=8) :: nconteq

  !common/params/
  integer :: nconteqn

  character(len=8) :: eqsym,eqdskalt,eqsource,eqmodel, fpsimodl

  !     ONETWO uses character*60 for eqdskin.
  character(len=256) :: eqdskin

  !common/readscal/ &
  real(c_double) :: atol, &
       ellptcty,eqpower,bsign, &
       povdelp, &
       rtol,rmag,rbox,rboxdst, &
       zbox
  integer :: methflag



  !*********************************************************************
  !     BEGIN arrays for LOWER HYBRID FAST WAVE and ECH Module.
  !*********************************************************************


  !common/params/ &
  integer :: nurftime

  character(len=8) :: urfdmp,iurfcoll(nmodsa),iurfl(nmodsa), &
       call_lh, &
       call_ech, &
       call_fw, &
       ech,fw,lh, &
       scaleurf,urfrstrt,urfwrray, &
       rftype(nmodsa)

  !common/readscal/ &
  integer :: ieqbrurf,urfncoef
  integer :: nmods,nbssltbl,nondamp,nrfstep2, &
       nrfpwr,nrfitr1,nrfitr2,nrfitr3
  real(c_double) :: urfmult
  integer :: nrdc

  !common/readvec/ &
  real(c_double) :: pwrscale(nmodsa),wdscale(nmodsa)
  integer :: nrfstep1(nmodsa), &
       nharms(nmodsa),nharm1(nmodsa),nrfspecies(nmodsa)

  !common/readvec/ &
  real(c_double) :: pwrscale1(nbctimea),urftime(nbctimea)

  !common/readvec/ &
  real(c_double) :: rdcscale(nrdca)
  integer :: nrdcspecies(nrdca)


  !-----------------------------------------------------------------------
  !     BEGIN variables for WP... modules for CQLP case
  !-----------------------------------------------------------------------

  character(len=8) :: special_calls,cqlpmod, &
       oldiag, &
       sbdry,scheck,ampfmod,eseswtch, &
       updown

  !      logical (lolz)
  character(len=8) :: nlrestrt,nlwritf

  logical :: nlotp1(noutpta),nlotp2(noutpta)
  logical :: nlotp3(noutpta),nlotp4(noutpta)

  !common /readscal/ &
  real(c_double) :: epsthet, elpar0
  integer :: lmidpln,lmidvel,laddbnd, &
       nchgdy,ngauss,nlagran, &
       nonavgf,nofavgf, &
       nontran, nofftran, nonelpr, noffelpr, &
       nummods,numixts
  real(c_double) :: ampferr
  integer :: nampfmax,nonampf

  !common /readvec/ &
  real(c_double) :: denpar(ntotala,0:lza+1)
  integer :: nkconro(ntotala)
  real(c_double) :: temppar(ntotala,0:lza+1)


  !.......................................................................
  !     Setup block for finite orbit width (FOW) calculations
  !.......................................................................
  character(len=8) ::  fow
  character(len=16) :: outorb 
  character(len=38) :: file_fow_plt ! for saving data on orbit to a file
  !common/fow_control/
  integer :: nmu,npfi
  integer :: nsteps_orb,nptsorb,i_orb_width,iorb2, &
       j0_ini,j0_end,inc_j0, i0_ini,i0_end,inc_i0, &
       j2_ini,j2_end,inc_j2, i2_ini,i2_end,inc_i2
  ! fow= 'enabled' or 'disabled' 
  ! outorb  ! 'detailed' or 'Not-detailed'
  ! (saving/not-saving data to a file for plotting)  
  ! nmu     ! grid sizes for ad.ivariant mu 
  ! npfi    ! and canonical momentum Pfi; 
  ! to setup COM->R lookup table.
  ! nsteps_orb ! Max.number of time steps for orbit integration.
  ! Also used to trace Pfi=const levels for COM->R table
  ! in order to find intersections with mu=const levels.
  ! nptsorb ! Number of points on a complete orbit 
  ! (ityp=0 "main" orbit)
  ! from which ityp=1 "secondary" orbits are launched.
  ! ityp=1 orbit is stopped when it reaches the midplane.
  ! (Note: secondary orbits are not traced usually, 
  ! see below, iorb2=0)
  ! i_orb_width ! 1 -> Normal finite-orbit-width calculations. 
  ! 0 -> V_drift_perp is set to 0 (ZOW approximation)
  ! iorb2  ! set to 1 to perform Runge-Kutta integration for tracing
  ! SECONDARY orbits to midplane; 0 - no RK tracing.
  ! This option (1) can be used for plotting orbits
  ! (together with outorb="detailed"),
  ! otherwise not needed.


  !.......................................................................
  !     variables that take on the values assigned to parameters.
  !.......................................................................

  !common /params/ &
  integer :: idim,iyjx,iyjxp1,iyp1jx, iyjx2, iyp1,jxp1, &
       lrors, &
       mxp1,mbet, &
       nonch,niong,nionm,ntotal


  !**********************************************************************
  !     Variables in common block diskx have as their last dimension lrza.
  !     Thus they are dimensioned with respect to the radial coordinate.
  !     If lrza=1  ==> single flux surface CQL run.
  !     Note that input variable lrz can be less than lrza. For some of
  !     the larger arrays we allocate space using lrz rather than
  !     dimensioning with lrza to save some space.
  !
  !     VECTORS
  !
  !**********************************************************************
  !
  !common /diskx/ &
  real(c_double) :: bmod0(lrza),btor0(lrza),bthr(lrza),btoru(lrza),bthr0(lrza), &
       consn(lrorsa),consn0(lrorsa),currt(lrza),currxj0(0:lrza), &
       currtp(lrza),currmtp(lrorsa),currmt(lrorsa),currxj(0:lrza), &
       currpar(lrza),curreq(lrza), &
       zreshin(lrza),zreskim(lrza),rovsc(lrza),rovsc_hi(lrza), &
       curtor(lrza),curpol(lrza),ccurtor(0:lrza),ccurpol(0:lrza), &
       deltapsi(lrza), &
       eps(lrza),etll0(lrza)
  integer :: itl_(lrorsa),itu_(lrorsa), &
       iy_(lrorsa),iyh_(lrorsa),iyjx_(lrorsa), &
       inew_(lrorsa),inewjx_(lrorsa),ieq_(lrorsa+1), &
       indxlr(0:lrorsa),indxls(0:lrorsa), &
       lorbit(lrza),lmdpln(0:lrza), &
       n_(lrorsa),nch(lrorsa), &
       nefiter_(lrza) ! counts iterations of el.field for each flux surface      &
  real(c_double) :: psimx(lrza),pibzmax(lrza),psidz(lrza), &
       qsafety(lrza), &
       r0geom(lrza), r0drdz(0:lrza),rgeom(lrza), zgeom(lrza), &
       rovs(lrza),rovsn(lrza),rovsloc(lrorsa)
  !common /diskx/ &
  real(c_double) :: sptzr(lrorsa),sgaint1(lrorsa),starnue(lrorsa), &
       thb(lrorsa),tauee(lrza),taueeh(lrorsa),time_(lrorsa), &
       twoint(lrorsa), &
       vthe(lrza), &
       xlbnd(lrza),xlndn0(lrza), &
       zmaxpsi(0:lrza),zmaxi(0:lrza), &
       zmaxpsii(0:lrza),zeff(lrza),zeff4(lrza), &
       vphipl(lrza), &
       srckotot(lrza),elecr(lrza), &
       denfl(lrza),denfl1(lrza),denfl2(lrza), &
       den_of_s(lza),den_of_s1(lza),den_of_s2(lza), &
       tauii(lrza),tau_neo(lrza),drr_gs(lrza), &
       rhol(lrza),rhol_pol(lrza)
  !common /diskx/ &
  real(c_double) :: taubi(lrza),tau_neo_b(lrza),drr_gs_b(lrza), &
       rhol_b(lrza),rhol_pol_b(lrza)

  !..................................................................
  !     2-D ARRAYS
  !..................................................................

  !common /diskx/ &
  real(c_double) :: energy(ntotala,lrza), vth(ntotala,lrza)


  !common /diskx/ &
  real(c_double) :: den_fsa(ngena,lrza),  &
       den_fsa_t0(ngena,lrza), reden_t0(ngena,lrza), &
       currm(ngena,lrorsa),curr(ngena,lrza), &
       hnis(ngena,lrza), &
       ratio(ngena,lrza), &
       wpar(ngena,lrza),wperp(ngena,lrza), &
       xlndn00(ngena,lrza),xlncur(ngena,lrza), &
       xlndn(ngena,lrza),energym(ngena,lrorsa), &
       xlndnr(ntotala,lrza),energyr(ntotala,lrza), &
       currr(ngena,lrza),xlndnv(ntotala,lrza), &
       energyv(ntotala,lrza),currv_(ngena,lrza),eoe0(ngena,lrza), &
       ucrit(ngena,lrza)
  integer :: jxcrit(ngena,lrza),      &
       jchang(ngena,lrorsa)
  real(c_double) :: denra(ngena,lrza),curra(ngena,lrza), &
       fdenra(ngena,lrza),fcurra(ngena,lrza),enn(lrza,npaproca)
  !common /diskx/ &
  real(c_double) :: alm(0:mbeta,lrza), &
       betta(0:mbeta,lrza), &
       entrintr(ngena,-1:15)

  !..................................................................
  !     Next follow variables and arrays not dimensioned by lrza.
  !..................................................................

  !..................................................................
  !     SCALARS.....
  !..................................................................

  character(len=8) :: outfile, prefix

  character(len=512) ::  t_

  !common &
  real(c_double) :: clight,charge,restmkev,symm, &
       cnorm,cnorm2,cnorm3,clite2,cnorm2i,cnormi, &
       ergtkev,eps0,elect,eratio,eovsz,eovedd, &
       em6,em8,em10,em12,em14,em37,em40,ep37,em90,ep90, &
       ep100,em100,em120,em300,four,fourpi,fc, &
       gszb,gsb,gsem,gsep,gszm,gszp,gseb,gszm2,gszp2,gszb2,gsla2, &
       gslb2,gsnm,half
  integer :: iyh,istate,itl,itu, &
       impcoef,imprf,impelec,impcah,irstart,impadi, &
       jxm1,jlwr,jmdl, &
       kgnande,kelecg,kelecm,kelec,kionn,kiong(ngena),kiongg(ngena), &
       kionm(nmaxa), &
       l_,lr_,indxlr_,indxls_,lmdpln_,ls_, &
       miyjx,ibeampon,ibeamponp

  !common &
  integer :: ipad1
  real(c_double) :: r0geomp, rgeomp, zgeomp, rgeom1, rgeom2
  integer :: niyjx,nccoef,ncentr,navnc,ncplt,n, &
       nflux,nframe,npltmx, ipad2
  real(c_double) :: one,one_, &
       pi,pio180,pio2,psir,proton, &
       rgyro,resist,resistn,rovsf,rmp0,rovscf,rbgn, &
       stopm, &
       tei,timet,three,third,sevenhun,two,twopi, &
       vparl,vpovth, &
       xconn,xmax
  integer :: ipxy,jpxy
  real(c_double) :: zero,zstar
  integer :: iplot,nplott,isave,nsavet,nirzplt, &
       nefiter ! counts iterations for el.field (efswtch="method5")
  !common &
  real(c_double) :: elecfldc,elecfldb,  & !V/cm, Added for ampfmod case
       it_ampf   !iteraction counter for ampfmod case

  !..................................................................
  !     VECTORS.....
  !..................................................................

  ! bug
  character(len=8) :: text(4)

  !common &
  integer :: ipad3, iota(0:6400)
  !BH111124     1  iota(0:750),realiota(0:750), &
  real(c_double) :: realiota(0:6400), tab(3)
  integer ::i1p(2),itab(3),ipad4, &
       imsh(20)
  real(c_double) :: frbuf(1024), &
       ws2d(2500), &
       work(tlfld1a)
  !common &
  real(c_double) :: ekev(ngena), &
       engain(ngena), &
       tnorm(ngena)
  !common &
  real(c_double) :: fions(ntotala), &
       gama0(ntotala), &
       tz1(lza),tz2(lza)

  !common &
  real(c_double) :: wkbc(3*nbctimea),iopbc(2)

  ! bug
  character(len=8) :: mplot(lrorsa)

  !common &
  real(c_double) ::  tplot(nplota),tsave(nsavea)

  !..................................................................
  !     TWO-DIMENSIONAL ARRAYS.....
  !..................................................................

  !common &
  real(c_double) :: gamt(ntotala,ntotala),gama(ntotala,ntotala), &
       satioz2(ntotala,ntotala),satiom(ntotala,ntotala), &
       sgain(8,ngena)

  !..................................................................
  !     Allocatable arrays allocated in subroutine ainalloc
  !..................................................................

  integer, pointer :: ix1(:),ix2(:),ix3(:),ix4(:),   &  !!!(0:mx) &
       ix5(:),ix6(:),ix7(:),ix8(:)      !!!(0:mx)

  real(c_double), pointer :: tom1(:),tom2(:),tom3(:),tom4(:)  !!!(0:mxp1)

  !BH_YP090809  WAS: choose(0:mx+2,0:mx+2),fctrl(0:mx+2)
  !BH_YP090809  First dimension of choose should be larger of 2*mx,mx+2
  real(c_double), pointer :: fctrl(:)     !!!(0:2*mx+2)
  real(c_double), pointer :: choose(:,:)  !!!(0:2*mx+2,0:mx+2)
  real(c_double), pointer :: cog(:,:)     !!!(0:mx,15)
  real(c_double), pointer :: pm(:,:)      !!!(0:mxp1,lrors)

  !common /dptr95/

  real(c_double), pointer :: f(:,:,:,:)  !f(0:iy+1,0:jx+1,ngen,lrors)
  !common /dptr95/ f
  real(c_double), pointer :: favg(:,:,:,:)
  !common /dptr95/ favg
  real(c_double), pointer :: fxsp(:,:,:,:)
  !common /dptr95/ fxsp
  real(c_double), pointer :: f_(:,:,:,:)
  !common /dptr95/ f_
  real(c_double), pointer :: spasou(:,:,:,:)
  !common /dptr95/ spasou
  real(c_double), pointer :: velsou(:,:,:,:)
  !common /dptr95/ velsou
  real(c_double), pointer :: velsou2(:,:,:,:)
  !common /dptr95/ velsou2
  real(c_double), pointer :: source(:,:,:,:)
  !common /dptr95/ source
  real(c_double), pointer :: gone(:,:,:,:)
  !common /dptr95/ gone
  real(c_double), pointer :: egylosa(:,:,:,:)
  !common /dptr95/ egylosa
  integer, pointer :: i0tran(:,:,:)
  !common /dptr95/ i0tran
  real(c_double), pointer :: cal(:,:,:,:)
  !common /dptr95/ cal
  real(c_double), pointer :: cbl(:,:,:,:)
  !common /dptr95/ cbl
  real(c_double), pointer :: ccl(:,:,:,:)
  !common /dptr95/ ccl
  real(c_double), pointer :: cdl(:,:,:,:)
  !common /dptr95/ cdl
  real(c_double), pointer :: cel(:,:,:,:)
  !common /dptr95/ cel
  real(c_double), pointer :: cfl(:,:,:,:)
  !common /dptr95/ cfl
  real(c_double), pointer :: eal(:,:,:,:,:)
  !common /dptr95/ eal
  real(c_double), pointer :: ebl(:,:,:,:,:)
  !common /dptr95/ ebl
  real(c_double), pointer :: scal(:,:)
  !common /dptr95/ scal
  real(c_double), pointer :: cet(:,:,:)
  !common /dptr95/ cet
  real(c_double), pointer :: cex(:,:,:)
  !common /dptr95/ cex
  real(c_double), pointer :: synca(:,:,:)
  !common /dptr95/ synca
  real(c_double), pointer :: syncd(:,:,:)
  !common /dptr95/ syncd
  real(c_double), pointer :: taulos(:,:,:)
  !common /dptr95/ taulos
  real(c_double), pointer :: elecfldn(:,:,:)
  !common /dptr95/ elecfldn
  real(c_double), pointer :: delecfld0n(:,:,:)
  !common /dptr95/ delecfld0n
  real(c_double), pointer :: elecn(:,:,:)
  !common /dptr95/ elecn
  real(c_double), pointer :: delecfld0(:,:)
  !common /dptr95/ delecfld0
  real(c_double), pointer :: psi0bar(:,:)
  !common /dptr95/ psi0bar
  real(c_double), pointer :: di(:,:,:,:)
  !common /dptr95/ di
  real(c_double), pointer :: dj(:,:,:,:)
  !common /dptr95/ dj
  real(c_double), pointer :: dym5(:,:)
  !common /dptr95/ dym5
  real(c_double), pointer :: dyp5(:,:)
  !common /dptr95/ dyp5
  real(c_double), pointer :: eyp5(:,:)
  !common /dptr95/ eyp5
  real(c_double), pointer :: eym5(:,:)
  !common /dptr95/ eym5
  real(c_double), pointer :: y(:,:)
  !common /dptr95/ y
  real(c_double), pointer :: dy(:,:),dyi(:,:)
  !common /dptr95/ dy,dyi
  real(c_double), pointer :: yptb(:,:)
  !common /dptr95/ yptb
  real(c_double), pointer :: coss(:,:)
  !common /dptr95/ coss
  real(c_double), pointer :: cynt2(:,:)
  !common /dptr95/ cynt2
  real(c_double), pointer :: batot(:,:)
  !common /dptr95/ batot
  integer, pointer :: lmax(:,:)
  !common /dptr95/ lmax
  real(c_double), pointer :: vpint(:,:)
  !common /dptr95/ vpint
  real(c_double), pointer :: psiiv(:,:)
  !common /dptr95/ psiiv
  real(c_double), pointer :: psiba(:,:)
  !common /dptr95/ psiba
  real(c_double), pointer :: psisq(:,:)
  !common /dptr95/ psisq
  real(c_double), pointer :: psicu(:,:)
  !common /dptr95/ psicu
  real(c_double), pointer :: psiqu(:,:)
  !common /dptr95/ psiqu
  real(c_double), pointer :: bavpd(:,:)
  !common /dptr95/ bavpd
  real(c_double), pointer :: bavdn(:,:)
  !common /dptr95/ bavdn
  real(c_double), pointer :: psiir(:,:)
  !common /dptr95/ psiir
  real(c_double), pointer :: vderb(:,:)
  !common /dptr95/ vderb
  real(c_double), pointer :: sinn(:,:)
  !common /dptr95/ sinn
  real(c_double), pointer :: tann(:,:)
  !common /dptr95/ tann
  real(c_double), pointer :: ymid(:,:)
  !common /dptr95/ ymid
  real(c_double), pointer :: tau(:,:)
  !common /dptr95/ tau
  real(c_double), pointer :: vptb(:,:)
  !common /dptr95/ vptb
  real(c_double), pointer :: zboun(:,:)
  !common /dptr95/ zboun
  integer, pointer :: idx(:,:)
  !common /dptr95/ idx
  integer, pointer :: imax(:,:)
  !common /dptr95/ imax
  real(c_double), pointer :: dz(:,:)
  !common /dptr95/ dz
  real(c_double), pointer :: pol(:,:)
  !common /dptr95/ pol
  real(c_double), pointer :: solrz(:,:)
  !common /dptr95/ solrz
  real(c_double), pointer :: solzz(:,:)
  !common /dptr95/ solzz
  real(c_double), pointer :: bpolz(:,:), btorz(:,:)  ! (lza,lrzmax)
  !common /dptr95/ bpolz, btorz

  ! YuP: [added Apr/2014] area and volume of a cell associated with each
  !                     (R,Z) point on flux surface, (R,Z)==(solrz,solzz)
  real(c_double), pointer :: ddarea(:,:), ddvol(:,:)  !  (lza,lrzmax)
  !common /dptr95/ ddarea, ddvol

  real(c_double), pointer :: thtab(:,:)
  !common /dptr95/ thtab
  real(c_double), pointer :: z(:,:)
  !common /dptr95/ z
  real(c_double), pointer :: zmid(:,:)
  !common /dptr95/ zmid
  real(c_double), pointer :: bbpsi(:,:)
  !common /dptr95/ bbpsi
  real(c_double), pointer :: consnp(:,:)
  !common /dptr95/ consnp
  real(c_double), pointer :: ptime(:,:)
  !common /dptr95/ ptime
  real(c_double), pointer :: sptzrp(:,:)
  !common /dptr95/ sptzrp
  real(c_double), pointer :: pefld(:,:)
  !common /dptr95/ pefld
  real(c_double), pointer :: rovsp(:,:)
  !common /dptr95/ rovsp
  real(c_double), pointer :: restp(:,:)
  !common /dptr95/ restp
  real(c_double), pointer :: restnp(:,:)
  !common /dptr95/ restnp
  real(c_double), pointer :: vpov(:,:)
  !common /dptr95/ vpov
  real(c_double), pointer :: es(:,:)
  !common /dptr95/ es
  real(c_double), pointer :: bpsi(:,:)
  !common /dptr95/ bpsi
  real(c_double), pointer :: d2bpsi(:,:)
  !common /dptr95/ d2bpsi
  real(c_double), pointer :: d2solrz(:,:)
  !common /dptr95/ d2solrz
  real(c_double), pointer :: d2solzz(:,:)
  !common /dptr95/ d2solzz
  real(c_double), pointer :: d2bpolz(:,:), d2btorz(:,:)
  !common /dptr95/ d2bpolz, d2btorz
  real(c_double), pointer :: d2thtpol(:,:)
  !common /dptr95/ d2thtpol
  real(c_double), pointer :: d2es(:,:)
  !common /dptr95/ d2es
  real(c_double), pointer :: thtpol(:,:)
  !common /dptr95/ thtpol
  real(c_double), pointer :: esfi(:,:)
  !common /dptr95/ esfi
  real(c_double), pointer :: psiesfi(:,:)
  !common /dptr95/ psiesfi
  real(c_double), pointer :: psifi(:,:)
  !common /dptr95/ psifi
  real(c_double), pointer :: espsifi(:,:)
  !common /dptr95/ espsifi
  real(c_double), pointer :: soupp(:,:)
  !common /dptr95/ soupp
  real(c_double), pointer :: waa(:,:,:)
  !common /dptr95/ waa
  real(c_double), pointer :: wbb(:,:,:)
  !common /dptr95/ wbb
  real(c_double), pointer :: cosz(:,:,:)
  !common /dptr95/ cosz
  real(c_double), pointer :: dtau(:,:,:)
  !common /dptr95/ dtau
  real(c_double), pointer :: sinz(:,:,:)
  !common /dptr95/ sinz
  real(c_double), pointer :: tanz(:,:,:)
  !common /dptr95/ tanz
  real(c_double), pointer :: yz(:,:,:)
  !common /dptr95/ yz
  real(c_double), pointer :: tot(:,:,:)
  !common /dptr95/ tot
  real(c_double), pointer :: vflux(:,:,:)
  !common /dptr95/ vflux
  real(c_double), pointer :: f_aveth(:,:,:,:)
  !common /dptr95/ f_aveth
  real(c_double), pointer :: sincosba(:,:,:)
  !common /dptr95/ sincosba
  real(c_double), pointer :: densz(:,:,:,:)
  !common /dptr95/ densz
  real(c_double), pointer :: ss(:,:,:,:)
  !common /dptr95/ ss
  real(c_double), pointer :: dcofleg(:,:,:,:)
  !common /dptr95/ dcofleg
  real(c_double), pointer :: dpcosz(:,:,:,:)
  !common /dptr95/ dpcosz
  real(c_double), pointer :: ssy(:,:,:,:)
  !common /dptr95/ ssy
  real(c_double), pointer :: ssyy(:,:,:,:)
  !common /dptr95/ ssyy
  real(c_double), pointer :: ssyi(:,:,:,:)
  !common /dptr95/ ssyi
  real(c_double), pointer :: ssyyy(:,:,:,:)
  !common /dptr95/ ssyyy
  real(c_double), pointer :: pcurr(:,:,:), pcurrm(:,:,:)
  !common /dptr95/ pcurr, pcurrm
  real(c_double), pointer :: pdens(:,:,:), pdenm(:,:,:)
  !common /dptr95/ pdens, pdenm
  real(c_double), pointer :: pengy(:,:,:), pengym(:,:,:)
  !common /dptr95/ pengy, pengym
  real(c_double), pointer :: pdenra(:,:)
  !common /dptr95/ pdenra
  real(c_double), pointer :: pcurra(:,:)
  !common /dptr95/ pcurra
  real(c_double), pointer :: pfdenra(:,:)
  !common /dptr95/ pfdenra
  real(c_double), pointer :: pfcurra(:,:)
  !common /dptr95/ pfcurra
  real(c_double), pointer :: pucrit(:,:)
  !common /dptr95/ pucrit
  real(c_double), pointer :: peoe0(:,:)
  !common /dptr95/ peoe0
  real(c_double), pointer :: psrc(:,:)
  !common /dptr95/ psrc
  real(c_double), pointer :: peoed(:,:)
  !common /dptr95/ peoed
  real(c_double), pointer :: cint2(:)
  !common /dptr95/ cint2
  real(c_double), pointer :: dx(:),dxi(:)
  !common /dptr95/ dx,dxi
  integer, pointer :: ifp(:)
  !common /dptr95/ ifp
  real(c_double), pointer :: sg(:)
  !common /dptr95/ sg
  real(c_double), pointer :: sgx(:)
  !common /dptr95/ sgx
  real(c_double), pointer :: sgxx(:)
  !common /dptr95/ sgxx
  real(c_double), pointer :: sh(:)
  !common /dptr95/ sh
  real(c_double), pointer :: shx(:)
  !common /dptr95/ shx
  real(c_double), pointer :: shxx(:)
  !common /dptr95/ shxx
  real(c_double), pointer :: shxxx(:)
  !common /dptr95/ shxxx
  real(c_double), pointer :: tam1(:)
  !common /dptr95/ tam1
  real(c_double), pointer :: tam2(:)
  !common /dptr95/ tam2
  real(c_double), pointer :: tam3(:)
  !common /dptr95/ tam3
  real(c_double), pointer :: tam4(:)
  !common /dptr95/ tam4
  real(c_double), pointer :: tam5(:)
  !common /dptr95/ tam5
  real(c_double), pointer :: tam6(:)
  !common /dptr95/ tam6
  real(c_double), pointer :: tam7(:)
  !common /dptr95/ tam7
  real(c_double), pointer :: tam8(:)
  !common /dptr95/ tam8
  real(c_double), pointer :: tam9(:)
  !common /dptr95/ tam9
  real(c_double), pointer :: tam10(:)
  !common /dptr95/ tam10
  real(c_double), pointer :: tam11(:)
  !common /dptr95/ tam11
  real(c_double), pointer :: tam12(:)
  !common /dptr95/ tam12
  real(c_double), pointer :: tam13(:)
  !common /dptr95/ tam13
  real(c_double), pointer :: tam14(:)
  !common /dptr95/ tam14
  real(c_double), pointer :: tam15(:)
  !common /dptr95/ tam15
  real(c_double), pointer :: tam16(:)
  !common /dptr95/ tam16
  real(c_double), pointer :: tam17(:)
  !common /dptr95/ tam17
  real(c_double), pointer :: tam18(:)
  !common /dptr95/ tam18
  real(c_double), pointer :: tam19(:)
  !common /dptr95/ tam19
  real(c_double), pointer :: tam20(:)
  !common /dptr95/ tam20
  real(c_double), pointer :: tam21(:)
  !common /dptr95/ tam21
  real(c_double), pointer :: tam22(:)
  !common /dptr95/ tam22
  real(c_double), pointer :: tam23(:)
  !common /dptr95/ tam23
  real(c_double), pointer :: tam24(:)
  !common /dptr95/ tam24
  real(c_double), pointer :: tam25(:)
  !common /dptr95/ tam25
  real(c_double), pointer :: tam26(:)
  !common /dptr95/ tam26
  real(c_double), pointer :: tam27(:)
  !common /dptr95/ tam27
  real(c_double), pointer :: tam28(:)
  !common /dptr95/ tam28
  real(c_double), pointer :: tam29(:)
  !common /dptr95/ tam29
  real(c_double), pointer :: tam30(:)
  !common /dptr95/ tam30
  real(c_double), pointer :: x(:)
  !common /dptr95/ x
  real(c_double), pointer :: xmidpt(:)
  !common /dptr95/ xmidpt
  real(c_double), pointer :: xi(:)
  !common /dptr95/ xi
  real(c_double), pointer :: xsq(:)
  !common /dptr95/ xsq
  real(c_double), pointer :: x3i(:)
  !common /dptr95/ x3i
  real(c_double), pointer :: x2i(:)
  !common /dptr95/ x2i
  real(c_double), pointer :: xcu(:)
  !common /dptr95/ xcu
  real(c_double), pointer :: xcenter(:)
  !common /dptr95/ xcenter
  real(c_double), pointer :: xcensq(:), xcent3(:)
  !common /dptr95/ xcensq, xcent3
  real(c_double), pointer :: uoc(:)
  !common /dptr95/ uoc
  real(c_double), pointer :: enerkev(:,:) !YuP[2018-01-08] added 2nd index (k)
  !common /dptr95/ enerkev
  real(c_double), pointer :: gamma(:)
  !common /dptr95/ gamma
  real(c_double), pointer :: gamsqr(:)
  !common /dptr95/ gamsqr
  real(c_double), pointer :: gamcub(:)
  !common /dptr95/ gamcub
  real(c_double), pointer :: gammi(:)
  !common /dptr95/ gammi
  real(c_double), pointer :: gamm2i(:)
  !common /dptr95/ gamm2i
  real(c_double), pointer :: gamm1(:)
  !common /dptr95/ gamm1
  real(c_double), pointer :: tcsgm1(:)
  !common /dptr95/ tcsgm1
  real(c_double), pointer :: gamefac(:)
  !common /dptr95/ gamefac
  integer, pointer :: ident(:)
  !common /dptr95/ ident
  real(c_double), pointer :: temc1(:)
  !common /dptr95/ temc1
  real(c_double), pointer :: temc2(:)
  !common /dptr95/ temc2
  real(c_double), pointer :: temc3(:)
  !common /dptr95/ temc3
  real(c_double), pointer :: temc4(:)
  !common /dptr95/ temc4
  integer, pointer :: itemc1(:)
  !common /dptr95/ itemc1
  integer, pointer :: itemc2(:)
  !common /dptr95/ itemc2
  integer, pointer :: l_lower(:)
  !common /dptr95/ l_lower
  integer, pointer :: lpt(:)
  !common /dptr95/ lpt
  real(c_double), pointer :: mun(:) !for real?
  !common /dptr95/ mun
  real(c_double), pointer :: fll(:)
  !common /dptr95/ fll
  real(c_double), pointer :: xpar(:)
  !common /dptr95/ xpar
  real(c_double), pointer :: rheads(:)
  !common /dptr95/ rheads
  real(c_double), pointer :: dfvlle(:)
  !common /dptr95/ dfvlle
  real(c_double), pointer :: dfvlli(:)
  !common /dptr95/ dfvlli
  real(c_double), pointer :: xperp(:)
  !common /dptr95/ xperp
  real(c_double), pointer :: xl(:)
  !common /dptr95/ xl
  integer, pointer :: jmaxxl(:)
  !common /dptr95/ jmaxxl
  real(c_double), pointer :: xlm(:)
  !common /dptr95/ xlm
  real(c_double), pointer :: dxl(:)
  !common /dptr95/ dxl
  real(c_double), pointer :: fl(:)
  !common /dptr95/ fl
  real(c_double), pointer :: fl1(:)
  !common /dptr95/ fl1
  real(c_double), pointer :: fl2(:)
  !common /dptr95/ fl2
  real(c_double), pointer :: ppars(:,:)
  !common /dptr95/ ppars
  real(c_double), pointer :: pprps(:,:)
  !common /dptr95/ pprps
  real(c_double), pointer :: faci(:,:)
  !common /dptr95/ faci
  real(c_double), pointer :: pparea(:,:)
  !common /dptr95/ pparea
  real(c_double), pointer :: wtfl0(:,:,:)
  !common /dptr95/ wtfl0
  real(c_double), pointer :: wtflm(:,:,:)
  !common /dptr95/ wtflm
  integer, pointer :: jflbin(:,:,:)
  !common /dptr95/ jflbin
  real(c_double), pointer :: xm(:,:)
  !common /dptr95/ xm
  real(c_double), pointer :: dbb(:,:)
  !common /dptr95/ dbb
  real(c_double), pointer :: dd(:,:)
  !common /dptr95/ dd
  real(c_double), pointer :: de(:,:)
  !common /dptr95/ de
  real(c_double), pointer :: df(:,:)
  !common /dptr95/ df
  real(c_double), pointer :: dff(:,:)
  !common /dptr95/ dff
  real(c_double), pointer :: cah(:,:)
  !common /dptr95/ cah
  real(c_double), pointer :: cthta(:,:)
  !common /dptr95/ cthta
  real(c_double), pointer :: gon(:,:)
  !common /dptr95/ gon
  real(c_double), pointer :: so(:,:)
  !common /dptr95/ so
  real(c_double), pointer :: currv(:,:,:)
  !common /dptr95/ currv
  real(c_double), pointer :: currvs(:,:)
  !common /dptr95/ currvs
  real(c_double), pointer :: pwrrf(:,:,:)
  !common /dptr95/ pwrrf
  real(c_double), pointer :: tal(:,:)
  !common /dptr95/ tal
  real(c_double), pointer :: tbl(:,:)
  !common /dptr95/ tbl
  real(c_double), pointer :: tfl(:,:)
  !common /dptr95/ tfl
  real(c_double), pointer :: pwrrfs(:,:,:)
  !common /dptr95/ pwrrfs
  real(c_double), pointer :: pleg(:,:)
  !common /dptr95/ pleg
  real(c_double), pointer :: feta(:,:)
  !common /dptr95/ feta
  real(c_double), pointer :: fetb(:,:)
  !common /dptr95/ fetb
  real(c_double), pointer :: wflux(:,:,:)
  !common /dptr95/ wflux
  !     NB:  rhs set up here for full 3d set of eqns (BH070525)
  real(c_double), pointer :: rhs(:)
  !common /dptr95/ rhs
  real(c_double), pointer :: sovt(:,:,:,:)
  !common /dptr95/ sovt
  real(c_double), pointer :: sigsxr(:,:,:,:)
  !common /dptr95/ sigsxr

  real(c_double), pointer :: pentr(:,:,:,:)  !!!(nonch,ngen,-1:15,lrors)
  !common /dptr95/ pentr

  real(c_double), pointer :: constp(:,:)  !!!(nonch,lrors)
  !common /dptr95/ constp

  real(c_double), pointer :: sigmtt(:,:),sigftt(:,:)  !!!(nonch,4)
  !common /dptr95/ sigmtt,sigftt

  real(c_double), pointer :: sgaint(:,:,:)  !!!(8,ngen,lrors)
  real(c_double), pointer :: entr(:,:,:)    !!!(ngen,-1:15,lrors)
  real(c_double), pointer :: xlndnz(:,:)    !!!(ngen+1,negyrga)
  real(c_double), pointer :: sounor(:,:,:,:)   !!!(ngen,nsoa,lz,lrz)
  !common /dptr95/ sgaint,entr,xlndnz,sounor


  !.......................................................................
  !*****arrays related to relativ=fully option
  !.......................................................................
  real(c_double), pointer :: gamman(:,:)
  !common /dptr95/ gamman
  real(c_double), pointer :: alphan(:,:)
  !common /dptr95/ alphan

  real(c_double), pointer :: asnha(:)
  !common /dptr95/ asnha
  integer, pointer :: item1(:)
  !common /dptr95/ item1
  integer, pointer :: item2(:)
  !common /dptr95/ item2
  integer, pointer :: item3(:)
  !common /dptr95/ item3
  integer, pointer :: item4(:)
  !common /dptr95/ item4
  integer, pointer :: item5(:)
  !common /dptr95/ item5
  integer, pointer :: item6(:)
  !common /dptr95/ item6
  real(c_double), pointer :: dxm5(:)
  !common /dptr95/ dxm5
  real(c_double), pointer :: exm5(:)
  !common /dptr95/ exm5
  real(c_double), pointer :: dxp5(:)
  !common /dptr95/ dxp5
  real(c_double), pointer :: exp5(:)
  !common /dptr95/ exp5
  real(c_double), pointer :: tamt1(:,:,:,:)
  !common /dptr95/ tamt1
  real(c_double), pointer :: tamt2(:,:,:,:)
  !common /dptr95/ tamt2
  real(c_double), pointer :: da(:,:)
  !common /dptr95/ da
  real(c_double), pointer :: db(:,:)
  !common /dptr95/ db
  real(c_double), pointer :: dc(:,:)
  !common /dptr95/ dc
  real(c_double), pointer :: ca(:,:)
  !common /dptr95/ ca
  real(c_double), pointer :: cb(:,:)
  !common /dptr95/ cb
  real(c_double), pointer :: cc(:,:)
  !common /dptr95/ cc
  real(c_double), pointer :: cd(:,:)
  !common /dptr95/ cd
  real(c_double), pointer :: ce(:,:)
  !common /dptr95/ ce
  real(c_double), pointer :: cf(:,:)
  !common /dptr95/ cf

  real(c_double), pointer :: tem1(:)
  !common /dptr95/ tem1
  real(c_double), pointer :: tem2(:)
  !common /dptr95/ tem2
  real(c_double), pointer :: tem3(:)
  !common /dptr95/ tem3
  real(c_double), pointer :: tem4(:)
  !common /dptr95/ tem4
  real(c_double), pointer :: tem5(:)
  !common /dptr95/ tem5
  real(c_double), pointer :: tem6(:)
  !common /dptr95/ tem6

  real(c_double), pointer :: egg(:,:)
  !common /dptr95/ egg
  real(c_double), pointer :: fgg(:,:)
  !common /dptr95/ fgg

  real(c_double), pointer :: xhead(:,:)
  !common /dptr95/ xhead
  real(c_double), pointer :: xtail(:,:)
  !common /dptr95/ xtail
  real(c_double), pointer :: ytail(:,:)
  !common /dptr95/ ytail
  real(c_double), pointer :: yhead(:,:)
  !common /dptr95/ yhead

  real(c_double), pointer :: fpn(:,:)
  !common /dptr95/ fpn

  real(c_double), pointer :: temp1(:,:)
  !common /dptr95/ temp1
  real(c_double), pointer :: temp2(:,:)
  !common /dptr95/ temp2
  real(c_double), pointer :: temp3(:,:)
  !common /dptr95/ temp3
  real(c_double), pointer :: temp4(:,:)
  !common /dptr95/ temp4
  real(c_double), pointer :: temp5(:,:)
  !common /dptr95/ temp5
  real(c_double), pointer :: temp6(:,:)
  !common /dptr95/ temp6

  real(c_double), pointer :: xllji(:,:)
  !common /dptr95/ xllji
  real(c_double), pointer :: xppji(:,:)
  !common /dptr95/ xppji

  !     Arrays used for first order orbit width calculations:
  real(c_double), pointer :: deltarho(:,:,:)
  !common /dptr95/ deltarho
  real(c_double), pointer :: deltarhop(:,:,:)
  !common /dptr95/ deltarhop
  real(c_double), pointer :: deltarz(:,:,:)
  !common /dptr95/ deltarz
  real(c_double), pointer :: r_delta(:)
  !common /dptr95/ r_delta
  real(c_double), pointer :: z_delta(:)
  !common /dptr95/ z_delta
  real(c_double), pointer :: t_delta(:)
  !common /dptr95/ t_delta
  real(c_double), pointer :: delta_bdb0(:,:)
  !common /dptr95/ delta_bdb0


  !*****************************************************************
  !     BEGIN arrays for analytic ion source (sou..) routines
  !*****************************************************************

  !common /diskx/ &
  real(c_double) :: bdre(lrza),bdrep(lrza), &
       sorpwt(lrza),sorpwti(0:lrza), &
       sorpw_nbii(1:ngena,0:lrza),sorpw_rfi(1:ngena,0:lrza), &
       xlncurt(lrza)

  !common /diskx/ &
  real(c_double) :: sorpw_rf(1:ngena,lrza),sorpw_nbi(1:ngena,lrza)

  !common /diskx/ &
  real(c_double) :: cosm1(ngena,nsoa,lrza),cosm2(ngena,nsoa,lrza), &
       sxllm1(ngena,nsoa,lrza),sxllm2(ngena,nsoa,lrza), &
       sxppm1(ngena,nsoa,lrza), &
       sxppm2(ngena,nsoa,lrza), &
       xem1(ngena,nsoa,lrza),xem2(ngena,nsoa,lrza), &
       zm1(ngena,nsoa,lrza),zm2(ngena,nsoa,lrza)


  !common /scalar/ &
  integer :: isounor


  !*****************************************************************
  !     BEGIN arrays for rf package..(rf...,vlh[B,...,vlf...) routines
  !*****************************************************************

  real(c_double), pointer :: cqlb(:,:,:,:),cqlc(:,:,:,:),cqle(:,:,:,:),cqlf(:,:,:,:)
  !common/qlcoef/cqlb,cqlc,cqle,cqlf 

  real(c_double), pointer :: bqlm(:,:)  ! (iy,jx)
  !common/qlcoef/ bqlm



  !****************************************************************
  !     BEGIN arrays for 3-d (td..) driver.
  !****************************************************************


  !..................................................................
  !     scalars used in CQL3D
  !..................................................................

  real(c_double) :: li
  !common/sc3d/ &
  real(c_double) :: dttr,dtreff,currtza,currtpza,conserv, &
       fom,fomp,fompla,fomtot,flxout, &
       eden,edenlavg,etemp,ethtemp,edntmp,pden,pdntmp,psynct, &
  !BH110314cBH070408(Not used except in fr routines):    1  smooth,
  !BH110314:  Restored, as needed in subroutine frsmooth, but needed
  !BH110314:  to change name, as conflict arises in tdreadf/tdwritef &
       smooth_, &
       toteqd,cursign,totcurza,total,total0,totcurtt,curxjtot
  integer :: ncount,iplt3d,ipacktp,n_d_rr

  !..................................................................
  !     all arrays used only in CQL3D
  !..................................................................

  ! maybe bug
  real(c_double) :: jparb(lrza),jparbt(lrza),jparbp(lrza)
  !common/ar3d/ 
  real(c_double) :: rrz(0:lrza), &
       tr(0:lrza),tr1(0:lrza),tr2(0:lrza), &
       tr3(0:lrza),tr4(0:lrza),tr5(0:lrza),drp5(0:lrza), &
       dpsi(0:lrza),dpsidrho(lrza)
  integer :: iytr(lrorsa)
  real(c_double) ::h_r(0:lrza), &
       area(0:lrza),equilpsp(0:lrza),equilpsi(0:lrza), &
       areamid(0:lrza),volmid(0:lrza), &
       psivalm(lrza),rpconz(lrza),rmconz(lrza),rpmconz(0:lrza), &
       bpolsqaz(0:lrza),aspin(lrza),trapfrac(lrza), &
       currz(ngena,lrza),currtpz(lrza),currtz(lrza), &
       currza(ngena),currtzi(0:lrza),currtpzi(0:lrza), &
       currmtz(lrorsa),currmtpz(lrorsa), &
       totcurzi(0:lrza),totcurz(lrza), &
       fpsiz2(0:lrza),fpsiz(lrza),ffpsiz(lrza)
  real(c_double) :: prestp(lrza),prest(lrza),d2prest(lrza),d2fpsiz(lrza), &
       d2ffpsiz(lrza), &
       bmdplne(lrza),d2bmdpl(lrza), &
       gkpwrz(ngena,lrza), &
       rfpwrz(ngena,lrza), &
       drrt(ngena),drt(ngena)
  !common/ar3d/ &
  real(c_double) :: dvol(lrza),darea(lrza), &
       psyncz(lrza),pegyz(ngena,lrza),pplossz(ngena,lrza), &
       wparzt(ngena),wperpzt(ngena),pegyt(ngena),pplosst(ngena), &
       rfpwrt(ngena), &
       gkpwrt(ngena),energyt(ntotala), &
       rz(0:lrza), &
       vfluxz(lrorsa),vol(0:lrza), &
       onovrpz(lrza,2), &
       tplt3d(nplota), &
       bscurma(2,2),bscurm(0:lrza,2,2),bscurmi(0:lrza,2,2)

  !common/sc3d/ &
  real(c_double) ::     sorpwtza

  !common /csxr/
  real(c_double) :: sxry(lrza,4),sang(lrza,4),spol(lrza,4)
  integer :: ibin(lrza,4)
  real(c_double) ::eflux(nena,nva),efluxt(nva),alphad(3),xs_(3), enk(nena),en_(nena)
  integer :: jval_(nena),inegsxr(nva),lensxr(nva)

  !common /csigma/
  integer :: mtab,msig,jxis
  real(c_double) :: elmin,delegy
  integer :: imaxwln(2,4),igenrl(2,4)
  real(c_double) :: sigm(4,lrorsa),sigf(4,lrorsa),sigmt(4),sigft(4), &
       fuspwrv(4,lrorsa),fuspwrvt(4),fuspwrm(4,lrorsa),fuspwrmt(4)

  real(c_double), pointer :: tamm1(:)
  !common /csigma/ tamm1  !(0:mmsv)

  real(c_double), pointer :: iind(:)
  !common /csigma/ iind  !(1:jx)

  !..............................................................
  !     Set up real(c_double), pointers
  !..............................................................

  real(c_double), pointer :: csv(:,:,:)
  !common /dptr95/ csv
  real(c_double), pointer :: svtab(:)
  !common /dptr95/ svtab


  !..............................................................
  !     Set up pointers to allocatable arrays for transport model.
  !     Space allocated in subroutine tdtraloc
  !..............................................................


  real(c_double), pointer :: frn_2(:,:,:,:)
  !common /dptr95/ frn_2
  real(c_double), pointer :: frn_1(:,:,:,:)
  !common /dptr95/ frn_1
  real(c_double), pointer :: frn(:,:,:,:)
  !common /dptr95/ frn
  real(c_double), pointer :: fvn_1(:,:,:,:)
  !common /dptr95/ fvn_1
  real(c_double), pointer :: fvn(:,:,:,:)
  !common /dptr95/ fvn
  real(c_double), pointer :: dl(:,:,:,:)
  !common /dptr95/ dl
  real(c_double), pointer :: d_rr(:,:,:,:)
  !common /dptr95/ d_rr
  real(c_double), pointer :: d_r(:,:,:,:)
  !common /dptr95/ d_r
  real(c_double), pointer :: f_lm(:,:,:)
  !common /dptr95/ f_lm
  real(c_double), pointer :: f_lp(:,:,:)
  !common /dptr95/ f_lp
  real(c_double), pointer :: f_up(:,:,:)
  !common /dptr95/ f_up
  real(c_double), pointer :: f_vtor(:,:,:,:)
  !common /dptr95/ f_vtor
  real(c_double), pointer :: cynt2_(:,:)
  !common /dptr95/ cynt2_
  real(c_double), pointer :: vpint_(:,:)
  !common /dptr95/ vpint_
  real(c_double), pointer :: vptb_(:,:)
  !common /dptr95/ vptb_
  real(c_double), pointer :: cosovb(:,:)
  !common /dptr95/ cosovb
  real(c_double), pointer :: bovcos(:,:)
  !common /dptr95/ bovcos
  real(c_double), pointer :: adv(:,:)
  !common /dptr95/ adv
  real(c_double), pointer :: dentarget(:)
  !common /dptr95/ dentarget
  real(c_double), pointer :: eg_(:,:,:)
  !common /dptr95/ eg_
  real(c_double), pointer :: fg_(:,:,:)
  !common /dptr95/ fg_



  !******************************************************************
  !     BEGIN arrays for EQUILIBRIUM MODEL (eq..) (NON-CIRCULAR CROSS
  !     SECTIONS).
  !******************************************************************


  !common /diskx/ &
  real(c_double) :: areacon(lrza), &
       bmidplne(lrza),bpolsqa(lrza), &
       eqdells(lrza),epsicon(lrza),erhocon(lrza), &
       fpsi(lrza),flxavgd(lrza), &
       psiovr(lrza),psiavg(2,lrza),onovrp(2,lrza),onovpsir3(lrza), &
       rpcon(lrza),rmcon(lrza),zpcon(lrza),zmcon(lrza), &
       volcon(lrza),fppsi(lrza),pppsi(lrza),es_bmax(lrza), &
       bpsi_max(lrza),bpsi_min(lrza),lbpsi_max(lrza), &
       z_bmax(lrza),bpsi_z_bmax(lrza), &
       dlpgpsii(lrza),dlpsii(lrza)
  integer :: lz_bmax(lrza),lbpsi_min(lrza)

  !common/params/
  integer :: nnz,nnr,nj12


  character(len=8) :: eqorb,eqcall

  !common &
  real(c_double) :: bpolsqlm
  integer :: imag,jmag
  integer :: nmag,nrc,nzc,nfp,nnv,iupdn
  real(c_double) ::psimag,psilim, &
       rmaxcon,rmincon,rhomax, &
       zmag,zmaxcon,zmincon,zshift

  !common &
  integer :: ibd(4)
  real(c_double) :: eqpsi(nconteqa),eqvol(nconteqa),eqfopsi(nconteqa), &
       q_(nconteqa),eqrho(nconteqa),d2eqrho(nconteqa),eqarea(nconteqa), &
       eqrpcon(nconteqa),eqrmcon(nconteqa), &
       eqzpcon(nconteqa),eqzmcon(nconteqa), &
       d2fpsiar(nnra), &
       ez(nnza),dummyaz(nnza), &
       fpsiar(nnra),ffpar(nnra),d2ffpar(nnra),qar(nnra),d2qar(nnra), &
       prar(nnra),d2prar(nnra),ppar(nnra),d2ppar(nnra),psiar(nnra), &
       er(nnra),dummyar(nnra), &
       wkepsi(nrz3p1a), &
       tlorb1(lfielda),tlorb2(lfielda)

  !common &
  real(c_double) :: epsi(nnra,nnza),epsirr(nnra,nnza), &
       epsizz(nnra,nnza),epsirz(nnra,nnza), &
       dummypsi(nnra,nnza),eqovrp(nconteqa,2)

  !common/output/ 
  integer :: lorbit_,ialign14
  real(c_double) :: rmcon_,rpcon_,zmcon_,zpcon_, &
       bthr_,btoru_,eqdells_,fpsi_,fppsi_,zmax_,btor0_,bthr0_, &
       es_bmax_,bpsi_max_,bpsi_min_
  integer :: lbpsi_max_,lbpsi_min_
  real(c_double) :: bmidplne_,solr_(lfielda),solz_(lfielda),es_(lfielda), &
       eqbpol_(lfielda),bpsi_(lfielda),thtpol_(lfielda), &
       eqdell_(lfielda)

  !..................................................................
  !     Allocatable arrays allocated in subroutine eqalloc
  !..................................................................

  real(c_double), pointer :: drpmconz(:)
  !common /dptr95/ drpmconz
  real(c_double), pointer :: eqdell(:,:)
  !common /dptr95/ eqdell
  real(c_double), pointer :: eqbpol(:,:)
  !common /dptr95/ eqbpol
  real(c_double), pointer :: solr(:,:)
  !common /dptr95/ solr
  real(c_double), pointer :: solz(:,:)
  !common /dptr95/ solz



  !*********************************************************************
  !     BEGIN arrays for LOWER HYBRID FAST WAVE and ECH Module.
  !*********************************************************************


  !common/params/ &
  integer :: jjx

  integer, pointer :: jbm1(:,:)
  !common jbm1
  integer, pointer :: jb0(:,:)
  !common jb0
  integer, pointer :: jbp1(:,:)
  !common jbp1

  ! bug
  character(len=8) :: irffile(nmodsa)

  !common &
  real(c_double) :: argmax, nurf
  !     1  lenj0, &
  integer :: mrf,mrfn,  &
       nrayn,nrayelts, & !-YuP 101122: added &
       irftype
  real(c_double) :: powray, &
       vnorm2,vnorm3,vnorm4, &
       dveps ! YuP[04-2016] for subr. urfb_add

  complex(c_double_complex),pointer :: cosz1(:),sinz1(:),sinz2(:)
  !common &
  real(c_double) :: bsslstp(nmodsa)
  integer  :: ncontrib(lrza)
  real(c_double) :: powrf(lrza,nmodsa),powrfc(lrza,nmodsa), &
       powrfl(lrza,nmodsa),powrft(lrza), &
       powurf(0:nmodsa),powurfc(0:nmodsa),powurfl(0:nmodsa), &
       powurfi(0:lrza,0:nmodsa)

  !common &
  real(c_double) :: freqcy(nmodsa),omega(nmodsa)
  integer :: nharm(nmodsa),nray(nmodsa)
  real(c_double) :: bsign1(nmodsa)
  integer :: krfn(nmodsa),irfn(nmodsa),irfm(nmodsa)

  !..................................................................
  !     Allocatable arrays allocated in subroutine urfalloc
  !..................................................................


  real(c_double), pointer :: urfb(:,:,:,:)
  !common /dptr95/ urfb
  real(c_double), pointer :: urfc(:,:,:,:)
  !common /dptr95/ urfc
  real(c_double), pointer :: cosmz(:,:,:)
  !common /dptr95/ cosmz
  real(c_double), pointer :: g_(:,:,:,:)
  !common /dptr95/ g_
  real(c_double), pointer :: alfag(:)
  !common /dptr95/ alfag
  real(c_double), pointer :: argmnt(:)
  !common /dptr95/ argmnt
  integer, pointer :: ilim1d(:)
  !common /dptr95/ ilim1d
  integer, pointer :: ilim2d(:)
  !common /dptr95/ ilim2d
  integer, pointer :: ilim1dd(:)
  !common /dptr95/ ilim1dd
  integer, pointer :: ilim2dd(:)
  !common /dptr95/ ilim2dd
  real(c_double), pointer :: sx(:)
  !common /dptr95/ sx
  real(c_double), pointer :: xmdx(:)
  !common /dptr95/ xmdx
  real(c_double), pointer :: thtf1(:)
  !common /dptr95/ thtf1
  real(c_double), pointer :: thtf2(:)
  !common /dptr95/ thtf2
  real(c_double), pointer :: alfi(:)
  !common /dptr95/ alfi
  real(c_double), pointer :: alfa(:)
  !common /dptr95/ alfa
  integer, pointer :: ilim1(:)
  !common /dptr95/ ilim1
  integer, pointer :: ilim2(:)
  !common /dptr95/ ilim2
  integer, pointer :: ifct1(:)
  !common /dptr95/ ifct1
  integer, pointer :: ifct2(:)
  !common /dptr95/ ifct2
  real(c_double), pointer :: urftmp(:)
  !common /dptr95/ urftmp
  real(c_double), pointer :: urfpwr(:,:,:)
  !common /dptr95/ urfpwr
  real(c_double), pointer :: urfpwrc(:,:,:)
  !common /dptr95/ urfpwrc
  real(c_double), pointer :: urfpwrl(:,:,:)
  !common /dptr95/ urfpwrl
  integer, pointer :: jminray(:,:,:)
  !common /dptr95/ jminray
  integer, pointer :: jmaxray(:,:,:)
  !common /dptr95/ jmaxray
  integer, pointer :: lloc(:,:,:)
  !common /dptr95/ lloc
  integer, pointer :: llray(:,:,:)
  !common /dptr95/ llray
  integer, pointer :: psiloc(:,:,:)
  !common /dptr95/ psiloc
  real(c_double), pointer :: scalurf(:,:,:)
  !common /dptr95/ scalurf
  complex(c_double_complex), pointer :: cwexde(:,:,:), cweyde(:,:,:), cwezde(:,:,:)
  !common /dptr95/ cwezde
  real(c_double), pointer :: delpwr(:,:,:)
  !common /dptr95/ delpwr
  real(c_double), pointer :: fluxn(:,:,:)
  !common /dptr95/ fluxn
  real(c_double), pointer :: seikon(:,:,:)
  !common /dptr95/ seikon
  real(c_double), pointer :: spsi(:,:,:)
  !common /dptr95/ spsi
  real(c_double), pointer :: sdpwr(:,:,:)
  !common /dptr95/ sdpwr
  real(c_double), pointer :: sbtot(:,:,:)
  !common /dptr95/ sbtot
  real(c_double), pointer :: sene(:,:,:)
  !common /dptr95/ sene
  real(c_double), pointer :: salphac(:,:,:)
  !common /dptr95/ salphac
  real(c_double), pointer :: salphal(:,:,:)
  !common /dptr95/ salphal
  real(c_double), pointer :: ws(:,:,:)
  !common /dptr95/ ws
  real(c_double), pointer :: wr(:,:,:)
  !common /dptr95/ wr
  real(c_double), pointer :: wz(:,:,:)
  !common /dptr95/ wz
  real(c_double), pointer :: wnpar(:,:,:)
  !common /dptr95/ wnpar
  real(c_double), pointer :: wdnpar(:,:,:)
  !common /dptr95/ wdnpar
  real(c_double), pointer :: wnper(:,:,:)
  !common /dptr95/ wnper
  real(c_double), pointer :: wphi(:,:,:)
  !common /dptr95/ wphi
  integer, pointer :: ilowp(:,:)
  !common /dptr95/ ilowp
  integer, pointer :: iupp(:,:)
  !common /dptr95/ iupp
  integer, pointer :: ifct1_(:,:)
  !common /dptr95/ ifct1_
  integer, pointer :: ifct2_(:,:)
  !common /dptr95/ ifct2_
  integer, pointer :: nrayelt(:,:)
  !common /dptr95/ nrayelt
  integer, pointer :: jslofas(:,:)
  !common /dptr95/ jslofas
  real(c_double), pointer :: nurefls(:,:)
  !common /dptr95/ nurefls
  real(c_double), pointer :: keiks(:,:)
  !common /dptr95/ keiks
  integer, pointer :: jpes(:,:)
  !common /dptr95/ jpes
  integer, pointer :: jpis(:,:)
  !common /dptr95/ jpis
  integer, pointer :: istarts(:,:)
  !common /dptr95/ istarts
  integer, pointer :: iprmt5(:,:)
  !common /dptr95/ iprmt5
  integer, pointer :: jhlfs(:,:)
  !common /dptr95/ jhlfs
  real(c_double), pointer :: sxxrt(:,:)
  !common /dptr95/ sxxrt
  real(c_double), pointer :: skpsi(:,:)
  !common /dptr95/ skpsi
  real(c_double), pointer :: skth(:,:)
  !common /dptr95/ skth
  real(c_double), pointer :: skphi(:,:)
  !common /dptr95/ skphi
  integer, pointer :: lrayelt(:,:)
  !common /dptr95/ lrayelt
  real(c_double), pointer :: delpwr0(:,:)
  !common /dptr95/ delpwr0
  real(c_double), pointer :: nrayelt0(:,:)
  !common /dptr95/ nrayelt0
  real(c_double), pointer :: truncd(:) ! 1:jx
  !common /dptr95/ truncd


  !..................................................................
  !     Allocatable arrays allocated in subroutine rdc_multi,
  !     used after subroutine execution.
  !     Here, we introduce f90 pointers, as they are easier
  !     to allocate.

  real(c_double), pointer :: rdcb(:,:,:,:)
  !common /dptr95/ rdcb
  real(c_double), pointer :: rdcc(:,:,:,:)
  !common /dptr95/ rdcc
  real(c_double), pointer :: rdce(:,:,:,:)
  !common /dptr95/ rdce
  real(c_double), pointer :: rdcf(:,:,:,:)
  !common /dptr95/ rdcf


  !..................................................................
  !     Allocatable arrays allocated in subroutine it3dalloc
  !     used after subroutine execution.
  !..................................................................

  !common /it3d/
  integer :: lapacki,lapackj,icsrij,icsrip,icsri2,krylov, &
       icsrikry,iwk_ilu
  !common / it3d/
  integer ::icsrijr,icsrijc

  real(c_double), pointer :: abd_lapack,a_csr,alu,w_ilu,rhs0,sol,vv
  integer, pointer :: ja_csr,ia_csr,jlu,ju,jw_ilu
  real(c_double), pointer :: ar_csr,ac_csr
  integer, pointer :: jar_csr,iar_csr,ipofi,jac_csr,iac_csr
  dimension abd_lapack(:,:),a_csr(:), alu(:),w_ilu(:),rhs0(:),sol(:),vv(:)
  dimension ja_csr(:),ia_csr(:),jlu(:),ju(:),jw_ilu(:)
  dimension ar_csr(:),ac_csr(:)
  dimension jar_csr(:),iar_csr(:),ipofi(:,:),jac_csr(:),iac_csr(:)
  !common /dptr95/ abd_lapack,a_csr,alu,w_ilu,rhs0,sol,vv
  !common /iptr95/ ja_csr,ia_csr,jlu,ju,jw_ilu
  !common /dptr95/ ar_csr,ac_csr
  !common /iptr95/ jar_csr,iar_csr,ipofi,jac_csr,iac_csr


  !..................................................................

  !common
  real(c_double) :: anecc(nrada),tekev(nrada),tikev(nint1a,nrada), &
       anicc(nint1a,nrada), &
       amass(nint1a),achrg(nint1a),elecf(nrada), &
       rho_(nrada),names(10),psiar_(nrada)
  integer :: nspc,ialign9
  real(c_double) :: fpsiar_(nrada),pary(nrada),ppary(nrada), &
       gpary(nrada),ztop_,zbot_,rleft,rright
  integer :: nx,nz,npsitm

  !-----------------------------------------------------------------------
  !     BEGIN variables for WP... modules for CQLP case
  !-----------------------------------------------------------------------

  character(len=8) :: analegco

  !common /wpscal/ &
  integer :: iymax, nsleft,nsrigt,numclas,numindx

  !common /wpvec/ &
  real(c_double) :: cofdfds(0:lsa1,2,ntrmdera), &
       enrgypa(ntotala,0:lsa1), &
       vthpar(ntotala,0:lsa1)
  integer :: lsbtopr(0:lsa1),lsprtob(0:lsa1),lpm1eff(0:lsa1,-1:+1)
  real(c_double) :: sz(0:lsa1),dsz(0:lsa1),dszm5(0:lsa1),dszp5(0:lsa1), &
       eszm5(0:lsa1),eszp5(0:lsa1), &
       psis(0:lsa1),psisp(0:lsa1),psipols(0:lsa1), &
       solrs(0:lsa1),solzs(0:lsa1), &
       elparol(0:lsa1),elparnw(0:lsa1), &
       flux1(0:lsa1),flux2(0:lsa1)

  !.......................................................................
  !     Arrays allocated in subroutine wpalloc for CQLP
  !.......................................................................

  integer, pointer :: l_upper(:)  !!! (1:iy)
  integer, pointer :: ilpm1ef(:,:,:)  !!! (0:iy+1,0:lsa1,-1:+1)

  real(c_double), pointer :: fnhalf(:,:,:,:)
  !common /dptr95/ fnhalf
  real(c_double), pointer :: fnp0(:,:,:,:)
  !common /dptr95/ fnp0
  real(c_double), pointer :: fnp1(:,:,:,:)
  !common /dptr95/ fnp1
  real(c_double), pointer :: dls(:,:,:,:)
  !common /dptr95/ dls
  real(c_double), pointer :: fh(:,:,:,:)
  !common /dptr95/ fh
  real(c_double), pointer :: fg(:,:,:,:)
  !common /dptr95/ fg
  real(c_double), pointer :: fedge(:,:,:,:)
  !common /dptr95/ fedge
  real(c_double), pointer :: rhspar(:,:,:)
  !common /dptr95/ rhspar
  real(c_double), pointer :: bndmats(:,:,:,:)
  !common /dptr95/ bndmats
  real(c_double), pointer :: wcqlb(:,:,:,:)
  !common /dptr95/ wcqlb
  real(c_double), pointer :: wcqlc(:,:,:,:)
  !common /dptr95/ wcqlc
  real(c_double), pointer :: wcqle(:,:,:,:)
  !common /dptr95/ wcqle
  real(c_double), pointer :: wcqlf(:,:,:,:)
  !common /dptr95/ wcqlf



  !.......................................................................
  !     Arrays allocated in ampfalloc
  !.......................................................................

  real(c_double), pointer :: ampfln(:)
  !common /dptr95/ ampfln
  real(c_double), pointer :: ampflh(:)
  !common /dptr95/ ampflh
  real(c_double), pointer :: ampflg(:)
  !common /dptr95/ ampflg
  real(c_double), pointer :: ampfa(:,:)
  !common /dptr95/ ampfa
  real(c_double), pointer :: ampfb(:,:)
  !common /dptr95/ ampfb
  real(c_double), pointer :: ampfaa(:,:)
  !common /dptr95/ ampfaa
  real(c_double), pointer :: ampfc(:)
  !common /dptr95/ ampfc
  real(c_double), pointer :: ampf2ebar(:)
  !common /dptr95/ ampf2ebar

  !.......................................................................
  !     Arrays for finite orbit width (FOW) calculations
  !.......................................................................

  !common/psiaxis/ psi_lim,psi_mag,R_axis,Z_axis ![cgs]

  real(c_double), pointer :: rcontr(:),zcontr(:),rlimiter(:),zlimiter(:)
  !common/limiter/
  integer :: ncontr
  integer :: nlimiter
  ! Setup by call equilib()

  !common/eqbox/
  real(c_double) :: ermin,ermax,ezmin,ezmax  ! Limits of equilibrium
  ! (R,Z)-grid (cm);
  ! Setup by call equilib()


  !BH170708: Removing FOW material from ZOW code      
  !$$$C---> Equilibrium B is calc-ed on (R,Z)-grid == (req(ir),zeq(iz))
  !$$$      common/Beq/ ireq,izeq,dreq,dzeq,odr,odz,req(nreqa),zeq(nzeqa), 
  !$$$     +       Beqr(nreqa,nzeqa),Beqz(nreqa,nzeqa),Beqphi(nreqa,nzeqa),
  !$$$     +       Beqmod(nreqa,nzeqa), psieq(nreqa,nzeqa) ! [cgs]
  !$$$C---> dreq=(ermax-ermin)/(nreqa-1); odr=1.d0/dreq
  !$$$C---> This block is used by subroutine gc_vel()
  !$$$C---> The fields at a given point along orbit are calculated  
  !$$$C---> by bilinear interpolation from four nearest nodes.
  !$$$ 
  !$$$      common/border/ iborder(nreqa,nzeqa),
  !$$$     +  Rplasma_min, Rplasma_max, 
  !$$$     +  Bplasma_min, Bplasma_max, 
  !$$$     +  PSIplasma_min, PSIplasma_max ! min/max values within border
  !$$$C---> iborder=0 in plasma; >0 at nodes representing wall.
  !$$$ 
  !$$$C---> Unit vector of eq.field (Beq) and its derivative-dep. functions:
  !$$$      common/BeqGrad/
  !$$$     +  bhri(nreqa,nzeqa),bhzi(nreqa,nzeqa),bhfi(nreqa,nzeqa), 
  !$$$     +  GRr(nreqa,nzeqa), GRz(nreqa,nzeqa), GRf(nreqa,nzeqa),
  !$$$     +  Ori(nreqa,nzeqa), Ozi(nreqa,nzeqa), Ofi(nreqa,nzeqa),
  !$$$     +  GV(nreqa,nzeqa), DRbbf(nreqa,nzeqa)
  !$$$C---> bh== B/|B|;   
  !$$$C---> GR== [Bxgrad|B|]/B^2
  !$$$C---> O == { [Bxgrad|B|]/B^2 + rot(B)/|B| - bhat(B.rotB)/B^2 }/|B| 
  !$$$C---> GV== dreq*{B.grad|B|}/|B| 
  !$$$
  !$$$
  !$$$
  !$$$      pointer cmu,cpfi_min,cpfi_max,dpfi,Rpfimax,cpfi,
  !$$$     +        com_rt1,com_rt2,com_rt3,com_rt4,
  !$$$     +        com_cp1,com_cp2,com_cp3,com_cp4
  !$$$      dimension cmu(:),      !(nmu) Normalized adiab.invariant mu (grid)
  !$$$     +          cpfi_min(:), !(jx) min of Canon. angular tor. momentum
  !$$$     +          cpfi_max(:), !(jx) max of Canon. angular tor. momentum
  !$$$     +          dpfi(:),     !(jx)
  !$$$     +          Rpfimax(:),  !(jx) R for the peak of Pfi on the midplane
  !$$$     +          cpfi(:,:),   !(jx,npfi) Canon. tor. momentum (grid)
  !$$$     +        com_rt1(:,:,:), !Rmidplane as a function of COM, 1st root
  !$$$     +        com_rt2(:,:,:), !Rmidplane as a function of COM, 2nd root
  !$$$     +        com_rt3(:,:,:), !Rmidplane as a function of COM, 3rd root
  !$$$     +        com_rt4(:,:,:), !Rmidplane as a function of COM, 4th root
  !$$$     +        com_cp1(:,:,:), !cosp as a function of COM, 1st root
  !$$$     +        com_cp2(:,:,:), !cosp as a function of COM, 2nd root
  !$$$     +        com_cp3(:,:,:), !cosp as a function of COM, 3rd root
  !$$$     +        com_cp4(:,:,:)  !cosp as a function of COM, 4th root
  !$$$
  !$$$      common/com/dmu,odmu,odpfi,dtcom,Ucom,cosps, cmu_min,cmu_max, cmu,
  !$$$     +           cpfi_min,cpfi_max,dpfi,Rpfimax, cpfi,
  !$$$     +           com_rt1,com_rt2,com_rt3,com_rt4,
  !$$$     +           com_cp1,com_cp2,com_cp3,com_cp4
  !$$$
  !$$$
  !$$$      REAL*8 mcq_drdt
  !$$$      common/orb/ nstp_orb,lostorb,not_complete,
  !$$$     +            cmuorb,vpar_ini,vprp_ini,mcq_drdt,renv,dtorb
  !$$$C---> nstp_orb= number of steps along orbit; found after orbit is traced
  !$$$C---> lostorb=1 if orbit is lost to walls,
  !$$$C---> not_complete=1 if not enough steps to complete the orbit,
  !$$$C---> cmuorb== 0.5(Vprpini^2/B)*(dtorb/dreq)^2, 
  !$$$C---> mcq_drdt== (dreq/dtorb)*mc/q
  !$$$C---> renv== dreq/dtorb
  !$$$C---> Used by subroutine gc_vel()
  !$$$
  !$$$      common/vdrift/ vdrift_r,vdrift_z,vdrift_phi 
  !$$$C---> Drift vel.*dtorb/dreq. Saved from gc_vel()
  !$$$      
  !$$$
  !$$$C---> Values at all steps along orbit:
  !$$$      pointer t_orb,R_orb,Z_orb,phi_orb,
  !$$$     +        upar_orb,uprp_orb,psi_orb,b_orb,bphi_orb ! (0:nsteps_orb)
  !$$$      dimension t_orb(:),     !-> [s]
  !$$$     &          R_orb(:),     !-> [cm]
  !$$$     &          Z_orb(:),     !-> [cm]
  !$$$     &          phi_orb(:),   !-> [rad]
  !$$$     &          upar_orb(:),
  !$$$     &          uprp_orb(:),  !-> [cm/s]
  !$$$     &          psi_orb(:),   !-> Psi/2pi  [cgs]
  !$$$     &          b_orb(:),bphi_orb(:)   !-> [cgs]
  !$$$      common/dptr95/ t_orb,R_orb,Z_orb,phi_orb,
  !$$$     +        upar_orb,uprp_orb,psi_orb,b_orb,bphi_orb
  !$$$
  !$$$
  !$$$C---> Values at selected points along MAIN orbit:
  !$$$      pointer t_orb0,R_orb0,Z_orb0,phi_orb0,
  !$$$     +        upar_orb0,uprp_orb0,psi_orb0,b_orb0,bphi_orb0 ! (1:nptsorb)
  !$$$      dimension t_orb0(:),     !-> [s]
  !$$$     &          R_orb0(:),     !-> [cm]
  !$$$     &          Z_orb0(:),     !-> [cm]
  !$$$     &          phi_orb0(:),   !-> [rad]
  !$$$     &          upar_orb0(:),
  !$$$     &          uprp_orb0(:),  !-> [cm/s]
  !$$$     &          psi_orb0(:),   !-> Psi/2pi  [cgs]
  !$$$     &          b_orb0(:),bphi_orb0(:)  !-> [cgs]
  !$$$      common/dptr95/ t_orb0,R_orb0,Z_orb0,phi_orb0,
  !$$$     +        upar_orb0,uprp_orb0,psi_orb0,b_orb0,bphi_orb0
  !$$$
  !$$$
  !$$$      common/orb0/qmc,Rorb0,Zorb0,borb0,psiorb0,cosp0,Rmid
  !$$$C---> '0' - values at the launching point of g.c. orbit [cgs]
  !$$$C---> cosp0 is cos(pitch-angle) at launching point.
  !$$$C---> Rmid is R at midplane for given orbit found by Newton iterations
  !.......................................................................
  !     variables transferred from freya
  !.......................................................................
  character(len=8) :: frmodp, fr_gyrop, beamplsep
  integer :: mfm1p
  real(c_double) :: beamponp, beampoffp  
  real(c_double) :: hibrzp(kz,ke,kb)  !kz=nconteqa+2, from param.h
  !common /freycomm/ &

  !These variables are set to frmod,fr_gyro,
  !beamplse,beampon,beampoff from frmod namelist.
  !The namelist is declared in frname.h and
  !passed to the comm.h related subroutines
  !as arguments of subroutine frnfreya
  !Similarly, hibrzp and mfmp1 are from freya 
  !routines through frnfreya arguments.
  !Purpose is communication with cql3d.
  !

  ! These are from tdtransp and friends, which can become modules
  ! on your time, I didn't feel like dealing with transp.h
  ! i've done the harder ones.
  !common/nob/ nobind
  character(len=8) :: nobind
  !
  !common /newt_norm/
  real(c_double) :: adv_norm(lrza),reden_norm(lrza)
  !COMMON /newtv/
  integer, parameter :: NP=300
  real(c_double) :: fvec(NP)
  integer :: newtv_nn



  save
contains
  ! lolz
end module comm_mod
