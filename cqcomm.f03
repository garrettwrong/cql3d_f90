!     cqcomm
! XXX need to convert name and name_decl before 
! XXX can get away from f77 here
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

module cqcomm_mod
  use iso_c_binding, only : c_float             !REAL*4
  use iso_c_binding, only : c_double            !REAL*8
  use iso_c_binding, only : c_double_complex    !COMPLEX*16
  use param_mod
  public
  !implicit none
  !      include 'name_decl.h'
  !     name_decl.h
  !
  !.......................................................................
  !     This file will hold all declarations of the namelist variables
  !     (given in name.h), type and dimensions.
  !.......................................................................



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

  character(len=256) :: lossfile

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
  integer :: ipad1,r0geomp, rgeomp, zgeomp, rgeom1, rgeom2, &
       niyjx,nccoef,ncentr,navnc,ncplt,n, &
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

  pointer :: ix1(:),ix2(:),ix3(:),ix4(:),   &  !!!(0:mx) &
       ix5(:),ix6(:),ix7(:),ix8(:)      !!!(0:mx)

  pointer :: tom1(:),tom2(:),tom3(:),tom4(:)  !!!(0:mxp1)

  !BH_YP090809  WAS: choose(0:mx+2,0:mx+2),fctrl(0:mx+2)
  !BH_YP090809  First dimension of choose should be larger of 2*mx,mx+2
  pointer :: fctrl(:)     !!!(0:2*mx+2)
  pointer :: choose(:,:)  !!!(0:2*mx+2,0:mx+2)
  pointer :: cog(:,:)     !!!(0:mx,15)
  pointer :: pm(:,:)      !!!(0:mxp1,lrors)

  !common /dptr95/
  integer :: ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8
  real(c_double) :: tom1,tom2,tom3,tom4,  &
       fctrl,choose,cog,pm

  pointer f
  dimension f(:,:,:,:)  !f(0:iy+1,0:jx+1,ngen,lrors)
  !common /dptr95/ f
  pointer favg
  dimension favg(:,:,:,:)
  !common /dptr95/ favg
  pointer fxsp
  dimension fxsp(:,:,:,:)
  !common /dptr95/ fxsp
  pointer f_
  dimension f_(:,:,:,:)
  !common /dptr95/ f_
  pointer spasou
  dimension spasou(:,:,:,:)
  !common /dptr95/ spasou
  pointer velsou
  dimension velsou(:,:,:,:)
  !common /dptr95/ velsou
  pointer velsou2
  dimension velsou2(:,:,:,:)
  !common /dptr95/ velsou2
  pointer source
  dimension source(:,:,:,:)
  !common /dptr95/ source
  pointer gone
  dimension gone(:,:,:,:)
  !common /dptr95/ gone
  pointer egylosa
  dimension egylosa(:,:,:,:)
  !common /dptr95/ egylosa
  pointer i0tran
  dimension i0tran(:,:,:)
  !common /dptr95/ i0tran
  pointer cal
  dimension cal(:,:,:,:)
  !common /dptr95/ cal
  pointer cbl
  dimension cbl(:,:,:,:)
  !common /dptr95/ cbl
  pointer ccl
  dimension ccl(:,:,:,:)
  !common /dptr95/ ccl
  pointer cdl
  dimension cdl(:,:,:,:)
  !common /dptr95/ cdl
  pointer cel
  dimension cel(:,:,:,:)
  !common /dptr95/ cel
  pointer cfl
  dimension cfl(:,:,:,:)
  !common /dptr95/ cfl
  pointer eal
  dimension eal(:,:,:,:,:)
  !common /dptr95/ eal
  pointer ebl
  dimension ebl(:,:,:,:,:)
  !common /dptr95/ ebl
  pointer scal
  dimension scal(:,:)
  !common /dptr95/ scal
  pointer cet
  dimension cet(:,:,:)
  !common /dptr95/ cet
  pointer cex
  dimension cex(:,:,:)
  !common /dptr95/ cex
  pointer synca
  dimension synca(:,:,:)
  !common /dptr95/ synca
  pointer syncd
  dimension syncd(:,:,:)
  !common /dptr95/ syncd
  pointer taulos
  dimension taulos(:,:,:)
  !common /dptr95/ taulos
  pointer elecfldn
  dimension elecfldn(:,:,:)
  !common /dptr95/ elecfldn
  pointer delecfld0n
  dimension delecfld0n(:,:,:)
  !common /dptr95/ delecfld0n
  pointer elecn
  dimension elecn(:,:,:)
  !common /dptr95/ elecn
  pointer delecfld0
  dimension delecfld0(:,:)
  !common /dptr95/ delecfld0
  pointer psi0bar
  dimension psi0bar(:,:)
  !common /dptr95/ psi0bar
  pointer di
  dimension di(:,:,:,:)
  !common /dptr95/ di
  pointer dj
  dimension dj(:,:,:,:)
  !common /dptr95/ dj
  pointer dym5
  dimension dym5(:,:)
  !common /dptr95/ dym5
  pointer dyp5
  dimension dyp5(:,:)
  !common /dptr95/ dyp5
  pointer eyp5
  dimension eyp5(:,:)
  !common /dptr95/ eyp5
  pointer eym5
  dimension eym5(:,:)
  !common /dptr95/ eym5
  pointer y
  dimension y(:,:)
  !common /dptr95/ y
  pointer dy,dyi
  dimension dy(:,:),dyi(:,:)
  !common /dptr95/ dy,dyi
  pointer yptb
  dimension yptb(:,:)
  !common /dptr95/ yptb
  pointer coss
  dimension coss(:,:)
  !common /dptr95/ coss
  pointer cynt2
  dimension cynt2(:,:)
  !common /dptr95/ cynt2
  pointer batot
  dimension batot(:,:)
  !common /dptr95/ batot
  pointer lmax
  dimension lmax(:,:)
  !common /dptr95/ lmax
  pointer vpint
  dimension vpint(:,:)
  !common /dptr95/ vpint
  pointer psiiv
  dimension psiiv(:,:)
  !common /dptr95/ psiiv
  pointer psiba
  dimension psiba(:,:)
  !common /dptr95/ psiba
  pointer psisq
  dimension psisq(:,:)
  !common /dptr95/ psisq
  pointer psicu
  dimension psicu(:,:)
  !common /dptr95/ psicu
  pointer psiqu
  dimension psiqu(:,:)
  !common /dptr95/ psiqu
  pointer bavpd
  dimension bavpd(:,:)
  !common /dptr95/ bavpd
  pointer bavdn
  dimension bavdn(:,:)
  !common /dptr95/ bavdn
  pointer psiir
  dimension psiir(:,:)
  !common /dptr95/ psiir
  pointer vderb
  dimension vderb(:,:)
  !common /dptr95/ vderb
  pointer sinn
  dimension sinn(:,:)
  !common /dptr95/ sinn
  pointer tann
  dimension tann(:,:)
  !common /dptr95/ tann
  pointer ymid
  dimension ymid(:,:)
  !common /dptr95/ ymid
  pointer tau
  dimension tau(:,:)
  !common /dptr95/ tau
  pointer vptb
  dimension vptb(:,:)
  !common /dptr95/ vptb
  pointer zboun
  dimension zboun(:,:)
  !common /dptr95/ zboun
  pointer idx
  dimension idx(:,:)
  !common /dptr95/ idx
  pointer imax
  dimension imax(:,:)
  !common /dptr95/ imax
  pointer dz
  dimension dz(:,:)
  !common /dptr95/ dz
  pointer pol
  dimension pol(:,:)
  !common /dptr95/ pol
  pointer solrz
  dimension solrz(:,:)
  !common /dptr95/ solrz
  pointer solzz
  dimension solzz(:,:)
  !common /dptr95/ solzz
  pointer bpolz, btorz ! Equil.field at pol.angle grid (lza)
  dimension bpolz(:,:), btorz(:,:)  ! (lza,lrzmax)
  !common /dptr95/ bpolz, btorz

  ! YuP: [added Apr/2014] area and volume of a cell associated with each
  !                     (R,Z) point on flux surface, (R,Z)==(solrz,solzz)
  pointer ddarea, ddvol
  dimension ddarea(:,:), ddvol(:,:)  !  (lza,lrzmax)
  !common /dptr95/ ddarea, ddvol

  pointer thtab
  dimension thtab(:,:)
  !common /dptr95/ thtab
  pointer z
  dimension z(:,:)
  !common /dptr95/ z
  pointer zmid
  dimension zmid(:,:)
  !common /dptr95/ zmid
  pointer bbpsi
  dimension bbpsi(:,:)
  !common /dptr95/ bbpsi
  pointer consnp
  dimension consnp(:,:)
  !common /dptr95/ consnp
  pointer ptime
  dimension ptime(:,:)
  !common /dptr95/ ptime
  pointer sptzrp
  dimension sptzrp(:,:)
  !common /dptr95/ sptzrp
  pointer pefld
  dimension pefld(:,:)
  !common /dptr95/ pefld
  pointer rovsp
  dimension rovsp(:,:)
  !common /dptr95/ rovsp
  pointer restp
  dimension restp(:,:)
  !common /dptr95/ restp
  pointer restnp
  dimension restnp(:,:)
  !common /dptr95/ restnp
  pointer vpov
  dimension vpov(:,:)
  !common /dptr95/ vpov
  pointer es
  dimension es(:,:)
  !common /dptr95/ es
  pointer bpsi
  dimension bpsi(:,:)
  !common /dptr95/ bpsi
  pointer d2bpsi
  dimension d2bpsi(:,:)
  !common /dptr95/ d2bpsi
  pointer d2solrz
  dimension d2solrz(:,:)
  !common /dptr95/ d2solrz
  pointer d2solzz
  dimension d2solzz(:,:)
  !common /dptr95/ d2solzz
  pointer d2bpolz, d2btorz
  dimension d2bpolz(:,:), d2btorz(:,:)
  !common /dptr95/ d2bpolz, d2btorz
  pointer d2thtpol
  dimension d2thtpol(:,:)
  !common /dptr95/ d2thtpol
  pointer d2es
  dimension d2es(:,:)
  !common /dptr95/ d2es
  pointer thtpol
  dimension thtpol(:,:)
  !common /dptr95/ thtpol
  pointer esfi
  dimension esfi(:,:)
  !common /dptr95/ esfi
  pointer psiesfi
  dimension psiesfi(:,:)
  !common /dptr95/ psiesfi
  pointer psifi
  dimension psifi(:,:)
  !common /dptr95/ psifi
  pointer espsifi
  dimension espsifi(:,:)
  !common /dptr95/ espsifi
  pointer soupp
  dimension soupp(:,:)
  !common /dptr95/ soupp
  pointer waa
  dimension waa(:,:,:)
  !common /dptr95/ waa
  pointer wbb
  dimension wbb(:,:,:)
  !common /dptr95/ wbb
  pointer cosz
  dimension cosz(:,:,:)
  !common /dptr95/ cosz
  pointer dtau
  dimension dtau(:,:,:)
  !common /dptr95/ dtau
  pointer sinz
  dimension sinz(:,:,:)
  !common /dptr95/ sinz
  pointer tanz
  dimension tanz(:,:,:)
  !common /dptr95/ tanz
  pointer yz
  dimension yz(:,:,:)
  !common /dptr95/ yz
  pointer tot
  dimension tot(:,:,:)
  !common /dptr95/ tot
  pointer vflux
  dimension vflux(:,:,:)
  !common /dptr95/ vflux
  pointer f_aveth
  dimension f_aveth(:,:,:,:)
  !common /dptr95/ f_aveth
  pointer sincosba
  dimension sincosba(:,:,:)
  !common /dptr95/ sincosba
  pointer densz
  dimension densz(:,:,:,:)
  !common /dptr95/ densz
  pointer ss
  dimension ss(:,:,:,:)
  !common /dptr95/ ss
  pointer dcofleg
  dimension dcofleg(:,:,:,:)
  !common /dptr95/ dcofleg
  pointer dpcosz
  dimension dpcosz(:,:,:,:)
  !common /dptr95/ dpcosz
  pointer ssy
  dimension ssy(:,:,:,:)
  !common /dptr95/ ssy
  pointer ssyy
  dimension ssyy(:,:,:,:)
  !common /dptr95/ ssyy
  pointer ssyi
  dimension ssyi(:,:,:,:)
  !common /dptr95/ ssyi
  pointer ssyyy
  dimension ssyyy(:,:,:,:)
  !common /dptr95/ ssyyy
  pointer pcurr, pcurrm
  dimension pcurr(:,:,:), pcurrm(:,:,:)
  !common /dptr95/ pcurr, pcurrm
  pointer pdens, pdenm
  dimension pdens(:,:,:), pdenm(:,:,:)
  !common /dptr95/ pdens, pdenm
  pointer pengy, pengym
  dimension pengy(:,:,:), pengym(:,:,:)
  !common /dptr95/ pengy, pengym
  pointer pdenra
  dimension pdenra(:,:)
  !common /dptr95/ pdenra
  pointer pcurra
  dimension pcurra(:,:)
  !common /dptr95/ pcurra
  pointer pfdenra
  dimension pfdenra(:,:)
  !common /dptr95/ pfdenra
  pointer pfcurra
  dimension pfcurra(:,:)
  !common /dptr95/ pfcurra
  pointer pucrit
  dimension pucrit(:,:)
  !common /dptr95/ pucrit
  pointer peoe0
  dimension peoe0(:,:)
  !common /dptr95/ peoe0
  pointer psrc
  dimension psrc(:,:)
  !common /dptr95/ psrc
  pointer peoed
  dimension peoed(:,:)
  !common /dptr95/ peoed
  pointer cint2
  dimension cint2(:)
  !common /dptr95/ cint2
  pointer dx,dxi
  dimension dx(:),dxi(:)
  !common /dptr95/ dx,dxi
  pointer ifp
  dimension ifp(:)
  !common /dptr95/ ifp
  pointer sg
  dimension sg(:)
  !common /dptr95/ sg
  pointer sgx
  dimension sgx(:)
  !common /dptr95/ sgx
  pointer sgxx
  dimension sgxx(:)
  !common /dptr95/ sgxx
  pointer sh
  dimension sh(:)
  !common /dptr95/ sh
  pointer shx
  dimension shx(:)
  !common /dptr95/ shx
  pointer shxx
  dimension shxx(:)
  !common /dptr95/ shxx
  pointer shxxx
  dimension shxxx(:)
  !common /dptr95/ shxxx
  pointer tam1
  dimension tam1(:)
  !common /dptr95/ tam1
  pointer tam2
  dimension tam2(:)
  !common /dptr95/ tam2
  pointer tam3
  dimension tam3(:)
  !common /dptr95/ tam3
  pointer tam4
  dimension tam4(:)
  !common /dptr95/ tam4
  pointer tam5
  dimension tam5(:)
  !common /dptr95/ tam5
  pointer tam6
  dimension tam6(:)
  !common /dptr95/ tam6
  pointer tam7
  dimension tam7(:)
  !common /dptr95/ tam7
  pointer tam8
  dimension tam8(:)
  !common /dptr95/ tam8
  pointer tam9
  dimension tam9(:)
  !common /dptr95/ tam9
  pointer tam10
  dimension tam10(:)
  !common /dptr95/ tam10
  pointer tam11
  dimension tam11(:)
  !common /dptr95/ tam11
  pointer tam12
  dimension tam12(:)
  !common /dptr95/ tam12
  pointer tam13
  dimension tam13(:)
  !common /dptr95/ tam13
  pointer tam14
  dimension tam14(:)
  !common /dptr95/ tam14
  pointer tam15
  dimension tam15(:)
  !common /dptr95/ tam15
  pointer tam16
  dimension tam16(:)
  !common /dptr95/ tam16
  pointer tam17
  dimension tam17(:)
  !common /dptr95/ tam17
  pointer tam18
  dimension tam18(:)
  !common /dptr95/ tam18
  pointer tam19
  dimension tam19(:)
  !common /dptr95/ tam19
  pointer tam20
  dimension tam20(:)
  !common /dptr95/ tam20
  pointer tam21
  dimension tam21(:)
  !common /dptr95/ tam21
  pointer tam22
  dimension tam22(:)
  !common /dptr95/ tam22
  pointer tam23
  dimension tam23(:)
  !common /dptr95/ tam23
  pointer tam24
  dimension tam24(:)
  !common /dptr95/ tam24
  pointer tam25
  dimension tam25(:)
  !common /dptr95/ tam25
  pointer tam26
  dimension tam26(:)
  !common /dptr95/ tam26
  pointer tam27
  dimension tam27(:)
  !common /dptr95/ tam27
  pointer tam28
  dimension tam28(:)
  !common /dptr95/ tam28
  pointer tam29
  dimension tam29(:)
  !common /dptr95/ tam29
  pointer tam30
  dimension tam30(:)
  !common /dptr95/ tam30
  pointer x
  dimension x(:)
  !common /dptr95/ x
  pointer xmidpt
  dimension xmidpt(:)
  !common /dptr95/ xmidpt
  pointer xi
  dimension xi(:)
  !common /dptr95/ xi
  pointer xsq
  dimension xsq(:)
  !common /dptr95/ xsq
  pointer x3i
  dimension x3i(:)
  !common /dptr95/ x3i
  pointer x2i
  dimension x2i(:)
  !common /dptr95/ x2i
  pointer xcu
  dimension xcu(:)
  !common /dptr95/ xcu
  pointer xcenter
  dimension xcenter(:)
  !common /dptr95/ xcenter
  pointer xcensq, xcent3
  dimension xcensq(:), xcent3(:)
  !common /dptr95/ xcensq, xcent3
  pointer uoc
  dimension uoc(:)
  !common /dptr95/ uoc
  pointer enerkev
  dimension enerkev(:,:) !YuP[2018-01-08] added 2nd index (k)
  !common /dptr95/ enerkev
  pointer gamma
  dimension gamma(:)
  !common /dptr95/ gamma
  pointer gamsqr
  dimension gamsqr(:)
  !common /dptr95/ gamsqr
  pointer gamcub
  dimension gamcub(:)
  !common /dptr95/ gamcub
  pointer gammi
  dimension gammi(:)
  !common /dptr95/ gammi
  pointer gamm2i
  dimension gamm2i(:)
  !common /dptr95/ gamm2i
  pointer gamm1
  dimension gamm1(:)
  !common /dptr95/ gamm1
  pointer tcsgm1
  dimension tcsgm1(:)
  !common /dptr95/ tcsgm1
  pointer gamefac
  dimension gamefac(:)
  !common /dptr95/ gamefac
  pointer ident
  dimension ident(:)
  !common /dptr95/ ident
  pointer temc1
  dimension temc1(:)
  !common /dptr95/ temc1
  pointer temc2
  dimension temc2(:)
  !common /dptr95/ temc2
  pointer temc3
  dimension temc3(:)
  !common /dptr95/ temc3
  pointer temc4
  dimension temc4(:)
  !common /dptr95/ temc4
  pointer itemc1
  dimension itemc1(:)
  !common /dptr95/ itemc1
  pointer itemc2
  dimension itemc2(:)
  !common /dptr95/ itemc2
  pointer l_lower
  dimension l_lower(:)
  !common /dptr95/ l_lower
  pointer lpt
  dimension lpt(:)
  !common /dptr95/ lpt
  pointer mun
  dimension mun(:)
  !common /dptr95/ mun
  pointer fll
  dimension fll(:)
  !common /dptr95/ fll
  pointer xpar
  dimension xpar(:)
  !common /dptr95/ xpar
  pointer rheads
  dimension rheads(:)
  !common /dptr95/ rheads
  pointer dfvlle
  dimension dfvlle(:)
  !common /dptr95/ dfvlle
  pointer dfvlli
  dimension dfvlli(:)
  !common /dptr95/ dfvlli
  pointer xperp
  dimension xperp(:)
  !common /dptr95/ xperp
  pointer xl
  dimension xl(:)
  !common /dptr95/ xl
  pointer jmaxxl
  dimension jmaxxl(:)
  !common /dptr95/ jmaxxl
  pointer xlm
  dimension xlm(:)
  !common /dptr95/ xlm
  pointer dxl
  dimension dxl(:)
  !common /dptr95/ dxl
  pointer fl
  dimension fl(:)
  !common /dptr95/ fl
  pointer fl1
  dimension fl1(:)
  !common /dptr95/ fl1
  pointer fl2
  dimension fl2(:)
  !common /dptr95/ fl2
  pointer ppars
  dimension ppars(:,:)
  !common /dptr95/ ppars
  pointer pprps
  dimension pprps(:,:)
  !common /dptr95/ pprps
  pointer faci
  dimension faci(:,:)
  !common /dptr95/ faci
  pointer pparea
  dimension pparea(:,:)
  !common /dptr95/ pparea
  pointer wtfl0
  dimension wtfl0(:,:,:)
  !common /dptr95/ wtfl0
  pointer wtflm
  dimension wtflm(:,:,:)
  !common /dptr95/ wtflm
  pointer jflbin
  dimension jflbin(:,:,:)
  !common /dptr95/ jflbin
  pointer xm
  dimension xm(:,:)
  !common /dptr95/ xm
  pointer dbb
  dimension dbb(:,:)
  !common /dptr95/ dbb
  pointer dd
  dimension dd(:,:)
  !common /dptr95/ dd
  pointer de
  dimension de(:,:)
  !common /dptr95/ de
  pointer df
  dimension df(:,:)
  !common /dptr95/ df
  pointer dff
  dimension dff(:,:)
  !common /dptr95/ dff
  pointer cah
  dimension cah(:,:)
  !common /dptr95/ cah
  pointer cthta
  dimension cthta(:,:)
  !common /dptr95/ cthta
  pointer gon
  dimension gon(:,:)
  !common /dptr95/ gon
  pointer so
  dimension so(:,:)
  !common /dptr95/ so
  pointer currv
  dimension currv(:,:,:)
  !common /dptr95/ currv
  pointer currvs
  dimension currvs(:,:)
  !common /dptr95/ currvs
  pointer pwrrf
  dimension pwrrf(:,:,:)
  !common /dptr95/ pwrrf
  pointer tal
  dimension tal(:,:)
  !common /dptr95/ tal
  pointer tbl
  dimension tbl(:,:)
  !common /dptr95/ tbl
  pointer tfl
  dimension tfl(:,:)
  !common /dptr95/ tfl
  pointer pwrrfs
  dimension pwrrfs(:,:,:)
  !common /dptr95/ pwrrfs
  pointer pleg
  dimension pleg(:,:)
  !common /dptr95/ pleg
  pointer feta
  dimension feta(:,:)
  !common /dptr95/ feta
  pointer fetb
  dimension fetb(:,:)
  !common /dptr95/ fetb
  pointer wflux
  dimension wflux(:,:,:)
  !common /dptr95/ wflux
  !     NB:  rhs set up here for full 3d set of eqns (BH070525)
  pointer rhs
  dimension rhs(:)
  !common /dptr95/ rhs
  pointer sovt
  dimension sovt(:,:,:,:)
  !common /dptr95/ sovt
  pointer sigsxr
  dimension sigsxr(:,:,:,:)
  !common /dptr95/ sigsxr

  pointer pentr
  dimension pentr(:,:,:,:)  !!!(nonch,ngen,-1:15,lrors)
  !common /dptr95/ pentr

  pointer constp
  dimension constp(:,:)  !!!(nonch,lrors)
  !common /dptr95/ constp

  pointer :: sigmtt(:,:),sigftt(:,:)  !!!(nonch,4)
  !common /dptr95/ sigmtt,sigftt

  pointer sgaint
  dimension sgaint(:,:,:)  !!!(8,ngen,lrors)
  pointer entr
  dimension entr(:,:,:)    !!!(ngen,-1:15,lrors)
  pointer xlndnz
  dimension xlndnz(:,:)    !!!(ngen+1,negyrga)
  pointer sounor
  dimension sounor(:,:,:,:)   !!!(ngen,nsoa,lz,lrz)
  !common /dptr95/ sgaint,entr,xlndnz,sounor


  !.......................................................................
  !*****arrays related to relativ=fully option
  !.......................................................................
  pointer gamman
  dimension gamman(:,:)
  !common /dptr95/ gamman
  pointer alphan
  dimension alphan(:,:)
  !common /dptr95/ alphan

  pointer asnha
  dimension asnha(:)
  !common /dptr95/ asnha
  pointer item1
  dimension item1(:)
  !common /dptr95/ item1
  pointer item2
  dimension item2(:)
  !common /dptr95/ item2
  pointer item3
  dimension item3(:)
  !common /dptr95/ item3
  pointer item4
  dimension item4(:)
  !common /dptr95/ item4
  pointer item5
  dimension item5(:)
  !common /dptr95/ item5
  pointer item6
  dimension item6(:)
  !common /dptr95/ item6
  pointer dxm5
  dimension dxm5(:)
  !common /dptr95/ dxm5
  pointer exm5
  dimension exm5(:)
  !common /dptr95/ exm5
  pointer dxp5
  dimension dxp5(:)
  !common /dptr95/ dxp5
  pointer exp5
  dimension exp5(:)
  !common /dptr95/ exp5
  pointer tamt1
  dimension tamt1(:,:,:,:)
  !common /dptr95/ tamt1
  pointer tamt2
  dimension tamt2(:,:,:,:)
  !common /dptr95/ tamt2
  pointer da
  dimension da(:,:)
  !common /dptr95/ da
  pointer db
  dimension db(:,:)
  !common /dptr95/ db
  pointer dc
  dimension dc(:,:)
  !common /dptr95/ dc
  pointer ca
  dimension ca(:,:)
  !common /dptr95/ ca
  pointer cb
  dimension cb(:,:)
  !common /dptr95/ cb
  pointer cc
  dimension cc(:,:)
  !common /dptr95/ cc
  pointer cd
  dimension cd(:,:)
  !common /dptr95/ cd
  pointer ce
  dimension ce(:,:)
  !common /dptr95/ ce
  pointer cf
  dimension cf(:,:)
  !common /dptr95/ cf

  pointer tem1
  dimension tem1(:)
  !common /dptr95/ tem1
  pointer tem2
  dimension tem2(:)
  !common /dptr95/ tem2
  pointer tem3
  dimension tem3(:)
  !common /dptr95/ tem3
  pointer tem4
  dimension tem4(:)
  !common /dptr95/ tem4
  pointer tem5
  dimension tem5(:)
  !common /dptr95/ tem5
  pointer tem6
  dimension tem6(:)
  !common /dptr95/ tem6

  pointer egg
  dimension egg(:,:)
  !common /dptr95/ egg
  pointer fgg
  dimension fgg(:,:)
  !common /dptr95/ fgg

  pointer xhead
  dimension xhead(:,:)
  !common /dptr95/ xhead
  pointer xtail
  dimension xtail(:,:)
  !common /dptr95/ xtail
  pointer ytail
  dimension ytail(:,:)
  !common /dptr95/ ytail
  pointer yhead
  dimension yhead(:,:)
  !common /dptr95/ yhead

  pointer fpn
  dimension fpn(:,:)
  !common /dptr95/ fpn

  pointer temp1
  dimension temp1(:,:)
  !common /dptr95/ temp1
  pointer temp2
  dimension temp2(:,:)
  !common /dptr95/ temp2
  pointer temp3
  dimension temp3(:,:)
  !common /dptr95/ temp3
  pointer temp4
  dimension temp4(:,:)
  !common /dptr95/ temp4
  pointer temp5
  dimension temp5(:,:)
  !common /dptr95/ temp5
  pointer temp6
  dimension temp6(:,:)
  !common /dptr95/ temp6

  pointer xllji
  dimension xllji(:,:)
  !common /dptr95/ xllji
  pointer xppji
  dimension xppji(:,:)
  !common /dptr95/ xppji

  !     Arrays used for first order orbit width calculations:
  pointer deltarho
  dimension deltarho(:,:,:)
  !common /dptr95/ deltarho
  pointer deltarhop
  dimension deltarhop(:,:,:)
  !common /dptr95/ deltarhop
  pointer deltarz
  dimension deltarz(:,:,:)
  !common /dptr95/ deltarz
  pointer r_delta
  dimension r_delta(:)
  !common /dptr95/ r_delta
  pointer z_delta
  dimension z_delta(:)
  !common /dptr95/ z_delta
  pointer t_delta
  dimension t_delta(:)
  !common /dptr95/ t_delta
  pointer delta_bdb0
  dimension delta_bdb0(:,:)
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

  pointer cqlb,cqlc,cqle,cqlf     ! (iy,jx,lrz,mrfn)
  dimension cqlb(:,:,:,:),cqlc(:,:,:,:),cqle(:,:,:,:),cqlf(:,:,:,:)
  !common/qlcoef/cqlb,cqlc,cqle,cqlf 

  pointer bqlm
  dimension bqlm(:,:)  ! (iy,jx)
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
  real(c_double) :: mun
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
  real(c_double) :: sxry(lrza,4),sang(lrza,4),spol(lrza,4), &
       ibin(lrza,4),eflux(nena,nva),efluxt(nva),alphad(3),xs_(3), &
       enk(nena),en_(nena)
  integer :: jval_(nena),inegsxr(nva),lensxr(nva)

  !common /csigma/
  integer :: mtab,msig,jxis
  real(c_double) :: elmin,delegy
  integer :: imaxwln(2,4),igenrl(2,4)
  real(c_double) :: sigm(4,lrorsa),sigf(4,lrorsa),sigmt(4),sigft(4), &
       fuspwrv(4,lrorsa),fuspwrvt(4),fuspwrm(4,lrorsa),fuspwrmt(4)

  pointer tamm1
  dimension tamm1(:)
  !common /csigma/ tamm1  !(0:mmsv)

  pointer iind
  dimension iind(:)
  !common /csigma/ iind  !(1:jx)

  !..............................................................
  !     Set up pointers for sigma-v
  !..............................................................

  pointer csv
  dimension csv(:,:,:)
  !common /dptr95/ csv
  pointer svtab
  dimension svtab(:)
  !common /dptr95/ svtab


  !..............................................................
  !     Set up pointers to allocatable arrays for transport model.
  !     Space allocated in subroutine tdtraloc
  !..............................................................


  pointer frn_2
  dimension frn_2(:,:,:,:)
  !common /dptr95/ frn_2
  pointer frn_1
  dimension frn_1(:,:,:,:)
  !common /dptr95/ frn_1
  pointer frn
  dimension frn(:,:,:,:)
  !common /dptr95/ frn
  pointer fvn_1
  dimension fvn_1(:,:,:,:)
  !common /dptr95/ fvn_1
  pointer fvn
  dimension fvn(:,:,:,:)
  !common /dptr95/ fvn
  pointer dl
  dimension dl(:,:,:,:)
  !common /dptr95/ dl
  pointer d_rr
  dimension d_rr(:,:,:,:)
  !common /dptr95/ d_rr
  pointer d_r
  dimension d_r(:,:,:,:)
  !common /dptr95/ d_r
  pointer f_lm
  dimension f_lm(:,:,:)
  !common /dptr95/ f_lm
  pointer f_lp
  dimension f_lp(:,:,:)
  !common /dptr95/ f_lp
  pointer f_up
  dimension f_up(:,:,:)
  !common /dptr95/ f_up
  pointer f_vtor
  dimension f_vtor(:,:,:,:)
  !common /dptr95/ f_vtor
  pointer cynt2_
  dimension cynt2_(:,:)
  !common /dptr95/ cynt2_
  pointer vpint_
  dimension vpint_(:,:)
  !common /dptr95/ vpint_
  pointer vptb_
  dimension vptb_(:,:)
  !common /dptr95/ vptb_
  pointer cosovb
  dimension cosovb(:,:)
  !common /dptr95/ cosovb
  pointer bovcos
  dimension bovcos(:,:)
  !common /dptr95/ bovcos
  pointer adv
  dimension adv(:,:)
  !common /dptr95/ adv
  pointer dentarget
  dimension dentarget(:)
  !common /dptr95/ dentarget
  pointer eg_
  dimension eg_(:,:,:)
  !common /dptr95/ eg_
  pointer fg_
  dimension fg_(:,:,:)
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

  pointer drpmconz
  dimension drpmconz(:)
  !common /dptr95/ drpmconz
  pointer eqdell
  dimension eqdell(:,:)
  !common /dptr95/ eqdell
  pointer eqbpol
  dimension eqbpol(:,:)
  !common /dptr95/ eqbpol
  pointer solr
  dimension solr(:,:)
  !common /dptr95/ solr
  pointer solz
  dimension solz(:,:)
  !common /dptr95/ solz



  !*********************************************************************
  !     BEGIN arrays for LOWER HYBRID FAST WAVE and ECH Module.
  !*********************************************************************


  !common/params/ &
  integer :: jjx, jbm1,jb0,jbp1

  pointer jbm1
  dimension jbm1(:,:)
  !common jbm1
  pointer jb0
  dimension jb0(:,:)
  !common jb0
  pointer jbp1
  dimension jbp1(:,:)
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

  complex(c_double_complex) :: cosz1,sinz1,sinz2
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


  complex(c_double_complex) :: cwexde,cweyde,cwezde

  pointer urfb
  dimension urfb(:,:,:,:)
  !common /dptr95/ urfb
  pointer urfc
  dimension urfc(:,:,:,:)
  !common /dptr95/ urfc
  pointer cosmz
  dimension cosmz(:,:,:)
  !common /dptr95/ cosmz
  pointer g_
  dimension g_(:,:,:,:)
  !common /dptr95/ g_
  pointer alfag
  dimension alfag(:)
  !common /dptr95/ alfag
  pointer argmnt
  dimension argmnt(:)
  !common /dptr95/ argmnt
  pointer ilim1d
  dimension ilim1d(:)
  !common /dptr95/ ilim1d
  pointer ilim2d
  dimension ilim2d(:)
  !common /dptr95/ ilim2d
  pointer ilim1dd
  dimension ilim1dd(:)
  !common /dptr95/ ilim1dd
  pointer ilim2dd
  dimension ilim2dd(:)
  !common /dptr95/ ilim2dd
  pointer sx
  dimension sx(:)
  !common /dptr95/ sx
  pointer xmdx
  dimension xmdx(:)
  !common /dptr95/ xmdx
  pointer cosz1
  dimension cosz1(:)
  !common /dptr95/ cosz1
  pointer sinz1
  dimension sinz1(:)
  !common /dptr95/ sinz1
  pointer sinz2
  dimension sinz2(:)
  !common /dptr95/ sinz2
  pointer thtf1
  dimension thtf1(:)
  !common /dptr95/ thtf1
  pointer thtf2
  dimension thtf2(:)
  !common /dptr95/ thtf2
  pointer alfi
  dimension alfi(:)
  !common /dptr95/ alfi
  pointer alfa
  dimension alfa(:)
  !common /dptr95/ alfa
  pointer ilim1
  dimension ilim1(:)
  !common /dptr95/ ilim1
  pointer ilim2
  dimension ilim2(:)
  !common /dptr95/ ilim2
  pointer ifct1
  dimension ifct1(:)
  !common /dptr95/ ifct1
  pointer ifct2
  dimension ifct2(:)
  !common /dptr95/ ifct2
  pointer urftmp
  dimension urftmp(:)
  !common /dptr95/ urftmp
  pointer urfpwr
  dimension urfpwr(:,:,:)
  !common /dptr95/ urfpwr
  pointer urfpwrc
  dimension urfpwrc(:,:,:)
  !common /dptr95/ urfpwrc
  pointer urfpwrl
  dimension urfpwrl(:,:,:)
  !common /dptr95/ urfpwrl
  pointer jminray
  dimension jminray(:,:,:)
  !common /dptr95/ jminray
  pointer jmaxray
  dimension jmaxray(:,:,:)
  !common /dptr95/ jmaxray
  pointer lloc
  dimension lloc(:,:,:)
  !common /dptr95/ lloc
  pointer llray
  dimension llray(:,:,:)
  !common /dptr95/ llray
  pointer psiloc
  dimension psiloc(:,:,:)
  !common /dptr95/ psiloc
  pointer scalurf
  dimension scalurf(:,:,:)
  !common /dptr95/ scalurf
  pointer cwexde
  dimension cwexde(:,:,:)
  !common /dptr95/ cwexde
  pointer cweyde
  dimension cweyde(:,:,:)
  !common /dptr95/ cweyde
  pointer cwezde
  dimension cwezde(:,:,:)
  !common /dptr95/ cwezde
  pointer delpwr
  dimension delpwr(:,:,:)
  !common /dptr95/ delpwr
  pointer fluxn
  dimension fluxn(:,:,:)
  !common /dptr95/ fluxn
  pointer seikon
  dimension seikon(:,:,:)
  !common /dptr95/ seikon
  pointer spsi
  dimension spsi(:,:,:)
  !common /dptr95/ spsi
  pointer sdpwr
  dimension sdpwr(:,:,:)
  !common /dptr95/ sdpwr
  pointer sbtot
  dimension sbtot(:,:,:)
  !common /dptr95/ sbtot
  pointer sene
  dimension sene(:,:,:)
  !common /dptr95/ sene
  pointer salphac
  dimension salphac(:,:,:)
  !common /dptr95/ salphac
  pointer salphal
  dimension salphal(:,:,:)
  !common /dptr95/ salphal
  pointer ws
  dimension ws(:,:,:)
  !common /dptr95/ ws
  pointer wr
  dimension wr(:,:,:)
  !common /dptr95/ wr
  pointer wz
  dimension wz(:,:,:)
  !common /dptr95/ wz
  pointer wnpar
  dimension wnpar(:,:,:)
  !common /dptr95/ wnpar
  pointer wdnpar
  dimension wdnpar(:,:,:)
  !common /dptr95/ wdnpar
  pointer wnper
  dimension wnper(:,:,:)
  !common /dptr95/ wnper
  pointer wphi
  dimension wphi(:,:,:)
  !common /dptr95/ wphi
  pointer ilowp
  dimension ilowp(:,:)
  !common /dptr95/ ilowp
  pointer iupp
  dimension iupp(:,:)
  !common /dptr95/ iupp
  pointer ifct1_
  dimension ifct1_(:,:)
  !common /dptr95/ ifct1_
  pointer ifct2_
  dimension ifct2_(:,:)
  !common /dptr95/ ifct2_
  pointer nrayelt
  dimension nrayelt(:,:)
  !common /dptr95/ nrayelt
  pointer jslofas
  dimension jslofas(:,:)
  !common /dptr95/ jslofas
  pointer nurefls
  dimension nurefls(:,:)
  !common /dptr95/ nurefls
  pointer keiks
  dimension keiks(:,:)
  !common /dptr95/ keiks
  pointer jpes
  dimension jpes(:,:)
  !common /dptr95/ jpes
  pointer jpis
  dimension jpis(:,:)
  !common /dptr95/ jpis
  pointer istarts
  dimension istarts(:,:)
  !common /dptr95/ istarts
  pointer iprmt5
  dimension iprmt5(:,:)
  !common /dptr95/ iprmt5
  pointer jhlfs
  dimension jhlfs(:,:)
  !common /dptr95/ jhlfs
  pointer sxxrt
  dimension sxxrt(:,:)
  !common /dptr95/ sxxrt
  pointer skpsi
  dimension skpsi(:,:)
  !common /dptr95/ skpsi
  pointer skth
  dimension skth(:,:)
  !common /dptr95/ skth
  pointer skphi
  dimension skphi(:,:)
  !common /dptr95/ skphi
  pointer lrayelt
  dimension lrayelt(:,:)
  !common /dptr95/ lrayelt
  pointer delpwr0
  dimension delpwr0(:,:)
  !common /dptr95/ delpwr0
  pointer nrayelt0
  dimension nrayelt0(:,:)
  !common /dptr95/ nrayelt0
  pointer truncd
  dimension truncd(:) ! 1:jx
  !common /dptr95/ truncd


  !..................................................................
  !     Allocatable arrays allocated in subroutine rdc_multi,
  !     used after subroutine execution.
  !     Here, we introduce f90 pointers, as they are easier
  !     to allocate.

  pointer rdcb
  dimension rdcb(:,:,:,:)
  !common /dptr95/ rdcb
  pointer rdcc
  dimension rdcc(:,:,:,:)
  !common /dptr95/ rdcc
  pointer rdce
  dimension rdce(:,:,:,:)
  !common /dptr95/ rdce
  pointer rdcf
  dimension rdcf(:,:,:,:)
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

  pointer abd_lapack,a_csr,alu,w_ilu,rhs0,sol,vv
  pointer ja_csr,ia_csr,jlu,ju,jw_ilu
  pointer ar_csr,ac_csr
  pointer jar_csr,iar_csr,ipofi,jac_csr,iac_csr
  dimension abd_lapack(:,:),a_csr(:), &
       alu(:),w_ilu(:),rhs0(:),sol(:),vv(:)
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

  pointer l_upper
  dimension l_upper(:)  !!! (1:iy)
  pointer ilpm1ef
  dimension ilpm1ef(:,:,:)  !!! (0:iy+1,0:lsa1,-1:+1)

  pointer fnhalf
  dimension fnhalf(:,:,:,:)
  !common /dptr95/ fnhalf
  pointer fnp0
  dimension fnp0(:,:,:,:)
  !common /dptr95/ fnp0
  pointer fnp1
  dimension fnp1(:,:,:,:)
  !common /dptr95/ fnp1
  pointer dls
  dimension dls(:,:,:,:)
  !common /dptr95/ dls
  pointer fh
  dimension fh(:,:,:,:)
  !common /dptr95/ fh
  pointer fg
  dimension fg(:,:,:,:)
  !common /dptr95/ fg
  pointer fedge
  dimension fedge(:,:,:,:)
  !common /dptr95/ fedge
  pointer rhspar
  dimension rhspar(:,:,:)
  !common /dptr95/ rhspar
  pointer bndmats
  dimension bndmats(:,:,:,:)
  !common /dptr95/ bndmats
  pointer wcqlb
  dimension wcqlb(:,:,:,:)
  !common /dptr95/ wcqlb
  pointer wcqlc
  dimension wcqlc(:,:,:,:)
  !common /dptr95/ wcqlc
  pointer wcqle
  dimension wcqle(:,:,:,:)
  !common /dptr95/ wcqle
  pointer wcqlf
  dimension wcqlf(:,:,:,:)
  !common /dptr95/ wcqlf



  !.......................................................................
  !     Arrays allocated in ampfalloc
  !.......................................................................

  pointer ampfln
  dimension ampfln(:)
  !common /dptr95/ ampfln
  pointer ampflh
  dimension ampflh(:)
  !common /dptr95/ ampflh
  pointer ampflg
  dimension ampflg(:)
  !common /dptr95/ ampflg
  pointer ampfa
  dimension ampfa(:,:)
  !common /dptr95/ ampfa
  pointer ampfb
  dimension ampfb(:,:)
  !common /dptr95/ ampfb
  pointer ampfaa
  dimension ampfaa(:,:)
  !common /dptr95/ ampfaa
  pointer ampfc
  dimension ampfc(:)
  !common /dptr95/ ampfc
  pointer ampf2ebar
  dimension ampf2ebar(:)
  !common /dptr95/ ampf2ebar

  !.......................................................................
  !     Arrays for finite orbit width (FOW) calculations
  !.......................................................................

  !common/psiaxis/ psi_lim,psi_mag,R_axis,Z_axis ![cgs]

  pointer rcontr,zcontr,rlimiter,zlimiter
  dimension rcontr(:),zcontr(:),rlimiter(:),zlimiter(:)
  !common/limiter/
  integer :: ncontr
  integer :: nlimiter
  real(c_double) :: rcontr, zcontr,  &  !-> Last closed flux surface &
       rlimiter, zlimiter !-> Limiter surface [cm]
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

  save
contains
  ! lolz
end module cqcomm_mod
