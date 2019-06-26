!     cqlcomm
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
!BH180527: Dimension of tem[1-6] modified to iyjx2l=max(iy+2,setup0%lrz)*(jx+2),
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
! # YuP[2019-05-31] name_decl.h is removed - not used anymore
!   [all declarations are here, in cqlcomm]
!.......................................................................

module cqlcomm_mod

  !---BEGIN USE

  use cqlconf_mod, only : setup0, eqsetup, rfsetup, trsetup, sousetup
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double_complex
  use param_mod

  !---END USE

  implicit none
  ! since we didn't like using derived types in the code,
  ! we can use pointers, but we can't associate them
  !    (well you can but is silently wrong)
  ! during decleration in f90/95.  I think f08 or something; too bad  I guess.
  logical, private :: initialized_cqlcomm = .FALSE.
  logical, private :: initialized_eq_pointers = .FALSE.
  logical, private :: initialized_rf_pointers = .FALSE.
  logical, private :: initialized_tr_pointers = .FALSE.
  logical, private :: initialized_sou_pointers = .FALSE.
  logical, private :: initialized_setup_pointers = .FALSE.
  private :: initialize_eq_pointers
  private :: initialize_rf_pointers
  private :: initialize_tr_pointers
  private :: initialize_sou_pointers
  private :: initialize_setup_pointers

  public

  !.......................................................................
  !     nml variables that take on the values assigned to parameters.
  !.......................................................................
  integer :: iy
  integer :: jx
  integer :: lz
  integer :: mx
  integer :: nbctime
  integer :: negyrg
  integer :: ngen
  integer :: nmax

  !.......................................................................
  !     SCALAR INPUT FOR NAMELIST SETUP...
  !.......................................................................

  character(len=8) :: chang
  character(len=8) :: eleccomp
  character(len=8) :: f4d_out
  character(len=8) :: tavg
  character(len=8) :: iactst
  character(len=8) :: ineg
  character(len=8) :: idrop
  character(len=8) :: idskf
  character(len=8) :: idskrf
  character(len=8) :: ichkpnt
  character(len=8) :: implct
  character(len=8) :: lbdry0
  character(len=8) :: locquas
  character(len=8) :: taunew
  character(len=8) :: machine
  character(len=8) :: meshy
  character(len=8) :: manymat
  character(len=8) :: netcdfvecal
  character(len=8) :: netcdfvecc
  character(len=8) :: netcdfvece
  character(len=8) :: netcdfvecrf
  character(len=8) :: netcdfvecs
  character(len=8) :: psimodel
  character(len=8) :: pltpowe
  character(len=8) :: pltend
  character(len=8) :: pltinput
  character(len=8) :: pltlim
  character(len=8) :: pltrdc
  character(len=8) :: pltrst
  character(len=8) :: plturfb
  character(len=8) :: pltvflu
  character(len=8) :: pltra
  character(len=8) :: pltfvs
  character(len=8) :: pltd
  character(len=8) :: pltprpp
  character(len=8) :: pltfofv
  character(len=8) :: pltlos
  character(len=8) :: pltdn
  character(len=8) :: pltvecal
  character(len=8) :: pltvecc
  character(len=8) :: pltvecrf
  character(len=8) :: pltvece
  character(len=8) :: pltstrm
  character(len=8) :: pltflux
  character(len=8) :: pltsig
  character(len=8) :: pltdnpos
  character(len=8) :: profpsi
  character(len=8) :: qsineut
  character(len=8) :: trapmod
  character(len=8) :: scatmod
  character(len=8) :: relativ
  character(len=8) :: sigmamod
  character(len=8) :: soln_method
  character(len=8) :: symtrap
  character(len=8) :: syncrad
  character(len=8) :: bremsrad
  character(len=8) :: tandem
  character(len=8) :: gamafac
  character(len=8) :: yreset

  character(len=256) :: netcdfnm
  character(len=8) :: izeff
  character(len=8) :: netcdfshort

  integer :: colmodl


  !common /readscal/ &
  real(c_double) :: btor
  real(c_double) :: bth
  real(c_double) :: contrmin
  real(c_double) :: constr
  real(c_double) :: deltabdb
  real(c_double) :: droptol
  real(c_double) :: dtr
  real(c_double) :: dtr0
  real(c_double) :: xsink
  real(c_double) :: esink
  real(c_double) :: ephicc
  real(c_double) :: esfac
  real(c_double) :: eoved
  real(c_double) :: enorm
  real(c_double) :: enorme
  real(c_double) :: enormi
  real(c_double) :: gsla
  real(c_double) :: gslb
  real(c_double) :: gamaset
  integer :: isigtst
  integer :: isigsgv1
  integer :: isigsgv2
  integer :: kenorm
  integer :: lfil
  integer :: nnspec

  complex(c_double_complex), pointer :: vlfemin(:) => null()
  complex(c_double_complex), pointer :: vlfeplus(:) => null()

  !common /readscal/ &
  integer :: nstop
  integer :: ncoef
  integer :: nchec
  integer :: ncont
  integer :: nrstrt
  integer :: nstps
  integer :: nfpld
  integer :: noncntrl
  integer :: nonel
  integer :: noffel
  integer :: nonvphi
  integer :: noffvphi
  integer :: nonloss
  integer :: noffloss
  integer :: numby
  real(c_double) :: pltmag
  real(c_double) :: xprpmax
  real(c_double) :: trapredc
  real(c_double) :: scatfrac
  real(c_double) :: radmaj
  real(c_double) :: radmin
  real(c_double) :: rmirror
  real(c_double) :: sigvi
  real(c_double) :: sigvcx
  real(c_double) :: brfac
  real(c_double) :: brfac1
  real(c_double) :: brfacgm3
  real(c_double) :: tfac
  real(c_double) :: tfacz
  real(c_double) :: temp_den
  real(c_double) :: veclnth
  real(c_double) :: vnorm

  real(c_double) :: npwr(0:ntotala)
  real(c_double) :: negy(ngena)
  real(c_double) :: ntorloss(ngena)
  real(c_double) :: mpwr(0:ntotala)
  real(c_double) :: megy(ngena)
  real(c_double) :: mtorloss(ngena)
  real(c_double) :: mpwrzeff
  real(c_double) :: npwrzeff
  real(c_double) :: mpwrvphi
  real(c_double) :: npwrvphi
  real(c_double) :: mpwrelec
  real(c_double) :: npwrelec
  real(c_double) :: mpwrxj
  real(c_double) :: npwrxj

  !common /readscal/ &
  real(c_double) :: xfac
  real(c_double) :: xpctlwr
  real(c_double) :: xpctmdl
  real(c_double) :: xlwr
  real(c_double) :: xmdl
  real(c_double) :: xsinkm
  real(c_double) :: ylower
  real(c_double) :: yupper
  real(c_double) :: elecscal
  real(c_double) :: enescal
  real(c_double) :: zeffscal
  real(c_double) :: vphiscal
  real(c_double) :: tescal
  real(c_double) :: tiscal
  real(c_double), pointer :: bctimescal => null()

  !..................................................................
  !     VECTOR DIMENSIONED NAMELIST SETUP COMMON BLOCK.....
  !..................................................................

  character(len=8), :: kpress(ntotala)
  character(len=8) :: kfield(ntotala)
  character(len=8) :: lbdry(ngena)
  character(len=8) :: lossmode(ngena)
  character(len=8) :: regy(ngena)
  character(len=8) :: torloss(ngena)
  character(len=8), pointer :: difus_type(:) => null()
  character(len=8), pointer :: difus_io(:) => null()

  character(len=256) :: lossfile(ngena)

  !common /readvec/ &
  integer :: isigmas(6)
  integer :: nplot(nplota)
  integer :: nsave(nsavea)
  real(c_double) :: bnumb(ntotala)
  real(c_double) :: fmass(ntotala)
  real(c_double) :: pltflux1(7)
  real(c_double) :: enloss(ngena)
  real(c_double) :: gamegy(ngena)
  real(c_double) :: paregy(ngena)
  real(c_double) :: peregy(ngena)
  real(c_double) :: pegy(ngena)
  real(c_double) :: zeffin(0:njenea)
  real(c_double) :: vphiplin(0:njenea)
  real(c_double) :: ennl(npaproca)
  real(c_double) :: ennb(npaproca)
  real(c_double) :: ennscal(npaproca)
  real(c_double) :: bctime(nbctimea)
  real(c_double) :: dtr1(ndtr1a)
  real(c_double) :: zeffc(nbctimea)
  real(c_double) :: zeffb(nbctimea)
  real(c_double) :: elecc(nbctimea)
  real(c_double) :: elecb(nbctimea)
  real(c_double) :: vphic(nbctimea)
  real(c_double) :: vphib(nbctimea)
  real(c_double) :: xjc(nbctimea)
  real(c_double) :: xjb(nbctimea)
  real(c_double) :: totcrt(nbctimea)
  real(c_double) :: tavg1(ntavga)
  real(c_double) :: tavg2(ntavga)
  integer :: nondtr1(ndtr1a)

  !..................................................................
  !     TWO DIMENSIONAL NAMELIST SETUP COMMON BLOCK.....
  !..................................................................

  !common /readarr/ &
  real(c_double) :: tauloss(3,ngena)
  character(len=8) :: kspeci(2,ntotala)
  real(c_double) :: fpld(10,ngena)
  real(c_double) :: redenc(nbctimea,ntotala)
  real(c_double) :: redenb(nbctimea,ntotala)
  real(c_double) :: tempc(nbctimea,ntotala)
  real(c_double) :: tempb(nbctimea,ntotala)



  !**********************************************************************
  !     Variables in common block diskx have as their last dimension lrza.
  !     Thus they are dimensioned with respect to the radial coordinate.
  !     If lrza=1  ==> single flux surface CQL run.
  !     Note that input variable setup0%lrz can be less than lrza. For some of
  !     the larger arrays we allocate space using setup0%lrz rather than
  !     dimensioning with lrza to save some space.
  !
  !     VECTORS
  !
  !**********************************************************************
  !
  !common /diskx/ &
  real(c_double) :: elecfld(0:lrza)
  real(c_double) :: tbnd(lrorsa)
  real(c_double) :: rovera(lrza)
  real(c_double) :: zmax(0:lrza)

  !..................................................................
  !     2-D ARRAYS
  !..................................................................

  !common /diskx/ &
  real(c_double) :: reden(ntotala,0:lrza)
  real(c_double) :: temp(ntotala,0:lrza)

  !common /diskx/ &
  real(c_double) :: tauegy(ngena,0:lrza)
  real(c_double) :: eparc(ngena,0:lrza)
  real(c_double) :: eperc(ngena,0:lrza)
  real(c_double) :: simpbfac

  real(c_double), pointer :: isoucof => null()
  real(c_double), pointer :: faccof => null()

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
  real(c_double), pointer :: mpwrsou(:) => null()
  real(c_double), pointer :: npwrsou(:) => null()
  real(c_double), pointer :: asor(:,:,:) => null()

  !common /params/ &
  integer, pointer :: nso => null()

  character(len=8), pointer :: pltso => null()
  character(len=8), pointer :: soucoord => null()
  character(len=8), pointer :: knockon => null()
  character(len=8), pointer :: komodel => null()
  character(len=8), pointer :: flemodel => null()

  !common /readscal/ &
  integer, pointer ::  nsou => null()
  integer, pointer :: nkorfn => null()
  integer, pointer :: nonko => null()
  integer, pointer :: noffko => null()
  real(c_double), pointer :: soffvte => null()
  real(c_double), pointer :: soffpr => null()
  real(c_double), pointer :: xlfac => null()
  real(c_double), pointer :: xlpctlwr => null()
  real(c_double), pointer :: xlpctmdl => null()
  real(c_double), pointer :: xllwr => null()
  real(c_double), pointer :: xlmdl => null()
  integer, pointer :: jfl => null()

  !common /readarr/ &
  integer, pointer :: nonso(:,:) => null()
  integer, pointer :: noffso(:,:) => null()
  real(c_double), pointer :: sellm1(:,:) => null()
  real(c_double), pointer :: sellm2(:,:) => null()
  real(c_double), pointer :: seppm1(:,:) => null()
  real(c_double), pointer :: seppm2(:,:) => null()
  real(c_double), pointer :: sem1(:,:) => null()
  real(c_double), pointer :: sem2(:,:) => null()
  real(c_double), pointer :: sthm1(:,:) => null()
  real(c_double), pointer :: scm2(:,:) => null()
  real(c_double), pointer :: szm1(:,:) => null()
  real(c_double), pointer :: szm2(:,:) => null()


  !*****************************************************************
  !     BEGIN arrays for rf package..(rf...,vlh[B,...,vlf...) routines
  !*****************************************************************
  character(len=8), pointer :: urfmod => null()
  character(len=8), pointer :: vlfmod => null()
  character(len=8), pointer :: vlfbes => null()
  character(len=8), pointer :: vlfnpvar => null()
  character(len=8), pointer :: vlhmod => null()
  character(len=8), pointer :: vprprop => null()
  character(len=8), pointer :: vlhplse => null()
  character(len=8), pointer :: vlhprprp(:) => null()
  character(len=8), pointer :: rfread => null()
  character(len=8), pointer :: rdcmod => null()
  character(len=8), pointer :: rdc_clipping => null()
  character(len=8), pointer :: rdc_netcdf => null()
  character(len=256), pointer :: rffile(:) => null()
  character(len=256), pointer :: rdcfile(:) => null()

  !common /params/ &
  integer, pointer :: nrf => null()
  !common/readscal/
  real(c_double), pointer :: vlhmodes => null()
  real(c_double), pointer :: vdalp => null()
  real(c_double), pointer :: vlh_karney => null()
  real(c_double), pointer :: vlhpon => null()
  real(c_double), pointer :: vlhpoff => null()
  real(c_double), pointer :: vlfmodes => null()
  real(c_double), pointer :: rdc_upar_sign => null()

  !common /readvec/ &
  integer, pointer :: nonrf(:) => null()
  integer, pointer :: noffrf(:) => null()
  real(c_double), pointer :: dlndau(:) => null()
  real(c_double), pointer :: vparmin(:) => null()
  real(c_double), pointer :: vparmax(:) => null()
  real(c_double), pointer :: vprpmin(:) => null()
  real(c_double), pointer :: vprpmax(:) => null()
  real(c_double), pointer :: vlhpolmn(:) => null()
  real(c_double), pointer :: vlhpolmx(:) => null()
  real(c_double), pointer :: vlffreq(:) => null()
  real(c_double), pointer :: vlfharms(:) => null()
  real(c_double), pointer :: vlfharm1(:) => null()
  real(c_double), pointer :: vlfnp(:) => null()
  real(c_double), pointer :: vlfdnp(:) => null()
  real(c_double), pointer :: vlfddnp(:) => null()
  real(c_double), pointer :: vlfpol(:) => null()
  real(c_double), pointer :: vlfdpol(:) => null()
  real(c_double), pointer :: vlfddpol(:) => null()
  real(c_double), pointer :: vlfnperp(:) => null()
  real(c_double), pointer :: vlfdnorm(:) => null()
  real(c_double), pointer :: vlfparmn(:) => null()
  real(c_double), pointer :: vlfparmx(:) => null()
  real(c_double), pointer :: vlfprpmn(:) => null()
  real(c_double), pointer :: vlfprpmx(:) => null()


  !..................................................................
  !     arrays in input
  !..................................................................

  !common/arr3d/ &
  real(c_double), pointer :: acoefne(:) => null()
  real(c_double), pointer :: acoefte(:) => null()
  real(c_double), pointer :: ryain(njenea) => null()
  real(c_double) :: elecin(njenea)
  real(c_double) :: enein(njenea,ntotala)
  real(c_double) :: tein(njenea)
  real(c_double) :: tiin(njenea)
  real(c_double) :: ennin(njenea,npaproca)
  integer :: irzplt(lrorsa)
  real(c_double) :: rya(0:lrza+1)
  real(c_double) :: thet1(nva),thet2(nva)
  real(c_double) :: thet1_npa(nva)
  real(c_double) :: thet2_npa(nva)
  integer :: nplt3d(nplota)

  !common/arr3d/ &
  real(c_double) :: enein_t(njenea,ntotala,nbctimea)
  real(c_double) :: tein_t(njenea,nbctimea)
  real(c_double) :: tiin_t(njenea,nbctimea)
  real(c_double) :: zeffin_t(njenea,nbctimea)
  real(c_double) :: elecin_t(njenea,nbctimea)
  real(c_double) :: xjin_t(njenea,nbctimea)
  real(c_double) :: vphiplin_t(njenea,nbctimea)
  real(c_double) :: ennin_t(njenea,nbctimea,npaproca) !neutrals,impurities,etc.

  !common/arr3d/ &
  real(c_double), pointer :: sellm1z(:,:,:) => null()
  real(c_double), pointer :: sellm2z(:,:,:) => null()
  real(c_double), pointer :: seppm1z(:,:,:) => null()
  real(c_double), pointer :: sem1z(:,:,:) => null()
  real(c_double), pointer :: sem2z(:,:,:) => null()
  real(c_double), pointer :: sthm1z(:,:,:) => null()
  real(c_double), pointer :: scm2z(:,:,:) => null()
  real(c_double), pointer :: szm1z(:,:,:) => null()
  real(c_double), pointer :: seppm2z(:,:,:) => null()
  real(c_double), pointer :: szm2z(:,:,:) => null()
  real(c_double), pointer :: asorz(:,:,:) => null()


  !****************************************************************
  !     BEGIN arrays for 3-d (td..) driver.
  !****************************************************************

  !common/params/ ndifus_io_t  !Max to be nbctimea
  integer, pointer :: ndifus_io_t => null()

  !..................................................................
  !     scalars in input
  !..................................................................

  character(len=8) :: bootst
  character(len=8) :: bootcalc
  character(len=8) :: bootupdt
  character(len=8) :: iprone
  character(len=8) :: iprote
  character(len=8) :: iproti
  character(len=8) :: iprozeff
  character(len=8) :: iprovphi
  character(len=8) :: iproelec
  character(len=8) :: ipronn
  character(len=8) :: iprocur
  character(len=8) :: tmdmeth
  character(len=8) :: partner
  character(len=8) :: plt3d
  character(len=8) :: pltvs
  character(len=8) :: radcoord
  character(len=8) :: rzset
  character(len=8) :: ndeltarho
  character(len=8) :: softxry
  character(len=8) :: npa_diag
  character(len=8) :: atten_npa
  character(len=8) :: efswtch
  character(len=8) :: efswtchn
  character(len=8) :: efiter
  character(len=8) :: efflag
  character(len=8), pointer :: pinch => null()
  character(len=8), pointer :: relaxtsp => null()
  character(len=8), pointer :: transp => null()
  character(len=8), pointer :: adimeth => null()

  character(len=8) :: npa_process(npaproca)

  character(len=256), pointer :: difus_io_file => null()

  !common /s3d/ &
  real(c_double), pointer :: difusr => null()
  real(c_double), pointer :: advectr => null()
  real(c_double) :: enmin
  real(c_double) :: enmax
  real(c_double) :: bootsign
  real(c_double) :: fds
  integer :: kfrsou
  integer :: mmsv
  integer :: msxr
  integer :: njene
  integer :: njte
  integer :: njti
  integer :: nonboot
  integer :: jhirsh
  integer :: nrskip
  integer :: nen
  integer :: nv
  integer :: nen_npa
  integer :: nv_npa
  integer :: npaproc
  integer :: nr_delta
  integer :: nz_delta
  integer :: nt_delta
  integer :: nr_f4d
  integer :: nz_f4d
  integer :: nv_f4d
  integer :: nt_f4d
  real(c_double) :: rfacz
  real(c_double) :: roveram
  real(c_double), pointer :: relaxden => null()
  real(c_double) :: enmin_npa
  real(c_double) :: enmax_npa
  real(c_double) :: fds_npa
  real(c_double) :: curr_edge
  real(c_double) :: efrelax
  real(c_double) :: efrelax1
  real(c_double) :: currerr
  integer, pointer :: nonadi

  !common/ar3d/ &
  real(c_double), pointer :: difus_rshape(:) => null()
  real(c_double), pointer :: difus_vshape(:) => null()
  real(c_double), pointer :: difin(:) => null()
  real(c_double) :: rd(nva)
  real(c_double) :: thetd(nva)
  real(c_double) :: x_sxr(nva)
  real(c_double) :: z_sxr(nva)
  real(c_double) :: rd_npa(nva)
  real(c_double) :: thetd_npa(nva)
  real(c_double) :: x_npa(nva)
  real(c_double) :: z_npa(nva)
  real(c_double), pointer :: difus_io_drrscale(:,:) => null()
  real(c_double), pointer :: difus_io_drscale(:,:) => null()
  real(c_double), pointer :: difus_io_t(:) => null()

  !******************************************************************
  !     BEGIN arrays for EQUILIBRIUM MODEL (eq..) (NON-CIRCULAR CROSS
  !     SECTIONS).
  !******************************************************************

  character(len=8), pointer :: nconteq => null()
  character(len=8), pointer :: eqmod => null()

  integer, pointer :: lfield => null()
  !common/params/
  integer, pointer :: nconteqn => null()

  character(len=8), pointer :: eqsym => null()
  character(len=8), pointer :: eqdskalt => null()
  character(len=8), pointer :: eqsource => null()
  character(len=8), pointer :: eqmodel=> null()
  character(len=8), pointer :: fpsimodl => null()

  !     ONETWO uses character*60 for eqdskin.
  character(len=256), pointer :: eqdskin => null()

  !common/readscal/ &
  real(c_double), pointer :: atol => null()
  real(c_double), pointer :: ellptcty => null()
  real(c_double), pointer :: eqpower => null()
  real(c_double), pointer :: bsign => null()
  real(c_double), pointer :: povdelp => null()
  real(c_double), pointer :: rtol => null()
  real(c_double), pointer :: rmag => null()
  real(c_double), pointer :: rbox => null()
  real(c_double), pointer :: rboxdst => null()
  real(c_double), pointer :: zbox => null()
  integer, pointer :: methflag => null()


  !*********************************************************************
  !     BEGIN arrays for LOWER HYBRID FAST WAVE and ECH Module.
  !*********************************************************************


  !common/params/ &
  integer, pointer :: nurftime => null()
  character(len=8), pointer :: urfdmp => null()
  character(len=8), pointer :: iurfcoll(:) => null()
  character(len=8), pointer :: iurfl(:) => null()
  character(len=8), pointer :: call_lh => null()
  character(len=8), pointer :: call_ech => null()
  character(len=8), pointer :: call_fw => null()
  character(len=8), pointer :: ech => null()
  character(len=8), pointer :: fw => null()
  character(len=8), pointer :: lh => null()
  character(len=8), pointer :: scaleurf => null()
  character(len=8), pointer :: urfrstrt => null()
  character(len=8), pointer :: urfwrray => null()
  character(len=8), pointer :: rftype(:) => null()

  !common/readscal/ &
  integer, pointer :: ieqbrurf => null()
  integer, pointer :: urfncoef => null()
  integer :: nmods ! xxx used in aindflt?
  integer, pointer :: nbssltbl => null()
  integer, pointer :: nondamp => null()
  integer, pointer :: nrfstep2 => null()
  integer, pointer :: nrfpwr => null()
  integer, pointer :: nrfitr1 => null()
  integer, pointer :: nrfitr2 => null()
  integer, pointer :: nrfitr3 => null()
  real(c_double), pointer :: urfmult => null()
  integer, pointer :: nrdc => null()

  !common/readvec/ &
  real(c_double), pointer :: pwrscale(:) => null()
  real(c_double), pointer :: wdscale(:) => null()
  integer, pointer :: nrfstep1(:) => null()
  integer, pointer :: nharms(:) => null()
  integer, pointer :: nharm1(:) => null()
  integer, pointer :: nrfspecies(:) => null()

  !common/readvec/ &
  real(c_double), pointer :: pwrscale1(:) => null()
  real(c_double), pointer :: urftime(:) => null()

  !common/readvec/ &
  real(c_double), pointer :: rdcscale(:) => null()
  integer, pointer :: nrdcspecies(:) => null()


  !-----------------------------------------------------------------------
  !     BEGIN variables for WP... modules for CQLP case
  !-----------------------------------------------------------------------

  character(len=8), pointer :: oldiag => null()
  character(len=8), pointer :: sbdry => null()
  character(len=8), pointer :: scheck => null()
  character(len=8), pointer :: ampfmod => null()
  character(len=8), pointer :: eseswtch => null()
  character(len=8), pointer :: updown => null()

  logical :: nlotp1(noutpta),nlotp2(noutpta)
  logical :: nlotp3(noutpta),nlotp4(noutpta)

  !common /readscal/ &
  real(c_double) :: epsthet
  real(c_double) :: elpar0
  integer :: lmidpln
  integer :: lmidvel
  integer :: laddbnd
  integer :: nchgdy
  integer :: ngauss
  integer :: nlagran
  integer :: nonavgf
  integer :: nofavgf
  integer :: nummods
  integer :: numixts
  integer, pointer :: nontran => null()
  integer, pointer :: nofftran => null()
  integer, pointer :: nonelpr => null()
  integer, pointer :: noffelpr => null()
  real(c_double), pointer :: ampferr => null()
  integer, pointer :: nampfmax => null()
  integer, pointer :: nonampf => null()

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
  integer :: nsteps_orb
  integer :: nptsorb
  integer :: i_orb_width
  integer :: iorb2
  integer :: j0_ini
  integer :: j0_end
  integer :: inc_j0
  integer :: i0_ini
  integer :: i0_end
  integer :: inc_i0
  integer :: j2_ini
  integer :: j2_end
  integer :: inc_j2
  integer :: i2_ini
  integer :: i2_end
  integer :: inc_i2
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
  integer :: idim
  integer :: iyjx
  integer :: iyjxp1
  integer :: iyp1jx
  integer :: iyjx2
  integer :: iyp1
  integer :: jxp1
  integer :: lrors
  integer :: mxp1
  integer :: mbet
  integer :: nonch
  integer :: niong
  integer :: nionm
  integer :: ntotal


  !**********************************************************************
  !     Variables in common block diskx have as their last dimension lrza.
  !     Thus they are dimensioned with respect to the radial coordinate.
  !     If lrza=1  ==> null()
  !     Note that input variable lrz can be less than lrza. For some of
  !     the larger arrays we allocate space using lrz rather than
  !     dimensioning with lrza to save some space.
  !
  !     VECTORS
  !
  !**********************************************************************
  !
  !common /diskx/ &
  real(c_double) :: bmod0(lrza)
  real(c_double) :: btor0(lrza)
  real(c_double) :: bthr(lrza)
  real(c_double) :: btoru(lrza)
  real(c_double) :: bthr0(lrza)
  real(c_double) :: consn(lrorsa)
  real(c_double) :: consn0(lrorsa)
  real(c_double) :: currt(lrza)
  real(c_double) :: currxj0(0:lrza)
  real(c_double) :: currtp(lrza)
  real(c_double) :: currmtp(lrorsa)
  real(c_double) :: currmt(lrorsa)
  real(c_double) :: currxj(0:lrza)
  real(c_double) :: currpar(lrza)
  real(c_double) :: curreq(lrza)
  real(c_double) :: zreshin(lrza)
  real(c_double) :: zreskim(lrza)
  real(c_double) :: rovsc(lrza)
  real(c_double) :: rovsc_hi(lrza)
  real(c_double) :: curtor(lrza)
  real(c_double) :: curpol(lrza)
  real(c_double) :: ccurtor(0:lrza)
  real(c_double) :: ccurpol(0:lrza)
  real(c_double) :: deltapsi(lrza)
  real(c_double) :: eps(lrza)
  real(c_double) :: etll0(lrza)
  integer :: itl_(lrorsa)
  real(c_double) :: itu_(lrorsa)
  real(c_double) :: iy_(lrorsa)
  real(c_double) :: iyh_(lrorsa)
  real(c_double) :: iyjx_(lrorsa)
  real(c_double) :: inew_(lrorsa)
  real(c_double) :: inewjx_(lrorsa)
  real(c_double) :: ieq_(lrorsa+1)
  real(c_double) :: indxlr(0:lrorsa)
  real(c_double) :: indxls(0:lrorsa)
  real(c_double) :: lorbit(lrza)
  real(c_double) :: lmdpln(0:lrza)
  real(c_double) :: n_(lrorsa)
  real(c_double) :: nch(lrorsa)
  real(c_double) :: nefiter_(lrza) ! counts iterations of el.field for each flux surface 
  real(c_double) :: psimx(lrza)
  real(c_double) :: pibzmax(lrza)
  real(c_double) :: psidz(lrza)
  real(c_double) :: qsafety(lrza)
  real(c_double) :: r0geom(lrza)
  real(c_double) :: r0drdz(0:lrza)
  real(c_double) :: rgeom(lrza)
  real(c_double) :: zgeom(lrza)
  real(c_double) :: rovs(lrza),rovsn(lrza),rovsloc(lrorsa)

  !common /diskx/ &
  real(c_double) :: sptzr(lrorsa)
  real(c_double) :: sgaint1(lrorsa)
  real(c_double) :: starnue(lrorsa)
  real(c_double) :: thb(lrorsa)
  real(c_double) :: tauee(lrza)
  real(c_double) :: taueeh(lrorsa)
  real(c_double) :: time_(lrorsa)
  real(c_double) :: twoint(lrorsa)
  real(c_double) :: vthe(lrza)
  real(c_double) :: xlbnd(lrza)
  real(c_double) :: xlndn0(lrza)
  real(c_double) :: zmaxpsi(0:lrza)
  real(c_double) :: zmaxi(0:lrza)
  real(c_double) :: zmaxpsii(0:lrza)
  real(c_double) :: zeff(lrza)
  real(c_double) :: zeff4(lrza)
  real(c_double) :: vphipl(lrza)
  real(c_double) :: srckotot(lrza)
  real(c_double) :: elecr(lrza)
  real(c_double) :: denfl(lrza)
  real(c_double) :: denfl1(lrza)
  real(c_double) :: denfl2(lrza)
  real(c_double) :: den_of_s(lza)
  real(c_double) :: den_of_s1(lza)
  real(c_double) :: den_of_s2(lza)
  real(c_double) :: tauii(lrza)
  real(c_double) :: tau_neo(lrza)
  real(c_double) :: drr_gs(lrza)
  real(c_double) :: rhol(lrza)
  real(c_double) :: rhol_pol(lrza)

  !common /diskx/ &
  real(c_double) :: taubi(lrza)
  real(c_double) :: tau_neo_b(lrza)
  real(c_double) :: drr_gs_b(lrza)
  real(c_double) :: rhol_b(lrza)
  real(c_double) :: rhol_pol_b(lrza)

  !..................................................................
  !     2-D ARRAYS
  !..................................................................

  !common /diskx/ &
  real(c_double) :: energy(ntotala,lrza)
  real(c_double) :: vth(ntotala,lrza)


  !common /diskx/ &
  real(c_double) :: den_fsa(ngena,lrza)
  real(c_double) :: den_fsa_t0(ngena,lrza)
  real(c_double) :: reden_t0(ngena,lrza)
  real(c_double) :: currm(ngena,lrorsa)
  real(c_double) :: curr(ngena,lrza)
  real(c_double) :: hnis(ngena,lrza)
  real(c_double) :: ratio(ngena,lrza)
  real(c_double) :: wpar(ngena,lrza)
  real(c_double) :: wperp(ngena,lrza)
  real(c_double) :: xlndn00(ngena,lrza)
  real(c_double) :: xlncur(ngena,lrza)
  real(c_double) :: xlndn(ngena,lrza)
  real(c_double) :: energym(ngena,lrorsa)
  real(c_double) :: xlndnr(ntotala,lrza)
  real(c_double) :: energyr(ntotala,lrza)
  real(c_double) :: currr(ngena,lrza)
  real(c_double) :: xlndnv(ntotala,lrza)
  real(c_double) :: energyv(ntotala,lrza)
  real(c_double) :: currv_(ngena,lrza)
  real(c_double) :: eoe0(ngena,lrza)
  real(c_double) :: ucrit(ngena,lrza)
  integer :: jxcrit(ngena,lrza),      &
       jchang(ngena,lrorsa)
  real(c_double) :: denra(ngena,lrza)
  real(c_double) :: curra(ngena,lrza)
  real(c_double) :: fdenra(ngena,lrza)
  real(c_double) :: fcurra(ngena,lrza)
  real(c_double) :: enn(lrza,npaproca)
  !common /diskx/ &
  real(c_double) :: alm(0:mbeta,lrza)
  real(c_double) :: betta(0:mbeta,lrza)
  real(c_double) :: entrintr(ngena,-1:15)

  !..................................................................
  !     Next follow variables and arrays not dimensioned by lrza.
  !..................................................................

  !..................................................................
  !     SCALARS.....
  !..................................................................

  character(len=8) :: outfile
  character(len=8) :: prefix

  character(len=512) ::  t_

  !common &
  real(c_double) :: clight
  real(c_double) :: charge
  real(c_double) :: restmkev
  real(c_double) :: symm
  real(c_double) :: cnorm
  real(c_double) :: cnorm2
  real(c_double) :: cnorm3
  real(c_double) :: clite2
  real(c_double) :: cnorm2i
  real(c_double) :: cnormi
  real(c_double) :: ergtkev
  real(c_double) :: eps0
  real(c_double) :: elect
  real(c_double) :: eratio
  real(c_double) :: eovsz
  real(c_double) :: eovedd
  real(c_double) :: em6
  real(c_double) :: em8
  real(c_double) :: em10
  real(c_double) :: em12
  real(c_double) :: em14
  real(c_double) :: em37
  real(c_double) :: em40
  real(c_double) :: ep37
  real(c_double) :: em90
  real(c_double) :: ep90
  real(c_double) :: em100
  real(c_double) :: em120
  real(c_double) :: em300
  real(c_double) :: four
  real(c_double) :: fourpi
  real(c_double) :: fc
  real(c_double) :: gszb
  real(c_double) :: gsb
  real(c_double) :: gsem
  real(c_double) :: gsep
  real(c_double) :: gszm
  real(c_double) :: gszp
  real(c_double) :: gseb
  real(c_double) :: gszm2
  real(c_double) :: gszp2
  real(c_double) :: gszb2
  real(c_double) :: gsla2
  real(c_double) :: gslb2
  real(c_double) :: gsnm
  real(c_double) :: half
  integer :: iyh
  integer :: istate
  integer :: itl
  integer :: itu
  integer :: impcoef
  integer :: imprf
  integer :: impelec
  integer :: impcah
  integer :: irstart
  integer :: impadi
  integer :: jxm1
  integer :: jlwr
  integer :: jmdl
  integer :: kgnande
  integer :: kelecg
  integer :: kelecm
  integer :: kelec
  integer :: kionn
  integer :: kiong(ngena)
  integer :: kiongg(ngena)
  integer :: kionm(nmaxa)
  integer :: l_
  integer :: lr_
  integer :: indxlr_
  integer :: indxls_
  integer :: lmdpln_
  integer :: ls_
  integer :: miyjx
  integer :: ibeampon
  integer :: ibeamponp

  !common &
  integer :: ipad1
  real(c_double) :: r0geomp
  real(c_double) :: rgeomp
  real(c_double) :: zgeomp
  real(c_double) :: rgeom1
  real(c_double) :: rgeom2
  integer :: niyjx
  integer :: nccoef
  integer :: ncentr
  integer :: navnc
  integer :: ncplt
  integer :: n
  integer :: nflux
  integer :: nframe
  integer :: npltmx
  real(c_double) :: ipad2
  real(c_double) :: one
  real(c_double) :: one_
  real(c_double) :: pi
  real(c_double) :: pio180
  real(c_double) :: pio2
  real(c_double) :: psir
  real(c_double) :: proton
  real(c_double) :: rgyro
  real(c_double) :: resist
  real(c_double) :: resistn
  real(c_double) :: rovsf
  real(c_double) :: rmp0
  real(c_double) :: rovscf
  real(c_double) :: rbgn
  real(c_double) :: stopm
  real(c_double) :: tei
  real(c_double) :: timet
  real(c_double) :: three
  real(c_double) :: third
  real(c_double) :: sevenhun
  real(c_double) :: two
  real(c_double) :: twopi
  real(c_double) :: vparl
  real(c_double) :: vpovth
  real(c_double) :: xconn
  real(c_double) :: xmax
  integer :: ipxy,jpxy
  real(c_double) :: zero
  real(c_double) :: zstar
  integer :: iplot
  integer :: nplott
  integer :: isave
  integer :: nsavet
  integer :: nirzplt
  integer :: nefiter ! counts iterations for el.field (efswtch="method5")
  !common &
  real(c_double) :: elecfldc
  real(c_double) :: elecfldb !V/cm, Added for ampfmod case
  integer :: it_ampf  !iteraction counter for ampfmod case !YuP[2019-04-24]

  !..................................................................
  !     VECTORS.....
  !..................................................................

  ! bug
  character(len=8) :: text(4)

  !common &
  integer :: ipad3, iota(0:6400)
  !BH111124     1  iota(0:750),realiota(0:750),  real(c_double) :: realiota(0:6400), tab(3)
  integer ::i1p(2),itab(3),ipad4,       imsh(20)
  real(c_double) :: frbuf(1024)
  real(c_double) :: ws2d(2500)
  real(c_double) :: work(tlfld1a)
  !common &
  real(c_double) :: ekev(ngena)
  real(c_double) :: engain(ngena)
  real(c_double) :: tnorm(ngena)
  !common &
  real(c_double) :: fions(ntotala)
  real(c_double) :: gama0(ntotala)
  real(c_double) :: tz1(lza)
  real(c_double) :: tz2(lza)

  !common &
  !!real(c_double) :: wkbc(3*nbctimea),iopbc(2) ! YuP[2019-04-24] Not used

  character(len=8) :: mplot(lrorsa)

  !common &
  real(c_double) :: tplot(nplota)
  real(c_double) :: tsave(nsavea)

  !..................................................................
  !     TWO-DIMENSIONAL ARRAYS.....
  !..................................................................

  !common &
  real(c_double) :: gamt(ntotala,ntotala)
  real(c_double) :: gama(ntotala,ntotala)
  real(c_double) :: satioz2(ntotala,ntotala)
  real(c_double) :: satiom(ntotala,ntotala)
  real(c_double) :: sgain(8,ngena)

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
  real(c_double), pointer :: bpolz(:,:), btorz(:,:)  ! (lza,setup0%lrzmax)
  !common /dptr95/ bpolz, btorz

  ! YuP: [added Apr/2014] area and volume of a cell associated with each
  !                     (R,Z) point on flux surface, (R,Z)==(solrz,solzz)
  real(c_double), pointer :: ddarea(:,:), ddvol(:,:)  !  (lza,setup0%lrzmax)
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
  real(c_double), pointer :: mun(:) ! real is correct
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
  real(c_double), pointer :: sounor(:,:,:,:)   !!!(ngen,nsoa,lz,setup0%lrz)
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
  real(c_double) :: bdre(lrza)
  real(c_double) :: bdrep(lrza)
  real(c_double) :: sorpwt(lrza)
  real(c_double) :: sorpwti(0:lrza)
  real(c_double) :: sorpw_nbii(1:ngena,0:lrza)
  real(c_double) :: sorpw_rfi(1:ngena,0:lrza)
  real(c_double) :: xlncurt(lrza)

  !common /diskx/ &
  real(c_double) :: sorpw_rf(1:ngena,lrza)
  real(c_double) :: sorpw_nbi(1:ngena,lrza)

  !common /diskx/ &
  real(c_double) :: cosm1(ngena,nsoa,lrza)
  real(c_double) :: cosm2(ngena,nsoa,lrza)
  real(c_double) :: sxllm1(ngena,nsoa,lrza)
  real(c_double) :: sxllm2(ngena,nsoa,lrza)
  real(c_double) :: sxppm1(ngena,nsoa,lrza)
  real(c_double) :: sxppm2(ngena,nsoa,lrza)
  real(c_double) :: xem1(ngena,nsoa,lrza)
  real(c_double) :: xem2(ngena,nsoa,lrza)
  real(c_double) :: zm1(ngena,nsoa,lrza)
  real(c_double) :: zm2(ngena,nsoa,lrza)


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

  real(c_double) :: li  ! real is correct
  !common/sc3d/ &
  real(c_double) :: dttr
  real(c_double) :: dtreff
  real(c_double) :: currtza
  real(c_double) :: currtpza
  real(c_double) :: conserv
  real(c_double) :: fom
  real(c_double) :: fomp
  real(c_double) :: fompla
  real(c_double) :: fomtot
  real(c_double) :: flxout
  real(c_double) :: eden
  real(c_double) :: edenlavg
  real(c_double) :: etemp
  real(c_double) :: ethtemp
  real(c_double) :: edntmp
  real(c_double) :: pden
  real(c_double) :: pdntmp
  real(c_double) :: psynct
  !BH110314cBH070408(Not used except in fr routines):    1  smooth,
  !BH110314:  Restored, as needed in subroutine frsmooth, but needed
  !BH110314:  to change name, as conflict arises in tdreadf/tdwritef &

  real(c_double) :: smooth_
  real(c_double) :: toteqd
  real(c_double) :: cursign
  real(c_double) :: totcurza
  real(c_double) :: total
  real(c_double) :: total0
  real(c_double) :: totcurtt
  real(c_double) :: curxjtot
  integer :: ncount
  integer :: iplt3d
  integer :: ipacktp
  integer :: n_d_rr

  !..................................................................
  !     all arrays used only in CQL3D
  !..................................................................

  real(c_double) :: jparb(lrza)
  real(c_double) :: jparbt(lrza)
  real(c_double) :: jparbp(lrza) ! real is correct
  !common/ar3d/
  real(c_double) :: rrz(0:lrza)
  real(c_double) :: tr(0:lrza)
  real(c_double) :: tr1(0:lrza)
  real(c_double) :: tr2(0:lrza)
  real(c_double) :: tr3(0:lrza)
  real(c_double) :: tr4(0:lrza)
  real(c_double) :: tr5(0:lrza)
  real(c_double) :: drp5(0:lrza)
  real(c_double) :: dpsi(0:lrza)
  real(c_double) :: dpsidrho(lrza)
  integer :: iytr(lrorsa)
  real(c_double) ::h_r(0:lrza)
  real(c_double) :: area(0:lrza)
  real(c_double) :: equilpsp(0:lrza)
  real(c_double) :: equilpsi(0:lrza)
  real(c_double) :: areamid(0:lrza)
  real(c_double) :: volmid(0:lrza)
  real(c_double) :: psivalm(lrza)
  real(c_double) :: rpconz(lrza)
  real(c_double) :: rmconz(lrza)
  real(c_double) :: rpmconz(0:lrza)
  real(c_double) :: bpolsqaz(0:lrza)
  real(c_double) :: aspin(lrza)
  real(c_double) :: trapfrac(lrza)
  real(c_double) :: currz(ngena
  real(c_double) :: lrza)
  real(c_double) :: currtpz(lrza)
  real(c_double) :: currtz(lrza)
  real(c_double) :: currza(ngena)
  real(c_double) :: currtzi(0:lrza)
  real(c_double) :: currtpzi(0:lrza)
  real(c_double) :: currmtz(lrorsa)
  real(c_double) :: currmtpz(lrorsa)
  real(c_double) :: totcurzi(0:lrza)
  real(c_double) :: totcurz(lrza)
  real(c_double) :: fpsiz2(0:lrza)
  real(c_double) :: fpsiz(lrza)
  real(c_double) :: ffpsiz(lrza)
  real(c_double) :: prestp(lrza)
  real(c_double) :: prest(lrza)
  real(c_double) :: d2prest(lrza)
  real(c_double) :: d2fpsiz(lrza)
  real(c_double) :: d2ffpsiz(lrza)
  real(c_double) :: bmdplne(lrza)
  real(c_double) :: d2bmdpl(lrza)
  real(c_double) :: gkpwrz(ngena,lrza)
  real(c_double) :: rfpwrz(ngena,lrza)
  real(c_double) :: drrt(ngena)
  real(c_double) :: drt(ngena)
  !common/ar3d/ &
  real(c_double) :: dvol(lrza)
  real(c_double) :: darea(lrza)
  real(c_double) :: psyncz(lrza)
  real(c_double) :: pegyz(ngena,lrza)
  real(c_double) :: pplossz(ngena,lrza)
  real(c_double) :: wparzt(ngena)
  real(c_double) :: wperpzt(ngena)
  real(c_double) :: pegyt(ngena)
  real(c_double) :: pplosst(ngena)
  real(c_double) :: rfpwrt(ngena)
  real(c_double) :: gkpwrt(ngena)
  real(c_double) :: energyt(ntotala)
  real(c_double) :: rz(0:lrza)
  real(c_double) :: vfluxz(lrorsa)
  real(c_double) :: vol(0:lrza)
  real(c_double) :: onovrpz(lrza,2)
  real(c_double) :: tplt3d(nplota)
  real(c_double) :: bscurma(2,2)
  real(c_double) :: bscurm(0:lrza,2,2)
  real(c_double) :: bscurmi(0:lrza,2,2)

  !common/sc3d/ &
  real(c_double) ::     sorpwtza

  !common /csxr/
  real(c_double) :: sxry(lrza,4)
  real(c_double) :: sang(lrza,4)
  real(c_double) :: spol(lrza,4)
  integer :: ibin(lrza,4)
  real(c_double) ::eflux(nena,nva)
  real(c_double) :: efluxt(nva)
  real(c_double) :: alphad(3)
  real(c_double) :: xs_(3)
  real(c_double) :: enk(nena)
  real(c_double) :: en_(nena)
  integer :: jval_(nena)
  integer :: inegsxr(nva)
  integer :: lensxr(nva)

  !common /csigma/
  integer :: mtab
  integer :: msig
  integer :: jxis
  real(c_double) :: elmin
  real(c_double) :: delegy
  integer :: imaxwln(2,4)
  real(c_double) :: igenrl(2,4)
  real(c_double) :: sigm(4,lrorsa)
  real(c_double) :: sigf(4,lrorsa)
  real(c_double) :: sigmt(4)
  real(c_double) :: sigft(4)
  real(c_double) :: fuspwrv(4,lrorsa)
  real(c_double) :: fuspwrvt(4)
  real(c_double) :: fuspwrm(4,lrorsa)
  real(c_double) :: fuspwrmt(4)

  real(c_double), pointer :: tamm1(:)
  !common /csigma/ tamm1  !(0:mmsv)

  integer, pointer :: iind(:) ! YuP[2019-04-24] integer
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
  real(c_double) :: areacon(lrza)
  real(c_double) :: bmidplne(lrza)
  real(c_double) :: bpolsqa(lrza)
  real(c_double) :: eqdells(lrza)
  real(c_double) :: epsicon(lrza)
  real(c_double) :: erhocon(lrza)
  real(c_double) :: fpsi(lrza)
  real(c_double) :: flxavgd(lrza)
  real(c_double) :: psiovr(lrza)
  real(c_double) :: psiavg(2,lrza)
  real(c_double) :: onovrp(2,lrza)
  real(c_double) :: onovpsir3(lrza)
  real(c_double) :: rpcon(lrza)
  real(c_double) :: rmcon(lrza)
  real(c_double) :: zpcon(lrza)
  real(c_double) :: zmcon(lrza)
  real(c_double) :: volcon(lrza)
  real(c_double) :: fppsi(lrza)
  real(c_double) :: pppsi(lrza)
  real(c_double) :: es_bmax(lrza)
  real(c_double) :: bpsi_max(lrza)
  real(c_double) :: bpsi_min(lrza)
  real(c_double) :: lbpsi_max(lrza)
  real(c_double) :: z_bmax(lrza)
  real(c_double) :: bpsi_z_bmax(lrza)
  real(c_double) :: dlpgpsii(lrza)
  real(c_double) :: dlpsii(lrza)
  integer :: lz_bmax(lrza)
  integer :: lbpsi_min(lrza)

  !common/params/
  integer :: nnz
  integer :: nnr
  integer :: nj12


  character(len=8) :: eqorb
  character(len=8) :: eqcall

  !common &
  real(c_double) :: bpolsqlm
  integer :: imag
  integer :: jmag
  integer :: nmag
  integer :: nrc
  integer :: nzc
  integer :: nfp
  integer :: nnv
  integer :: iupdn
  real(c_double) ::psimag
  real(c_double) :: psilim
  real(c_double) :: rmaxcon
  real(c_double) :: rmincon
  real(c_double) :: rhomax
  real(c_double) :: zmag
  real(c_double) :: zmaxcon
  real(c_double) :: zmincon
  real(c_double) :: zshift

  !common &
  integer :: ibd(4)
  real(c_double) :: eqpsi(nconteqa)
  real(c_double) :: eqvol(nconteqa)
  real(c_double) :: eqfopsi(nconteqa)
  real(c_double) :: q_(nconteqa)
  real(c_double) :: eqrho(nconteqa)
  real(c_double) :: d2eqrho(nconteqa)
  real(c_double) :: eqarea(nconteqa)
  real(c_double) :: eqrpcon(nconteqa)
  real(c_double) :: eqrmcon(nconteqa)
  real(c_double) :: eqzpcon(nconteqa)
  real(c_double) :: eqzmcon(nconteqa)
  real(c_double) :: d2fpsiar(nnra)
  real(c_double) :: ez(nnza)
  real(c_double) :: dummyaz(nnza)
  real(c_double) :: fpsiar(nnra)
  real(c_double) :: ffpar(nnra)
  real(c_double) :: d2ffpar(nnra)
  real(c_double) :: qar(nnra)
  real(c_double) :: d2qar(nnra)
  real(c_double) :: prar(nnra)
  real(c_double) :: d2prar(nnra)
  real(c_double) :: ppar(nnra)
  real(c_double) :: d2ppar(nnra)
  real(c_double) :: psiar(nnra)
  real(c_double) :: er(nnra)
  real(c_double) :: dummyar(nnra)
  real(c_double) :: wkepsi(nrz3p1a)
  real(c_double) :: tlorb1(lfielda)
  real(c_double) :: tlorb2(lfielda)

  !common &
  real(c_double) :: epsi(nnra,nnza)
  real(c_double) :: epsirr(nnra,nnza)
  real(c_double) :: epsizz(nnra,nnza)
  real(c_double) :: epsirz(nnra,nnza)
  real(c_double) :: dummypsi(nnra,nnza)
  real(c_double) :: eqovrp(nconteqa,2)

  !common/output/
  integer :: lorbit_
  integer :: ialign14
  real(c_double) :: rmcon_
  real(c_double) :: rpcon_
  real(c_double) :: zmcon_
  real(c_double) :: zpcon_
  real(c_double) :: bthr_
  real(c_double) :: btoru_
  real(c_double) :: eqdells_
  real(c_double) :: fpsi_
  real(c_double) :: fppsi_
  real(c_double) :: zmax_
  real(c_double) :: btor0_
  real(c_double) :: bthr0_
  real(c_double) :: es_bmax_
  real(c_double) :: bpsi_max_
  real(c_double) :: bpsi_min_
  integer :: lbpsi_max_
  real(c_double) :: lbpsi_min_
  real(c_double) :: bmidplne_
  real(c_double) :: solr_(lfielda)
  real(c_double) :: solz_(lfielda)
  real(c_double) :: es_(lfielda)
  real(c_double) :: eqbpol_(lfielda)
  real(c_double) :: bpsi_(lfielda)
  real(c_double) :: thtpol_(lfielda)
  real(c_double) :: eqdell_(lfielda)

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

  real(c_double), pointer :: jbm1(:,:) !YuP[2019-04-24] real is correct
  !common jbm1
  real(c_double), pointer :: jb0(:,:)  !YuP[2019-04-24] real is correct
  !common jb0
  real(c_double), pointer :: jbp1(:,:) !YuP[2019-04-24] real is correct
  !common jbp1

  character(len=8) :: irffile(nmodsa)

  !common &
  real(c_double) :: argmax
  integer :: nurf ! YuP[2019-04-24] integer
  !     1  lenj0,  integer :: mrf,mrfn,  &
  integer :: nrayn
  integer :: nrayelts
  integer :: !-YuP 101122: added
  integer ::  irftype
  real(c_double) :: powray
  real(c_double) :: vnorm2
  real(c_double) :: vnorm3
  real(c_double) :: vnorm4
  real(c_double) :: dveps ! YuP[04-2016] for subr. urfb_add

  complex(c_double_complex), pointer :: cosz1(:),sinz1(:),sinz2(:)
  !common &
  real(c_double) :: bsslstp(nmodsa)
  integer  :: ncontrib(lrza)
  real(c_double) :: powrf(lrza,nmodsa)
  real(c_double) :: powrfc(lrza,nmodsa)
  real(c_double) :: powrfl(lrza,nmodsa)
  real(c_double) :: powrft(lrza)
  real(c_double) :: powurf(0:nmodsa)
  real(c_double) :: powurfc(0:nmodsa)
  real(c_double) :: powurfl(0:nmodsa)
  real(c_double) :: powurfi(0:lrza,0:nmodsa)

  !common &
  real(c_double) :: freqcy(nmodsa)
  real(c_double) :: omega(nmodsa)
  integer :: nharm(nmodsa)
  integer :: nray(nmodsa)
  real(c_double) :: bsign1(nmodsa)
  integer :: krfn(nmodsa)
  integer :: irfn(nmodsa)
  integer :: irfm(nmodsa)

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
  real(c_double), pointer :: psiloc(:,:,:) !YuP[2019-04-24] real
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
  integer, pointer :: nurefls(:,:) !YuP[2019-04-24] integer
  !common /dptr95/ nurefls
  integer, pointer :: keiks(:,:) !YuP[2019-04-24] integer
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
  integer, pointer :: nrayelt0(:,:) !YuP[2019-04-24] integer
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
  integer :: lapacki
  integer :: lapackj
  integer :: icsrij
  integer :: icsrip
  integer :: icsri2
  integer :: krylov
  integer :: icsrikry
  integer :: iwk_ilu
  !common / it3d/
  integer ::icsrijr
  integer :: icsrijc

  real(c_double), pointer :: abd_lapack
  real(c_double), pointer :: a_csr
  real(c_double), pointer :: alu
  real(c_double), pointer :: w_ilu
  real(c_double), pointer :: rhs0
  real(c_double), pointer :: sol
  real(c_double), pointer :: vv
  integer, pointer :: ja_csr
  integer, pointer :: ia_csr
  integer, pointer :: jlu
  integer, pointer :: ju
  integer, pointer :: jw_ilu
  real(c_double), pointer :: ar_csr
  real(c_double), pointer :: ac_csr
  integer, pointer :: jar_csr
  integer, pointer :: iar_csr
  integer, pointer :: ipofi
  integer, pointer :: jac_csr
  integer, pointer :: iac_csr
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
  real(c_double) :: anecc(nrada)
  real(c_double) :: tekev(nrada)
  real(c_double) :: tikev(nint1a,nrada)
  real(c_double) :: anicc(nint1a,nrada)
  real(c_double) :: amass(nint1a)
  real(c_double) :: achrg(nint1a)
  real(c_double) :: elecf(nrada)
  real(c_double) :: rho_(nrada)
  real(c_double) :: psiar_(nrada)
  character(len=80) :: names(10) !YuP[2019-04-24]
  integer :: nspc
  integer :: ialign9
  real(c_double) :: fpsiar_(nrada)
  real(c_double) :: pary(nrada)
  real(c_double) :: ppary(nrada)
  real(c_double) :: gpary(nrada)
  real(c_double) :: ztop_
  real(c_double) :: zbot_
  real(c_double) :: rleft
  real(c_double) :: rright
  integer :: nx
  integer :: nz
  integer :: npsitm

  !-----------------------------------------------------------------------
  !     BEGIN variables for WP... modules for CQLP case
  !-----------------------------------------------------------------------

  character(len=8) :: analegco

  !common /wpscal/ &
  integer :: iymax
  integer :: nsleft
  integer :: nsrigt
  integer :: numclas
  integer :: numindx

  !common /wpvec/ &
  real(c_double) :: cofdfds(0:lsa1,2,ntrmdera)
  real(c_double) :: enrgypa(ntotala,0:lsa1)
  real(c_double) :: vthpar(ntotala,0:lsa1)
  integer :: lsbtopr(0:lsa1)
  integer :: lsprtob(0:lsa1)
  integer :: lpm1eff(0:lsa1,-1:+1)
  real(c_double) :: sz(0:lsa1)
  real(c_double) :: dsz(0:lsa1)
  real(c_double) :: dszm5(0:lsa1)
  real(c_double) :: dszp5(0:lsa1)
  real(c_double) :: eszm5(0:lsa1)
  real(c_double) :: eszp5(0:lsa1)
  real(c_double) :: psis(0:lsa1)
  real(c_double) :: psisp(0:lsa1)
  real(c_double) :: psipols(0:lsa1)
  real(c_double) :: solrs(0:lsa1)
  real(c_double) :: solzs(0:lsa1)
  real(c_double) :: elparol(0:lsa1)
  real(c_double) :: elparnw(0:lsa1)
  real(c_double) :: flux1(0:lsa1)
  real(c_double) :: flux2(0:lsa1)

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

  real(c_double), pointer :: rcontr(:)
  real(c_double), pointer :: zcontr(:)
  real(c_double), pointer :: rlimiter(:)
  real(c_double), pointer :: zlimiter(:)
  !common/limiter/
  integer :: ncontr
  integer :: nlimiter
  ! Setup by call equilib()

  !common/eqbox/
  real(c_double) :: ermin
  real(c_double) :: ermax
  real(c_double) :: ezmin
  real(c_double) :: ezmax  ! Limits of equilibrium
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
  !$$$      real(c_double) mcq_drdt
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
  !$$$C---> '0' - values at the launching point of g.c.orbit [cgs]
  !$$$C---> cosp0 is cos(pitch-angle) at launching point.
  !$$$C---> Rmid is R at midplane for given orbit found by Newton iterations
  !.......................................................................
  !     variables transferred from freya
  !.......................................................................
  character(len=8) :: frmodp
  character(len=8) :: fr_gyrop
  character(len=8) :: beamplsep
  integer :: mfm1p
  real(c_double) :: beamponp
  real(c_double) :: beampoffp
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
  real(c_double) :: adv_norm(lrza)
  real(c_double) :: reden_norm(lrza)
  !COMMON /newtv/
  integer, parameter :: NP=300
  real(c_double) :: fvec(NP)
  integer :: newtv_nn

  save
contains

  subroutine initialize_rf_pointers()
    if(initialized_rf_pointers) call abort
    vlfemin => rfsetup%vlfemin
    vlfeplus => rfsetup%vlfeplus
    urfmod => rfsetup%urfmod
    vlfmod => rfsetup%vlfmod
    vlfbes => rfsetup%vlfbes
    vlfnpvar => rfsetup%vlfnpvar
    vlhmod => rfsetup%vlhmod
    vprprop => rfsetup%vprprop
    vlhplse => rfsetup%vlhplse
    vlhprprp => rfsetup%vlhprprp
    rfread => rfsetup%rfread
    rdcmod => rfsetup%rdcmod
    rdc_clipping => rfsetup%rdc_clipping
    rdc_netcdf => rfsetup%rdc_netcdf
    rffile => rfsetup%rffile
    rdcfile => rfsetup%rdcfile
    nrf => rfsetup%nrf
    vlhmodes => rfsetup%vlhmodes
    vdalp => rfsetup%vdalp
    vlh_karney => rfsetup%vlh_karney
    vlhpon => rfsetup%vlhpon
    vlhpoff => rfsetup%vlhpoff
    vlfmodes => rfsetup%vlfmodes
    rdc_upar_sign => rfsetup%rdc_upar_sign
    nonrf => rfsetup%nonrf
    noffrf => rfsetup%noffrf
    dlndau => rfsetup%dlndau
    vparmin => rfsetup%vparmin
    vparmax => rfsetup%vparmax
    vprpmin => rfsetup%vprpmin
    vprpmax => rfsetup%vprpmax
    vlhpolmn => rfsetup%vlhpolmn
    vlhpolmx => rfsetup%vlhpolmx
    vlffreq => rfsetup%vlffreq
    vlfharms => rfsetup%vlfharms
    vlfharm1 => rfsetup%vlfharm1
    vlfnp => rfsetup%vlfnp
    vlfdnp => rfsetup%vlfdnp
    vlfddnp => rfsetup%vlfddnp
    vlfpol => rfsetup%vlfpol
    vlfdpol => rfsetup%vlfdpol
    vlfddpol => rfsetup%vlfddpol
    vlfnperp => rfsetup%vlfnperp
    vlfdnorm => rfsetup%vlfdnorm
    vlfparmn => rfsetup%vlfparmn
    vlfparmx => rfsetup%vlfparmx
    vlfprpmn => rfsetup%vlfprpmn
    vlfprpmx => rfsetup%vlfprpmx
    nurftime => rfsetup%nurftime
    urfdmp => rfsetup%urfdmp
    iurfcoll => rfsetup%iurfcoll
    iurfl => rfsetup%iurfl
    call_lh => rfsetup%call_lh
    call_ech => rfsetup%call_ech
    call_fw => rfsetup%call_fw
    ech => rfsetup%ech
    fw => rfsetup%fw
    lh => rfsetup%lh
    scaleurf => rfsetup%scaleurf
    urfrstrt => rfsetup%urfrstrt
    urfwrray => rfsetup%urfwrray
    rftype => rfsetup%rftype
    ieqbrurf => rfsetup%ieqbrurf
    urfncoef => rfsetup%urfncoef
    nbssltbl => rfsetup%nbssltbl
    nondamp => rfsetup%nondamp
    nrfstep2 => rfsetup%nrfstep2
    nrfpwr => rfsetup%nrfpwr
    nrfitr1 => rfsetup%nrfitr1
    nrfitr2 => rfsetup%nrfitr2
    nrfitr3 => rfsetup%nrfitr3
    urfmult => rfsetup%urfmult
    nrdc => rfsetup%nrdc
    pwrscale => rfsetup%pwrscale
    wdscale => rfsetup%wdscale
    nrfstep1 => rfsetup%nrfstep1
    nharms => rfsetup%nharms
    nharm1 => rfsetup%nharm1
    nrfspecies => rfsetup%nrfspecies
    pwrscale1 => rfsetup%pwrscale1
    urftime => rfsetup%urftime
    rdcscale => rfsetup%rdcscale
    nrdcspecies => rfsetup%nrdcspecies

    initialized_rf_pointers = .TRUE.
  end subroutine initialize_rf_pointers

  subroutine initialize_eq_pointers
    if(initialized_eq_pointers) call abort
    nconteq => eqsetup%nconteq
    eqmod => eqsetup%eqmod
    lfield => eqsetup%lfield
    nconteqn => eqsetup%nconteqn
    eqsym => eqsetup%eqsym
    eqdskalt => eqsetup%eqdskalt
    eqsource => eqsetup%eqsource
    eqmodel=> eqsetup%eqmodel
    fpsimodl => eqsetup%fpsimodl
    eqdskin => eqsetup%eqdskin
    atol => eqsetup%atol
    ellptcty => eqsetup%ellptcty
    eqpower => eqsetup%eqpower
    bsign => eqsetup%bsign
    povdelp => eqsetup%povdelp
    rtol => eqsetup%rtol
    rmag => eqsetup%rmag
    rbox => eqsetup%rbox
    rboxdst => eqsetup%rboxdst
    zbox => eqsetup%zbox
    methflag => eqsetup%methflag

    initialized_eq_pointers = .TRUE.
  end subroutine initialize_eq_pointers

  subroutine initialize_tr_pointers
    if(initialized_tr_pointers) call abort
    advectr => trsetup%advectr
    difus_type => trsetup%difus_type
    difusr => trsetup%difusr
    difus_rshape => trsetup%difus_rshape
    difus_vshape => trsetup%difus_vshape
    difin => trsetup%difin
    difus_io => trsetup%difus_io
    difus_io_file => trsetup%difus_io_file
    difus_io_drrscale => trsetup%difus_io_drrscale
    difus_io_drscale => trsetup%difus_io_drscale
    difus_io_t => trsetup%difus_io_t
    pinch => trsetup%pinch
    relaxden => trsetup%relaxden
    relaxtsp => trsetup%relaxtsp
    transp => trsetup%transp
    adimeth => trsetup%adimeth
    nonadi => trsetup%nonadi
    nontran => trsetup%nontran
    nofftran => trsetup%nofftran
    nonelpr => trsetup%nonelpr
    noffelpr => trsetup%noffelpr
    ndifus_io_t => trsetup%ndifus_io_t
    initialized_tr_pointers = .TRUE.
  end subroutine initialize_tr_pointers

  subroutine initialize_sou_pointers
    if(initialized_sou_pointers) call abort
    asorz => sousetup%asorz
    asor => sousetup%asor
    flemodel => sousetup%flemodel
    nonso => sousetup%nonso
    noffso => sousetup%noffso
    nso => sousetup%nso
    nsou => sousetup%nsou
    pltso => sousetup%pltso
    mpwrsou => sousetup%mpwrsou
    npwrsou => sousetup%npwrsou
    scm2z => sousetup%scm2z
    szm1z => sousetup%szm1z
    scm2 => sousetup%scm2
    sellm1 => sousetup%sellm1
    sellm2 => sousetup%sellm2
    seppm1 => sousetup%seppm1
    sellm1z => sousetup%sellm1z
    sellm2z => sousetup%sellm2z
    seppm2 => sousetup%seppm2
    sem1 => sousetup%sem1
    sem2 => sousetup%sem2
    seppm1z => sousetup%seppm1z
    sem1z => sousetup%sem1z
    sem2z => sousetup%sem2z
    sthm1z => sousetup%sthm1z
    seppm2z => sousetup%seppm2z
    soucoord => sousetup%soucoord
    knockon => sousetup%knockon
    komodel => sousetup%komodel
    nkorfn => sousetup%nkorfn
    nonko => sousetup%nonko
    noffko => sousetup%noffko
    soffvte => sousetup%soffvte
    soffpr => sousetup%soffpr
    isoucof => sousetup%isoucof
    faccof => sousetup%faccof
    jfl => sousetup%jfl
    xlfac => sousetup%xlfac
    xlpctlwr => sousetup%xlpctlwr
    xlpctmdl => sousetup%xlpctmdl
    xllwr => sousetup%xllwr
    xlmdl => sousetup%xlmdl
    szm2z => sousetup%szm2z
    sthm1 => sousetup%sthm1
    Szm1 => sousetup%Szm1
    szm2 => sousetup%szm2
    initialized_sou_pointers = .FALSE.
  end subroutine initialize_sou_pointers

  subroutine initialize_setup_pointers
    if(initialized_setup_pointers) call abort

    initialized_setup_pointers = .TRUE.
  end subroutine initialize_setup_pointers

  subroutine initialize_cqlcomm
    if(initialized_cqlcomm) call abort
    call initialize_eq_pointers
    call initialize_rf_pointers
    call initialize_tr_pointers
    call initialize_sou_pointers
    call initialize_setup_pointers
    initialized_cqlcomm = .TRUE.
  end subroutine initialize_cqlcomm

end module cqlcomm_mod
