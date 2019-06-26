module cqlconf_mod
  use param_mod, only : lrorsa, lfielda
  use param_mod, only : ep100, nmodsa, ngena, nrdca, nbctimea
  use param_mod, only : njenea
  use param_mod, only : nsoa, lrza
  use param_mod, only : ntotala, npaproca, lza, nmaxa, nplota, nsavea, nva, nena
  use param_mod, only : ndtr1a, negyrga, ntavga, noutpta
  use param_mod, only : pi
  use iso_c_binding, only : c_double, c_double_complex
  implicit none

  private
  integer :: ll
  logical, save :: nml_file_open = .FALSE.
  integer, save :: nml_fd = -1
  real(c_double), parameter :: drya=1.d0/DBLE(lrza)
  public nml_close

  public get_setup0_from_nml
  public print_setup0
  public set_setup0

  public get_eqsetup_from_nml
  public print_eqsetup
  public set_eqsetup

  public get_rfsetup_from_nml
  public print_rfsetup
  public set_rfsetup

  public get_trsetup_from_nml
  public print_trsetup
  public set_trsetup

  public get_sousetup_from_nml
  public print_sousetup
  public set_sousetup

  public get_setup_from_nml
  public print_setup
  public set_setup

  public print_all_conf_nml


  type, public ::  setup0_t
     character(len=256) :: mnemonic = "default_output"
     integer :: ioutput(2) = (/ 6, 0 /)
     character(len=8) :: iuser = "unset"
     character(len=8) :: ibox(3) = (/ "unset", "unset", "unset" /)
     character(len=8) :: noplots = "disabled"
     integer :: lnwidth = 3
     character(len=8) :: nmlstout = "trnscrib"
     character(len=8) :: special_calls = "enabled"
     !     Avoid some special calls, if setup0%special_calls=disabled in
     !     &setup0 namelist.
     !     [System calls not supported on some machines.]
     !     [Could use this for other system dependent branching, in future....]

     !     main parameters (used to allocate memory and determine model)
     character(len=8) :: cqlpmod = "disabled"
     integer :: lrz = 0
     character(len=8) :: lrzdiff = "disabled"
     integer :: lrzmax = 0
     integer :: lrindx(0:lrorsa) =  (/ (ll, ll=0,lrorsa) /)
     integer :: ls = 0
     integer :: lsmax = 0
     character(len=8) :: lsdiff = "disabled"
     integer :: lsindx(0:lrorsa) =  (/ (ll, ll=0,lrorsa) /)
     character(len=8) :: nlrestrt = "disabled"
     character(len=8) :: nlwritf = "disabled"
  end type setup0_t

  type, public :: eqsetup_t
     real(c_double) :: atol = 1.e-8
     real(c_double) :: ellptcty = 0.
     character(len=8) :: eqdskalt = "disabled"
     character(len=256) :: eqdskin = "eqdsk"
     character(len=8) :: eqmodel = "power"
     real(c_double) :: eqpower = 1.
     character(len=8) :: eqsym = "average"
     character(len=8) :: eqsource = "eqdsk"
     real(c_double) :: bsign = +1.0
     character(len=8) :: eqmod = "disabled"
     character(len=8) :: fpsimodl = "constant"
     integer :: lfield = lfielda
     integer :: methflag = 10
     character(len=8) :: nconteq = "disabled"
     integer :: nconteqn = 50
     real(c_double) :: povdelp = .2
     real(c_double) :: rtol = 1.e-8
     real(c_double) :: rmag = 100.
     real(c_double) :: rbox = 100.
     real(c_double) :: rboxdst = 50.
     real(c_double) :: zbox = 100.
  end type eqsetup_t

  type, public :: rfsetup_t
     !..................................................................
     !     Variables determine which, if any,  ray tracing code is called.
     !..................................................................
     character(len=8) :: call_lh = "disabled"
     character(len=8) :: call_ech = "disabled"
     character(len=8) :: call_fw = "disabled"
     !..................................................................
     !     ieqbrurf designates the source of equilibrium data to be used by
     !     xbram.  Appropriate values are: ieqbrurf=0 to use Brambilla
     !     ieqbrurf=0, to use Brambilla analytic "equilibria",
     !     =3, to use standard eqdsk input.
     !     =4, to use extended eqdsk, including density, and
     !     temperature profiles,...[output by GA ONETWO transport code].
     !     If eqsource="tsc", ieqbrurf is reset to = 4.
     !..................................................................
     integer :: ieqbrurf = 1
     !..................................................................
     !     Step interval at which to recalc diffusion coefficients,
     !     as a fraction of ncoef.
     !..................................................................
     integer :: urfncoef = 1.d0
     real(c_double)  :: dlndau(1:nmodsa) = 1.
     !..................................................................
     !     lh,ech and fw determine which wave modes are utilized
     !     Alternatively(060314), this data can be input through rftype()
     !..................................................................
     character(len=8) :: lh = "disabled"
     character(len=8) :: ech = "disabled"
     character(len=8) :: fw = "disabled"
     character(len=8) :: rftype(1:nmodsa) = "notset"
     character(len=256) :: rffile(1:nmodsa) = "notset"
     character(len=256) :: rdcfile(1:nrdca) = (/ (merge("du0u0_input", "notset     ", ll==1), ll=1,nrdca) /)
     character(len=8) :: rfread = "text"
     !..................................................................
     !     nharms.gt.0, then calculate damping for harmonics
     !     nharm1 to nharm1+(nharms-1), rather than according to nharm
     !     in the rayop file.
     !     This option is only viable for one rf mode, i.e., only one
     !     of lh, ech, or fw may be enabled.
     !     060214: Extended to multi-rf-wave types.
     !..................................................................
     integer :: nharms(1:nmodsa) = 0
     integer :: nharm1(1:nmodsa) = 0
     integer :: nrfspecies(1:nmodsa) = 1
     !..................................................................
     !     iurfcoll (iurfl) indicate use of collisional (additional linear)
     !     absorption coeff. Passed in ray data files.
     !..................................................................
     character(len=8) :: iurfcoll(1:nmodsa) = "disabled"
     character(len=8) :: iurfl(1:nmodsa) = "disabled"
     !     number of elements in Bessel table at setup time
     integer :: nbssltbl = 10000
     !     damping on at step nondamp
     integer :: nondamp = 0
     !..................................................................
     !     maximum number of ray elements per ray for first call to R.T code
     !     [Not presently implemented, BH060314]
     !..................................................................
     integer :: nrfstep1(1:nmodsa) = 100
     !..................................................................
     !     number of ray elements added at each step in addition of new data.
     !     [Not presently implemented, BH060314]
     !.................................................................
     integer :: nrfstep2 = 50
     !     number of steps in power ramp-up
     integer :: nrfpwr = 3
     !     number of iterations at full power, for first ray elements
     integer :: nrfitr1 = 1
     !     number of iterations after additional ray data
     integer :: nrfitr2 = 1
     !     number of steps adding new ray data
     integer :: nrfitr3 = 2
     integer :: nonrf(1:ngena) = 0
     integer :: noffrf(1:ngena) = 10000
     integer :: nrf = 0
     !..................................................................
     !     scaleurf="enabled" rescale contribution to B so that a particular
     !     ray does not "overdamp" on a given flux surface volume.
     !..................................................................
     character(len=8) :: scaleurf = "enabled"
     !..................................................................
     !     scaling of power and nparallel- width in rays
     !..................................................................
     real(c_double) :: pwrscale(1:nmodsa) = 1.d0
     real(c_double) :: wdscale(1:nmodsa) = 1.d0
     character(len=8) :: urfrstrt = "disabled"
     character(len=8) :: urfwrray = "disabled"
     !..................................................................
     !     Additional time-dependent scaling of power
     !..................................................................
     integer :: nurftime = 0
     real(c_double) :: urftime(nbctimea) = 0.0
     real(c_double) :: pwrscale1(1:nbctimea) = 1.d0
     !..................................................................
     !     urfdmp="secondd" means utilize the "second derivative" damping
     !     diagnostic in computing the damping due to each ray as
     !     it moves through the plasma. If urfdmp .ne. "secondd" the
     !     diagnostic uses an integration by parts technique to compute
     !     the damping. We highly recommend "secondd" because convergence
     !     is the agreement between the sum of the damping of all rays
     !     and the absorbed power as computed from dF/dt. This latter
     !     diagnostic utilizes the "second derivative" approach so
     !     consistency demands "secondd" for the rays.
     !..................................................................
     character(len=8) :: urfdmp = "secondd"
     real(c_double) :: urfmult = 1.0
     character(len=8) :: urfmod = "disabled"
     character(len=8) :: vlhmod = "disabled"
     real(c_double) :: vlhmodes = 1.
     real(c_double) :: vparmin(1:nmodsa) = 1.
     real(c_double) :: vparmax(1:nmodsa) = 1.
     character(len=8) :: vprprop = "disabled"
     real(c_double) :: vdalp = .03
     real(c_double) :: vlh_karney = 0.
     character(len=8) :: vlhprprp(1:nmodsa) = "parallel"
     character(len=8) :: vlhplse = "disabled"
     real(c_double) :: vlhpon = .1
     real(c_double) :: vlhpoff = .11
     real(c_double) :: vprpmin(1:nmodsa) = 0.
     real(c_double) :: vprpmax(1:nmodsa) = ep100
     real(c_double) :: vlhpolmn(1:nmodsa) = 0.
     real(c_double) :: vlhpolmx(1:nmodsa) = 100.
     character(len=8) :: vlfmod = "disabled"
     real(c_double) :: vlfmodes = 1.
     real(c_double) :: vlffreq(1:nmodsa) = .8e9
     real(c_double) :: vlfnp(1:nmodsa) = 5.
     real(c_double) :: vlfdnp(1:nmodsa) = .2
     real(c_double) :: vlfddnp(1:nmodsa) = .1
     complex(c_double_complex) :: vlfeplus(1:nmodsa) = (0.,0.)
     complex(c_double_complex) :: vlfemin(1:nmodsa) = (0.,0.)
     real(c_double) :: vlfpol(1:nmodsa) = 0.
     real(c_double) :: vlfdpol(1:nmodsa) = 360.
     real(c_double) :: vlfddpol(1:nmodsa) = 20.
     character(len=8) :: vlfbes = "enabled"
     character(len=8) :: vlfnpvar = "1/r"
     real(c_double) :: vlfharms(1:nmodsa) = 1.
     real(c_double) :: vlfharm1(1:nmodsa) = 0.
     real(c_double) :: vlfnperp(1:nmodsa) = 5.
     real(c_double) :: vlfdnorm(1:nmodsa) = 10.
     real(c_double) :: vlfparmn(1:nmodsa) = -ep100
     real(c_double) :: vlfparmx(1:nmodsa) = +ep100
     real(c_double) :: vlfprpmn(1:nmodsa) = 0.
     real(c_double) :: vlfprpmx(1:nmodsa) = ep100
     real(c_double) :: rdc_upar_sign = +1
     integer :: nrdc = 1
     character(len=8) :: rdcmod = "disabled"
     character(len=8) :: rdc_clipping = "disabled"
     integer :: nrdcspecies(1:nrdca) = (/ (merge(1, 0, ll==1), ll=1,nrdca) /)
     real(c_double) :: rdcscale(1:nrdca) = 1.d0
     ! this may have been undefined in original code
     character(len=8) :: rdc_netcdf = "disabled"
  end type rfsetup_t

  type, public :: trsetup_t
     real(c_double) :: advectr = 1.0
     character(len=8) :: difus_type(1:ngena) = "specify"
     real(c_double) :: difusr = 1.d4
     real(c_double) :: difus_rshape(8) = (/ 1.0, 3.0, 3.0, 1.0, -1.0, 0.0, 0.0, 0.0 /)
     real(c_double) :: difus_vshape(4) = 0.
     real(c_double) :: difin(1:njenea) = 0.
     character(len=8) :: difus_io(1:ngena) = "disabled"
     character(len=256) :: difus_io_file = "drrin.nc"
     real(c_double) :: difus_io_drrscale(1:nbctimea, 1:ngena) = 1.
     real(c_double) :: difus_io_drscale(1:nbctimea, 1:ngena) = 1.
     real(c_double) :: difus_io_t(1:nbctimea) = 0.
     character(len=8) :: pinch = "disabled"
     real(c_double) :: relaxden = 1.
     character(len=8) :: relaxtsp = "disabled"
     character(len=8) :: transp= " disabled"
     character(len=8) :: adimeth = "disabled"
     integer :: nonadi = 5
     integer :: nontran = 0
     integer :: nofftran = 10000
     integer :: nonelpr = 10000
     integer :: noffelpr = 0
     integer :: ndifus_io_t = 0
  end type trsetup_t

  type, public :: sousetup_t
     real(c_double) :: asorz(1:ngena,1:nsoa,0:lrza) = 0
     real(c_double) :: asor(1:ngena,1:nsoa,1:lrza) = 0
     character(len=8) :: flemodel = "fsa"
     integer :: nonso(1:ngena,1:nsoa) = 100000
     integer :: noffso(1:ngena,1:nsoa) = 100000
     integer :: nso = 0
     integer :: nsou = 1
     character(len=8) :: pltso = "enabled"
     real(c_double) :: mpwrsou(0:ngena) = 3.
     real(c_double) :: npwrsou(0:ngena) = 2.
     character(len=8) :: soucoord = "disabled"
     character(len=8) :: knockon = "disabled"
     character(len=8) :: komodel = "mr"
     integer :: nkorfn = 1
     integer :: nonko = 10000
     integer :: noffko = 10000
     real(c_double) :: soffvte = 3.
     real(c_double) :: soffpr = 0.5
     real(c_double) :: isoucof = 0 ! hrmmm, maybe should be int?
     real(c_double) :: faccof = 1.e0
     integer :: jfl = 151  !!!min(201,jx)
     real(c_double) :: xlfac = 1.
     real(c_double) :: xlpctlwr = .1
     real(c_double) :: xlpctmdl = .4
     real(c_double) :: xllwr = 1./43.
     real(c_double) :: xlmdl = .25
     ! the following will be initilized to non trivial values in the setter
     ! we should probably init these to nans, but fortran makes that a chore
     real(c_double) :: scm2z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: scm2(1:ngena,1:nsoa) = 0.
     real(c_double) :: sellm1(1:ngena,1:nsoa) = 1.
     real(c_double) :: sellm2(1:ngena,1:nsoa) = 1.
     real(c_double) :: seppm1(1:ngena,1:nsoa) = 0.
     real(c_double) :: sellm1z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: sellm2z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: seppm2(1:ngena,1:nsoa) = 1.
     real(c_double) :: sem1(ngena,nsoa) = 0.
     real(c_double) :: sem2(ngena,nsoa) = 0.
     real(c_double) :: seppm1z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: sem1z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: sem2z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: sthm1z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: seppm2z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: szm1z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: szm2z(ngena,nsoa,0:lrza) = 0.
     real(c_double) :: sthm1(ngena,nsoa) = 0.
     real(c_double) :: szm1(ngena,nsoa) = 0.
     real(c_double) :: szm2(ngena,nsoa) = 1.
  end type sousetup_t

  type, public ::  setup_t
     ! forgive me, i'm just a code custodian
     !..................................................................
     !     acoef's specify ASDEX exponentail profiles
     !     (ti profiles given by te profile, for "asdex" option).
     !..................................................................
     real(c_double) :: acoefne(4) = (/ -1.87, -0.57, -24.78, -181.38 /)
     real(c_double) :: acoefte(4) = (/ 7.51, -13.45, 6.21, -125.64 /)
     character(len=8) :: ampfmod = "disabled"
     integer :: nampfmax = 2
     integer :: nonampf = 0
     real(c_double) :: ampferr = 1.d-3
     real(c_double) :: bctimescal = 1.d0
     real(c_double) :: bnumb(ntotala) = 1.
     real(c_double) :: bth = 1000.
     real(c_double) :: btor = 10000.
     character(len=8) :: bootst = "disabled"   ! analytic (Hinton and Haseltine)
     ! bootstrap current
     character(len=8) :: bootcalc = "disabled" !computational bootstrap current off.
     character(len=8) :: bootupdt = "disabled" !updating 0th order distn for bs radial derv.
     real(c_double) :: bootsign = +1.0
     integer :: nonboot = 2           !turn on computational bootstrap at n=nonboot.
     integer :: jhirsh = 0
     real(c_double) :: contrmin = 1.d-12
     real(c_double) :: constr = 1.d-3
     character(len=8) :: chang = "enabled"
     integer :: colmodl = 1
     real(c_double) :: deltabdb = 0.
     real(c_double) :: denpar(ntotala,0:lza+1) =  1.
     real(c_double) :: droptol = 0.001d0
     real(c_double) :: dtr = 5.d0
     real(c_double) :: dtr1(ndtr1a) = 0.
     real(c_double) :: eegy(negyrga,2,ngena,lrza) = 0.
     real(c_double) :: eparc(ngena,0:lrza) = -1.
     real(c_double) :: eperc(ngena,0:lrza) = -1.
     real(c_double) :: simpbfac = 1.d0
     !unused? epar
     !unused? eper
     real(c_double) :: elecfld(0:lrza) = 0.0
     real(c_double) :: elpar0 = 0.
     real(c_double) :: enorm = 200.
     real(c_double) :: enorme = 200. ! (re)set to enom at runtime, sort of dangerous
     real(c_double) :: enormi = 200. ! (re)set to enom at runtime, sort of dangerous
     character(len=8) :: eleccomp = "enabled"
     real(c_double) :: elecin(njenea) = 0.0d0
     real(c_double) :: elecin_t(njenea,nbctimea) = 0.0d0
     real(c_double) :: elecscal = 1.
     real(c_double) :: enein(njenea,ntotala) = 0.0d0
     real(c_double) :: enein_t(njenea,ntotala,nbctimea) = 0.0d0
     ! xxx this var was undefined in original code, setting 0
     real(c_double) :: ennin_t(njenea,nbctimea,npaproca) = 0.0d0 !neutrals,impurities,etc.
     real(c_double) :: ennin(njenea,npaproca) = 0.0d0
     real(c_double) :: enescal = 1.
     real(c_double) :: enloss(ngena) = 200.
     real(c_double) :: epsthet = 0.1
     real(c_double) :: enmin = 5.
     real(c_double) :: enmax = 50.
     real(c_double) :: ennb(npaproca) = 1.e10
     real(c_double) :: ennl(npaproca) = 5.0
     real(c_double) :: ennscal(npaproca) = 1.0
     real(c_double) :: enmin_npa = 5.
     real(c_double) :: enmax_npa = 50.
     character(len=8) :: eseswtch = "disabled"
     real(c_double) :: xsink = 0.
     real(c_double) :: esink = 0.
     real(c_double) :: ephicc = 0.
     !unused real(c_double) :: esfac
     real(c_double) :: eoved = -.01
     real(c_double) :: fds = 0.2
     real(c_double) :: fds_npa = 0.2
     real(c_double) :: fmass(ntotala) = 1.e-29
     character(len=8) :: f4d_out = "disabled"
     character(len=8) :: tavg = "disabled"
     real(c_double) :: tavg1(ntavga) = 0.
     real(c_double) :: tavg2(ntavga) = 0.
     real(c_double) :: gsla = 270.
     real(c_double) :: gslb = 35.
     real(c_double) :: gamaset = 0.
     character(len=8) :: gamafac = "disabled"
     real(c_double) :: gamegy(ngena) = 0.
     character(len=8) :: iactst = "disabled"
     character(len=8) :: ineg = "disabled"
     character(len=8) :: idskf = "disabled"
     character(len=8) :: idskrf = "disabled"
     character(len=8) :: ichkpnt = "disabled"
     character(len=8) :: implct = "enabled"
     !..................................................................
     !     Profile options are "parabola", "splines", and "asdex":
     !..................................................................
     character(len=8) :: iprone = "parabola"
     character(len=8) :: iprote = "parabola"
     character(len=8) :: iproti = "parabola"
     character(len=8) :: iprozeff = "disabled"
     character(len=8) :: iprovphi = "disabled"
     character(len=8) :: iproelec = "parabola"
     character(len=8) :: ipronn = "disabled"
     character(len=8) :: iprocur = "parabola"
     character(len=8) :: tmdmeth = "method1"
     integer :: isigmas(6) = 0
     integer :: isigtst = 1
     integer :: isigsgv1 = 0
     integer :: isigsgv2 = 0
     real(c_double) :: pltflux1(7) = (/ 1., 1., 1., 1., 1., 1., 0. /)
     integer :: irzplt(lrorsa) = 0
     character(len=8) :: izeff = "backgrnd"
     ! unused ioutime, yes indeed
     integer :: iy = 200 ! default value; will be over-written by cqlinput value
     integer :: jx = 300 ! default value; will be over-written by cqlinput value
     integer :: kenorm = 1
     integer :: lfil = 30
     integer :: kfrsou = 0
     character(len=8) :: kpress(ntotala) = "enabled"
     character(len=8) :: kfield(ntotala) = "enabled"
     character(len=8) :: kspeci(2,ntotala) = " "
     ! this one is a pain to init
     real(c_double) :: fpld(10,ngena) = reshape( &
          (/ ( merge(pi, merge(1.d10, 0.d0 , mod(ll,10)==8), &
          mod(ll,10)==0), ll=1,10*ngena) /), &
          shape=(/10 ,ngena/))
     integer :: lmidpln = 1
     character(len=8) :: locquas = "disabled"
     character(len=8) :: lbdry(ngena) = "conserv"
     character(len=8) :: lbdry0 = "enabled"
     character(len=256) :: lossfile(ngena) = "./prompt_loss.txt"
     character(len=8) :: lossmode(ngena) = "disabled"
     integer :: lmidvel = 0
     integer :: laddbnd = 1
     integer :: lz = lza
     character(len=8) :: machine = "toroidal"
     character(len=8) :: meshy = "free"
     character(len=8) :: manymat = "disabled"
     character(len=256) :: netcdfnm = "disabled"
     character(len=8) :: netcdfshort = "disabled"
     character(len=8) :: netcdfvecal = "disabled"
     character(len=8) :: netcdfvecc = "disabled"
     character(len=8) :: netcdfvece = "disabled"
     character(len=8) :: netcdfvecrf = "disabled"
     character(len=8) :: netcdfvecs = "all"
     integer :: nnspec = 1
     real(c_double) :: mpwr(0:ntotala) = (/ ( merge(1, 3, ll<4), ll=0,ntotala) /)
     real(c_double) :: megy(ngena) = 0.0
     real(c_double) :: mtorloss(ngena) = 0.0
     integer :: mmsv = 3 ! xxx reset by mx at run time, seems dangerous
     integer :: msxr = 3 ! xxx reset by mx at run time, seems dangerous
     integer :: mx = 3   ! default value; will be over-written by cqlinput value
     !     if nchgdy = 1 adapt dy(ith)
     integer :: nchgdy = 0
     integer :: ngauss = 0
     integer :: nlagran = 4
     logical :: nlotp1(noutpta) = (/ (merge(.false., .true., ll==4), ll=1,noutpta) /) 
     logical :: nlotp2(noutpta) = .false.
     logical :: nlotp3(noutpta) = .false.
     logical :: nlotp4(noutpta) = .false.
     integer :: nmax = nmaxa
     integer :: ngen = ngena
     integer :: nkconro(ntotala) = (/ (merge(ll, 0, ll<3), ll=1,ntotala) /)
     integer :: nplt3d(nplota) = -10000
     integer :: nrskip = 10
     integer :: nen = nena
     integer :: nv = 1
     integer :: nv_npa = 1
     integer :: nen_npa = nena
     integer :: npaproc = 1
     character(len=8) :: npa_process(npaproca) = (/ (merge("cxh     ", "notset  ", ll==1), ll=1,npaproca) /)
     integer :: nr_delta = 65
     integer :: nz_delta = 65
     integer :: nt_delta = 80 !Needs to be even
     ! For saving f4d== f(R,Z,u,theta) distribution:
     integer :: nr_f4d = 20
     integer :: nz_f4d = 21
     integer :: nv_f4d = 20
     integer :: nt_f4d = 20
     real(c_double) :: npwr(0:ntotala) = 2.
     real(c_double) :: negy(ngena) = 0.0
     real(c_double) :: ntorloss(ngena) = 0.0
     integer :: njene = 0
     integer :: njte = 0 ! reset to njene at runtime, seems dangerous
     integer :: njti = 0 ! reset to njene at runtime, seems dangerous
     integer :: nstop = 5
     integer :: nondtr1(ndtr1a) = -1
     integer :: nplot(nplota) = -10000
     integer :: nsave(nsavea) = -10000
     integer :: ncoef = 1
     integer :: nchec = 1
     integer :: ncont = 25
     integer :: nrstrt = 1
     integer :: nstps = 100
     !     if ngauss.ge.1 => analegco=disabled
     !     good numbers are nlagran=4 and ngauss=4 or 6
     !     max. nlagran allowed: 15
     integer :: nfpld = 0
     integer :: noncntrl = 0
     integer :: nonel = 0
     integer :: noffel = 10000
     integer :: nonvphi = 10000
     integer :: noffvphi = 10000
     integer :: nonavgf = 5
     integer :: nofavgf = 10
     integer :: nonloss = 0
     integer :: noffloss = 10000
     integer :: nummods = 1
     !     if numixts= 1=>forw/back; -1=>back/forw for numindx=2
     integer :: numixts = 1
     integer :: numby = 20
     integer :: negyrg = 0
     !     old way of integrating dens,cur in diaggnde
     character(len=8) :: oldiag = "enabled"
     character(len=8) :: plt3d = "enabled"
     character(len=8) :: pltfvs = "disabled"
     character(len=8) :: partner = "disabled"
     real(c_double) :: paregy(ngena) = 0.0
     real(c_double) :: peregy(ngena) = 0.0
     real(c_double) :: pegy(ngena) = 0.0
     real(c_double) :: zeffin(0:njenea) = 1.
     real(c_double) :: zeffin_t(njenea,nbctimea) = 0.0d0
     real(c_double) :: zeffscal = 1.
     real(c_double) :: vphiplin(0:njenea) = 0.0d0
     real(c_double) :: vphiplin_t(njenea,nbctimea) = 0.0d0
     real(c_double) :: vphiscal = 1.
     character(len=8) :: pltdn = "disabled"
     character(len=8) :: pltvecal = "disabled"
     character(len=8) :: pltvecc = "disabled"
     character(len=8) :: pltvecrf = "disabled"
     character(len=8) :: pltvece = "disabled"
     character(len=8) :: pltstrm = "disabled"
     character(len=8) :: pltflux = "disabled"
     real(c_double) :: pltmag = 1.
     character(len=8) :: pltsig = "enabled"
     character(len=8) :: pltlim = "disabled"
     ! pltlimm use seems sparse, not in comm...
     real(c_double) :: pltlimm = 0.
     character(len=8) :: pltrst = "enabled"
     character(len=8) :: plturfb = "enabled" !YuP[2018-02-07] New: 'color' for color contour plots
     character(len=8) :: pltvflu = "disabled"
     character(len=8) :: pltra = "disabled"
     character(len=8) :: pltvs = "rho"
     character(len=8) :: pltd = "enabled"
     character(len=8) :: pltprpp = "disabled"
     character(len=8) :: pltfofv = "disabled"
     character(len=8) :: pltlos = "disabled"
     character(len=8) :: profpsi = "disabled"
     character(len=8) :: psimodel = "axitorus"
     character(len=8) :: pltpowe = "disabled"
     character(len=8) :: pltend = "enabled"
     character(len=8) :: pltinput = "enabled"
     !XXXX pltrdc was found uninitialized... keeping that way XXX
     character(len=8) :: pltrdc
     !pltview was unused
     character(len=8) :: qsineut = "disabled"
     character(len=8) :: trapmod = "disabled"
     real(c_double) :: trapredc = 0.
     character(len=8) :: scatmod = "disabled"
     real(c_double) :: scatfrac = 1.
     real(c_double) :: ryain(njenea) = 0.0d0
     real(c_double) :: radmaj = 100.
     real(c_double) :: radmin = 50.
     real(c_double) :: rmirror = 2.
     character(len=8) :: relativ = "enabled"
     real(c_double) :: reden(ntotala,0:lrza) = 1.
     character(len=8) :: regy(ngena) = "disabled"
     real(c_double) :: rfacz = .7
     character(len=8) :: rzset = "disabled"
     real(c_double) :: rd(nva) = (/ (merge(100.d0, 0.d0, ll==1), ll=1,nva) /)
     real(c_double) :: roveram = 1.e-6
     real(c_double) :: rovera(lrza) = .1
     real(c_double) :: rya(0:lrza+1) = (/ ( min(1., (ll - min(dble(ll),0.5))*drya), ll= 0, lrza+1) /)
     character(len=8) :: radcoord = "sqtorflx"
     character(len=8) :: sbdry = "bounded"
     character(len=8) :: scheck = "enabled"
     character(len=8) :: ndeltarho = "disabled"
     character(len=8) :: softxry = "disabled"
     character(len=8) :: npa_diag = "disabled"
     character(len=8) :: symtrap = "enabled"
     character(len=8) :: syncrad = "disabled"
     character(len=8) :: bremsrad = "disabled"
     real(c_double) :: brfac = 0.
     real(c_double) :: brfac1 = 0.
     real(c_double) :: brfacgm3 = 1.0
     character(len=8) :: sigmamod = "disabled"
     real(c_double) :: sigvcx = 0.
     real(c_double) :: sigvi = 0.
     character(len=8) :: soln_method = "direct"
     real(c_double) :: tauegy(ngena,0:lrza) = -1.
     character(len=8) :: taunew = "disabled"
     real(c_double) :: tescal = 1.
     real(c_double) :: tein(njenea) = 0.
     real(c_double) :: tein_t(njenea,nbctimea) = 0.
     real(c_double) :: tiscal = 1.
     real(c_double) :: tiin_t(njenea,nbctimea) = 0.
     real(c_double) :: tiin(njenea) = 0.
     real(c_double) :: tauloss(3,ngena) = 0.
     real(c_double) :: temp(ntotala,0:lrza) = 1.
     real(c_double) :: temppar(ntotala,0:lza+1) = 1.
     real(c_double) :: tfac = 1.
     real(c_double) :: tfacz = 1.
     real(c_double) :: tbnd(lrorsa) = (/ (merge( 0.002, 0., ll==1), ll=1,lrorsa ) /)
     character(len=8) :: tandem = "disabled"
     real(c_double) :: thetd(nva) = 0.0
     character(len=8) :: torloss(ngena) = "disabled"
     real(c_double) :: thet1(nva) = 90.
     real(c_double) :: thet2(nva) = 180.
     real(c_double) :: thet1_npa(nva) = 90.
     real(c_double) :: thet2_npa(nva) = 180.
     real(c_double) :: x_sxr(nva) = 0.
     real(c_double) :: z_sxr(nva) = 0.
     real(c_double) :: rd_npa(nva) = 100.d0
     real(c_double) :: thetd_npa(nva) = 0.
     real(c_double) :: x_npa(nva) = 0.
     real(c_double) :: z_npa(nva) = 0.
     character(len=8) :: atten_npa = "enabled"
     character(len=8) :: updown = "symmetry"
     real(c_double) :: veclnth = 1.0
     real(c_double) :: vnorm = 4.e10   !  Usually set through enorm
     real(c_double) :: xfac = 1.
     real(c_double) :: xpctlwr = .1
     real(c_double) :: xpctmdl = .4
     real(c_double) :: xlwr = 1./43.
     real(c_double) :: xmdl = .25
     real(c_double) :: xsinkm = 1.
     real(c_double) :: xprpmax = 1.
     ! original code, these were defined in aindflt1 at runtime, dangerous, why nml
     integer :: ipxy,jpxy
     character(len=8) :: yreset = "disabled"
     real(c_double) :: ylower = 1.22
     real(c_double) :: yupper = 1.275
     real(c_double) :: mpwrzeff = 1.
     real(c_double) :: npwrzeff = 1.
     real(c_double) :: mpwrvphi = 0.
     real(c_double) :: npwrvphi = 2.
     real(c_double) :: mpwrxj = 1.
     real(c_double) :: npwrxj = 1.
     real(c_double) :: npwrelec = 1.
     real(c_double) :: mpwrelec = 1.
     real(c_double) :: redenc(nbctimea,ntotala) = 0.d0
     real(c_double) :: redenb(nbctimea,ntotala) = 0.d0
     real(c_double) :: temp_den = 0.0
     real(c_double) :: tempc(nbctimea,ntotala) = 0.d0
     real(c_double) :: tempb(nbctimea,ntotala) = 0.d0
     real(c_double) :: zeffc(nbctimea) = 0.d0
     real(c_double) :: zeffb(nbctimea) = 0.d0
     real(c_double) :: elecc(nbctimea) = 0.d0
     real(c_double) :: elecb(nbctimea) = 0.d0
     real(c_double) :: vphic(nbctimea) = 0.d0
     real(c_double) :: vphib(nbctimea) = 0.d0
     real(c_double) :: xjc(nbctimea) = 0.d0
     real(c_double) :: xjb(nbctimea) = 0.d0
     real(c_double) :: xjin_t(njenea,nbctimea) = 0.0d0
     real(c_double) :: totcrt(nbctimea) = 0.d0
     character(len=8) :: efswtch = "method1"
     character(len=8) :: efswtchn = "disabled"
     character(len=8) :: efiter = "enabled"
     character(len=8) :: efflag = 'toroidal'
     real(c_double) :: curr_edge = 0.
     real(c_double) :: efrelax = 0.5
     real(c_double) :: efrelax1 = 0.5 !0.8
     real(c_double) :: currerr = 0.1 !0.1
     real(c_double) :: bctime(nbctimea) = (/ (DBLE(ll-1), ll = 1, nbctimea) /)
     integer :: nbctime = 0
     real(c_double) :: zmax(0:lrza) = 1000.
     character(len=8) :: fow = "disabled" ! "disabled" is to use ZOW model as the main model in CQL3D
     ! from here, remaining vars in compx version of nml are unused in this code
  end type setup_t


  ! rest of cql3d will access this
  type(setup0_t), public, save :: setup0
  type(eqsetup_t), public, target, save :: eqsetup
  type(rfsetup_t), public, target, save :: rfsetup
  type(trsetup_t), public, target, save :: trsetup
  type(sousetup_t), public, target, save :: sousetup
  type(setup_t), public, target, save :: setup
  type(setup_t) :: setup_default

contains

  subroutine maybe_nml_open(nml_file)
    character(len=*), intent(in) :: nml_file
    if (.not. nml_file_open) then
       print *, "Opening nml_file: ", nml_file
       nml_fd = newunit()
       open(unit=nml_fd, file=nml_file, status='old')
       nml_file_open = .TRUE.
    end if
  end subroutine maybe_nml_open

  subroutine nml_close()
    if (nml_file_open) then
       close(nml_fd)
       nml_file_open = .FALSE.
    end if
  end subroutine nml_close

  subroutine get_setup0_from_nml(nml_file, close_nml_file, debug_print)
    implicit none
    character(len=*), intent(in) :: nml_file
    logical, intent(in), optional :: close_nml_file
    logical, intent(in), optional :: debug_print
    ! local
    type(setup0_t) :: setup0_

    ! make private local variables to read in the namelist

    character(len=256) :: mnemonic
    integer :: ioutput(2)
    character(len=8) :: iuser
    character(len=8) :: ibox(3)
    character(len=8) :: noplots
    integer :: lnwidth
    character(len=8) :: nmlstout
    character(len=8) :: special_calls
    character(len=8) :: cqlpmod
    integer :: lrz
    character(len=8) :: lrzdiff
    integer :: lrzmax
    integer :: lrindx(0:lrorsa)
    integer :: ls
    integer :: lsmax
    character(len=8) :: lsdiff
    integer :: lsindx(0:lrorsa)
    character(len=8) :: nlrestrt
    character(len=8) :: nlwritf

    ! state the namelist, with associated vars

    namelist/setup0/ mnemonic,ioutput,iuser,ibox,noplots,lnwidth, &
         nmlstout,special_calls,cqlpmod,lrz,lrzdiff,lrzmax,lrindx, &
         ls,lsmax,lsdiff,lsindx,nlrestrt,nlwritf

    ! copy defaults to local vars

    mnemonic = setup0_%mnemonic
    ioutput = setup0_%ioutput
    iuser = setup0_%iuser
    ibox = setup0_%ibox
    noplots = setup0_%noplots
    lnwidth = setup0_%lnwidth
    nmlstout = setup0_%nmlstout
    special_calls = setup0_%special_calls
    cqlpmod = setup0_%cqlpmod
    lrz = setup0_%lrz
    lrzdiff = setup0_%lrzdiff
    lrzmax = setup0_%lrzmax
    lrindx = setup0_%lrindx
    ls = setup0_%ls
    lsmax = setup0_%lsmax
    lsdiff = setup0_%lsdiff
    lsindx = setup0_%lsindx
    nlrestrt = setup0_%nlrestrt
    nlwritf = setup0_%nlwritf

    ! read the nml, which will write into the local vars

    call maybe_nml_open(nml_file)
    read(nml_fd, setup0)

    ! external codes can call this, which packs the setup0 derived type.
    call set_setup0(mnemonic,ioutput,iuser,ibox,noplots,lnwidth, &
         nmlstout,special_calls,cqlpmod,lrz,lrzdiff,lrzmax,lrindx, &
         ls,lsmax,lsdiff,lsindx,nlrestrt,nlwritf, debug_print)

    ! we optionally close the nml file.
    if (present(close_nml_file)) then
       if(close_nml_file) then
          call nml_close()
       end if
    endif

  end subroutine get_setup0_from_nml

  subroutine set_setup0(mnemonic,ioutput,iuser,ibox,noplots,lnwidth, &
       nmlstout,special_calls,cqlpmod,lrz,lrzdiff,lrzmax,lrindx, &
       ls,lsmax,lsdiff,lsindx,nlrestrt,nlwritf, &
       debug_print)
    logical, intent(in), optional :: debug_print
    !
    character(len=256), intent(in), optional :: mnemonic
    integer, intent(in), optional :: ioutput(2)
    character(len=8), intent(in), optional :: iuser
    character(len=8), intent(in), optional :: ibox(3)
    character(len=8), intent(in), optional :: noplots
    integer, intent(in), optional :: lnwidth
    character(len=8), intent(in), optional :: nmlstout
    character(len=8), intent(in), optional :: special_calls
    character(len=8), intent(in), optional :: cqlpmod
    integer, intent(in), optional :: lrz
    character(len=8), intent(in), optional :: lrzdiff
    integer, intent(in), optional :: lrzmax
    integer, intent(in), optional :: lrindx(0:lrorsa)
    integer, intent(in), optional :: ls
    integer, intent(in), optional :: lsmax
    character(len=8), intent(in), optional :: lsdiff
    integer, intent(in), optional :: lsindx(0:lrorsa)
    character(len=8), intent(in), optional :: nlrestrt
    character(len=8), intent(in), optional :: nlwritf

    ! All this code should do is override the defaults
    ! in setup0 with optional args.

    if (present(mnemonic)) then
       setup0%mnemonic = mnemonic
    end if
    if (present(ioutput)) then
       setup0%ioutput = ioutput
    end if
    if (present(iuser)) then
       setup0%iuser = iuser
    end if
    if (present(ibox)) then
       setup0%ibox = ibox
    end if
    if (present(noplots)) then
       setup0%noplots = noplots
    end if
    if (present(lnwidth)) then
       setup0%lnwidth = lnwidth
    end if
    if (present(nmlstout)) then
       setup0%nmlstout = nmlstout
    end if
    if (present(special_calls)) then
       setup0%special_calls = special_calls
    end if
    if (present(cqlpmod)) then
       setup0%cqlpmod = cqlpmod
    end if
    if (present(lrz)) then
       setup0%lrz = lrz
    else
       stop 'setup0%lrz is required'
    end if
    if (present(lrzdiff)) then
       setup0%lrzdiff = lrzdiff
    end if
    if (present(lrzmax)) then
       setup0%lrzmax = lrzmax
    end if
    if (present(lrindx)) then
       setup0%lrindx = lrindx
    end if
    if (present(ls)) then
       setup0%ls = ls
    end if
    if (present(lsmax)) then
       setup0%lsmax = lsmax
    end if
    if (present(lsdiff)) then
       setup0%lsdiff = lsdiff
    end if
    if (present(lsindx)) then
       setup0%lsindx = lsindx
    end if
    if (present(nlrestrt)) then
       setup0%nlrestrt = nlrestrt
    end if
    if (present(nlwritf)) then
       setup0%nlwritf = nlwritf
    end if

    if ( present(debug_print)) then
       if (debug_print) call print_setup0
    end if

  end subroutine set_setup0

  subroutine print_setup0()
    namelist /setup0_nml/ setup0
    WRITE(*, *) "!----  BEGIN SETUP0 DUMP"
    WRITE(*, nml = setup0_nml)
    WRITE(*, *)  "!----  END SETUP0 DUMP"
  end subroutine print_setup0


  subroutine get_eqsetup_from_nml(nml_file, close_nml_file, debug_print)
    implicit none
    character(len=*), intent(in) :: nml_file
    logical, intent(in), optional :: close_nml_file
    logical, intent(in), optional :: debug_print
    ! local
    type(eqsetup_t) :: eqsetup_

    ! make local vars and assign defaults from type
    real(c_double) :: atol
    real(c_double) :: ellptcty
    character(len=8) :: eqdskalt
    character(len=256) :: eqdskin
    character(len=8) :: eqmodel
    real(c_double) :: eqpower
    character(len=8) :: eqsym
    character(len=8) :: eqsource
    real(c_double) :: bsign
    character(len=8) :: eqmod
    character(len=8) :: fpsimodl
    integer :: lfield
    integer :: methflag
    character(len=8) :: nconteq
    integer :: nconteqn
    real(c_double) :: povdelp
    real(c_double) :: rtol
    real(c_double) :: rmag
    real(c_double) :: rbox
    real(c_double) :: rboxdst
    real(c_double) :: zbox

    ! state the namelist, with associated vars

    namelist/eqsetup/ &
         atol, &
         ellptcty,eqmodel,eqpower,eqsource,eqdskin,bsign,eqmod, &
         eqsym,eqdskalt, &
         fpsimodl, &
         lfield, &
         methflag, &
         nconteq,nconteqn, &
         povdelp, &
         rtol,rmag,rbox,rboxdst, &
         zbox

    ! copy defaults to local vars

    atol = eqsetup_%atol
    ellptcty = eqsetup_%ellptcty
    eqdskalt = eqsetup_%eqdskalt
    eqdskin = eqsetup_%eqdskin
    eqmodel = eqsetup_%eqmodel
    eqpower = eqsetup_%eqpower
    eqsym = eqsetup_%eqsym
    eqsource = eqsetup_%eqsource
    bsign = eqsetup_%bsign
    eqmod = eqsetup_%eqmod
    fpsimodl = eqsetup_%fpsimodl
    lfield = eqsetup_%lfield
    methflag = eqsetup_%methflag
    nconteq = eqsetup_%nconteq
    nconteqn = eqsetup_%nconteqn
    povdelp = eqsetup_%povdelp
    rtol = eqsetup_%rtol
    rmag = eqsetup_%rmag
    rbox = eqsetup_%rbox
    rboxdst = eqsetup_%rboxdst
    zbox = eqsetup_%zbox

    call maybe_nml_open(nml_file)
    read(nml_fd, eqsetup)

    call set_eqsetup(atol, ellptcty, eqdskalt, eqdskin, eqmodel, eqpower, &
         eqsym, eqsource, bsign, eqmod, fpsimodl, lfield, methflag, nconteq, &
         nconteqn, povdelp, rtol, rmag, rbox, rboxdst, zbox, debug_print)

    if (present(close_nml_file)) then
       if(close_nml_file) then
          call nml_close()
       end if
    endif

  end subroutine get_eqsetup_from_nml

  subroutine set_eqsetup(atol, ellptcty, eqdskalt, eqdskin, eqmodel, eqpower, &
       eqsym, eqsource, bsign, eqmod, fpsimodl, lfield, methflag, nconteq, &
       nconteqn, povdelp, rtol, rmag, rbox, rboxdst, zbox, debug_print)
    implicit none
    logical, intent(in), optional :: debug_print
    !
    real(c_double), optional, intent(in) :: atol
    real(c_double), optional, intent(in) :: ellptcty
    character(len=8), optional, intent(in) :: eqdskalt
    character(len=256), optional, intent(in) :: eqdskin
    character(len=8), optional, intent(in) :: eqmodel
    real(c_double), optional, intent(in) :: eqpower
    character(len=8), optional, intent(in) :: eqsym
    character(len=8), optional, intent(in) :: eqsource
    real(c_double), optional, intent(in) :: bsign
    character(len=8), optional, intent(in) :: eqmod
    character(len=8), optional, intent(in) :: fpsimodl
    integer, optional, intent(in) :: lfield
    integer, optional, intent(in) :: methflag
    character(len=8), optional, intent(in) :: nconteq
    integer, optional, intent(in) :: nconteqn
    real(c_double), optional, intent(in) :: povdelp
    real(c_double), optional, intent(in) :: rtol
    real(c_double), optional, intent(in) :: rmag
    real(c_double), optional, intent(in) :: rbox
    real(c_double), optional, intent(in) :: rboxdst
    real(c_double), optional, intent(in) :: zbox

    if(present(atol)) eqsetup%atol = atol
    if(present(ellptcty)) eqsetup%ellptcty = ellptcty
    if(present(eqdskalt)) eqsetup%eqdskalt = eqdskalt
    if(present(eqdskin)) eqsetup%eqdskin = eqdskin
    if(present(eqmodel)) eqsetup%eqmodel = eqmodel
    if(present(eqpower)) eqsetup%eqpower = eqpower
    if(present(eqsym)) eqsetup%eqsym = eqsym
    if(present(eqsource)) eqsetup%eqsource = eqsource
    if(present(bsign)) eqsetup%bsign = bsign
    if(present(eqmod)) eqsetup%eqmod = eqmod
    if(present(fpsimodl)) eqsetup%fpsimodl = fpsimodl
    if(present(lfield)) eqsetup%lfield = lfield
    if(present(methflag)) eqsetup%methflag = methflag
    if(present(nconteq)) eqsetup%nconteq = nconteq
    if(present(nconteqn)) eqsetup%nconteqn = nconteqn
    if(present(povdelp)) eqsetup%povdelp = povdelp
    if(present(rtol)) eqsetup%rtol = rtol
    if(present(rmag)) eqsetup%rmag = rmag
    if(present(rbox)) eqsetup%rbox = rbox
    if(present(rboxdst)) eqsetup%rboxdst = rboxdst
    if(present(zbox)) eqsetup%zbox = zbox

    if ( present(debug_print)) then
       if (debug_print) call print_eqsetup
    end if

  end subroutine set_eqsetup


  subroutine print_eqsetup()
    namelist /eqsetup_nml/ eqsetup
    WRITE(*, *) "!----  BEGIN EQSETUP DUMP"
    WRITE(*, nml = eqsetup_nml)
    WRITE(*, *)  "!----  END EQSETUP DUMP"
  end subroutine print_eqsetup


  subroutine get_rfsetup_from_nml(nml_file, close_nml_file, debug_print)
    implicit none
    character(len=*), intent(in) :: nml_file
    logical, intent(in), optional :: close_nml_file
    logical, intent(in), optional :: debug_print
    ! local
    type(rfsetup_t) :: rfsetup_
    character(len=8) :: call_lh
    character(len=8) :: call_ech
    character(len=8) :: call_fw
    integer :: ieqbrurf
    integer :: urfncoef
    real(c_double)  :: dlndau(nmodsa)
    character(len=8) :: lh
    character(len=8) :: ech
    character(len=8) :: fw
    character(len=8) :: rftype(nmodsa)
    character(len=256) :: rffile(nmodsa)
    character(len=256) :: rdcfile(nrdca)
    character(len=8) :: rfread
    integer :: nharms(nmodsa)
    integer :: nharm1(nmodsa)
    integer :: nrfspecies(nmodsa)
    character(len=8) :: iurfcoll(nmodsa)
    character(len=8) :: iurfl(nmodsa)
    integer :: nbssltbl
    integer :: nondamp
    integer :: nrfstep1(nmodsa)
    integer :: nrfstep2
    integer :: nrfpwr
    integer :: nrfitr1
    integer :: nrfitr2
    integer :: nrfitr3
    integer :: nonrf(ngena)
    integer :: noffrf(ngena)
    integer :: nrf
    character(len=8) :: scaleurf
    real(c_double) :: pwrscale(nmodsa)
    real(c_double) :: wdscale(nmodsa)
    character(len=8) :: urfrstrt
    character(len=8) :: urfwrray
    integer :: nurftime
    real(c_double) :: urftime(nbctimea)
    real(c_double) :: pwrscale1(nbctimea)
    character(len=8) :: urfdmp
    real(c_double) :: urfmult
    character(len=8) :: urfmod
    character(len=8) :: vlhmod
    real(c_double) :: vlhmodes
    real(c_double) :: vparmin(nmodsa)
    real(c_double) :: vparmax(nmodsa)
    character(len=8) :: vprprop
    real(c_double) :: vdalp
    real(c_double) :: vlh_karney
    character(len=8) :: vlhprprp(nmodsa)
    character(len=8) :: vlhplse
    real(c_double) :: vlhpon
    real(c_double) :: vlhpoff
    real(c_double) :: vprpmin(nmodsa)
    real(c_double) :: vprpmax(nmodsa)
    real(c_double) :: vlhpolmn(nmodsa)
    real(c_double) :: vlhpolmx(nmodsa)
    character(len=8) :: vlfmod
    real(c_double) :: vlfmodes
    real(c_double) :: vlffreq(nmodsa)
    real(c_double) :: vlfnp(nmodsa)
    real(c_double) :: vlfdnp(nmodsa)
    real(c_double) :: vlfddnp(nmodsa)
    complex(c_double_complex) :: vlfeplus(nmodsa)
    complex(c_double_complex) :: vlfemin(nmodsa)
    real(c_double) :: vlfpol(nmodsa)
    real(c_double) :: vlfdpol(nmodsa)
    real(c_double) :: vlfddpol(nmodsa)
    character(len=8) :: vlfbes
    character(len=8) :: vlfnpvar
    real(c_double) :: vlfharms(nmodsa)
    real(c_double) :: vlfharm1(nmodsa)
    real(c_double) :: vlfnperp(nmodsa)
    real(c_double) :: vlfdnorm(nmodsa)
    real(c_double) :: vlfparmn(nmodsa)
    real(c_double) :: vlfparmx(nmodsa)
    real(c_double) :: vlfprpmn(nmodsa)
    real(c_double) :: vlfprpmx(nmodsa)
    real(c_double) :: rdc_upar_sign
    integer :: nrdc
    character(len=8) :: rdcmod
    character(len=8) :: rdc_clipping
    integer :: nrdcspecies(nrdca)
    real(c_double) :: rdcscale(nrdca)
    character(len=8) :: rdc_netcdf

    namelist/rfsetup/ &
         call_lh,call_ech,call_fw,ieqbrurf,urfncoef, &
         dlndau, &
         lh, &
         ech, &
         fw, &
         rftype, &
         rffile,rdcfile,rfread, &
         nharms,nharm1,nrfspecies,iurfcoll,iurfl, &
         nbssltbl,nondamp,nrfstep1,nrfstep2, &
         nrfpwr,nrfitr1,nrfitr2,nrfitr3, &
         nonrf,noffrf,nrf, &
         scaleurf,pwrscale,wdscale,urfrstrt,urfwrray, &
         nurftime,urftime,pwrscale1, &
         urfdmp,urfmult,urfmod, &
         vlhmod,vlhmodes,vparmin,vparmax,vprprop,vdalp,vlh_karney, &
         vlhprprp,vlhplse,vlhpon,vlhpoff,vprpmin,vprpmax,vlhpolmn, &
         vlhpolmx,vlfmod,vlfmodes,vlffreq,vlfnp,vlfdnp,vlfddnp, &
         vlfeplus,vlfemin,vlfpol,vlfdpol,vlfddpol,vlfbes,vlfnpvar, &
         vlfharms,vlfharm1,vlfnperp,vlfdnorm, &
         vlfparmn,vlfparmx,vlfprpmn,vlfprpmx,rdc_upar_sign,nrdc, &
         rdcmod,rdc_clipping,nrdcspecies,rdcscale,rdc_netcdf

    ! copy defaults to local vars

    call_lh = rfsetup_%call_lh
    call_ech = rfsetup_%call_ech
    call_fw = rfsetup_%call_fw
    ieqbrurf = rfsetup_%ieqbrurf
    urfncoef = rfsetup_%urfncoef
    dlndau = rfsetup_%dlndau
    lh = rfsetup_%lh
    ech = rfsetup_%ech
    fw = rfsetup_%fw
    rftype = rfsetup_%rftype
    rffile = rfsetup_%rffile
    rdcfile = rfsetup_%rdcfile
    rfread = rfsetup_%rfread
    nharms = rfsetup_%nharms
    nharm1 = rfsetup_%nharm1
    nrfspecies = rfsetup_%nrfspecies
    iurfcoll = rfsetup_%iurfcoll
    iurfl = rfsetup_%iurfl
    nbssltbl = rfsetup_%nbssltbl
    nondamp = rfsetup_%nondamp
    nrfstep1 = rfsetup_%nrfstep1
    nrfstep2 = rfsetup_%nrfstep2
    nrfpwr = rfsetup_%nrfpwr
    nrfitr1 = rfsetup_%nrfitr1
    nrfitr2 = rfsetup_%nrfitr2
    nrfitr3 = rfsetup_%nrfitr3
    nonrf = rfsetup_%nonrf
    noffrf = rfsetup_%noffrf
    nrf = rfsetup_%nrf
    scaleurf = rfsetup_%scaleurf
    pwrscale = rfsetup_%pwrscale
    wdscale = rfsetup_%wdscale
    urfrstrt = rfsetup_%urfrstrt
    urfwrray = rfsetup_%urfwrray
    nurftime = rfsetup_%nurftime
    urftime = rfsetup_%urftime
    pwrscale1 = rfsetup_%pwrscale1
    urfdmp = rfsetup_%urfdmp
    urfmult = rfsetup_%urfmult
    urfmod = rfsetup_%urfmod
    vlhmod = rfsetup_%vlhmod
    vlhmodes = rfsetup_%vlhmodes
    vparmin = rfsetup_%vparmin
    vparmax = rfsetup_%vparmax
    vprprop = rfsetup_%vprprop
    vdalp = rfsetup_%vdalp
    vlh_karney = rfsetup_%vlh_karney
    vlhprprp = rfsetup_%vlhprprp
    vlhplse = rfsetup_%vlhplse
    vlhpon = rfsetup_%vlhpon
    vlhpoff = rfsetup_%vlhpoff
    vprpmin = rfsetup_%vprpmin
    vprpmax = rfsetup_%vprpmax
    vlhpolmn = rfsetup_%vlhpolmn
    vlhpolmx = rfsetup_%vlhpolmx
    vlfmod = rfsetup_%vlfmod
    vlfmodes = rfsetup_%vlfmodes
    vlffreq = rfsetup_%vlffreq
    vlfnp = rfsetup_%vlfnp
    vlfdnp = rfsetup_%vlfdnp
    vlfddnp = rfsetup_%vlfddnp
    vlfeplus = rfsetup_%vlfeplus
    vlfemin = rfsetup_%vlfemin
    vlfpol = rfsetup_%vlfpol
    vlfdpol = rfsetup_%vlfdpol
    vlfddpol = rfsetup_%vlfddpol
    vlfbes = rfsetup_%vlfbes
    vlfnpvar = rfsetup_%vlfnpvar
    vlfharms = rfsetup_%vlfharms
    vlfharm1 = rfsetup_%vlfharm1
    vlfnperp = rfsetup_%vlfnperp
    vlfdnorm = rfsetup_%vlfdnorm
    vlfparmn = rfsetup_%vlfparmn
    vlfparmx = rfsetup_%vlfparmx
    vlfprpmn = rfsetup_%vlfprpmn
    vlfprpmx = rfsetup_%vlfprpmx
    rdc_upar_sign = rfsetup_%rdc_upar_sign
    nrdc = rfsetup_%nrdc
    rdcmod = rfsetup_%rdcmod
    rdc_clipping = rfsetup_%rdc_clipping
    nrdcspecies = rfsetup_%nrdcspecies
    rdcscale = rfsetup_%rdcscale
    rdc_netcdf = rfsetup_%rdc_netcdf

    ! read the nml, which will write into the local vars

    call maybe_nml_open(nml_file)
    read(nml_fd, rfsetup)

    ! external codes can call this, which packs the setup0 derived type.
    call set_rfsetup(call_lh, call_ech, call_fw, ieqbrurf, urfncoef, &
         dlndau, lh, ech, fw, rftype, rffile, rdcfile, rfread, nharms, nharm1, &
         nrfspecies, iurfcoll, iurfl, nbssltbl, nondamp, nrfstep1, nrfstep2, &
         nrfpwr, nrfitr1, nrfitr2, nrfitr3, nonrf, noffrf, nrf, scaleurf, &
         pwrscale, wdscale, urfrstrt, urfwrray, nurftime, urftime, pwrscale1, &
         urfdmp, urfmult, urfmod, vlhmod, vlhmodes, vparmin, vparmax, vprprop, &
         vdalp, vlh_karney, vlhprprp, vlhplse, vlhpon, vlhpoff, vprpmin, &
         vprpmax, vlhpolmn, vlhpolmx, vlfmod, vlfmodes, vlffreq, vlfnp, &
         vlfdnp, vlfddnp, vlfeplus, vlfemin, vlfpol, vlfdpol, vlfddpol, vlfbes, &
         vlfnpvar, vlfharms, vlfharm1, vlfnperp, vlfdnorm, vlfparmn, vlfparmx, &
         vlfprpmn, vlfprpmx, rdc_upar_sign, nrdc, rdcmod, rdc_clipping, &
         nrdcspecies, rdcscale, rdc_netcdf, debug_print)

    ! we optionally close the nml file.
    if (present(close_nml_file)) then
       if(close_nml_file) then
          call nml_close()
       end if
    endif

  end subroutine get_rfsetup_from_nml


  subroutine set_rfsetup(call_lh, call_ech, call_fw, ieqbrurf, urfncoef, &
       dlndau, lh, ech, fw, rftype, rffile, rdcfile, rfread, nharms, nharm1, &
       nrfspecies, iurfcoll, iurfl, nbssltbl, nondamp, nrfstep1, nrfstep2, &
       nrfpwr, nrfitr1, nrfitr2, nrfitr3, nonrf, noffrf, nrf, scaleurf, &
       pwrscale, wdscale, urfrstrt, urfwrray, nurftime, urftime, pwrscale1, &
       urfdmp, urfmult, urfmod, vlhmod, vlhmodes, vparmin, vparmax, vprprop, &
       vdalp, vlh_karney, vlhprprp, vlhplse, vlhpon, vlhpoff, vprpmin, &
       vprpmax, vlhpolmn, vlhpolmx, vlfmod, vlfmodes, vlffreq, vlfnp, &
       vlfdnp, vlfddnp, vlfeplus, vlfemin, vlfpol, vlfdpol, vlfddpol, vlfbes, &
       vlfnpvar, vlfharms, vlfharm1, vlfnperp, vlfdnorm, vlfparmn, vlfparmx, &
       vlfprpmn, vlfprpmx, rdc_upar_sign, nrdc, rdcmod, rdc_clipping, &
       nrdcspecies, rdcscale, rdc_netcdf, debug_print)
    implicit none
    logical, intent(in), optional :: debug_print
    !
    character(len=8), intent(in), optional :: call_lh
    character(len=8), intent(in), optional :: call_ech
    character(len=8), intent(in), optional :: call_fw
    integer, intent(in), optional :: ieqbrurf
    integer, intent(in), optional :: urfncoef
    real(c_double) , intent(in), optional :: dlndau(nmodsa)
    character(len=8), intent(in), optional :: lh
    character(len=8), intent(in), optional :: ech
    character(len=8), intent(in), optional :: fw
    character(len=8), intent(in), optional :: rftype(nmodsa)
    character(len=256), intent(in), optional :: rffile(nmodsa)
    character(len=256), intent(in), optional :: rdcfile(nrdca)
    character(len=8), intent(in), optional :: rfread
    integer, intent(in), optional :: nharms(nmodsa)
    integer, intent(in), optional :: nharm1(nmodsa)
    integer, intent(in), optional :: nrfspecies(nmodsa)
    character(len=8), intent(in), optional :: iurfcoll(nmodsa)
    character(len=8), intent(in), optional :: iurfl(nmodsa)
    integer, intent(in), optional :: nbssltbl
    integer, intent(in), optional :: nondamp
    integer, intent(in), optional ::nrfstep1(nmodsa)
    integer, intent(in), optional :: nrfstep2
    integer, intent(in), optional :: nrfpwr
    integer, intent(in), optional :: nrfitr1
    integer, intent(in), optional :: nrfitr2
    integer, intent(in), optional :: nrfitr3
    integer, intent(in), optional :: nonrf(ngena)
    integer, intent(in), optional :: noffrf(ngena)
    integer, intent(in), optional :: nrf
    character(len=8), intent(in), optional :: scaleurf
    real(c_double), intent(in), optional :: pwrscale(nmodsa)
    real(c_double), intent(in), optional :: wdscale(nmodsa)
    character(len=8), intent(in), optional :: urfrstrt
    character(len=8), intent(in), optional :: urfwrray
    integer, intent(in), optional :: nurftime
    real(c_double), intent(in), optional :: urftime(nbctimea)
    real(c_double), intent(in), optional :: pwrscale1(nbctimea)
    character(len=8), intent(in), optional :: urfdmp
    real(c_double), intent(in), optional :: urfmult
    character(len=8), intent(in), optional :: urfmod
    character(len=8), intent(in), optional :: vlhmod
    real(c_double), intent(in), optional :: vlhmodes
    real(c_double), intent(in), optional :: vparmin(nmodsa)
    real(c_double), intent(in), optional :: vparmax(nmodsa)
    character(len=8), intent(in), optional :: vprprop
    real(c_double), intent(in), optional :: vdalp
    real(c_double), intent(in), optional :: vlh_karney
    character(len=8), intent(in), optional :: vlhprprp(nmodsa)
    character(len=8), intent(in), optional :: vlhplse
    real(c_double), intent(in), optional :: vlhpon
    real(c_double), intent(in), optional :: vlhpoff
    real(c_double), intent(in), optional :: vprpmin(nmodsa)
    real(c_double), intent(in), optional :: vprpmax(nmodsa)
    real(c_double), intent(in), optional :: vlhpolmn(nmodsa)
    real(c_double), intent(in), optional :: vlhpolmx(nmodsa)
    character(len=8), intent(in), optional :: vlfmod
    real(c_double), intent(in), optional :: vlfmodes
    real(c_double), intent(in), optional :: vlffreq(nmodsa)
    real(c_double), intent(in), optional :: vlfnp(nmodsa)
    real(c_double), intent(in), optional :: vlfdnp(nmodsa)
    real(c_double), intent(in), optional :: vlfddnp(nmodsa)
    complex(c_double_complex), intent(in), optional :: vlfeplus(nmodsa)
    complex(c_double_complex), intent(in), optional :: vlfemin(nmodsa)
    real(c_double), intent(in), optional :: vlfpol(nmodsa)
    real(c_double), intent(in), optional :: vlfdpol(nmodsa)
    real(c_double), intent(in), optional :: vlfddpol(nmodsa)
    character(len=8), intent(in), optional :: vlfbes
    character(len=8), intent(in), optional :: vlfnpvar
    real(c_double), intent(in), optional :: vlfharms(nmodsa)
    real(c_double), intent(in), optional :: vlfharm1(nmodsa)
    real(c_double), intent(in), optional :: vlfnperp(nmodsa)
    real(c_double), intent(in), optional :: vlfdnorm(nmodsa)
    real(c_double), intent(in), optional :: vlfparmn(nmodsa)
    real(c_double), intent(in), optional :: vlfparmx(nmodsa)
    real(c_double), intent(in), optional :: vlfprpmn(nmodsa)
    real(c_double), intent(in), optional :: vlfprpmx(nmodsa)
    real(c_double), intent(in), optional :: rdc_upar_sign
    integer, intent(in), optional :: nrdc
    character(len=8), intent(in), optional :: rdcmod
    character(len=8), intent(in), optional :: rdc_clipping
    integer, intent(in), optional :: nrdcspecies(nrdca)
    real(c_double), intent(in), optional :: rdcscale(nrdca)
    character(len=8), intent(in), optional :: rdc_netcdf


    if(present(call_lh)) rfsetup%call_lh = call_lh
    if(present(call_ech)) rfsetup%call_ech = call_ech
    if(present(call_fw)) rfsetup%call_fw = call_fw
    if(present(ieqbrurf)) rfsetup%ieqbrurf = ieqbrurf
    if(present(urfncoef)) rfsetup%urfncoef = urfncoef
    if(present(dlndau)) rfsetup%dlndau = dlndau
    if(present(lh)) rfsetup%lh = lh
    if(present(ech)) rfsetup%ech = ech
    if(present(fw)) rfsetup%fw = fw
    if(present(rftype)) rfsetup%rftype = rftype
    if(present(rffile)) rfsetup%rffile = rffile
    if(present(rdcfile)) rfsetup%rdcfile = rdcfile
    if(present(rfread)) rfsetup%rfread = rfread
    if(present(nharms)) rfsetup%nharms = nharms
    if(present(nharm1)) rfsetup%nharm1 = nharm1
    if(present(nrfspecies)) rfsetup%nrfspecies = nrfspecies
    if(present(iurfcoll)) rfsetup%iurfcoll = iurfcoll
    if(present(iurfl)) rfsetup%iurfl = iurfl
    if(present(nbssltbl)) rfsetup%nbssltbl = nbssltbl
    if(present(nondamp)) rfsetup%nondamp = nondamp
    if(present(nrfstep1)) rfsetup%nrfstep1 = nrfstep1
    if(present(nrfstep2)) rfsetup%nrfstep2 = nrfstep2
    if(present(nrfpwr)) rfsetup%nrfpwr = nrfpwr
    if(present(nrfitr1)) rfsetup%nrfitr1 = nrfitr1
    if(present(nrfitr2)) rfsetup%nrfitr2 = nrfitr2
    if(present(nrfitr3)) rfsetup%nrfitr3 = nrfitr3
    if(present(nonrf)) rfsetup%nonrf = nonrf
    if(present(noffrf)) rfsetup%noffrf = noffrf
    if(present(nrf)) rfsetup%nrf = nrf
    if(present(scaleurf)) rfsetup%scaleurf = scaleurf
    if(present(pwrscale)) rfsetup%pwrscale = pwrscale
    if(present(wdscale)) rfsetup%wdscale = wdscale
    if(present(urfrstrt)) rfsetup%urfrstrt = urfrstrt
    if(present(urfwrray)) rfsetup%urfwrray = urfwrray
    if(present(nurftime)) rfsetup%nurftime = nurftime
    if(present(urftime)) rfsetup%urftime = urftime
    if(present(pwrscale1)) rfsetup%pwrscale1 = pwrscale1
    if(present(urfdmp)) rfsetup%urfdmp = urfdmp
    if(present(urfmult)) rfsetup%urfmult = urfmult
    if(present(urfmod)) rfsetup%urfmod = urfmod
    if(present(vlhmod)) rfsetup%vlhmod = vlhmod
    if(present(vlhmodes)) rfsetup%vlhmodes = vlhmodes
    if(present(vparmin)) rfsetup%vparmin = vparmin
    if(present(vparmax)) rfsetup%vparmax = vparmax
    if(present(vprprop)) rfsetup%vprprop = vprprop
    if(present(vdalp)) rfsetup%vdalp = vdalp
    if(present(vlh_karney)) rfsetup%vlh_karney = vlh_karney
    if(present(vlhprprp)) rfsetup%vlhprprp = vlhprprp
    if(present(vlhplse)) rfsetup%vlhplse = vlhplse
    if(present(vlhpon)) rfsetup%vlhpon = vlhpon
    if(present(vlhpoff)) rfsetup%vlhpoff = vlhpoff
    if(present(vprpmin)) rfsetup%vprpmin = vprpmin
    if(present(vprpmax)) rfsetup%vprpmax = vprpmax
    if(present(vlhpolmn)) rfsetup%vlhpolmn = vlhpolmn
    if(present(vlhpolmx)) rfsetup%vlhpolmx = vlhpolmx
    if(present(vlfmod)) rfsetup%vlfmod = vlfmod
    if(present(vlfmodes)) rfsetup%vlfmodes = vlfmodes
    if(present(vlffreq)) rfsetup%vlffreq = vlffreq
    if(present(vlfnp)) rfsetup%vlfnp = vlfnp
    if(present(vlfdnp)) rfsetup%vlfdnp = vlfdnp
    if(present(vlfddnp)) rfsetup%vlfddnp = vlfddnp
    if(present(vlfeplus)) rfsetup%vlfeplus = vlfeplus
    if(present(vlfemin)) rfsetup%vlfemin = vlfemin
    if(present(vlfpol)) rfsetup%vlfpol = vlfpol
    if(present(vlfdpol)) rfsetup%vlfdpol = vlfdpol
    if(present(vlfddpol)) rfsetup%vlfddpol = vlfddpol
    if(present(vlfbes)) rfsetup%vlfbes = vlfbes
    if(present(vlfnpvar)) rfsetup%vlfnpvar = vlfnpvar
    if(present(vlfharms)) rfsetup%vlfharms = vlfharms
    if(present(vlfharm1)) rfsetup%vlfharm1 = vlfharm1
    if(present(vlfnperp)) rfsetup%vlfnperp = vlfnperp
    if(present(vlfdnorm)) rfsetup%vlfdnorm = vlfdnorm
    if(present(vlfparmn)) rfsetup%vlfparmn = vlfparmn
    if(present(vlfparmx)) rfsetup%vlfparmx = vlfparmx
    if(present(vlfprpmn)) rfsetup%vlfprpmn = vlfprpmn
    if(present(vlfprpmx)) rfsetup%vlfprpmx = vlfprpmx
    if(present(rdc_upar_sign)) rfsetup%rdc_upar_sign = rdc_upar_sign
    if(present(nrdc)) rfsetup%nrdc = nrdc
    if(present(rdcmod)) rfsetup%rdcmod = rdcmod
    if(present(rdc_clipping)) rfsetup%rdc_clipping = rdc_clipping
    if(present(nrdcspecies)) rfsetup%nrdcspecies = nrdcspecies
    if(present(rdcscale)) rfsetup%rdcscale = rdcscale
    if(present(rdc_netcdf)) rfsetup%rdc_netcdf = rdc_netcdf

    if ( present(debug_print)) then
       if (debug_print) call print_rfsetup
    end if

  end subroutine set_rfsetup

  subroutine print_rfsetup()
    namelist /rfsetup_nml/ rfsetup
    WRITE(*, *) "!----  BEGIN RFSETUP DUMP"
    WRITE(*, nml = rfsetup_nml)
    WRITE(*, *)  "!----  END RFSETUP DUMP"
  end subroutine print_rfsetup


  integer function newunit(unit)
    ! Thanks fortran wiki !
    ! This is a simple function to search for an available unit.
    ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
    ! The UNIT value is returned by the function, and also by the optional
    ! argument. This allows the function to be used directly in an OPEN
    ! statement, and optionally save the result in a local variable.
    ! If no units are available, -1 is returned.
    integer, intent(out), optional :: unit
    ! local
    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: opened
    integer :: lun
    ! begin
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
       inquire(unit=lun,opened=opened)
       if (.not. opened) then
          newunit=lun
          exit
       end if
    end do
    if (present(unit)) unit=newunit
  end function newunit

  subroutine get_trsetup_from_nml(nml_file, close_nml_file, debug_print)
    implicit none
    character(len=*), intent(in) :: nml_file
    logical, intent(in), optional :: close_nml_file
    logical, intent(in), optional :: debug_print
    ! local
    type(trsetup_t) :: trsetup_
    !
    real(c_double) :: advectr
    character(len=8) :: difus_type(1:ngena)
    real(c_double) :: difusr
    real(c_double) :: difus_rshape(8)
    real(c_double) :: difus_vshape(4)
    real(c_double) :: difin(1:njenea)
    character(len=8) :: difus_io(1:ngena)
    character(len=256) :: difus_io_file
    real(c_double) :: difus_io_drrscale(1:nbctimea, 1:ngena)
    real(c_double) :: difus_io_drscale(1:nbctimea, 1:ngena)
    real(c_double) :: difus_io_t(1:nbctimea)
    character(len=8) :: pinch
    real(c_double) :: relaxden
    character(len=8) :: relaxtsp
    character(len=8) :: transp= " disabled"
    character(len=8) :: adimeth
    integer :: nonadi
    integer :: nontran
    integer :: nofftran
    integer :: nonelpr
    integer :: noffelpr
    integer :: ndifus_io_t

    namelist/trsetup/ &
         advectr, &
         difus_type,difusr,difus_rshape,difus_vshape,difin, &
         difus_io,difus_io_file, &
         difus_io_drrscale,difus_io_drscale,difus_io_t, &
         pinch, &
         relaxden,relaxtsp, &
         transp,adimeth,nonadi, &
         nontran,nofftran,nonelpr,noffelpr,ndifus_io_t

    ! copy defaults to local vars

    advectr = trsetup_%advectr
    difus_type = trsetup_%difus_type
    difusr = trsetup_%difusr
    difus_rshape = trsetup_%difus_rshape
    difus_vshape = trsetup_%difus_vshape
    difin = trsetup_%difin
    difus_io = trsetup_%difus_io
    difus_io_file = trsetup_%difus_io_file
    difus_io_drrscale = trsetup_%difus_io_drrscale
    difus_io_drscale = trsetup_%difus_io_drscale
    difus_io_t = trsetup_%difus_io_t
    pinch = trsetup_%pinch
    relaxden = trsetup_%relaxden
    relaxtsp = trsetup_%relaxtsp
    transp = trsetup_%transp
    adimeth = trsetup_%adimeth
    nonadi = trsetup_%nonadi
    nontran = trsetup_%nontran
    nofftran = trsetup_%nofftran
    nonelpr = trsetup_%nonelpr
    noffelpr = trsetup_%noffelpr
    ndifus_io_t = trsetup_%ndifus_io_t


    ! read the nml, which will write into the local vars

    call maybe_nml_open(nml_file)
    read(nml_fd, trsetup)

    ! external codes can call this, which packs the setup0 derived type.
    call set_trsetup(advectr, &
         difus_type,difusr,difus_rshape,difus_vshape,difin, &
         difus_io,difus_io_file, &
         difus_io_drrscale,difus_io_drscale,difus_io_t, &
         pinch, &
         relaxden,relaxtsp, &
         transp,adimeth,nonadi, &
         nontran,nofftran,nonelpr,noffelpr,ndifus_io_t, &
         debug_print)

    ! we optionally close the nml file.
    if (present(close_nml_file)) then
       if(close_nml_file) then
          call nml_close()
       end if
    endif
  end subroutine get_trsetup_from_nml

  subroutine set_trsetup(advectr, &
       difus_type,difusr,difus_rshape,difus_vshape,difin, &
       difus_io,difus_io_file, &
       difus_io_drrscale,difus_io_drscale,difus_io_t, &
       pinch, &
       relaxden,relaxtsp, &
       transp,adimeth,nonadi, &
       nontran,nofftran,nonelpr,noffelpr,ndifus_io_t, &
       debug_print)
    logical, intent(in), optional :: debug_print
    !
    real(c_double), intent(in), optional :: advectr
    character(len=8), intent(in), optional :: difus_type(1:ngena)
    real(c_double), intent(in), optional :: difusr
    real(c_double), intent(in), optional :: difus_rshape(8)
    real(c_double), intent(in), optional :: difus_vshape(4)
    real(c_double), intent(in), optional :: difin(1:njenea)
    character(len=8), intent(in), optional :: difus_io(1:ngena)
    character(len=256), intent(in), optional :: difus_io_file
    real(c_double), intent(in), optional :: difus_io_drrscale(1:nbctimea, 1:ngena)
    real(c_double), intent(in), optional :: difus_io_drscale(1:nbctimea, 1:ngena)
    real(c_double), intent(in), optional :: difus_io_t(1:nbctimea)
    character(len=8), intent(in), optional :: pinch
    real(c_double), intent(in), optional :: relaxden
    character(len=8), intent(in), optional :: relaxtsp
    character(len=8), intent(in), optional :: transp
    character(len=8), intent(in), optional :: adimeth
    integer, intent(in), optional :: nonadi
    integer, intent(in), optional :: nontran
    integer, intent(in), optional :: nofftran
    integer, intent(in), optional :: nonelpr
    integer, intent(in), optional :: noffelpr
    integer, intent(in), optional :: ndifus_io_t

    if(present(advectr)) trsetup%advectr = advectr
    if(present(difus_type)) trsetup%difus_type = difus_type
    if(present(difusr)) trsetup%difusr = difusr
    if(present(difus_rshape)) trsetup%difus_rshape = difus_rshape
    if(present(difus_vshape)) trsetup%difus_vshape = difus_vshape
    if(present(difin)) trsetup%difin = difin
    if(present(difus_io)) trsetup%difus_io = difus_io
    if(present(difus_io_file)) trsetup%difus_io_file = difus_io_file
    if(present(difus_io_drrscale)) trsetup%difus_io_drrscale = difus_io_drrscale
    if(present(difus_io_drscale)) trsetup%difus_io_drscale = difus_io_drscale
    if(present(difus_io_t)) trsetup%difus_io_t = difus_io_t
    if(present(pinch)) trsetup%pinch = pinch
    if(present(relaxden)) trsetup%relaxden = relaxden
    if(present(relaxtsp)) trsetup%relaxtsp = relaxtsp
    if(present(transp)) trsetup%transp = transp
    if(present(adimeth)) trsetup%adimeth = adimeth
    if(present(nonadi)) trsetup%nonadi = nonadi
    if(present(nontran)) trsetup%nontran = nontran
    if(present(nofftran)) trsetup%nofftran = nofftran
    if(present(nonelpr)) trsetup%nonelpr = nonelpr
    if(present(noffelpr)) trsetup%noffelpr = noffelpr
    if(present(ndifus_io_t)) trsetup%ndifus_io_t = ndifus_io_t

    if ( present(debug_print)) then
       if (debug_print) call print_trsetup
    end if

  end subroutine set_trsetup

  subroutine print_trsetup()
    namelist /trsetup_nml/ trsetup
    WRITE(*, *) "!----  BEGIN TRSETUP DUMP"
    WRITE(*, nml = trsetup_nml)
    WRITE(*, *)  "!----  END TRSETUP DUMP"
  end subroutine print_trsetup

  subroutine get_sousetup_from_nml(nml_file, close_nml_file, debug_print)
    implicit none
    character(len=*), intent(in) :: nml_file
    logical, intent(in), optional :: close_nml_file
    logical, intent(in), optional :: debug_print
    ! local
    type(sousetup_t) :: sousetup_
    !
    real(c_double) :: asorz(1:ngena,1:nsoa,0:lrza)
    real(c_double) :: asor(1:ngena,1:nsoa,1:lrza)
    character(len=8) :: flemodel
    integer :: nonso(1:ngena,1:nsoa)
    integer :: noffso(1:ngena,1:nsoa)
    integer :: nso
    integer :: nsou
    character(len=8) :: pltso
    real(c_double) :: mpwrsou(0:ngena)
    real(c_double) :: npwrsou(0:ngena)
    character(len=8) :: soucoord
    character(len=8) :: knockon
    character(len=8) :: komodel
    integer :: nkorfn
    integer :: nonko
    integer :: noffko
    real(c_double) :: soffvte
    real(c_double) :: soffpr
    real(c_double) :: isoucof
    real(c_double) :: faccof
    integer :: jfl
    real(c_double) :: xlfac
    real(c_double) :: xlpctlwr
    real(c_double) :: xlpctmdl
    real(c_double) :: xllwr
    real(c_double) :: xlmdl
    ! the following will be initilized to non trivial values in the setter
    ! we should probably init these to nans, but fortran makes that a chore
    real(c_double) :: scm2z(ngena,nsoa,0:lrza)
    real(c_double) :: scm2(1:ngena,1:nsoa)
    real(c_double) :: sellm1(1:ngena,1:nsoa)
    real(c_double) :: sellm2(1:ngena,1:nsoa)
    real(c_double) :: seppm1(1:ngena,1:nsoa)
    real(c_double) :: sellm1z(ngena,nsoa,0:lrza)
    real(c_double) :: sellm2z(ngena,nsoa,0:lrza)
    real(c_double) :: seppm2(1:ngena,1:nsoa)
    real(c_double) :: sem1(ngena,nsoa)
    real(c_double) :: sem2(ngena,nsoa)
    real(c_double) :: seppm1z(ngena,nsoa,0:lrza)
    real(c_double) :: sem1z(ngena,nsoa,0:lrza)
    real(c_double) :: sem2z(ngena,nsoa,0:lrza)
    real(c_double) :: sthm1z(ngena,nsoa,0:lrza)
    real(c_double) :: seppm2z(ngena,nsoa,0:lrza)
    real(c_double) :: szm1z(ngena,nsoa,0:lrza)
    real(c_double) :: szm2z(ngena,nsoa,0:lrza)
    real(c_double) :: sthm1(ngena,nsoa)
    real(c_double) :: szm1(ngena,nsoa)
    real(c_double) :: szm2(ngena,nsoa)

    !..................................................................
    !     Namelist (sousetup) for "sou" simple source routines.
    !..................................................................

    namelist/sousetup/ &
         asorz,asor,flemodel, &
         nonso,noffso,nso,nsou, &
         pltso,mpwrsou,npwrsou, &
         scm2z,szm1z,scm2,sellm1,sellm2,seppm1, &
         sellm1z,sellm2z,seppm2,sem1,sem2, &
         seppm1z,sem1z,sem2z,sthm1z, &
         seppm2z,soucoord,knockon,komodel,nkorfn,nonko,noffko,soffvte, &
         soffpr,isoucof,faccof,jfl,xlfac,xlpctlwr,xlpctmdl,xllwr,xlmdl, &
         szm2z,sthm1,szm1,szm2

    ! copy defaults to local vars
    asorz = sousetup_%asorz
    asor = sousetup_%asor
    flemodel = sousetup_%flemodel
    nonso = sousetup_%nonso
    noffso = sousetup_%noffso
    nso = sousetup_%nso
    nsou = sousetup_%nsou
    pltso = sousetup_%pltso
    mpwrsou = sousetup_%mpwrsou
    npwrsou = sousetup_%npwrsou
    scm2z = sousetup_%scm2z
    szm1z = sousetup_%szm1z
    scm2 = sousetup_%scm2
    sellm1 = sousetup_%sellm1
    sellm2 = sousetup_%sellm2
    seppm1 = sousetup_%seppm1
    sellm1z = sousetup_%sellm1z
    sellm2z = sousetup_%sellm2z
    seppm2 = sousetup_%seppm2
    sem1 = sousetup_%sem1
    sem2 = sousetup_%sem2
    seppm1z = sousetup_%seppm1z
    sem1z = sousetup_%sem1z
    sem2z = sousetup_%sem2z
    sthm1z = sousetup_%sthm1z
    seppm2z = sousetup_%seppm2z
    soucoord = sousetup_%soucoord
    knockon = sousetup_%knockon
    komodel = sousetup_%komodel
    nkorfn = sousetup_%nkorfn
    nonko = sousetup_%nonko
    noffko = sousetup_%noffko
    soffvte = sousetup_%soffvte
    soffpr = sousetup_%soffpr
    isoucof = sousetup_%isoucof
    faccof = sousetup_%faccof
    jfl = sousetup_%jfl
    xlfac = sousetup_%xlfac
    xlpctlwr = sousetup_%xlpctlwr
    xlpctmdl = sousetup_%xlpctmdl
    xllwr = sousetup_%xllwr
    xlmdl = sousetup_%xlmdl
    szm1z = sousetup_%szm2z
    szm2z = sousetup_%szm2z
    sthm1 = sousetup_%sthm1
    szm1 = sousetup_%szm1
    szm2 = sousetup_%szm2

    ! read the nml, which will write into the local vars

    call maybe_nml_open(nml_file)
    read(nml_fd, sousetup)

    ! external codes can call this, which packs the setup0 derived type.
    call set_sousetup(asorz,asor,flemodel, &
         nonso,noffso,nso,nsou, &
         pltso,mpwrsou,npwrsou, &
         scm2z,szm1z,scm2,sellm1,sellm2,seppm1, &
         sellm1z,sellm2z,seppm2,sem1,sem2, &
         seppm1z,sem1z,sem2z,sthm1z, &
         seppm2z,soucoord,knockon,komodel,nkorfn,nonko,noffko,soffvte, &
         soffpr,isoucof,faccof,jfl,xlfac,xlpctlwr,xlpctmdl,xllwr,xlmdl, &
         szm2z,sthm1,szm1,szm2, &
         debug_print)

    ! we optionally close the nml file.
    if (present(close_nml_file)) then
       if(close_nml_file) then
          call nml_close()
       end if
    endif
  end subroutine get_sousetup_from_nml

  subroutine sou_init_loop(lhs, rhs)
    real(c_double), intent(inout) :: lhs(:,:,:)
    real(c_double), intent(inout) :: rhs(:,:)
    integer :: k, m, ll
    ! many of the variables were assigned with this loop in aindflt
    do k=1,ngena
       do m=1,nsoa
          do ll=0,lrza
             lhs(k,m,ll)=rhs(k,m)
          enddo
       enddo
    enddo
  end subroutine sou_init_loop

  subroutine set_sousetup(asorz,asor,flemodel, &
       nonso,noffso,nso,nsou, &
       pltso,mpwrsou,npwrsou, &
       scm2z,szm1z,scm2,sellm1,sellm2,seppm1, &
       sellm1z,sellm2z,seppm2,sem1,sem2, &
       seppm1z,sem1z,sem2z,sthm1z, &
       seppm2z,soucoord,knockon,komodel,nkorfn,nonko,noffko,soffvte, &
       soffpr,isoucof,faccof,jfl,xlfac,xlpctlwr,xlpctmdl,xllwr,xlmdl, &
       szm2z,sthm1,szm1,szm2, &
       debug_print)

    logical, intent(in), optional :: debug_print
    integer :: k, m, ll
    !
    real(c_double), intent(in), optional :: asorz(1:ngena,1:nsoa,0:lrza)
    real(c_double), intent(in), optional :: asor(1:ngena,1:nsoa,1:lrza)
    character(len=8), intent(in), optional :: flemodel
    integer, intent(in), optional :: nonso(1:ngena,1:nsoa)
    integer, intent(in), optional :: noffso(1:ngena,1:nsoa)
    integer, intent(in), optional :: nso
    integer, intent(in), optional :: nsou
    character(len=8), intent(in), optional :: pltso
    real(c_double), intent(in), optional :: mpwrsou(0:ngena)
    real(c_double), intent(in), optional :: npwrsou(0:ngena)
    character(len=8), intent(in), optional :: soucoord
    character(len=8), intent(in), optional :: knockon
    character(len=8), intent(in), optional :: komodel
    integer, intent(in), optional :: nkorfn
    integer, intent(in), optional :: nonko
    integer, intent(in), optional :: noffko
    real(c_double), intent(in), optional :: soffvte
    real(c_double), intent(in), optional :: soffpr
    real(c_double), intent(in), optional :: isoucof
    real(c_double), intent(in), optional :: faccof
    integer, intent(in), optional :: jfl
    real(c_double), intent(in), optional :: xlfac
    real(c_double), intent(in), optional :: xlpctlwr
    real(c_double), intent(in), optional :: xlpctmdl
    real(c_double), intent(in), optional :: xllwr
    real(c_double), intent(in), optional :: xlmdl
    real(c_double), intent(in), optional :: scm2z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: scm2(1:ngena,1:nsoa)
    real(c_double), intent(in), optional :: sellm1(1:ngena,1:nsoa)
    real(c_double), intent(in), optional :: sellm2(1:ngena,1:nsoa)
    real(c_double), intent(in), optional :: seppm1(1:ngena,1:nsoa)
    real(c_double), intent(in), optional :: sellm1z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: sellm2z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: seppm2(1:ngena,1:nsoa)
    real(c_double), intent(in), optional :: sem1(ngena,nsoa)
    real(c_double), intent(in), optional :: sem2(ngena,nsoa)
    real(c_double), intent(in), optional :: seppm1z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: sem1z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: sem2z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: sthm1z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: seppm2z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: szm1z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: szm2z(ngena,nsoa,0:lrza)
    real(c_double), intent(in), optional :: sthm1(ngena,nsoa)
    real(c_double), intent(in), optional :: szm1(ngena,nsoa)
    real(c_double), intent(in), optional :: szm2(ngena,nsoa)


    if(present(asor)) then
       sousetup%asor = asor
    else
       !BH080125  usage of asor.
       !BH080125         asor(1,1,ll)=.25e+13
       !BH080125         asor(1,2,ll)=3.25e+1
    end if
    if(present(asorz)) then
       sousetup%asorz = asorz
    else
       do ll=1,lrza
          sousetup%asorz(k,m,ll)=asor(k,m,ll)
       enddo
    end if
    if(present(flemodel)) sousetup%flemodel = flemodel
    if(present(nonso)) sousetup%nonso = nonso
    if(present(noffso)) sousetup%noffso = noffso
    if(present(nso)) sousetup%nso = nso
    if(present(nsou)) sousetup%nsou = nsou
    if(present(pltso)) sousetup%pltso = pltso
    if(present(mpwrsou)) then
       sousetup%mpwrsou = mpwrsou
    else
       sousetup%mpwrsou(0)=1.
    end if
    if(present(npwrsou)) then
       sousetup%npwrsou = npwrsou
    else
       sousetup%npwrsou(0)=2.
    end if

    if(present(scm2)) then
       sousetup%scm2 = scm2
    else
       sousetup%scm2(1,1)=.001
       sousetup%scm2(1,2)=10000.
    end if
    if(present(scm2z)) then
       sousetup%scm2z = scm2z
    else
       call sou_init_loop(sousetup%scm2z, sousetup%scm2)
    end if

    if(present(szm1)) then
       sousetup%szm1 = szm1
    else
       sousetup%szm1(1,1)=0.
       sousetup%szm1(1,2)=0.
    end if

    if(present(szm1z)) then
       sousetup%szm1z = szm1z
    else
       call sou_init_loop(sousetup%szm1z, sousetup%szm1)
    end if

    if(present(szm2)) then
       sousetup%szm2 = szm2
    else
       sousetup%szm2(1,1)=1.e+5
       sousetup%szm2(1,2)=1.e+5
    end if
    if(present(szm2z)) then
       sousetup%szm2z = szm2z
    else
       call sou_init_loop(sousetup%szm2z, sousetup%szm2)
    end if

    if(present(sellm1)) then
       sousetup%sellm1 = sellm1
    else
       sousetup%sellm1(1,1)=1.
       sousetup%sellm1(1,2)=1.
    end if
    if(present(sellm1z)) then
       sousetup%sellm1z = sellm1z
    else
       call sou_init_loop(sousetup%sellm1z, sousetup%sellm1)
    end if

    if(present(sellm2)) then
       sousetup%sellm2 = sellm2
    else
       sousetup%sellm2(1,1)=1.
       sousetup%sellm2(1,2)=1.
    end if

    if(present(sellm2z)) then
       sousetup%sellm2z = sellm2z
    else
       call sou_init_loop(sousetup%sellm2z, sousetup%sellm2)
    end if

    if(present(seppm1)) then
       sousetup%seppm1 = seppm1
    else
       sousetup%seppm1(1,1)=1.
       sousetup%seppm1(1,2)=1.
    end if
    if(present(seppm1z)) then
       sousetup%seppm1z = seppm1z
    else
       call sou_init_loop(sousetup%seppm1z, sousetup%seppm1)
    end if

    if(present(seppm2)) then
       sousetup%seppm2 = seppm2
    else
       sousetup%seppm2(1,1)=1.
       sousetup%seppm2(1,2)=1.
    end if
    if(present(seppm2z)) then
       sousetup%seppm2z = seppm2z
    else
       call sou_init_loop(sousetup%seppm2z, sousetup%seppm2)
    end if

    if(present(sem1)) then
       sousetup%sem1 = sem1
    else
       sousetup%sem1(1,1)=1600.
       sousetup%sem1(1,2)=0.
    end if
    if(present(sem1z)) then
       sousetup%sem1z = sem1z
    else
       call sou_init_loop(sousetup%sem1z, sousetup%sem1)
    end if


    if(present(sem2)) then
       sousetup%sem2 = sem2
    else
       sousetup%sem2(1,1)=.5
       sousetup%sem2(1,2)=25.
    end if
    if(present(sem2z)) then
       sousetup%sem2z = sem2z
    else
       call sou_init_loop(sousetup%sem2z, sousetup%sem2)
    end if

    if(present(sthm1)) then
       sousetup%sthm1 = sthm1
    else
       sousetup%sthm1(1,1)=5.
       sousetup%sthm1(1,2)=0.
    end if
    if(present(sthm1z)) then
       sousetup%sthm1z = sthm1z
    else
       call sou_init_loop(sousetup%sthm1z, sousetup%sthm1)
    end if

    if(present(soucoord)) sousetup%soucoord = soucoord
    if(present(knockon)) sousetup%knockon = knockon
    if(present(komodel)) sousetup%komodel = komodel
    if(present(nkorfn)) sousetup%nkorfn = nkorfn
    if(present(nonko)) sousetup%nonko = nonko
    if(present(noffko)) sousetup%noffko = noffko
    if(present(soffvte)) sousetup%soffvte = soffvte
    if(present(soffpr)) sousetup%soffpr = soffpr
    if(present(isoucof)) sousetup%isoucof = isoucof
    if(present(faccof)) sousetup%faccof = faccof
    if(present(xlfac)) sousetup%xlfac = xlfac
    if(present(xlpctlwr)) sousetup%xlpctlwr = xlpctlwr
    if(present(xlpctmdl)) sousetup%xlpctmdl = xlpctmdl
    if(present(xllwr)) sousetup%xllwr = xllwr
    if(present(xlmdl)) sousetup%xlmdl = xlmdl

    if(present(jfl)) sousetup%jfl = jfl
    if (mod(jfl,2).eq.0) then
       print *, "WARNING, jfl needed to be odd because of jpxyh=(jfl+1)/2 in pltprppr"
       print *, "   Adjusting with sousetup%jfl=sousetup%jfl-1"
       sousetup%jfl=sousetup%jfl-1
    end if

    if ( present(debug_print)) then
       if (debug_print) call print_sousetup
    end if

  end subroutine set_sousetup

  subroutine print_sousetup()
    namelist /sousetup_nml/ sousetup
    WRITE(*, *) "!----  BEGIN SOUSETUP DUMP"
    WRITE(*, nml = sousetup_nml)
    WRITE(*, *)  "!----  END SOUSETUP DUMP"
  end subroutine print_sousetup

  subroutine get_setup_from_nml(nml_file, close_nml_file, debug_print)
    implicit none
    character(len=*), intent(in) :: nml_file
    logical, intent(in), optional :: close_nml_file
    logical, intent(in), optional :: debug_print
    ! local

    ! make private local variables to read in the namelist
    real(c_double) :: acoefne(4)
    real(c_double) :: acoefte(4)
    character(len=8) :: ampfmod
    integer :: nampfmax
    integer :: nonampf
    real(c_double) :: ampferr
    real(c_double) :: bctimescal
    real(c_double) :: bnumb(ntotala)
    real(c_double) :: bth
    real(c_double) :: btor
    character(len=8) :: bootst
    character(len=8) :: bootcalc
    character(len=8) :: bootupdt
    real(c_double) :: bootsign
    integer :: nonboot
    integer :: jhirsh
    real(c_double) :: contrmin
    real(c_double) :: constr
    character(len=8) :: chang
    integer :: colmodl
    real(c_double) :: deltabdb
    real(c_double) :: denpar(ntotala,0:lza+1)
    real(c_double) :: droptol
    real(c_double) :: dtr
    real(c_double) :: dtr1(ndtr1a)
    real(c_double) :: eegy(negyrga,2,ngena,lrza)
    real(c_double) :: eparc(ngena,0:lrza)
    real(c_double) :: eperc(ngena,0:lrza)
    real(c_double) :: simpbfac
    real(c_double) :: elecfld(0:lrza)
    real(c_double) :: elpar0
    real(c_double) :: enorm
    real(c_double) :: enorme
    real(c_double) :: enormi
    character(len=8) :: eleccomp
    real(c_double) :: elecin(njenea)
    real(c_double) :: elecin_t(njenea,nbctimea)
    real(c_double) :: elecscal
    real(c_double) :: enein(njenea,ntotala)
    real(c_double) :: enein_t(njenea,ntotala,nbctimea)
    real(c_double) :: ennin_t(njenea,nbctimea,npaproca)
    real(c_double) :: ennin(njenea,npaproca)
    real(c_double) :: enescal
    real(c_double) :: enloss(ngena)
    real(c_double) :: epsthet
    real(c_double) :: enmin
    real(c_double) :: enmax
    real(c_double) :: ennb(npaproca)
    real(c_double) :: ennl(npaproca)
    real(c_double) :: ennscal(npaproca)
    real(c_double) :: enmin_npa
    real(c_double) :: enmax_npa
    character(len=8) :: eseswtch
    real(c_double) :: xsink
    real(c_double) :: esink
    real(c_double) :: ephicc
    real(c_double) :: eoved
    real(c_double) :: fds
    real(c_double) :: fds_npa
    real(c_double) :: fmass(ntotala)
    character(len=8) :: f4d_out
    character(len=8) :: tavg
    real(c_double) :: tavg1(ntavga)
    real(c_double) :: tavg2(ntavga)
    real(c_double) :: gsla
    real(c_double) :: gslb
    real(c_double) :: gamaset
    character(len=8) :: gamafac
    real(c_double) :: gamegy(ngena)
    character(len=8) :: iactst
    character(len=8) :: ineg
    character(len=8) :: idskf
    character(len=8) :: idskrf
    character(len=8) :: ichkpnt
    character(len=8) :: implct
    character(len=8) :: iprone
    character(len=8) :: iprote
    character(len=8) :: iproti
    character(len=8) :: iprozeff
    character(len=8) :: iprovphi
    character(len=8) :: iproelec
    character(len=8) :: ipronn
    character(len=8) :: iprocur
    character(len=8) :: tmdmeth
    integer :: isigmas(6)
    integer :: isigtst
    integer :: isigsgv1
    integer :: isigsgv2
    real(c_double) :: pltflux1(7)
    integer :: irzplt(lrorsa)
    character(len=8) :: izeff
    integer :: iy
    integer :: jx
    integer :: kenorm
    integer :: lfil
    integer :: kfrsou
    character(len=8) :: kpress(ntotala)
    character(len=8) :: kfield(ntotala)
    character(len=8) :: kspeci(2,ntotala)
    real(c_double) :: fpld(10,ngena)
    integer :: lmidpln
    character(len=8) :: locquas
    character(len=8) :: lbdry(ngena)
    character(len=8) :: lbdry0
    character(len=256) :: lossfile(ngena)
    character(len=8) :: lossmode(ngena)
    integer :: lmidvel
    integer :: laddbnd
    integer :: lz
    character(len=8) :: machine
    character(len=8) :: meshy
    character(len=8) :: manymat
    character(len=256) :: netcdfnm
    character(len=8) :: netcdfshort
    character(len=8) :: netcdfvecal
    character(len=8) :: netcdfvecc
    character(len=8) :: netcdfvece
    character(len=8) :: netcdfvecrf
    character(len=8) :: netcdfvecs
    integer :: nnspec
    real(c_double) :: mpwr(0:ntotala) ! non trivial
    real(c_double) :: megy(ngena)
    real(c_double) :: mtorloss(ngena)
    integer :: mmsv
    integer :: msxr
    integer :: mx
    integer :: nchgdy
    integer :: ngauss
    integer :: nlagran
    logical :: nlotp1(noutpta)
    logical :: nlotp2(noutpta)
    logical :: nlotp3(noutpta)
    logical :: nlotp4(noutpta)
    integer :: nmax
    integer :: ngen
    integer :: nkconro(ntotala)
    integer :: nplt3d(nplota)
    integer :: nrskip
    integer :: nen
    integer :: nv
    integer :: nv_npa
    integer :: nen_npa
    integer :: npaproc
    character(len=8) :: npa_process(npaproca)
    integer :: nr_delta
    integer :: nz_delta
    integer :: nt_delta
    integer :: nr_f4d
    integer :: nz_f4d
    integer :: nv_f4d
    integer :: nt_f4d
    real(c_double) :: npwr(0:ntotala)
    real(c_double) :: negy(ngena)
    real(c_double) :: ntorloss(ngena)
    integer :: njene
    integer :: njte
    integer :: njti
    integer :: nstop
    integer :: nondtr1(ndtr1a)
    integer :: nplot(nplota)
    integer :: nsave(nsavea)
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
    integer :: nonavgf
    integer :: nofavgf
    integer :: nonloss
    integer :: noffloss
    integer :: nummods
    integer :: numixts
    integer :: numby
    integer :: negyrg
    character(len=8) :: oldiag
    character(len=8) :: plt3d
    character(len=8) :: pltfvs
    character(len=8) :: partner
    real(c_double) :: paregy(ngena)
    real(c_double) :: peregy(ngena)
    real(c_double) :: pegy(ngena)
    real(c_double) :: zeffin(0:njenea)
    real(c_double) :: zeffin_t(njenea,nbctimea)
    real(c_double) :: zeffscal
    real(c_double) :: vphiplin(0:njenea)
    real(c_double) :: vphiplin_t(njenea,nbctimea)
    real(c_double) :: vphiscal
    character(len=8) :: pltdn
    character(len=8) :: pltvecal
    character(len=8) :: pltvecc
    character(len=8) :: pltvecrf
    character(len=8) :: pltvece
    character(len=8) :: pltstrm
    character(len=8) :: pltflux
    real(c_double) :: pltmag
    character(len=8) :: pltsig
    character(len=8) :: pltlim
    real(c_double) :: pltlimm
    character(len=8) :: pltrst
    character(len=8) :: plturfb
    character(len=8) :: pltvflu
    character(len=8) :: pltra
    character(len=8) :: pltvs
    character(len=8) :: pltd
    character(len=8) :: pltprpp
    character(len=8) :: pltfofv
    character(len=8) :: pltlos
    character(len=8) :: profpsi
    character(len=8) :: psimodel
    character(len=8) :: pltpowe
    character(len=8) :: pltend
    character(len=8) :: pltinput
    character(len=8) :: pltrdc
    character(len=8) :: qsineut
    character(len=8) :: trapmod
    real(c_double) :: trapredc
    character(len=8) :: scatmod
    real(c_double) :: scatfrac
    real(c_double) :: ryain(njenea)
    real(c_double) :: radmaj
    real(c_double) :: radmin
    real(c_double) :: rmirror
    character(len=8) :: relativ
    real(c_double) :: reden(ntotala,0:lrza)
    character(len=8) :: regy(ngena)
    real(c_double) :: rfacz
    character(len=8) :: rzset
    real(c_double) :: rd(nva)
    real(c_double) :: roveram
    real(c_double) :: rovera(lrza)
    real(c_double) :: rya(0:lrza+1)
    character(len=8) :: radcoord
    character(len=8) :: sbdry
    character(len=8) :: scheck
    character(len=8) :: ndeltarho
    character(len=8) :: softxry
    character(len=8) :: npa_diag
    character(len=8) :: symtrap
    character(len=8) :: syncrad
    character(len=8) :: bremsrad
    real(c_double) :: brfac
    real(c_double) :: brfac1
    real(c_double) :: brfacgm3
    character(len=8) :: sigmamod
    real(c_double) :: sigvcx
    real(c_double) :: sigvi
    character(len=8) :: soln_method
    real(c_double) :: tauegy(ngena,0:lrza)
    character(len=8) :: taunew
    real(c_double) :: tescal
    real(c_double) :: tein(njenea)
    real(c_double) :: tein_t(njenea,nbctimea)
    real(c_double) :: tiscal
    real(c_double) :: tiin_t(njenea,nbctimea)
    real(c_double) :: tiin(njenea)
    real(c_double) :: tauloss(3,ngena)
    real(c_double) :: temp(ntotala,0:lrza)
    real(c_double) :: temppar(ntotala,0:lza+1)
    real(c_double) :: tfac
    real(c_double) :: tfacz
    real(c_double) :: tbnd(lrorsa)
    character(len=8) :: tandem
    real(c_double) :: thetd(nva)
    character(len=8) :: torloss(ngena)
    real(c_double) :: thet1(nva)
    real(c_double) :: thet2(nva)
    real(c_double) :: thet1_npa(nva)
    real(c_double) :: thet2_npa(nva)
    real(c_double) :: x_sxr(nva)
    real(c_double) :: z_sxr(nva)
    real(c_double) :: rd_npa(nva)
    real(c_double) :: thetd_npa(nva)
    real(c_double) :: x_npa(nva)
    real(c_double) :: z_npa(nva)
    character(len=8) :: atten_npa
    character(len=8) :: updown
    real(c_double) :: veclnth
    real(c_double) :: vnorm
    real(c_double) :: xfac
    real(c_double) :: xpctlwr
    real(c_double) :: xpctmdl
    real(c_double) :: xlwr
    real(c_double) :: xmdl
    real(c_double) :: xsinkm
    real(c_double) :: xprpmax
    integer :: ipxy,jpxy
    character(len=8) :: yreset
    real(c_double) :: ylower
    real(c_double) :: yupper
    real(c_double) :: mpwrzeff
    real(c_double) :: npwrzeff
    real(c_double) :: mpwrvphi
    real(c_double) :: npwrvphi
    real(c_double) :: mpwrxj
    real(c_double) :: npwrxj
    real(c_double) :: npwrelec
    real(c_double) :: mpwrelec
    real(c_double) :: redenc(nbctimea,ntotala)
    real(c_double) :: redenb(nbctimea,ntotala)
    real(c_double) :: temp_den
    real(c_double) :: tempc(nbctimea,ntotala)
    real(c_double) :: tempb(nbctimea,ntotala)
    real(c_double) :: zeffc(nbctimea)
    real(c_double) :: zeffb(nbctimea)
    real(c_double) :: elecc(nbctimea)
    real(c_double) :: elecb(nbctimea)
    real(c_double) :: vphic(nbctimea)
    real(c_double) :: vphib(nbctimea)
    real(c_double) :: xjc(nbctimea)
    real(c_double) :: xjb(nbctimea)
    real(c_double) :: xjin_t(njenea,nbctimea)
    real(c_double) :: totcrt(nbctimea)
    character(len=8) :: efswtch
    character(len=8) :: efswtchn
    character(len=8) :: efiter
    character(len=8) :: efflag
    real(c_double) :: curr_edge
    real(c_double) :: efrelax
    real(c_double) :: efrelax1
    real(c_double) :: currerr
    real(c_double) :: bctime(nbctimea)
    integer :: nbctime
    real(c_double) :: zmax(0:lrza)
    character(len=8) :: fow


    ! state the namelist, with associated vars
    namelist/setup/ &
         acoefne,acoefte, &
         ampfmod,nampfmax,nonampf,ampferr,bctimescal, &
         bnumb,btor,bth,bootst,bootcalc,bootupdt,bootsign,nonboot,jhirsh, &
         contrmin,constr,chang,colmodl, &
         deltabdb,denpar,droptol,dtr,dtr1, &
         eegy,eparc,eperc,simpbfac, &
         elecfld,elpar0,enorm,enorme,enormi,eleccomp, &
         elecin,elecin_t,elecscal,enein,enein_t,ennin_t, &
         enescal,enloss,epsthet, &
         enmin,enmax,ennb,ennin,ennl,ennscal,enmin_npa,enmax_npa, &
         eseswtch,xsink,esink,ephicc,eoved, &
         fds,fds_npa,fmass,f4d_out, &
         tavg,tavg1,tavg2, &
         gsla,gslb,gamaset,gamafac,gamegy, &
         iactst,ineg,idskf,idskrf,ichkpnt,implct, &
         iprone,iprote,iproti,iprozeff,iprovphi,iproelec,ipronn,iprocur, &
         tmdmeth,isigmas,isigtst,isigsgv1,isigsgv2, &
         pltflux1, &
         irzplt,izeff,&
         iy, &
         jx, &
         kenorm,lfil,kfrsou,kpress,kfield,kspeci,fpld, &
         lmidpln,locquas,lbdry,lbdry0,lossfile,lossmode,lmidvel,laddbnd, &
         lz, &
         machine,meshy,manymat,netcdfnm,netcdfshort, &
         netcdfvecal,netcdfvecc,netcdfvece,netcdfvecrf,netcdfvecs, &
         nnspec,mpwr,megy,mtorloss,mmsv,msxr,mx, &
         nchgdy,ngauss,nlagran, &
         nlotp1,nlotp2,nlotp3,nlotp4, &
         nmax,ngen,nkconro,nplt3d,nrskip,nen,nv,nen_npa,nv_npa, &
         npaproc,npa_process, &
         nr_delta,nz_delta,nt_delta, &
         nr_f4d,nz_f4d,nv_f4d,nt_f4d, &
         npwr,negy,ntorloss,njene,njte,njti, &
         nstop,nondtr1,nplot,nsave,ncoef,nchec,ncont,nrstrt,nstps,nfpld, &
         noncntrl,nonel,noffel,nonvphi,noffvphi,nonavgf,nofavgf, &
         nonloss,noffloss,nummods,numixts, &
         numby,negyrg, &
         oldiag, &
         plt3d,pltvs,partner,paregy,peregy,pegy, &
         zeffin,zeffin_t,zeffscal,vphiplin,vphiplin_t,vphiscal, &
         pltdn,pltvecal,pltvecc,pltvecrf,pltvece, &
         pltstrm,pltflux,pltmag,pltsig,pltlim,pltlimm, &
         pltrst,plturfb,pltvflu,pltra,pltfvs,pltd,pltprpp,pltfofv,pltlos, &
         profpsi, &
         psimodel,pltpowe,pltend,pltinput, pltrdc, &
         qsineut,trapmod,trapredc,scatmod,scatfrac, &
         ryain,radmaj,radmin,rmirror,relativ, &
         reden,regy,rfacz,rzset,rd,roveram, &
         rovera,rya,radcoord, &
         sbdry,scheck,ndeltarho,softxry,npa_diag,symtrap,syncrad, &
         bremsrad,brfac,brfac1,brfacgm3,sigmamod,sigvcx,sigvi, &
         soln_method,tauegy,taunew,tein,tein_t,tescal,tiin,tiin_t,tiscal, &
         tauloss,temp,temppar, &
         tfac,tfacz,tbnd,tandem, &
         thetd,torloss,thet1,thet2,x_sxr,z_sxr, &
         rd_npa,thetd_npa,x_npa,z_npa,thet1_npa,thet2_npa,atten_npa, &
         updown, &
         veclnth,vnorm, &
         xfac,xpctlwr,xpctmdl,xlwr,xmdl,xsinkm, &
         xprpmax,ipxy,jpxy, &
         yreset,ylower,yupper, &
         mpwrzeff,npwrzeff,mpwrvphi,npwrvphi,mpwrxj,npwrxj, &
         mpwrelec,npwrelec, &
         redenc,redenb,temp_den,tempc,tempb,zeffc,zeffb,elecc,elecb, &
         vphic,vphib,xjc,xjb,xjin_t,totcrt,efswtch,efswtchn, &
         efiter,efflag,curr_edge,efrelax,efrelax1,currerr, &
         bctime,nbctime, &
         zmax, &
         fow

    ! copy defaults to local vars
    acoefne = setup_default%acoefne
    acoefte = setup_default%acoefte
    ampfmod = setup_default%ampfmod
    nampfmax = setup_default%nampfmax
    nonampf = setup_default%nonampf
    ampferr = setup_default%ampferr
    bctimescal = setup_default%bctimescal
    bnumb = setup_default%bnumb
    btor = setup_default%btor
    bth = setup_default%bth
    bootst = setup_default%bootst
    bootcalc = setup_default%bootcalc
    bootupdt = setup_default%bootupdt
    bootsign = setup_default%bootsign
    nonboot = setup_default%nonboot
    jhirsh = setup_default%jhirsh
    contrmin = setup_default%contrmin
    constr = setup_default%constr
    chang = setup_default%chang
    colmodl = setup_default%colmodl
    deltabdb = setup_default%deltabdb
    denpar = setup_default%denpar
    droptol = setup_default%droptol
    dtr = setup_default%dtr
    dtr1 = setup_default%dtr1
    eegy = setup_default%eegy
    eparc = setup_default%eparc
    eperc = setup_default%eperc
    simpbfac = setup_default%simpbfac
    elecfld = setup_default%elecfld
    elpar0 = setup_default%elpar0
    enorm = setup_default%enorm
    enorme = setup_default%enorme
    enormi = setup_default%enormi
    eleccomp = setup_default%eleccomp
    elecin = setup_default%elecin
    elecin_t = setup_default%elecin_t
    elecscal = setup_default%elecscal
    enein = setup_default%enein
    enein_t = setup_default%enein_t
    ennin_t = setup_default%ennin_t
    enescal = setup_default%enescal
    enloss = setup_default%enloss
    epsthet = setup_default%epsthet
    enmin = setup_default%enmin
    enmax = setup_default%enmax
    ennb = setup_default%ennb
    ennin = setup_default%ennin
    ennl = setup_default%ennl
    ennscal = setup_default%ennscal
    enmin_npa = setup_default%enmin_npa
    enmax_npa = setup_default%enmax_npa
    eseswtch = setup_default%eseswtch
    xsink = setup_default%xsink
    esink = setup_default%esink
    ephicc = setup_default%ephicc
    eoved = setup_default%eoved
    fds = setup_default%fds
    fds_npa = setup_default%fds_npa
    fmass = setup_default%fmass
    f4d_out = setup_default%f4d_out
    tavg = setup_default%tavg
    tavg1 = setup_default%tavg1
    tavg2 = setup_default%tavg2
    gsla = setup_default%gsla
    gslb = setup_default%gslb
    gamaset = setup_default%gamaset
    gamafac = setup_default%gamafac
    gamegy = setup_default%gamegy
    iactst = setup_default%iactst
    ineg = setup_default%ineg
    idskf = setup_default%idskf
    idskrf = setup_default%idskrf
    ichkpnt = setup_default%ichkpnt
    implct = setup_default%implct
    iprone = setup_default%iprone
    iprote = setup_default%iprote
    iproti = setup_default%iproti
    iprozeff = setup_default%iprozeff
    iprovphi = setup_default%iprovphi
    iproelec = setup_default%iproelec
    ipronn = setup_default%ipronn
    iprocur = setup_default%iprocur
    tmdmeth = setup_default%tmdmeth
    isigmas = setup_default%isigmas
    isigtst = setup_default%isigtst
    isigsgv1 = setup_default%isigsgv1
    isigsgv2 = setup_default%isigsgv2
    pltflux1 = setup_default%pltflux1
    irzplt = setup_default%irzplt
    izeff = setup_default%izeff
    iy = setup_default%iy
    jx = setup_default%jx
    kenorm = setup_default%kenorm
    lfil = setup_default%lfil
    kfrsou = setup_default%kfrsou
    kpress = setup_default%kpress
    kfield = setup_default%kfield
    kspeci = setup_default%kspeci
    fpld = setup_default%fpld
    lmidpln = setup_default%lmidpln
    locquas = setup_default%locquas
    lbdry = setup_default%lbdry
    lbdry0 = setup_default%lbdry0
    lossfile = setup_default%lossfile
    lossmode = setup_default%lossmode
    lmidvel = setup_default%lmidvel
    laddbnd = setup_default%laddbnd
    lz = setup_default%lz
    machine = setup_default%machine
    meshy = setup_default%meshy
    manymat = setup_default%manymat
    netcdfnm = setup_default%netcdfnm
    netcdfshort = setup_default%netcdfshort
    netcdfvecal = setup_default%netcdfvecal
    netcdfvecc = setup_default%netcdfvecc
    netcdfvece = setup_default%netcdfvece
    netcdfvecrf = setup_default%netcdfvecrf
    netcdfvecs = setup_default%netcdfvecs
    nnspec = setup_default%nnspec
    mpwr = setup_default%mpwr
    megy = setup_default%megy
    mtorloss = setup_default%mtorloss
    mmsv = setup_default%mmsv
    msxr = setup_default%msxr
    mx = setup_default%mx
    nchgdy = setup_default%nchgdy
    ngauss = setup_default%ngauss
    nlagran = setup_default%nlagran
    nlotp1 = setup_default%nlotp1
    nlotp2 = setup_default%nlotp2
    nlotp3 = setup_default%nlotp3
    nlotp4 = setup_default%nlotp4
    nmax = setup_default%nmax
    ngen = setup_default%ngen
    nkconro = setup_default%nkconro
    nplt3d = setup_default%nplt3d
    nrskip = setup_default%nrskip
    nen = setup_default%nen
    nv = setup_default%nv
    nen_npa = setup_default%nen_npa
    nv_npa = setup_default%nv_npa
    npaproc = setup_default%npaproc
    npa_process = setup_default%npa_process
    nr_delta = setup_default%nr_delta
    nz_delta = setup_default%nz_delta
    nt_delta = setup_default%nt_delta
    nr_f4d = setup_default%nr_f4d
    nz_f4d = setup_default%nz_f4d
    nv_f4d = setup_default%nv_f4d
    nt_f4d = setup_default%nt_f4d
    npwr = setup_default%npwr
    negy = setup_default%negy
    ntorloss = setup_default%ntorloss
    njene = setup_default%njene
    njte = setup_default%njte
    njti = setup_default%njti
    nstop = setup_default%nstop
    nondtr1 = setup_default%nondtr1
    nplot = setup_default%nplot
    nsave = setup_default%nsave
    ncoef = setup_default%ncoef
    nchec = setup_default%nchec
    ncont = setup_default%ncont
    nrstrt = setup_default%nrstrt
    nstps = setup_default%nstps
    nfpld = setup_default%nfpld
    noncntrl = setup_default%noncntrl
    nonel = setup_default%nonel
    noffel = setup_default%noffel
    nonvphi = setup_default%nonvphi
    noffvphi = setup_default%noffvphi
    nonavgf = setup_default%nonavgf
    nofavgf = setup_default%nofavgf
    nonloss = setup_default%nonloss
    noffloss = setup_default%noffloss
    nummods = setup_default%nummods
    numixts = setup_default%numixts
    numby = setup_default%numby
    negyrg = setup_default%negyrg
    oldiag = setup_default%oldiag
    plt3d = setup_default%plt3d
    pltvs = setup_default%pltvs
    partner = setup_default%partner
    paregy = setup_default%paregy
    peregy = setup_default%peregy
    pegy = setup_default%pegy
    zeffin = setup_default%zeffin
    zeffin_t = setup_default%zeffin_t
    zeffscal = setup_default%zeffscal
    vphiplin = setup_default%vphiplin
    vphiplin_t = setup_default%vphiplin_t
    vphiscal = setup_default%vphiscal
    pltdn = setup_default%pltdn
    pltvecal = setup_default%pltvecal
    pltvecc = setup_default%pltvecc
    pltvecrf = setup_default%pltvecrf
    pltvece = setup_default%pltvece
    pltstrm = setup_default%pltstrm
    pltflux = setup_default%pltflux
    pltmag = setup_default%pltmag
    pltsig = setup_default%pltsig
    pltlim = setup_default%pltlim
    pltlimm = setup_default%pltlimm
    pltrst = setup_default%pltrst
    plturfb = setup_default%plturfb
    pltvflu = setup_default%pltvflu
    pltra = setup_default%pltra
    pltfvs = setup_default%pltfvs
    pltd = setup_default%pltd
    pltprpp = setup_default%pltprpp
    pltfofv = setup_default%pltfofv
    pltlos = setup_default%pltlos
    profpsi = setup_default%profpsi
    psimodel = setup_default%psimodel
    pltpowe = setup_default%pltpowe
    pltend = setup_default%pltend
    pltinput = setup_default%pltinput
    pltrdc = setup_default%pltrdc
    qsineut = setup_default%qsineut
    trapmod = setup_default%trapmod
    trapredc = setup_default%trapredc
    scatmod = setup_default%scatmod
    scatfrac = setup_default%scatfrac
    ryain = setup_default%ryain
    radmaj = setup_default%radmaj
    radmin = setup_default%radmin
    rmirror = setup_default%rmirror
    relativ = setup_default%relativ
    reden = setup_default%reden
    regy = setup_default%regy
    rfacz = setup_default%rfacz
    rzset = setup_default%rzset
    rd = setup_default%rd
    roveram = setup_default%roveram
    rovera = setup_default%rovera
    rya = setup_default%rya
    radcoord = setup_default%radcoord
    sbdry = setup_default%sbdry
    scheck = setup_default%scheck
    ndeltarho = setup_default%ndeltarho
    softxry = setup_default%softxry
    npa_diag = setup_default%npa_diag
    symtrap = setup_default%symtrap
    syncrad = setup_default%syncrad
    bremsrad = setup_default%bremsrad
    brfac = setup_default%brfac
    brfac1 = setup_default%brfac1
    brfacgm3 = setup_default%brfacgm3
    sigmamod = setup_default%sigmamod
    sigvcx = setup_default%sigvcx
    sigvi = setup_default%sigvi
    soln_method = setup_default%soln_method
    tauegy = setup_default%tauegy
    taunew = setup_default%taunew
    tein = setup_default%tein
    tein_t = setup_default%tein_t
    tescal = setup_default%tescal
    tiin = setup_default%tiin
    tiin_t = setup_default%tiin_t
    tiscal = setup_default%tiscal
    tauloss = setup_default%tauloss
    temp = setup_default%temp
    temppar = setup_default%temppar
    tfac = setup_default%tfac
    tfacz = setup_default%tfacz
    tbnd = setup_default%tbnd
    tandem = setup_default%tandem
    thetd = setup_default%thetd
    torloss = setup_default%torloss
    thet1 = setup_default%thet1
    thet2 = setup_default%thet2
    x_sxr = setup_default%x_sxr
    z_sxr = setup_default%z_sxr
    rd_npa = setup_default%rd_npa
    thetd_npa = setup_default%thetd_npa
    x_npa = setup_default%x_npa
    z_npa = setup_default%z_npa
    thet1_npa = setup_default%thet1_npa
    thet2_npa = setup_default%thet2_npa
    atten_npa = setup_default%atten_npa
    updown = setup_default%updown
    veclnth = setup_default%veclnth
    vnorm = setup_default%vnorm
    xfac = setup_default%xfac
    xpctlwr = setup_default%xpctlwr
    xpctmdl = setup_default%xpctmdl
    xlwr = setup_default%xlwr
    xmdl = setup_default%xmdl
    xsinkm = setup_default%xsinkm
    xprpmax = setup_default%xprpmax
    ipxy = setup_default%ipxy
    jpxy = setup_default%jpxy
    yreset = setup_default%yreset
    ylower = setup_default%ylower
    yupper = setup_default%yupper
    mpwrzeff = setup_default%mpwrzeff
    npwrzeff = setup_default%npwrzeff
    mpwrvphi = setup_default%mpwrvphi
    npwrvphi = setup_default%npwrvphi
    mpwrxj = setup_default%mpwrxj
    npwrxj = setup_default%npwrxj
    mpwrelec = setup_default%mpwrelec
    npwrelec = setup_default%npwrelec
    redenc = setup_default%redenc
    redenb = setup_default%redenb
    temp_den = setup_default%temp_den
    tempc = setup_default%tempc
    tempb = setup_default%tempb
    zeffc = setup_default%zeffc
    zeffb = setup_default%zeffb
    elecc = setup_default%elecc
    elecb = setup_default%elecb
    vphic = setup_default%vphic
    vphib = setup_default%vphib
    xjc = setup_default%xjc
    xjb = setup_default%xjb
    xjin_t = setup_default%xjin_t
    totcrt = setup_default%totcrt
    efswtch = setup_default%efswtch
    efswtchn = setup_default%efswtchn
    efiter = setup_default%efiter
    efflag = setup_default%efflag
    curr_edge = setup_default%curr_edge
    efrelax = setup_default%efrelax
    efrelax1 = setup_default%efrelax1
    currerr = setup_default%currerr
    bctime = setup_default%bctime
    nbctime = setup_default%nbctime
    zmax = setup_default%zmax
    fow = setup_default%fow

    ! read the nml, which will write into the local vars

    call maybe_nml_open(nml_file)
    read(nml_fd, setup)

    ! external codes can call this, which packs the setup0 derived type.
    call set_setup(        acoefne,acoefte, &
         ampfmod,nampfmax,nonampf,ampferr,bctimescal, &
         bnumb,btor,bth,bootst,bootcalc,bootupdt,bootsign,nonboot,jhirsh, &
         contrmin,constr,chang,colmodl, &
         deltabdb,denpar,droptol,dtr,dtr1, &
         eegy,eparc,eperc,simpbfac, &
         elecfld,elpar0,enorm,enorme,enormi,eleccomp, &
         elecin,elecin_t,elecscal,enein,enein_t,ennin_t, &
         enescal,enloss,epsthet, &
         enmin,enmax,ennb,ennin,ennl,ennscal,enmin_npa,enmax_npa, &
         eseswtch,xsink,esink,ephicc,eoved, &
         fds,fds_npa,fmass,f4d_out, &
         tavg,tavg1,tavg2, &
         gsla,gslb,gamaset,gamafac,gamegy, &
         iactst,ineg,idskf,idskrf,ichkpnt,implct, &
         iprone,iprote,iproti,iprozeff,iprovphi,iproelec,ipronn,iprocur, &
         tmdmeth,isigmas,isigtst,isigsgv1,isigsgv2, &
         pltflux1, &
         irzplt,izeff,&
         iy, &
         jx, &
         kenorm,lfil,kfrsou,kpress,kfield,kspeci,fpld, &
         lmidpln,locquas,lbdry,lbdry0,lossfile,lossmode,lmidvel,laddbnd, &
         lz, &
         machine,meshy,manymat,netcdfnm,netcdfshort, &
         netcdfvecal,netcdfvecc,netcdfvece,netcdfvecrf,netcdfvecs, &
         nnspec,mpwr,megy,mtorloss,mmsv,msxr,mx, &
         nchgdy,ngauss,nlagran, &
         nlotp1,nlotp2,nlotp3,nlotp4, &
         nmax,ngen,nkconro,nplt3d,nrskip,nen,nv,nen_npa,nv_npa, &
         npaproc,npa_process, &
         nr_delta,nz_delta,nt_delta, &
         nr_f4d,nz_f4d,nv_f4d,nt_f4d, &
         npwr,negy,ntorloss,njene,njte,njti, &
         nstop,nondtr1,nplot,nsave,ncoef,nchec,ncont,nrstrt,nstps,nfpld, &
         noncntrl,nonel,noffel,nonvphi,noffvphi,nonavgf,nofavgf, &
         nonloss,noffloss,nummods,numixts, &
         numby,negyrg, &
         oldiag, &
         plt3d,pltvs,partner,paregy,peregy,pegy, &
         zeffin,zeffin_t,zeffscal,vphiplin,vphiplin_t,vphiscal, &
         pltdn,pltvecal,pltvecc,pltvecrf,pltvece, &
         pltstrm,pltflux,pltmag,pltsig,pltlim,pltlimm, &
         pltrst,plturfb,pltvflu,pltra,pltfvs,pltd,pltprpp,pltfofv,pltlos, &
         profpsi, &
         psimodel,pltpowe,pltend,pltinput, pltrdc, &
         qsineut,trapmod,trapredc,scatmod,scatfrac, &
         ryain,radmaj,radmin,rmirror,relativ, &
         reden,regy,rfacz,rzset,rd,roveram, &
         rovera,rya,radcoord, &
         sbdry,scheck,ndeltarho,softxry,npa_diag,symtrap,syncrad, &
         bremsrad,brfac,brfac1,brfacgm3,sigmamod,sigvcx,sigvi, &
         soln_method,tauegy,taunew,tein,tein_t,tescal,tiin,tiin_t,tiscal, &
         tauloss,temp,temppar, &
         tfac,tfacz,tbnd,tandem, &
         thetd,torloss,thet1,thet2,x_sxr,z_sxr, &
         rd_npa,thetd_npa,x_npa,z_npa,thet1_npa,thet2_npa,atten_npa, &
         updown, &
         veclnth,vnorm, &
         xfac,xpctlwr,xpctmdl,xlwr,xmdl,xsinkm, &
         xprpmax,ipxy,jpxy, &
         yreset,ylower,yupper, &
         mpwrzeff,npwrzeff,mpwrvphi,npwrvphi,mpwrxj,npwrxj, &
         mpwrelec,npwrelec, &
         redenc,redenb,temp_den,tempc,tempb,zeffc,zeffb,elecc,elecb, &
         vphic,vphib,xjc,xjb,xjin_t,totcrt,efswtch,efswtchn, &
         efiter,efflag,curr_edge,efrelax,efrelax1,currerr, &
         bctime,nbctime, &
         zmax, &
         fow, &
         debug_print)

    ! we optionally close the nml file.
    if (present(close_nml_file)) then
       if(close_nml_file) then
          call nml_close()
       end if
    endif

  end subroutine get_setup_from_nml

  subroutine set_setup(        acoefne,acoefte, &
       ampfmod,nampfmax,nonampf,ampferr,bctimescal, &
       bnumb,btor,bth,bootst,bootcalc,bootupdt,bootsign,nonboot,jhirsh, &
       contrmin,constr,chang,colmodl, &
       deltabdb,denpar,droptol,dtr,dtr1, &
       eegy,eparc,eperc,simpbfac,&
       elecfld,elpar0,enorm,enorme,enormi,eleccomp, &
       elecin,elecin_t,elecscal,enein,enein_t,ennin_t, &
       enescal,enloss,epsthet, &
       enmin,enmax,ennb,ennin,ennl,ennscal,enmin_npa,enmax_npa, &
       eseswtch,xsink,esink,ephicc,eoved, &
       fds,fds_npa,fmass,f4d_out, &
       tavg,tavg1,tavg2, &
       gsla,gslb,gamaset,gamafac,gamegy, &
       iactst,ineg,idskf,idskrf,ichkpnt,implct, &
       iprone,iprote,iproti,iprozeff,iprovphi,iproelec,ipronn,iprocur, &
       tmdmeth,isigmas,isigtst,isigsgv1,isigsgv2, &
       pltflux1, &
       irzplt,izeff,&
       iy, &
       jx, &
       kenorm,lfil,kfrsou,kpress,kfield,kspeci,fpld, &
       lmidpln,locquas,lbdry,lbdry0,lossfile,lossmode,lmidvel,laddbnd, &
       lz, &
       machine,meshy,manymat,netcdfnm,netcdfshort, &
       netcdfvecal,netcdfvecc,netcdfvece,netcdfvecrf,netcdfvecs, &
       nnspec,mpwr,megy,mtorloss,mmsv,msxr,mx, &
       nchgdy,ngauss,nlagran, &
       nlotp1,nlotp2,nlotp3,nlotp4, &
       nmax,ngen,nkconro,nplt3d,nrskip,nen,nv,nen_npa,nv_npa, &
       npaproc,npa_process, &
       nr_delta,nz_delta,nt_delta, &
       nr_f4d,nz_f4d,nv_f4d,nt_f4d, &
       npwr,negy,ntorloss,njene,njte,njti, &
       nstop,nondtr1,nplot,nsave,ncoef,nchec,ncont,nrstrt,nstps,nfpld, &
       noncntrl,nonel,noffel,nonvphi,noffvphi,nonavgf,nofavgf, &
       nonloss,noffloss,nummods,numixts, &
       numby,negyrg, &
       oldiag, &
       plt3d,pltvs,partner,paregy,peregy,pegy, &
       zeffin,zeffin_t,zeffscal,vphiplin,vphiplin_t,vphiscal, &
       pltdn,pltvecal,pltvecc,pltvecrf,pltvece, &
       pltstrm,pltflux,pltmag,pltsig,pltlim,pltlimm, &
       pltrst,plturfb,pltvflu,pltra,pltfvs,pltd,pltprpp,pltfofv,pltlos, &
       profpsi, &
       psimodel,pltpowe,pltend,pltinput, pltrdc, &
       qsineut,trapmod,trapredc,scatmod,scatfrac, &
       ryain,radmaj,radmin,rmirror,relativ, &
       reden,regy,rfacz,rzset,rd,roveram, &
       rovera,rya,radcoord, &
       sbdry,scheck,ndeltarho,softxry,npa_diag,symtrap,syncrad, &
       bremsrad,brfac,brfac1,brfacgm3,sigmamod,sigvcx,sigvi, &
       soln_method,tauegy,taunew,tein,tein_t,tescal,tiin,tiin_t,tiscal, &
       tauloss,temp,temppar, &
       tfac,tfacz,tbnd,tandem, &
       thetd,torloss,thet1,thet2,x_sxr,z_sxr, &
       rd_npa,thetd_npa,x_npa,z_npa,thet1_npa,thet2_npa,atten_npa, &
       updown, &
       veclnth,vnorm, &
       xfac,xpctlwr,xpctmdl,xlwr,xmdl,xsinkm, &
       xprpmax,ipxy,jpxy, &
       yreset,ylower,yupper, &
       mpwrzeff,npwrzeff,mpwrvphi,npwrvphi,mpwrxj,npwrxj, &
       mpwrelec,npwrelec, &
       redenc,redenb,temp_den,tempc,tempb,zeffc,zeffb,elecc,elecb, &
       vphic,vphib,xjc,xjb,xjin_t,totcrt,efswtch,efswtchn, &
       efiter,efflag,curr_edge,efrelax,efrelax1,currerr, &
       bctime,nbctime, &
       zmax, &
       fow, &
       debug_print)
    logical, intent(in), optional :: debug_print
    !
    real(c_double), intent(in), optional :: acoefne(4)
    real(c_double), intent(in), optional :: acoefte(4)
    character(len=8), intent(in), optional :: ampfmod
    integer, intent(in), optional :: nampfmax
    integer, intent(in), optional :: nonampf
    real(c_double), intent(in), optional :: ampferr
    real(c_double), intent(in), optional :: bctimescal
    real(c_double), intent(in), optional :: bnumb(ntotala)
    real(c_double), intent(in), optional :: bth
    real(c_double), intent(in), optional :: btor
    character(len=8), intent(in), optional :: bootst
    character(len=8), intent(in), optional :: bootcalc
    character(len=8), intent(in), optional :: bootupdt
    real(c_double), intent(in), optional :: bootsign
    integer, intent(in), optional :: nonboot
    integer, intent(in), optional :: jhirsh
    real(c_double), intent(in), optional :: contrmin
    real(c_double), intent(in), optional :: constr
    character(len=8), intent(in), optional :: chang
    integer, intent(in), optional :: colmodl
    real(c_double), intent(in), optional :: deltabdb
    real(c_double), intent(in), optional :: denpar(ntotala,0:lza+1)
    real(c_double), intent(in), optional :: droptol
    real(c_double), intent(in), optional :: dtr
    real(c_double), intent(in), optional :: dtr1(ndtr1a)
    real(c_double), intent(in), optional :: eegy(negyrga,2,ngena,lrza)
    real(c_double), intent(in), optional :: eparc(ngena,0:lrza)
    real(c_double), intent(in), optional :: eperc(ngena,0:lrza)
    real(c_double), intent(in), optional :: simpbfac
    real(c_double), intent(in), optional :: elecfld(0:lrza)
    real(c_double), intent(in), optional :: elpar0
    real(c_double), intent(in), optional :: enorm
    real(c_double), intent(in), optional :: enorme
    real(c_double), intent(in), optional :: enormi
    character(len=8), intent(in), optional :: eleccomp
    real(c_double), intent(in), optional :: elecin(njenea)
    real(c_double), intent(in), optional :: elecin_t(njenea,nbctimea)
    real(c_double), intent(in), optional :: elecscal
    real(c_double), intent(in), optional :: enein(njenea,ntotala)
    real(c_double), intent(in), optional :: enein_t(njenea,ntotala,nbctimea)
    real(c_double), intent(in), optional :: ennin_t(njenea,nbctimea,npaproca)
    real(c_double), intent(in), optional :: ennin(njenea,npaproca)
    real(c_double), intent(in), optional :: enescal
    real(c_double), intent(in), optional :: enloss(ngena)
    real(c_double), intent(in), optional :: epsthet
    real(c_double), intent(in), optional :: enmin
    real(c_double), intent(in), optional :: enmax
    real(c_double), intent(in), optional :: ennb(npaproca)
    real(c_double), intent(in), optional :: ennl(npaproca)
    real(c_double), intent(in), optional :: ennscal(npaproca)
    real(c_double), intent(in), optional :: enmin_npa
    real(c_double), intent(in), optional :: enmax_npa
    character(len=8), intent(in), optional :: eseswtch
    real(c_double), intent(in), optional :: xsink
    real(c_double), intent(in), optional :: esink
    real(c_double), intent(in), optional :: ephicc
    real(c_double), intent(in), optional :: eoved
    real(c_double), intent(in), optional :: fds
    real(c_double), intent(in), optional :: fds_npa
    real(c_double), intent(in), optional :: fmass(ntotala)
    character(len=8), intent(in), optional :: f4d_out
    character(len=8), intent(in), optional :: tavg
    real(c_double), intent(in), optional :: tavg1(ntavga)
    real(c_double), intent(in), optional :: tavg2(ntavga)
    real(c_double), intent(in), optional :: gsla
    real(c_double), intent(in), optional :: gslb
    real(c_double), intent(in), optional :: gamaset
    character(len=8), intent(in), optional :: gamafac
    real(c_double), intent(in), optional :: gamegy(ngena)
    character(len=8), intent(in), optional :: iactst
    character(len=8), intent(in), optional :: ineg
    character(len=8), intent(in), optional :: idskf
    character(len=8), intent(in), optional :: idskrf
    character(len=8), intent(in), optional :: ichkpnt
    character(len=8), intent(in), optional :: implct
    character(len=8), intent(in), optional :: iprone
    character(len=8), intent(in), optional :: iprote
    character(len=8), intent(in), optional :: iproti
    character(len=8), intent(in), optional :: iprozeff
    character(len=8), intent(in), optional :: iprovphi
    character(len=8), intent(in), optional :: iproelec
    character(len=8), intent(in), optional :: ipronn
    character(len=8), intent(in), optional :: iprocur
    character(len=8), intent(in), optional :: tmdmeth
    integer, intent(in), optional :: isigmas(6)
    integer, intent(in), optional :: isigtst
    integer, intent(in), optional :: isigsgv1
    integer, intent(in), optional :: isigsgv2
    real(c_double), intent(in), optional :: pltflux1(7)
    integer, intent(in), optional :: irzplt(lrorsa)
    character(len=8), intent(in), optional :: izeff
    integer, intent(in), optional :: iy
    integer, intent(in), optional :: jx
    integer, intent(in), optional :: kenorm
    integer, intent(in), optional :: lfil
    integer, intent(in), optional :: kfrsou
    character(len=8), intent(in), optional :: kpress(ntotala)
    character(len=8), intent(in), optional :: kfield(ntotala)
    character(len=8), intent(in), optional :: kspeci(2,ntotala)
    real(c_double), intent(in), optional :: fpld(10,ngena)
    integer, intent(in), optional :: lmidpln
    character(len=8), intent(in), optional :: locquas
    character(len=8), intent(in), optional :: lbdry(ngena)
    character(len=8), intent(in), optional :: lbdry0
    character(len=256), intent(in), optional :: lossfile(ngena)
    character(len=8), intent(in), optional :: lossmode(ngena)
    integer, intent(in), optional :: lmidvel
    integer, intent(in), optional :: laddbnd
    integer, intent(in), optional :: lz
    character(len=8), intent(in), optional :: machine
    character(len=8), intent(in), optional :: meshy
    character(len=8), intent(in), optional :: manymat
    character(len=256), intent(in), optional :: netcdfnm
    character(len=8), intent(in), optional :: netcdfshort
    character(len=8), intent(in), optional :: netcdfvecal
    character(len=8), intent(in), optional :: netcdfvecc
    character(len=8), intent(in), optional :: netcdfvece
    character(len=8), intent(in), optional :: netcdfvecrf
    character(len=8), intent(in), optional :: netcdfvecs
    integer, intent(in), optional :: nnspec
    real(c_double), intent(in), optional :: mpwr(0:ntotala) ! non trivial
    real(c_double), intent(in), optional :: megy(ngena)
    real(c_double), intent(in), optional :: mtorloss(ngena)
    integer, intent(in), optional :: mmsv
    integer, intent(in), optional :: msxr
    integer, intent(in), optional :: mx
    integer, intent(in), optional :: nchgdy
    integer, intent(in), optional :: ngauss
    integer, intent(in), optional :: nlagran
    logical, intent(in), optional :: nlotp1(noutpta)
    logical, intent(in), optional :: nlotp2(noutpta)
    logical, intent(in), optional :: nlotp3(noutpta)
    logical, intent(in), optional :: nlotp4(noutpta)
    integer, intent(in), optional :: nmax
    integer, intent(in), optional :: ngen
    integer, intent(in), optional :: nkconro(ntotala)
    integer, intent(in), optional :: nplt3d(nplota)
    integer, intent(in), optional :: nrskip
    integer, intent(in), optional :: nen
    integer, intent(in), optional :: nv
    integer, intent(in), optional :: nv_npa
    integer, intent(in), optional :: nen_npa
    integer, intent(in), optional :: npaproc
    character(len=8), intent(in), optional :: npa_process(npaproca)
    integer, intent(in), optional :: nr_delta
    integer, intent(in), optional :: nz_delta
    integer, intent(in), optional :: nt_delta
    integer, intent(in), optional :: nr_f4d
    integer, intent(in), optional :: nz_f4d
    integer, intent(in), optional :: nv_f4d
    integer, intent(in), optional :: nt_f4d
    real(c_double), intent(in), optional :: npwr(0:ntotala)
    real(c_double), intent(in), optional :: negy(ngena)
    real(c_double), intent(in), optional :: ntorloss(ngena)
    integer, intent(in), optional :: njene
    integer, intent(in), optional :: njte
    integer, intent(in), optional :: njti
    integer, intent(in), optional :: nstop
    integer, intent(in), optional :: nondtr1(ndtr1a)
    integer, intent(in), optional :: nplot(nplota)
    integer, intent(in), optional :: nsave(nsavea)
    integer, intent(in), optional :: ncoef
    integer, intent(in), optional :: nchec
    integer, intent(in), optional :: ncont
    integer, intent(in), optional :: nrstrt
    integer, intent(in), optional :: nstps
    integer, intent(in), optional :: nfpld
    integer, intent(in), optional :: noncntrl
    integer, intent(in), optional :: nonel
    integer, intent(in), optional :: noffel
    integer, intent(in), optional :: nonvphi
    integer, intent(in), optional :: noffvphi
    integer, intent(in), optional :: nonavgf
    integer, intent(in), optional :: nofavgf
    integer, intent(in), optional :: nonloss
    integer, intent(in), optional :: noffloss
    integer, intent(in), optional :: nummods
    integer, intent(in), optional :: numixts
    integer, intent(in), optional :: numby
    integer, intent(in), optional :: negyrg
    character(len=8), intent(in), optional :: oldiag
    character(len=8), intent(in), optional :: plt3d
    character(len=8), intent(in), optional :: pltfvs
    character(len=8), intent(in), optional :: partner
    real(c_double), intent(in), optional :: paregy(ngena)
    real(c_double), intent(in), optional :: peregy(ngena)
    real(c_double), intent(in), optional :: pegy(ngena)
    real(c_double), intent(in), optional :: zeffin(0:njenea)
    real(c_double), intent(in), optional :: zeffin_t(njenea,nbctimea)
    real(c_double), intent(in), optional :: zeffscal
    real(c_double), intent(in), optional :: vphiplin(0:njenea)
    real(c_double), intent(in), optional :: vphiplin_t(njenea,nbctimea)
    real(c_double), intent(in), optional :: vphiscal
    character(len=8), intent(in), optional :: pltdn
    character(len=8), intent(in), optional :: pltvecal
    character(len=8), intent(in), optional :: pltvecc
    character(len=8), intent(in), optional :: pltvecrf
    character(len=8), intent(in), optional :: pltvece
    character(len=8), intent(in), optional :: pltstrm
    character(len=8), intent(in), optional :: pltflux
    real(c_double), intent(in), optional :: pltmag
    character(len=8), intent(in), optional :: pltsig
    character(len=8), intent(in), optional :: pltlim
    real(c_double), intent(in), optional :: pltlimm
    character(len=8), intent(in), optional :: pltrst
    character(len=8), intent(in), optional :: plturfb
    character(len=8), intent(in), optional :: pltvflu
    character(len=8), intent(in), optional :: pltra
    character(len=8), intent(in), optional :: pltvs
    character(len=8), intent(in), optional :: pltd
    character(len=8), intent(in), optional :: pltprpp
    character(len=8), intent(in), optional :: pltfofv
    character(len=8), intent(in), optional :: pltlos
    character(len=8), intent(in), optional :: profpsi
    character(len=8), intent(in), optional :: psimodel
    character(len=8), intent(in), optional :: pltpowe
    character(len=8), intent(in), optional :: pltend
    character(len=8), intent(in), optional :: pltinput
    character(len=8), intent(in), optional :: pltrdc
    character(len=8), intent(in), optional :: qsineut
    character(len=8), intent(in), optional :: trapmod
    real(c_double), intent(in), optional :: trapredc
    character(len=8), intent(in), optional :: scatmod
    real(c_double), intent(in), optional :: scatfrac
    real(c_double), intent(in), optional :: ryain(njenea)
    real(c_double), intent(in), optional :: radmaj
    real(c_double), intent(in), optional :: radmin
    real(c_double), intent(in), optional :: rmirror
    character(len=8), intent(in), optional :: relativ
    real(c_double), intent(in), optional :: reden(ntotala,0:lrza)
    character(len=8), intent(in), optional :: regy(ngena)
    real(c_double), intent(in), optional :: rfacz
    character(len=8), intent(in), optional :: rzset
    real(c_double), intent(in), optional :: rd(nva)
    real(c_double), intent(in), optional :: roveram
    real(c_double), intent(in), optional :: rovera(lrza)
    real(c_double), intent(in), optional :: rya(0:lrza+1)
    character(len=8), intent(in), optional :: radcoord
    character(len=8), intent(in), optional :: sbdry
    character(len=8), intent(in), optional :: scheck
    character(len=8), intent(in), optional :: ndeltarho
    character(len=8), intent(in), optional :: softxry
    character(len=8), intent(in), optional :: npa_diag
    character(len=8), intent(in), optional :: symtrap
    character(len=8), intent(in), optional :: syncrad
    character(len=8), intent(in), optional :: bremsrad
    real(c_double), intent(in), optional :: brfac
    real(c_double), intent(in), optional :: brfac1
    real(c_double), intent(in), optional :: brfacgm3
    character(len=8), intent(in), optional :: sigmamod
    real(c_double), intent(in), optional :: sigvcx
    real(c_double), intent(in), optional :: sigvi
    character(len=8), intent(in), optional :: soln_method
    real(c_double), intent(in), optional :: tauegy(ngena,0:lrza)
    character(len=8), intent(in), optional :: taunew
    real(c_double), intent(in), optional :: tescal
    real(c_double), intent(in), optional :: tein(njenea)
    real(c_double), intent(in), optional :: tein_t(njenea,nbctimea)
    real(c_double), intent(in), optional :: tiscal
    real(c_double), intent(in), optional :: tiin_t(njenea,nbctimea)
    real(c_double), intent(in), optional :: tiin(njenea)
    real(c_double), intent(in), optional :: tauloss(3,ngena)
    real(c_double), intent(in), optional :: temp(ntotala,0:lrza)
    real(c_double), intent(in), optional :: temppar(ntotala,0:lza+1)
    real(c_double), intent(in), optional :: tfac
    real(c_double), intent(in), optional :: tfacz
    real(c_double), intent(in), optional :: tbnd(lrorsa)
    character(len=8), intent(in), optional :: tandem
    real(c_double), intent(in), optional :: thetd(nva)
    character(len=8), intent(in), optional :: torloss(ngena)
    real(c_double), intent(in), optional :: thet1(nva)
    real(c_double), intent(in), optional :: thet2(nva)
    real(c_double), intent(in), optional :: thet1_npa(nva)
    real(c_double), intent(in), optional :: thet2_npa(nva)
    real(c_double), intent(in), optional :: x_sxr(nva)
    real(c_double), intent(in), optional :: z_sxr(nva)
    real(c_double), intent(in), optional :: rd_npa(nva)
    real(c_double), intent(in), optional :: thetd_npa(nva)
    real(c_double), intent(in), optional :: x_npa(nva)
    real(c_double), intent(in), optional :: z_npa(nva)
    character(len=8), intent(in), optional :: atten_npa
    character(len=8), intent(in), optional :: updown
    real(c_double), intent(in), optional :: veclnth
    real(c_double), intent(in), optional :: vnorm
    real(c_double), intent(in), optional :: xfac
    real(c_double), intent(in), optional :: xpctlwr
    real(c_double), intent(in), optional :: xpctmdl
    real(c_double), intent(in), optional :: xlwr
    real(c_double), intent(in), optional :: xmdl
    real(c_double), intent(in), optional :: xsinkm
    real(c_double), intent(in), optional :: xprpmax
    integer, intent(in), optional :: ipxy
    integer, intent(in), optional :: jpxy
    character(len=8), intent(in), optional :: yreset
    real(c_double), intent(in), optional :: ylower
    real(c_double), intent(in), optional :: yupper
    real(c_double), intent(in), optional :: mpwrzeff
    real(c_double), intent(in), optional :: npwrzeff
    real(c_double), intent(in), optional :: mpwrvphi
    real(c_double), intent(in), optional :: npwrvphi
    real(c_double), intent(in), optional :: mpwrxj
    real(c_double), intent(in), optional :: npwrxj
    real(c_double), intent(in), optional :: npwrelec
    real(c_double), intent(in), optional :: mpwrelec
    real(c_double), intent(in), optional :: redenc(nbctimea,ntotala)
    real(c_double), intent(in), optional :: redenb(nbctimea,ntotala)
    real(c_double), intent(in), optional :: temp_den
    real(c_double), intent(in), optional :: tempc(nbctimea,ntotala)
    real(c_double), intent(in), optional :: tempb(nbctimea,ntotala)
    real(c_double), intent(in), optional :: zeffc(nbctimea)
    real(c_double), intent(in), optional :: zeffb(nbctimea)
    real(c_double), intent(in), optional :: elecc(nbctimea)
    real(c_double), intent(in), optional :: elecb(nbctimea)
    real(c_double), intent(in), optional :: vphic(nbctimea)
    real(c_double), intent(in), optional :: vphib(nbctimea)
    real(c_double), intent(in), optional :: xjc(nbctimea)
    real(c_double), intent(in), optional :: xjb(nbctimea)
    real(c_double), intent(in), optional :: xjin_t(njenea,nbctimea)
    real(c_double), intent(in), optional :: totcrt(nbctimea)
    character(len=8), intent(in), optional :: efswtch
    character(len=8), intent(in), optional :: efswtchn
    character(len=8), intent(in), optional :: efiter
    character(len=8), intent(in), optional :: efflag
    real(c_double), intent(in), optional :: curr_edge
    real(c_double), intent(in), optional :: efrelax
    real(c_double), intent(in), optional :: efrelax1
    real(c_double), intent(in), optional :: currerr
    real(c_double), intent(in), optional :: bctime(nbctimea)
    integer, intent(in), optional :: nbctime
    real(c_double), intent(in), optional :: zmax(0:lrza)
    character(len=8), intent(in), optional :: fow
    ! local scratch
    integer :: i

    ! All this code should do is override the defaults
    ! in setup0 with optional args, or setup non trivial defaults

    if(present(acoefne)) setup%acoefne = acoefne
    if(present(acoefte)) setup%acoefte = acoefte
    if(present(ampfmod)) setup%ampfmod = ampfmod
    if(present(nampfmax)) setup%nampfmax = nampfmax
    if(present(nonampf)) setup%nonampf = nonampf
    if(present(ampferr)) setup%ampferr = ampferr
    if(present(bctimescal)) setup%bctimescal = bctimescal
    if(present(bnumb)) setup%bnumb = bnumb
    if(present(btor)) setup%btor = btor
    if(present(bth)) setup%bth = bth
    if(present(bootst)) setup%bootst = bootst
    if(present(bootcalc)) setup%bootcalc = bootcalc
    if(present(bootupdt)) setup%bootupdt = bootupdt
    if(present(bootsign)) setup%bootsign = bootsign
    if(present(nonboot)) setup%nonboot = nonboot
    if(present(jhirsh)) setup%jhirsh = jhirsh
    if(present(contrmin)) setup%contrmin = contrmin
    if(present(constr)) setup%constr = constr
    if(present(chang)) setup%chang = chang
    if(present(colmodl)) setup%colmodl = colmodl
    if(present(deltabdb)) setup%deltabdb = deltabdb
    if(present(denpar)) setup%denpar = denpar
    if(present(droptol)) setup%droptol = droptol
    if(present(dtr)) setup%dtr = dtr
    if(present(dtr1)) setup%dtr1 = dtr1
    if(present(eegy)) setup%eegy = eegy
    if(present(eparc)) setup%eparc = eparc
    if(present(eperc)) setup%eperc = eperc
    if(present(simpbfac)) setup%simpbfac = simpbfac
    if(present(elecfld)) setup%elecfld = elecfld
    if(present(elpar0)) setup%elpar0 = elpar0
    if(present(enorm)) setup%enorm = enorm
    if(present(enorme) .and. setup%enorme .ne. enorme) then
       setup%enorme = enorme
    else
       setup%enorme = setup%enorm   ! (re)set to enom at runtime, sort of dangerous
    end if
    if(present(enormi) .and. setup%enormi .ne. enormi) then
       setup%enormi = enormi
    else
       setup%enormi = setup%enorm  ! (re)set to enom at runtime, sort of dangerous
    end if
    if(present(eleccomp)) setup%eleccomp = eleccomp
    if(present(elecin)) setup%elecin = elecin
    if(present(elecin_t)) setup%elecin_t = elecin_t
    if(present(elecscal)) setup%elecscal = elecscal
    if(present(enein)) setup%enein = enein
    if(present(enein_t)) setup%enein_t = enein_t
    if(present(ennin_t)) setup%ennin_t = ennin_t
    if(present(enescal)) setup%enescal = enescal
    if(present(enloss)) setup%enloss = enloss
    if(present(epsthet)) setup%epsthet = epsthet
    if(present(enmin)) setup%enmin = enmin
    if(present(enmax)) setup%enmax = enmax
    if(present(ennb)) setup%ennb = ennb
    if(present(ennin)) setup%ennin = ennin
    if(present(ennl)) setup%ennl = ennl
    if(present(ennscal)) setup%ennscal = ennscal
    if(present(enmin_npa)) setup%enmin_npa = enmin_npa
    if(present(enmax_npa)) setup%enmax_npa = enmax_npa
    if(present(eseswtch)) setup%eseswtch = eseswtch
    if(present(xsink)) setup%xsink = xsink
    if(present(esink)) setup%esink = esink
    if(present(ephicc)) setup%ephicc = ephicc
    if(present(eoved)) setup%eoved = eoved
    if(present(fds)) setup%fds = fds
    if(present(fds_npa)) setup%fds_npa = fds_npa
    if(present(fmass)) setup%fmass = fmass
    if(present(f4d_out)) setup%f4d_out = f4d_out
    if(present(tavg)) setup%tavg = tavg
    if(present(tavg1)) setup%tavg1 = tavg1
    if(present(tavg2)) setup%tavg2 = tavg2
    if(present(gsla)) setup%gsla = gsla
    if(present(gslb)) setup%gslb = gslb
    if(present(gamaset)) setup%gamaset = gamaset
    if(present(gamafac)) setup%gamafac = gamafac
    if(present(gamegy)) setup%gamegy = gamegy
    if(present(iactst)) setup%iactst = iactst
    if(present(ineg)) setup%ineg = ineg
    if(present(idskf)) setup%idskf = idskf
    if(present(idskrf)) setup%idskrf = idskrf
    if(present(ichkpnt)) setup%ichkpnt = ichkpnt
    if(present(implct)) setup%implct = implct
    if(present(iprone)) setup%iprone = iprone
    if(present(iprote)) setup%iprote = iprote
    if(present(iproti)) setup%iproti = iproti
    if(present(iprozeff)) setup%iprozeff = iprozeff
    if(present(iprovphi)) setup%iprovphi = iprovphi
    if(present(iproelec)) setup%iproelec = iproelec
    if(present(ipronn)) setup%ipronn = ipronn
    if(present(iprocur)) setup%iprocur = iprocur
    if(present(tmdmeth)) setup%tmdmeth = tmdmeth
    if(present(isigmas)) setup%isigmas = isigmas
    if(present(isigtst)) setup%isigtst = isigtst
    if(present(isigsgv1)) setup%isigsgv1 = isigsgv1
    if(present(isigsgv2)) setup%isigsgv2 = isigsgv2
    if(present(pltflux1)) setup%pltflux1 = pltflux1
    if(present(irzplt)) setup%irzplt = irzplt
    if(present(izeff)) setup%izeff = izeff
    if(present(iy)) setup%iy = iy
    if(present(jx)) setup%jx = jx
    if(present(kenorm)) setup%kenorm = kenorm
    if(present(lfil)) setup%lfil = lfil
    if(present(kfrsou)) setup%kfrsou = kfrsou
    if(present(kpress)) setup%kpress = kpress
    if(present(kfield)) setup%kfield = kfield
    if(present(kspeci)) setup%kspeci = kspeci
    if(present(fpld)) setup%fpld = fpld
    if(present(lmidpln)) setup%lmidpln = lmidpln
    if(present(locquas)) setup%locquas = locquas
    if(present(lbdry)) setup%lbdry = lbdry
    if(present(lbdry0)) setup%lbdry0 = lbdry0
    if(present(lossfile)) setup%lossfile = lossfile
    if(present(lossmode)) setup%lossmode = lossmode
    if(present(lmidvel)) setup%lmidvel = lmidvel
    if(present(laddbnd)) setup%laddbnd = laddbnd
    if(present(lz)) setup%lz = lz
    if(present(machine)) setup%machine = machine
    if(present(meshy)) setup%meshy = meshy
    if(present(manymat)) setup%manymat = manymat
    if(present(netcdfnm)) setup%netcdfnm = netcdfnm
    if(present(netcdfshort)) setup%netcdfshort = netcdfshort
    if(present(netcdfvecal)) setup%netcdfvecal = netcdfvecal
    if(present(netcdfvecc)) setup%netcdfvecc = netcdfvecc
    if(present(netcdfvece)) setup%netcdfvece = netcdfvece
    if(present(netcdfvecrf)) setup%netcdfvecrf = netcdfvecrf
    if(present(netcdfvecs)) setup%netcdfvecs = netcdfvecs
    if(present(nnspec)) setup%nnspec = nnspec
    if(present(mpwr)) setup%mpwr = mpwr
    if(present(megy)) setup%megy = megy
    if(present(mtorloss)) setup%mtorloss = mtorloss
    if(present(mx)) setup%mx = mx
    if(present(mmsv) .and. mmsv .ne. setup_default%mmsv) then
       setup%mmsv = mmsv
    else
       setup%mmsv = setup%mx ! xxx reset by mx at run time, seems dangerous
    end if
    if(present(msxr) .and. msxr .ne. setup_default%msxr) then
       setup%msxr = msxr
    else
       setup%msxr = setup%mx ! xxx reset by mx at run time, seems dangerous
    end if
    if(present(nchgdy)) setup%nchgdy = nchgdy
    if(present(ngauss)) setup%ngauss = ngauss
    if(present(nlagran)) setup%nlagran = nlagran
    if(present(nlotp1)) setup%nlotp1 = nlotp1
    if(present(nlotp2)) setup%nlotp2 = nlotp2
    if(present(nlotp3)) setup%nlotp3 = nlotp3
    if(present(nlotp4)) setup%nlotp4 = nlotp4
    if(present(nmax)) setup%nmax = nmax
    if(present(ngen)) setup%ngen = ngen
    if(present(nkconro)) setup%nkconro = nkconro
    if(present(nplt3d)) setup%nplt3d = nplt3d
    if(present(nrskip)) setup%nrskip = nrskip
    if(present(nen)) setup%nen = nen
    if(present(nv)) setup%nv = nv
    if(present(nen_npa)) setup%nen_npa = nen_npa
    if(present(nv_npa)) setup%nv_npa = nv_npa
    if(present(npaproc)) setup%npaproc = npaproc
    if(present(npa_process)) setup%npa_process = npa_process
    if(present(nr_delta)) setup%nr_delta = nr_delta
    if(present(nz_delta)) setup%nz_delta = nz_delta
    if(present(nt_delta)) setup%nt_delta = nt_delta
    if(present(nr_f4d)) setup%nr_f4d = nr_f4d
    if(present(nz_f4d)) setup%nz_f4d = nz_f4d
    if(present(nv_f4d)) setup%nv_f4d = nv_f4d
    if(present(nt_f4d)) setup%nt_f4d = nt_f4d
    if(present(npwr)) setup%npwr = npwr
    if(present(negy)) setup%negy = negy
    if(present(ntorloss)) setup%ntorloss = ntorloss
    if(present(njene)) setup%njene = njene
    if(present(njte) .and. njte .ne. setup_default%njte) then
       setup%njte = njte   ! reset to njene at runtime, seems dangerous
    else
       setup%njte = setup%njene
    end if
    if(present(njti) .and. njti .ne. setup_default%njti) then
       setup%njti = njti  ! reset to njene at runtime, seems dangerous
    else
       setup%njti = setup%njene
    end if
    if(present(nstop)) setup%nstop = nstop
    if(present(nondtr1)) setup%nondtr1 = nondtr1
    if(present(nplot)) setup%nplot = nplot
    if(present(nsave)) setup%nsave = nsave
    if(present(ncoef)) setup%ncoef = ncoef
    if(present(nchec)) setup%nchec = nchec
    if(present(ncont)) setup%ncont = ncont
    if(present(nrstrt)) setup%nrstrt = nrstrt
    if(present(nstps)) setup%nstps = nstps
    if(present(nfpld)) setup%nfpld = nfpld
    if(present(noncntrl)) setup%noncntrl = noncntrl
    if(present(nonel)) setup%nonel = nonel
    if(present(noffel)) setup%noffel = noffel
    if(present(nonvphi)) setup%nonvphi = nonvphi
    if(present(noffvphi)) setup%noffvphi = noffvphi
    if(present(nonavgf)) setup%nonavgf = nonavgf
    if(present(nofavgf)) setup%nofavgf = nofavgf
    if(present(nonloss)) setup%nonloss = nonloss
    if(present(noffloss)) setup%noffloss = noffloss
    if(present(nummods)) setup%nummods = nummods
    if(present(numixts)) setup%numixts = numixts
    if(present(numby)) setup%numby = numby
    if(present(negyrg)) setup%negyrg = negyrg
    if(present(oldiag)) setup%oldiag = oldiag
    if(present(plt3d)) setup%plt3d = plt3d
    if(present(pltvs)) setup%pltvs = pltvs
    if(present(partner)) setup%partner = partner
    if(present(paregy)) setup%paregy = paregy
    if(present(peregy)) setup%peregy = peregy
    if(present(pegy)) setup%pegy = pegy
    if(present(zeffin)) setup%zeffin = zeffin
    if(present(zeffin_t)) setup%zeffin_t = zeffin_t
    if(present(zeffscal)) setup%zeffscal = zeffscal
    if(present(vphiplin)) setup%vphiplin = vphiplin
    if(present(vphiplin_t)) setup%vphiplin_t = vphiplin_t
    if(present(vphiscal)) setup%vphiscal = vphiscal
    if(present(pltdn)) setup%pltdn = pltdn
    if(present(pltvecal)) setup%pltvecal = pltvecal
    if(present(pltvecc)) setup%pltvecc = pltvecc
    if(present(pltvecrf)) setup%pltvecrf = pltvecrf
    if(present(pltvece)) setup%pltvece = pltvece
    if(present(pltstrm)) setup%pltstrm = pltstrm
    if(present(pltflux)) setup%pltflux = pltflux
    if(present(pltmag)) setup%pltmag = pltmag
    if(present(pltsig)) setup%pltsig = pltsig
    if(present(pltlim)) setup%pltlim = pltlim
    if(present(pltlimm)) setup%pltlimm = pltlimm
    if(present(pltrst)) setup%pltrst = pltrst
    if(present(plturfb)) setup%plturfb = plturfb
    if(present(pltvflu)) setup%pltvflu = pltvflu
    if(present(pltra)) setup%pltra = pltra
    if(present(pltfvs)) setup%pltfvs = pltfvs
    if(present(pltd)) setup%pltd = pltd
    if(present(pltprpp)) setup%pltprpp = pltprpp
    if(present(pltfofv)) setup%pltfofv = pltfofv
    if(present(pltlos)) setup%pltlos = pltlos
    if(present(profpsi)) setup%profpsi = profpsi
    if(present(psimodel)) setup%psimodel = psimodel
    if(present(pltpowe)) setup%pltpowe = pltpowe
    if(present(pltend)) setup%pltend = pltend
    if(present(pltinput)) setup%pltinput = pltinput
    if(present(pltrdc)) setup%pltrdc = pltrdc
    if(present(qsineut)) setup%qsineut = qsineut
    if(present(trapmod)) setup%trapmod = trapmod
    if(present(trapredc)) setup%trapredc = trapredc
    if(present(scatmod)) setup%scatmod = scatmod
    if(present(scatfrac)) setup%scatfrac = scatfrac
    if(present(ryain)) setup%ryain = ryain
    if(present(radmaj)) setup%radmaj = radmaj
    if(present(radmin)) setup%radmin = radmin
    if(present(rmirror)) setup%rmirror = rmirror
    if(present(relativ)) setup%relativ = relativ
    if(present(reden)) setup%reden = reden
    if(present(regy)) setup%regy = regy
    if(present(rfacz)) setup%rfacz = rfacz
    if(present(rzset)) setup%rzset = rzset
    if(present(rd)) setup%rd = rd
    if(present(roveram)) setup%roveram = roveram
    if(present(rovera)) setup%rovera = rovera
    if(present(rya)) setup%rya = rya
    if(present(radcoord)) setup%radcoord = radcoord
    if(present(sbdry)) setup%sbdry = sbdry
    if(present(scheck)) setup%scheck = scheck
    if(present(ndeltarho)) setup%ndeltarho = ndeltarho
    if(present(softxry)) setup%softxry = softxry
    if(present(npa_diag)) setup%npa_diag = npa_diag
    if(present(symtrap)) setup%symtrap = symtrap
    if(present(syncrad)) setup%syncrad = syncrad
    if(present(bremsrad)) setup%bremsrad = bremsrad
    if(present(brfac)) setup%brfac = brfac
    if(present(brfac1)) setup%brfac1 = brfac1
    if(present(brfacgm3)) setup%brfacgm3 = brfacgm3
    if(present(sigmamod)) setup%sigmamod = sigmamod
    if(present(sigvcx)) setup%sigvcx = sigvcx
    if(present(sigvi)) setup%sigvi = sigvi
    if(present(soln_method)) setup%soln_method = soln_method
    if(present(tauegy)) setup%tauegy = tauegy
    if(present(taunew)) setup%taunew = taunew
    if(present(tein)) setup%tein = tein
    if(present(tein_t)) setup%tein_t = tein_t
    if(present(tescal)) setup%tescal = tescal
    if(present(tiin)) setup%tiin = tiin
    if(present(tiin_t)) setup%tiin_t = tiin_t
    if(present(tiscal)) setup%tiscal = tiscal
    if(present(tauloss)) setup%tauloss = tauloss
    if(present(temp)) setup%temp = temp
    if(present(temppar)) setup%temppar = temppar
    if(present(tfac)) setup%tfac = tfac
    if(present(tfacz)) setup%tfacz = tfacz
    if(present(tbnd)) setup%tbnd = tbnd
    if(present(tandem)) setup%tandem = tandem
    if(present(thetd)) setup%thetd = thetd
    if(present(torloss)) setup%torloss = torloss
    if(present(thet1)) setup%thet1 = thet1
    if(present(thet2)) setup%thet2 = thet2
    if(present(x_sxr)) setup%x_sxr = x_sxr
    if(present(z_sxr)) setup%z_sxr = z_sxr
    if(present(rd_npa)) setup%rd_npa = rd_npa
    if(present(thetd_npa)) setup%thetd_npa = thetd_npa
    if(present(x_npa)) setup%x_npa = x_npa
    if(present(z_npa)) setup%z_npa = z_npa
    if(present(thet1_npa)) setup%thet1_npa = thet1_npa
    if(present(thet2_npa)) setup%thet2_npa = thet2_npa
    if(present(atten_npa)) setup%atten_npa = atten_npa
    if(present(updown)) setup%updown = updown
    if(present(veclnth)) setup%veclnth = veclnth
    if(present(vnorm)) setup%vnorm = vnorm
    if(present(xfac)) setup%xfac = xfac
    if(present(xpctlwr)) setup%xpctlwr = xpctlwr
    if(present(xpctmdl)) setup%xpctmdl = xpctmdl
    if(present(xlwr)) setup%xlwr = xlwr
    if(present(xmdl)) setup%xmdl = xmdl
    if(present(xsinkm)) setup%xsinkm = xsinkm
    if(present(xprpmax)) setup%xprpmax = xprpmax
    if(present(ipxy)) setup%ipxy = ipxy
    if(present(jpxy)) setup%jpxy = jpxy
    if(present(yreset)) setup%yreset = yreset
    if(present(ylower)) setup%ylower = ylower
    if(present(yupper)) setup%yupper = yupper
    if(present(mpwrzeff)) setup%mpwrzeff = mpwrzeff
    if(present(npwrzeff)) setup%npwrzeff = npwrzeff
    if(present(mpwrvphi)) setup%mpwrvphi = mpwrvphi
    if(present(npwrvphi)) setup%npwrvphi = npwrvphi
    if(present(mpwrxj)) setup%mpwrxj = mpwrxj
    if(present(npwrxj)) setup%npwrxj = npwrxj
    if(present(mpwrelec)) setup%mpwrelec = mpwrelec
    if(present(npwrelec)) setup%npwrelec = npwrelec
    if(present(redenc)) setup%redenc = redenc
    if(present(redenb)) setup%redenb = redenb
    if(present(temp_den)) setup%temp_den = temp_den
    if(present(tempc)) setup%tempc = tempc
    if(present(tempb)) setup%tempb = tempb
    if(present(zeffc)) setup%zeffc = zeffc
    if(present(zeffb)) setup%zeffb = zeffb
    if(present(elecc)) setup%elecc = elecc
    if(present(elecb)) setup%elecb = elecb
    if(present(vphic)) setup%vphic = vphic
    if(present(vphib)) setup%vphib = vphib
    if(present(xjc)) setup%xjc = xjc
    if(present(xjb)) setup%xjb = xjb
    if(present(xjin_t)) setup%xjin_t = xjin_t
    if(present(totcrt)) setup%totcrt = totcrt
    if(present(efswtch)) setup%efswtch = efswtch
    if(present(efswtchn)) setup%efswtchn = efswtchn
    if(present(efiter)) setup%efiter = efiter
    if(present(efflag)) setup%efflag = efflag
    if(present(curr_edge)) setup%curr_edge = curr_edge
    if(present(efrelax)) setup%efrelax = efrelax
    if(present(efrelax1)) setup%efrelax1 = efrelax1
    if(present(currerr)) setup%currerr = currerr
    if(present(bctime)) setup%bctime = bctime
    if(present(nbctime)) setup%nbctime = nbctime
    if(present(zmax)) setup%zmax = zmax
    if(present(fow)) setup%fow = fow

    if ( present(debug_print)) then
       if (debug_print) call print_setup
    end if

  end subroutine set_setup

  subroutine print_setup()
    namelist /setup_nml/ setup
    WRITE(*, *) "!----  BEGIN SETUP DUMP"
    WRITE(*, nml = setup_nml)
    WRITE(*, *)  "!----  END SETUP DUMP"
  end subroutine print_setup

  subroutine print_all_conf_nml
    call print_setup0
    call print_setup
    call print_trsetup
    call print_sousetup
    call print_eqsetup
    call print_rfsetup
  end subroutine print_all_conf_nml

  ! function default_fpld(dim1, dim2)
  !   integer, intent(in) :: dim1, dim2
  !   real(c_double), dimension(dim1, dim2) :: default_fpld
  !   default_fpld(1:6,:)=0.
  !   default_fpld(7,:)=0.
  !   default_fpld(8,:)=1.e10
  !   default_fpld(9,:)=0.
  !   default_fpld(10,:)=pi
  ! end function default_fpld


end module cqlconf_mod
