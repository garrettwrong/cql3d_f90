module cqlconf_mod
  use param_mod, only : lrorsa, lfielda
  use param_mod, only : ep100, nmodsa, ngena, nrdca, nbctimea
  use param_mod, only : njenea
  use iso_c_binding, only : c_double, c_double_complex
  implicit none

  private
  integer :: ll
  logical, save :: nml_file_open = .FALSE.
  integer, save :: nml_fd = -1
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
     character(len=256) :: rdcfile(1:nrdca) = "notset"
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
     integer :: nrdcspecies(1:nrdca) = 0
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


  ! rest of cql3d will access this
  type(setup0_t), public, save :: setup0
  type(eqsetup_t), public, target, save :: eqsetup
  type(rfsetup_t), public, target, save :: rfsetup
  type(trsetup_t), public, target, save :: trsetup

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
    if(present(rdcfile)) then
       rfsetup%rdcfile = rdcfile
    else
       ! set non trivial default first value:
       rfsetup%rdcfile(1) = "du0u0_input"
    end if
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
    if(present(nrdcspecies)) then
       rfsetup%nrdcspecies = nrdcspecies
    else
       ! set non trivial default first value:
       rfsetup%nrdcspecies(1) = 1
    end if
    if(present(rdcscale)) rfsetup%rdcscale = rdcscale
    if(present(rdc_netcdf)) rfsetup%rdc_netcdf = rdc_netcdf

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

  end subroutine set_trsetup

  subroutine print_trsetup()
    namelist /trsetup_nml/ trsetup
    WRITE(*, *) "!----  BEGIN TRSETUP DUMP"
    WRITE(*, nml = trsetup_nml)
    WRITE(*, *)  "!----  END TRSETUP DUMP"
  end subroutine print_trsetup

  end module cqlconf_mod
