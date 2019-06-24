module cqlconf_mod
  use param_mod, only : lrorsa, lfielda
  use iso_c_binding, only : c_double
  implicit none

  private
  integer :: ll
  logical, save :: nml_file_open = .FALSE.
  integer, save:: nml_fd = -1
  public nml_close

  public get_setup0_from_nml
  public print_setup0
  public set_setup0

  public get_eqsetup_from_nml
  public print_eqsetup
  public set_eqsetup

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

  ! rest of cql3d will access this
  type(setup0_t), public, save :: setup0
  type(eqsetup_t), public, target, save :: eqsetup

  !..................................................................
  !     NAMELIST (SETUP0) DECLARATION FOR INPUT
  !..................................................................

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
    logical, intent(in), optional :: debug_print

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

end module cqlconf_mod
