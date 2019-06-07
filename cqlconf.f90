module cqlconf_mod
  use param_mod, only : lrorsa
  implicit none

  private
  logical :: nml_file_open = .FALSE.
  integer, save:: nml_fd = -1
  public get_setup0_from_nml
  public set_setup0
  public nml_close

  type setup0_t
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
  end type setup0_t

  ! rest of cql3d will access this
  type(setup0_t), target, public :: setup0
  type(setup0_t), pointer, private :: setup0_=>setup0! hack
  !..................................................................
  !     NAMELIST (SETUP0) DECLARATION FOR INPUT
  !..................................................................

contains

  subroutine maybe_nml_open(nml_file)
    character(len=*), intent(in) :: nml_file
    if (nml_file_open .eqv. .FALSE.) then
       nml_fd = newunit()
       open(unit=nml_fd, file=nml_file, status='old')
       nml_file_open = .TRUE.
    end if
  end subroutine maybe_nml_open

  subroutine nml_close()
    if (nml_file_open .eqv. .TRUE.) then
       close(nml_fd)
       nml_file_open = .FALSE.
    end if
  end subroutine nml_close

  subroutine get_setup0_from_nml(nml_file, close_nml_file)
    implicit none
    character(len=*), intent(in) :: nml_file
    logical, intent(in), optional :: close_nml_file

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

    namelist/setup0/ mnemonic,ioutput,iuser,ibox,noplots,lnwidth, &
         nmlstout,special_calls,cqlpmod,lrz,lrzdiff,lrzmax,lrindx, &
         ls,lsmax,lsdiff,lsindx,nlrestrt,nlwritf

    ! for now, we take thse from setup0,
    ! because there is existing defaulting code earlier.
    ! Ideally that code would live with the config setters and getters,
    ! and this goes away.
    ! To my knowledge fortran does not provide a better facility
    mnemonic = setup0_%mnemonic
    ioutput(2) = setup0_%ioutput(2)
    iuser = setup0_%iuser
    ibox(3) = setup0_%ibox(3)
    noplots = setup0_%noplots
    lnwidth = setup0_%lnwidth
    nmlstout = setup0_%nmlstout
    special_calls = setup0_%special_calls
    cqlpmod = setup0_%cqlpmod
    lrz = setup0_%lrz
    lrzdiff = setup0_%lrzdiff
    lrzmax = setup0_%lrzmax
    lrindx(0:lrorsa) = setup0_%lrindx(0:lrorsa)
    ls = setup0_%ls
    lsmax = setup0_%lsmax
    lsdiff = setup0_%lsdiff
    lsindx(0:lrorsa) = setup0_%lsindx(0:lrorsa)
    nlrestrt = setup0_%nlrestrt
    nlwritf = setup0_%nlwritf


    call maybe_nml_open(nml_file)
    read(nml_fd, setup0)


    ! external codes can call this
    call set_setup0(mnemonic,ioutput,iuser,ibox,noplots,lnwidth, &
         nmlstout,special_calls,cqlpmod,lrz,lrzdiff,lrzmax,lrindx, &
         ls,lsmax,lsdiff,lsindx,nlrestrt,nlwritf)

    if (present(close_nml_file)) then
       if(close_nml_file .eqv. .TRUE.) then
          call nml_close()
       end if
    endif

  end subroutine get_setup0_from_nml

  subroutine set_setup0(mnemonic,ioutput,iuser,ibox,noplots,lnwidth, &
         nmlstout,special_calls,cqlpmod,lrz,lrzdiff,lrzmax,lrindx, &
         ls,lsmax,lsdiff,lsindx,nlrestrt,nlwritf)
    character(len=256), optional :: mnemonic
    integer, optional :: ioutput(2)
    character(len=8), optional :: iuser
    character(len=8), optional :: ibox(3)
    character(len=8), optional :: noplots
    integer, optional :: lnwidth
    character(len=8), optional :: nmlstout
    character(len=8), optional :: special_calls
    character(len=8), optional :: cqlpmod
    integer, optional :: lrz
    character(len=8), optional :: lrzdiff
    integer, optional :: lrzmax
    integer, optional :: lrindx(0:lrorsa)
    integer, optional :: ls
    integer, optional :: lsmax
    character(len=8), optional :: lsdiff
    integer, optional :: lsindx(0:lrorsa)
    character(len=8), optional :: nlrestrt
    character(len=8), optional :: nlwritf

    ! Note, for all of these, we should collect old cql3d
    ! default logic sprinkled around the code, and condense
    ! it here

    if (present(mnemonic)) then
       setup0%mnemonic = mnemonic
    else
       setup0%mnemonic = 'default_output.nc'
    end if
    if (present(ioutput)) then
       setup0%ioutput = ioutput
    else
       ! setup0%ioutput = ! default(s) somewhere deep in code
    end if
    if (present(iuser)) then
       setup0%iuser = iuser
    else
       setup0%iuser = 'default_user'
    end if
    if (present(ibox)) then
       setup0%ibox = ibox
    else
       ! setup0%ibox =! default(s) somewhere deep in code
    end if
    if (present(noplots)) then
       noplots = noplots
    else
       !setup0%noplots = im not convinced this works the way I thought
    end if
    if (present(lnwidth)) then
       setup0%lnwidth = lnwidth
    else
       ! setup0%lnwidth =! default(s) somewhere deep in code
    end if
    if (present(nmlstout)) then
       setup0%nmlstout = nmlstout
    else
       setup0%nmlstout = 'disabled'
    end if
    if (present(special_calls)) then
       setup0%special_calls = special_calls
    else
       setup0%special_calls = 'disabled'
    end if
    if (present(cqlpmod)) then
       setup0%cqlpmod = cqlpmod
    else
       setup0%cqlpmod = 'disabled'
    end if
    if (present(lrz)) then
       setup0%lrz = lrz
    else
       print *, 'setup0%lrz is required'
       call abort
    end if
    if (present(lrzdiff)) then
       setup0%lrzdiff = lrzdiff
    else
       !setup0%lrzdiff = ! default(s) somewhere deep in code
    end if
    if (present(lrzmax)) then
       setup0%lrzmax = lrzmax
    else
       !setup0%lrzmax = ! default(s) somewhere deep in code
    end if
    if (present(lrindx)) then
       setup0%lrindx = lrindx
    else
       !setup0%lrindx = ! default(s) somewhere deep in code
    end if
    if (present(ls)) then
       setup0%ls = ls
    else
       !setup0%ls = ! default(s) somewhere deep in code
    end if
    if (present(lsmax)) then
       setup0%lsmax = lsmax
    else
       !setup0%lsmax = ! default(s) somewhere deep in code
    end if
    if (present(lsdiff)) then
       setup0%lsdiff = lsdiff
    else
       !setup0%lsdiff = ! default(s) somewhere deep in code
    end if
    if (present(lsindx)) then
       setup0%lsindx = lsindx
    else
       !setup0%lsindx = ! default(s) somewhere deep in code
    end if
    if (present(nlrestrt)) then
       setup0%nlrestrt = nlrestrt
    else
       ! this and nlwritf relate to dumping distribution
       ! function in nc files for restarts
       ! we may want for diagnostic debug etc
       !setup0%nlrestrt =! default(s) somewhere deep in code
    end if
    if (present(nlwritf)) then
       setup0%nlwritf = nlwritf
    else
       !setup0%nlwritf = ! default(s) somewhere deep in code
    end if
  end subroutine set_setup0

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
