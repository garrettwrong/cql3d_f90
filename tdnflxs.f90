module tdnflxs_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use cqlcomm_mod
  use diagwrng_mod, only : diagwrng
  use param_mod
  !---END USE

contains

  subroutine tdnflxs(ll)
    !implicit integer (i-n), real(c_double) (a-h,o-z)
    implicit none
    save

    integer :: ll
    integer :: iyy
    !.......................................................................
    !     This routine is used when the code jumps from calculations on
    !     one flux surface to another. It sets a few mesh and time variables
    !.......................................................................
    !

    !..................................................................
    !     l_  defines the current spatial variable mesh label
    !     (radial or parallel).
    !     lr_ defines the effective index of the current flux surface.
    !     ls_ defines the effective index of the current orbit along B.
    !     lmdpln_ defines the orbit index of midplane (l_=lmdpln_)
    !     on current lr_.
    !..................................................................

    !
    l_=ll
    lr_=lrindx(l_)
    indxlr_=indxlr(lr_)
    lmdpln_=lmdpln(indxlr_)
    ls_=lsindx(l_)
    indxls_=indxls(ls_)

    itl=itl_(l_)
    itu=itu_(l_)
    n=n_(l_)
    iyy=iy_(l_) !-YuP-101215: Don't use iy=; it's in common /params/
    ! Don't let overwrite the cqlinput value!
    iyh=iyh_(l_)
    iyjx=iyjx_(l_)
    timet=time_(l_)

    !.......................................................................
    !     check iy mesh
    !.......................................................................

    if (n .gt. 0) then
       if (iyh .ne. iyy/2) call diagwrng(20)
       if (itu .ne. iyy-itl+1) call diagwrng(21)
    endif

    return
  end subroutine tdnflxs

end module tdnflxs_mod
