module impnorm_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use iso_c_binding, only : c_double

  !---END USE

contains

  subroutine impnorm(xnorm,a,rhs,nn)
    !implicit integer (i-n), real(c_double) (a-h,o-z)
    implicit none

    !..................................................................
    !     This routine normalizes the matrix a so that the maximum
    !     coefficient for each equation is of order 1.
    !..................................................................
    real(c_double) :: a(nn)
    real(c_double) :: xnorm
    real(c_double) :: rhs
    integer :: i
    integer :: nn


    xnorm=0.d0
    do i=1,nn
       xnorm=xnorm+dabs(a(i))
    end do
    if (xnorm.gt.0.d0) then
       rhs=rhs/xnorm
       do i=1,nn
          a(i)=a(i)/xnorm
       end do
    endif

    return
  end subroutine impnorm

end module impnorm_mod
