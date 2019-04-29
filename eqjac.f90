module eqjac_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use zcunix_mod, only : terp1
  use zcunix_mod, only : terp2

  !---END USE

!
!

contains

      subroutine eqjac(neq,t,yv,ml,mu,pd,nrowpd)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      dimension yv(neq),pd(nrowpd,neq)
!..................................................................
!     This routine is used with the O.D.E. solver LSODE. It
!     provides the Jacobian of the system of equations.
!..................................................................

!..................................................................
!     Provide the various derivatives of epsi (first and second
!     order).
!..................................................................

      dpsidr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,1,0)
      dpsidz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,0,1)
      d2psidrr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,2,0)
      d2psidrz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,1,1)
      d2psidzz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,0,2)

!..................................................................
!     Compute B and it's first order derivatives.
!..................................................................

      alpha=(dpsidr**2+dpsidz**2+(fpsi_)**2)
      sqrtalp=sqrt(alpha)
      dalphadr=2.*(dpsidz*d2psidrz+dpsidr*d2psidrr)
      dalphadz=2.*(dpsidz*d2psidzz+dpsidr*d2psidrz)

!..................................................................
!     Compute the 4 elements of the Jacobian.
!..................................................................

      pd(1,1)=(d2psidrz-dpsidz*dalphadr/alpha)/alpha
      pd(1,2)=(d2psidzz-dpsidz/alpha*dalphadz)/alpha
      pd(2,1)=(-d2psidrr+dpsidr/alpha*dalphadr)/alpha
      pd(2,2)=(-d2psidrz+dpsidr/alpha*dalphadz)/alpha
      return
      end
end module eqjac_mod
