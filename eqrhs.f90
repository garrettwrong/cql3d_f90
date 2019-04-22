module eqrhs_mod

!
!

contains

      subroutine eqrhs(neq,t,yv,ydot)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

      dimension yv(2), ydot(2)
!..................................................................
!     This routine is used with O.D.E. solver LSODE. It provides
!     the r.h.s. of the set of coupled O.D.E.'s .
!..................................................................

      dpsidr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,1,0)
      dpsidz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz, &
        epsirz,nnra,0,1)
      rbval=sqrt(dpsidr**2+dpsidz**2+(fpsi_)**2)
!
!     Calc Br and Bz for integration along the  poloidal flux surface.
!     But, omit cursign factor, so start in pos-Z dirn.
!
      ydot(1)=dpsidz/(rbval)
      ydot(2)=-dpsidr/(rbval)
      return
      end
end module eqrhs_mod
