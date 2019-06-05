module wpcthta_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine wpcthta
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine calculates the coefficient cthta(i,j) encounting for
!     the mu*grad_parallel B force, when cqlpmod=enabled
!..................................................................

!.......................................................................

      do 100 j=1,jx
        if (mod(nummods,10).le.4 .or. lmidvel.eq.0) then
          ztra=-0.5*x(j)*vnorm*psisp(l_)/psis(l_)
        else
          ztra=-0.5*x(j)*vnorm*(psisp(l_)+psisp(l_+1)) &
            /(psis(l_)+psis(l_+1))
        endif
        do 110 i=1,iy
          cthta(i,j)=ztra*sinn(i,l_)*dyi(i,l_)
 110    continue
 100  continue

      return
      end subroutine wpcthta

end module wpcthta_mod
