module tdsxray_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use tdsxr0_mod, only : tdsxr0

  !---END USE

!
!

contains

      subroutine tdsxray(icall,iplotsxr)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     sets up call to soft-x-ray analyzer
!..................................................................

      character*8 icall,iplotsxr

      do 1 l=1,lrzmax
        tr1(l)=reden(kelec,l)
 1    continue
!BH081106:  In some radcoord cases, rrz contains normalized radial
!BH081106:  coord data, and is not suitable for eqmod.ne."enabled"
!BH081106:  circ plasma model, or the eqmod.eq."enabled" eqdsk
!BH081106:  equilibria.
      if (eqmod.ne.'enabled') then
         if (radcoord.ne.'sqtorflx') &
            write(*,*)'tdsxray: WARNING, check our radial coord further'
         call tdsxr0(rrz,tr1(1),icall,iplotsxr)
      else
         call tdsxr0(rpmconz,tr1(1),icall,iplotsxr)
      endif
      return
      end
end module tdsxray_mod
