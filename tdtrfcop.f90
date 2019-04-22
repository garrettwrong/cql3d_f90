module tdtrfcop_mod

!
!

contains

      subroutine tdtrfcop(kopt)
      use param_mod
      use comm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real*8 (a-h,o-z)

!..............................................................
!     copy a distribution function to another dist. function
!     Used for different diagnostics.
!
!     CQL3D:
!     kopt= 1: f -> f_ and frn -> f, (before first call to diaggnde)
!     2: fvn -> f (before second call to diaggnde)
!     3: f_ -> f (before third call to diaggnde)
!     CQLP:
!     kopt= 1: f -> f_
!     2: fvn -> f
!     3: f_ -> f (=> 1 and 3 give same f)
!..............................................................


      if (cqlpmod .ne. "enabled") then
        if (kopt .eq. 1) then
          call dcopy(iyjx2*ngen,f(0:iyjx2*ngen-1,0,1,l_),1, &
                f_(0:iyjx2*ngen-1,0,1,l_),1)
          call dcopy(iyjx2*ngen,frn(0:iyjx2*ngen-1,0,1,l_),1, &
               f(0:iyjx2*ngen-1,0,1,l_),1)
        else if (kopt .eq. 2) then
          call dcopy(iyjx2*ngen,fvn(0:iyjx2*ngen-1,0,1,l_),1, &
                f(0:iyjx2*ngen-1,0,1,l_),1)
        else if (kopt .eq. 3) then
          call dcopy(iyjx2*ngen,f_(0:iyjx2*ngen-1,0,1,l_),1, &
                f(0:iyjx2*ngen-1,0,1,l_),1)
        endif
      else
        if (kopt .eq. 1) then
          call dcopy(iyjx2*ngen,f(0:iyjx2*ngen-1,0,1,l_),1, &
                f_(0:iyjx2*ngen-1,0,1,l_),1)
          call dcopy(iyjx2*ngen,fnp1(0:iyjx2*ngen-1,0,1,l_),1, &
               f(0:iyjx2*ngen-1,0,1,l_),1)
        else if (kopt .eq. 2) then
          call dcopy(iyjx2*ngen,fnhalf(0:iyjx2*ngen-1,0,1,l_),1, &
                f(0:iyjx2*ngen-1,0,1,l_),1)
        else if (kopt .eq. 3) then
          call dcopy(iyjx2*ngen,f_(0:iyjx2*ngen-1,0,1,l_),1, &
                f(0:iyjx2*ngen-1,0,1,l_),1)
        endif
      endif
      return
      end
end module tdtrfcop_mod