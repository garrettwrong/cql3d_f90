c
c
      subroutine tdtrfcop(kopt)
      use param_mod
      use cqcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     copy a distribution function to another dist. function
c     Used for different diagnostics.
c
c     CQL3D:
c     kopt= 1: f -> f_ and frn -> f, (before first call to diaggnde)
c     2: fvn -> f (before second call to diaggnde)
c     3: f_ -> f (before third call to diaggnde)
c     CQLP:
c     kopt= 1: f -> f_  
c     2: fvn -> f
c     3: f_ -> f (=> 1 and 3 give same f)
c..............................................................


      if (cqlpmod .ne. "enabled") then
        if (kopt .eq. 1) then
          call dcopy(iyjx2*ngen,f(0:iyjx2*ngen-1,0,1,l_),1,
     +          f_(0:iyjx2*ngen-1,0,1,l_),1)
          call dcopy(iyjx2*ngen,frn(0:iyjx2*ngen-1,0,1,l_),1,
     +         f(0:iyjx2*ngen-1,0,1,l_),1)
        else if (kopt .eq. 2) then
          call dcopy(iyjx2*ngen,fvn(0:iyjx2*ngen-1,0,1,l_),1,
     +          f(0:iyjx2*ngen-1,0,1,l_),1)
        else if (kopt .eq. 3) then
          call dcopy(iyjx2*ngen,f_(0:iyjx2*ngen-1,0,1,l_),1,
     +          f(0:iyjx2*ngen-1,0,1,l_),1)
        endif
      else
        if (kopt .eq. 1) then
          call dcopy(iyjx2*ngen,f(0:iyjx2*ngen-1,0,1,l_),1,
     +          f_(0:iyjx2*ngen-1,0,1,l_),1)
          call dcopy(iyjx2*ngen,fnp1(0:iyjx2*ngen-1,0,1,l_),1,
     +         f(0:iyjx2*ngen-1,0,1,l_),1)
        else if (kopt .eq. 2) then
          call dcopy(iyjx2*ngen,fnhalf(0:iyjx2*ngen-1,0,1,l_),1,
     +          f(0:iyjx2*ngen-1,0,1,l_),1)
        else if (kopt .eq. 3) then
          call dcopy(iyjx2*ngen,f_(0:iyjx2*ngen-1,0,1,l_),1,
     +          f(0:iyjx2*ngen-1,0,1,l_),1)
        endif
      endif
      return
      end
