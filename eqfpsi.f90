module eqfpsi_mod

  !---BEGIN USE

  use eqwrng_mod, only : eqwrng
  use zcunix_mod, only : terp1

  !---END USE

!
!

contains

      subroutine eqfpsi(psval,fpsi__,fppsi__)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!..................................................................
!     This routine provides f(psi) to model the toroidal
!     magnetic field. For cases that eqsource="ellipse"
!     the f is ad-hoc and is determined through the namelist
!     model, fpsimodel. In the case that eqsource="filename", then
!     a file exists on disk which provides f and the equilibrium
!     psi. As of 9/21/88 filename=eqdsk or topeol.
!     Also provided is the derivative df/dpsi, fppsi.
!..................................................................

      if (eqsource.eq."ellipse") then
        if (fpsimodl.eq."constant") then
          fpsi__=btor*radmaj
          fppsi__=0.
        else
          call eqwrng(7)
        endif
      else
        itab(1)=1
        itab(2)=1
        itab(3)=0
        call terp1(nfp,psiar,fpsiar,d2fpsiar,psval,1,tab,itab)
        fpsi__=tab(1)
        fppsi__=tab(2)
      endif
      return
      end


!
!
      subroutine eqppsi(psval,ppsi__,pppsi__)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!..................................................................
!     This routine provides p(psi) to model the plasma
!     pressure. For cases that eqsource="ellipse"
!     the p i zero. In the case that eqsource="filename", then
!     a file exists on disk which provides f and the equilibrium
!     psi. As of 9/21/88 filename=eqdsk or topeol.
!     Also provided is the derivative dp/dpsi, pppsi.
!..................................................................

      if (eqsource.eq."ellipse") then
         ppsi__=0.
         pppsi__=0.
      else
        itab(1)=1
        itab(2)=1
        itab(3)=0
        call terp1(nfp,psiar,prar,d2prar,psval,1,tab,itab)
        ppsi__=tab(1)
        pppsi__=tab(2)
      endif
      return
      end
end module eqfpsi_mod
