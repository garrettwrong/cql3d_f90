module eqonovrp_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use eqflxavg_mod, only : eqflxavg
  use eqorbit_mod, only : eqorbit

  !---END USE

!
!

contains

      subroutine eqonovrp(epsicon_,onovrp1,onovrp2)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)


!..................................................................
!     This routine computes <1./R**2>, <1./R> for flux surface psival
!..................................................................

      if (eqorb.eq."enabled") then
        call eqorbit(epsicon_)
      endif
      do 20 ipower=1,2
        do 10 l=1,lorbit_
          tlorb1(l)=1./(solr_(l)**ipower)
 10     continue
        call eqflxavg(epsicon_,tlorb1,onovrs,flxavgd_)
        if (ipower.eq.1) onovrp1=onovrs
        if (ipower.eq.2) onovrp2=onovrs
 20   continue
      return
      end subroutine eqonovrp


end module eqonovrp_mod
