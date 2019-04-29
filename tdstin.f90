module tdstin_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine tdstin
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!


!..................................................................
!     This routine initializes flux surface dependent data so the
!     2-d code initialization routine ainitial can be called.
!..................................................................

      do 2 k=1,ngen
        do 3 m=1,nso
          sellm1(k,m)=sellm1z(k,m,lr_)
          seppm1(k,m)=seppm1z(k,m,lr_)
          sellm2(k,m)=sellm2z(k,m,lr_)
          seppm2(k,m)=seppm2z(k,m,lr_)
          sem1(k,m)=sem1z(k,m,lr_)
          sem2(k,m)=sem2z(k,m,lr_)
          sthm1(k,m)=sthm1z(k,m,lr_)
          scm2(k,m)=scm2z(k,m,lr_)
          szm1(k,m)=szm1z(k,m,lr_)
          szm2(k,m)=szm2z(k,m,lr_)
          asor(k,m,lr_)=asorz(k,m,lr_)
 3      continue
 2    continue

      return
      end
end module tdstin_mod
