!
!
      subroutine frsmooth(k,curnorm)
      use param_mod
      use cqlcomm_mod
      !lrzmax =  setup0%lrzmax
      !mnemonic = setup0%mnemonic
      use cqlconf_mod, only : setup0

      implicit none
      integer k ! input
      real(c_double) :: curnorm ! output
      integer l,ll ! local
      real(c_double) :: smth ! local

      ! YuP[2019-06-12] Getting lrz from setup0 type:
      integer :: lrz
      lrz = setup0%lrz


      if (smooth_ .lt. .005) return

!..................................................................
!     This routine smooths the source deposition current. The weighting
!     function is a Gaussian with a s.deviation (smooth_) which is
!     normally some small fraction (~.05-.1) of the normalized toroidal
!     coordinate radius radmin. This allows the resulting current
!     profile to be far smoother that what is seen with the model turned
!     off (smooth_<.005). This in turn allows for much less noisy
!     p' and ff' input profiles to the MHD equilibrium code.
!     If utilized, array tr3 and scalar curnorm are altered.
!..................................................................


!..................................................................
! Compute the average source: compute the total source current,curnorm
!..................................................................

      curnorm=0.
      smth=smooth_**2*radmin**2
      do 5 l=1,lrz
        tr5(l)=0.
        if (asorz(k,1,l).gt.em120) tr5(l)=1./asorz(k,1,l)
 5    continue
      do 10 l=1,lrz
        tr(l)=0.
        tr1(l)=0.
        do 20 ll=1,lrz
          tr(l)=exp(-(rz(ll)-rz(l))**2/smth)*dvol(l)+tr(l)
          tr1(l)=exp(-(rz(ll)-rz(l))**2/smth)*asorz(k,1,ll)*dvol(l)
     &      +tr1(l)
 20     continue
        tr4(l)=tr1(l)/tr(l)
        curnorm=curnorm+dvol(l)*tr4(l)
        tr3(l)=tr4(l)*tr5(l)
 10   continue
      return
      end subroutine frsmooth
