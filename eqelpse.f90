module eqelpse_mod

  !---BEGIN USE

  use eqfn_mod, only : eqfn

  !---END USE


!
!

contains

      subroutine eqelpse
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save
!
      character*8 ifirst


!..................................................................
!     This routine generates an ad-hoc psi (an "solution" to
!     the equilibrium calculation)  with elliptical contours of psi
!     centered at rmag with ellipticity = ellptcty. The funtional form
!     of psi(E,lr_) where E**2=z**2+(r-rmag)**2/(1.-ellptcty**2) is
!     arbitrary and is handled by function eqfn. The psi function
!     is normalized so that the poloidal field at r=rmag, z=radmin
!     is bth.
!..................................................................

!..................................................................
!     Create the z,r meshes..
!..................................................................

      data ifirst /"first"/
      if (lr_.ne.lrzmax) return
      if (ifirst.eq."first") then
        zst=zbox/(nnz-1)
        rst=rbox/(nnr-1)
        er(1)=rboxdst
        ez(1)=-zbox*.5
        do 10 nn=2,nnr
          er(nn)=er(nn-1)+rst
 10     continue
        do 11 nn=2,nnz
          ez(nn)=ez(nn-1)+zst
 11     continue
        zmag=zero
        ifirst="notfirst"
      endif

!..................................................................
!     Set up the psi array. First determine the normalization
!     constant from the poloidal field constraint.
!..................................................................

      fhp=eqfn(radmin*1.0001,one)
      fhm=eqfn(radmin*.9999,one)
      deriv=(fhp-fhm)/(radmin*2.e-4)
      scalfct=-bth*rmag/deriv

!..................................................................
!     scalfct is the factor, now for the psi array..
!..................................................................

      do 20 ir=1,nnr
        do 21 iz=1,nnz
          epsi(ir,iz)=sqrt(ez(iz)**2+(er(ir)-rmag)**2/(1.-ellptcty**2))
 21     continue
 20   continue
      do 30 ir=1,nnr
        do 31 iz=1,nnz
          epsi(ir,iz)=eqfn(epsi(ir,iz),scalfct)
 31     continue
 30   continue
      return
      end
end module eqelpse_mod
