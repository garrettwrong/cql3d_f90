module coefefld_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!
!
!

contains

      subroutine coefefld
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     compute coefficients required to represent parallel
!     electric field (flux surface averaged) due to ohmic current
!..................................................................


      call bcast(cex(1:iyjx,1,l_),zero,iyjx)
      call bcast(cet(1:iyjx,1,l_),zero,iyjx)
!     Division by 300. below is conversion to cgs: 300 volts/statvolt.
!%OS  coefld=-radmaj*fpsi(lr_)*onovrp(2,lr_)*flxavgd(lr_)
      if (cqlpmod .ne. "enabled") then
         if (efflag .eq. "toroidal") then
            coefld=-rmag*fpsi(lr_)*onovrp(2,lr_)*flxavgd(lr_) &
                 *charge/300./vnorm
!BH000926 Adding option for specification electric field parallel to B.
         elseif (efflag .eq. "parallel") then
            coefld=-r0drdz(lr_)*charge/300./vnorm
         endif
      endif
!%OS  if (cqlpmod .eq. "enabled") coefld=-radmaj*fpsi(lr_)/solrs(l_)**2
      if (cqlpmod .eq. "enabled") then
        if (mod(nummods,10).le.4 .or. transp.ne."enabled" .or. &
          lmidvel.eq.0) then
          coefld=-rmag*fpsi(lr_)/solrs(l_)**2 &
            /psis(l_)/bmidplne(lr_)*charge/300./vnorm
        else
          coefld=-rmag*fpsi(lr_)/0.125/(solrs(l_)+solrs(l_+1))**2 &
            /(psis(l_)+psis(l_+1))/bmidplne(lr_)*charge/300./vnorm
        endif
      endif

!%OS
!%OS  coefld=-rmag*fpsi(lr_)/solrs(1)**2
!%OS  +                         /psis(1)/bmidplne(lr_)*charge/300./vnorm
!%OS

      iend=itl-1
      if (cqlpmod.eq."enabled" .and. symtrap.ne."enabled") iend=iyh

      do 40 i=1,iend
        ii=iy+1-i
        do 50 j=1,jx
          cex(i,j,l_)=coefld*coss(i,l_)*xsq(j)
          cex(ii,j,l_)=-cex(i,j,l_)
          cet(i,j,l_)=-coefld*sinn(i,l_)**2*x(j)
          cet(ii,j,l_)=cet(i,j,l_)
 50     continue
 40   continue
      return
      end
end module coefefld_mod
