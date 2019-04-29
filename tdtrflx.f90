module tdtrflx_mod

  !---BEGIN USE

  use bcast_mod, only : bcast

  !---END USE


!
!

contains

      subroutine tdtrflx
      use param_mod
      use comm_mod
      use r8subs_mod, only : cvmgt

      implicit integer (i-n), real*8 (a-h,o-z)

!..............................................................
!     This routine checks the accuracy of the radial time advancement
!     by integrating the flux at the outside. This should equal
!     the change in the number of particles in the plasma over
!     the last time step.
!..............................................................


      include 'trans.h'

      call bcast(fxsp,zero,iyjx2*ngen*lrors)  !BH080428, why?
      flxout=0.
      if (soln_method.ne.'it3drv') then
      do 10 k=1,ngen
        do 20 j=1,jx
          do 30 i=1,iytr(lrors)
            l=lrors-1
            id=idx(i,l)
            if (l.ne.lpt(i)) then
              flxout=cosovb(id,l)*(h_r(l)*bovcos(id,l) &
                *sfu(i,j,k,l))*cynt2_(id,l)*cint2(j)*4.*pi**2*radmaj &
                +flxout
            else if (l.eq.lpt(i)) then
              id=idx(i,l)
              i_=iytr(lrors)+1-i
              i_d=idx(i_,l)
              flxout=flxout+(.5*cosovb(id,l)*h_r(l)*bovcos(id,l) &
                *sfu(i,j,k,l)+ &
                .5*cosovb(i_d,l)*h_r(l)*bovcos(i_d,l)*sfu(i_,j,k,l)) &
                *cynt2_(i_d,l)*cint2(j)*4.*pi**2 &
                *radmaj
            endif
 30       continue
 20     continue
 10   continue

      else   !Using sfup rather than sfu [BH:probably can make coding
             !more efficient].

      do k=1,ngen
        do j=1,jx
          do i=1,iytr(lrors)
            l=lrors-1
            id=idx(i,l)
            if (l.ne.lpt(i)) then
              flxout=cosovb(id,l)*(h_r(l)*bovcos(id,l) &
                *sfup(i,j,k,l))*cynt2_(id,l)*cint2(j)*4.*pi**2*radmaj &
                +flxout
            else if (l.eq.lpt(i)) then
              id=idx(i,l)
              i_=iytr(lrors)+1-i
              i_d=idx(i_,l)
              flxout=flxout+(.5*cosovb(id,l)*h_r(l)*bovcos(id,l) &
                *sfup(i,j,k,l)+ &
                .5*cosovb(i_d,l)*h_r(l)*bovcos(i_d,l)*sfup(i_,j,k,l)) &
                *cynt2_(i_d,l)*cint2(j)*4.*pi**2 &
                *radmaj
            endif
          enddo
        enddo
      enddo
      endif

      flxout=flxout*dttr
      return
      end
end module tdtrflx_mod
