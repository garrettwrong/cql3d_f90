module esefld_mod

  !---BEGIN USE

  use iso_c_binding, only : c_double
  use param_mod

  !---END USE

!
!

contains

      subroutine efld_cd(dz,ls,vnorm,flux1,flux2,elparnw,flux0)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     Calculates the electric field (V/cm) and flux0 (1/(cm**2*sec)),
!     to maintain constant flux (current)
!     with no external electric field, such as in a RF driven torus.
!     (The internal electric field is generated by charge separation.)
!     Uses fluxes from the h and g functions.
!.......................................................................


!..... input: ls,dz,flux0,flux1,flux2

      dimension dz(0:ls+1),flux1(0:ls+1),flux2(0:ls+1)

!.....input/output:
      dimension elparnw(0:ls+1)

!.......................................................................
!
!.......................................................................


      sum1=0.d0
      sum2=0.d0
      do kk=1,ls
         sum1=sum1 + flux1(kk)*dz(kk)/flux2(kk)
         sum2=sum2 + dz(kk)/flux2(kk)
      enddo

      flux0=vnorm*sum1/sum2

      do kk=1,ls
         delparnw=(flux0-flux1(kk)*vnorm)/(flux2(kk)*vnorm)
         elparnw(kk)=elparnw(kk) +300.d0*delparnw
      enddo

      return
      end
!
!
      real(c_double) function fluxpar &
           (kopt,x,coss,cynt2,cint2,f,iy,jx)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Compute the parallel flux ---
!     parallel flux (particles/sec per cm^2)=fluxpar*vnorm,
!     where vnorm is the velocity normalization in cm/sec,
!     i.e., fluxpar has the units of 1/cm^3.
!.................................................................
!     kopt=1 the flux from the total function f
!     kopt=2 the flux from half of f  at v_par >0
!     kopt=3 the flux from half of f  at v_par <0
!..................................................................

      real(c_double) x(jx),coss(iy),cynt2(iy),cint2(jx), &
                       f(0:iy+1,0:jx+1)

      fluxpar=0.d0

      if (kopt.eq.1) then
        do 20 i=1,iy
          dum=cynt2(i)*coss(i)
          do 10 j=1,jx
            fluxpar=fluxpar+f(i,j)*cint2(j)*x(j)*dum
 10        continue
 20      continue
      endif

      if (kopt.eq.2) then
        ihy=iy/2
        do i=1,ihy
          dum=cynt2(i)*coss(i)
          do j=1,jx
            fluxpar=fluxpar+f(i,j)*cint2(j)*x(j)*dum
           enddo
         enddo
      endif

      if (kopt.eq.3) then
        ihy=iy/2
        do i=ihy+1,iy
          dum=cynt2(i)*coss(i)
          do j=1,jx
            fluxpar=fluxpar+f(i,j)*cint2(j)*x(j)*dum
           enddo
         enddo
      endif

      return
      end
end module esefld_mod
