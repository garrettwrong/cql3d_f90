module tdboothi_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use zcunix_mod, only : coeff1
  use zcunix_mod, only : terp1

  !---END USE

!
!

contains

      subroutine tdboothi
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      parameter(itlrza=3*lrza+1)

!     Calculation of Hirshman bootstrap current expression,
!     according to Phys. Fluids 31, p3150-2 (1988),
!     for arbitrary aspect ratio, banana regime, single
!     electron and ion species.
!     The Hirshman expressions are adapted for multiple
!     species in CQL3D by using Zeff and a combined
!     Ti and ion pressure piofr
!
!     bscurm(0:lrza,2,2) refers to bootstrap current density as
!     a function of radius for electrons:ions, thermal:nonthermal.
!     bscurmi(0:lrza,2,2) are radially integrated cumulative values,
!     and bscurma(2,2) are total radially integrated values.
!     (amps/cm**2) and amps.


      dimension workk(itlrza)

!     dimensions for "temp of r", "pres of r", refer to
!     tempofr(radius, e:i, thermal:nonthermal, 0th:1st:2nd deriv):
      dimension tempofr(lrza,2,2,3),presofr(lrza,2,2,3)

!..................................................................
!     This version currently is written for eqmod="enabled" -
!     the circular analogue will be added later. The formula for
!     the current is MKS. So all quantities will be put into MKS
!     for the evaluation of the current which will remain in
!     Amps/cm**2
!..................................................................
!
! --- statement functions (x = trapped particle fraction, z = charge #,
! --- see equations 11a,11b,11c,and equation for D(X) ):
!     (Snatched from ONETWO, 990823).
      a31(xx,zz)  = ((zz+2.21)*zz+0.754)+xx*((zz+1.243)*zz+0.348)
      a32(zz)    =  0.884+2.074*zz
      alphai(xx) = -1.172/(1.0 + 0.462*xx)
!BobH990823      df(xx,zz)   = (1.414+zz)*zz+(xx*((zz+1.243)*zz+3.48)
      ddf(xx,zz)   = (1.414+zz)*zz+(xx*((zz+1.243)*zz+0.348) &
               + ((2.0*zz+2.657)*zz+0.754))*xx
!     (End-Snatch).
!
! --- Statement functions for Sauter, Angioni, Lin-Liu in banana limit
!     (O. Sauter, C. Angioni, and Y.R. Lin-Liu, "Neoclassical
!      conductivity and bootstrap current formulas for general
!      axisymmetric equilibria and arbitrary collisionality regime",
!      Phys. of Plasmas 6, 2834 (1999)).
      sa31(xx,zz)=xx+xx*(1.4+xx*(-1.9+xx*(0.3+0.2*xx)))/(1.+zz)

      sa32(xx,zz)=(0.05+0.62*zz)*(xx-xx**4)/zz/(1.+0.44*zz) &
                  +(xx**2-xx**4-1.2*(xx**3-xx**4))/(1.+0.22*zz) &
                  -(0.56+1.93*zz)*(xx-xx**4)/zz/(1.+0.44*zz) &
                  +4.95*(xx**2-xx**4-0.55*(xx**3-xx**4))/(1.+2.48*zz)
!    +            +1.2*xx**4/(1.+0.5*zz)
!    +            -1.2*xx**4/(1.+0.5*zz)
      salpha(xx)=-1.17*(1-xx)/(1.-0.22*xx-0.19*xx**2)


      if (eqmod.ne. "enabled") return

!..................................................................
!     Will use splines to determine derivatives of densities and
!     temperatures vs. rho (=rz)
!..................................................................


      call bcast(bscurm,zero,(lrza+1)*4)
      call bcast(bscurmi,zero,(lrza+1)*4)
      call bcast(bscurma,zero,4)

      if (bootst.ne."enabled" .and. jhirsh.ne.88) return


      call bcast(tempofr,zero,lrza*12)
      call bcast(presofr,zero,lrza*12)

      i1p(1)=4
      i1p(2)=4

!..................................................................
!     Thermal electron logrithmic derivatives
!     Find profiles.
!     Use coeff1/terp1 to obtain radial derivative.
!..................................................................

      if (kelecm.ne.0) then

      do k=ngen+1,ntotal
      if (k.eq.kelecm) then
        do l=1,lrzmax
           tempofr(l,1,1,1)=temp(k,l)
           presofr(l,1,1,1)=reden(k,l)*temp(k,l)
           tr(l)=rz(l)*.01
        enddo
      endif
      enddo

      call coeff1(lrzmax,tr(1),tempofr(1,1,1,1),tempofr(1,1,1,3), &
           i1p,1,workk)
      call coeff1(lrzmax,tr(1),presofr(1,1,1,1),presofr(1,1,1,3), &
           i1p,1,workk)

      itab(1)=0
      itab(2)=1
      itab(3)=0
      do l=1,lrzmax
         call terp1(lrzmax,tr(1),tempofr(1,1,1,1),tempofr(1,1,1,3), &
              tr(l),1,tab,itab)
         tempofr(l,1,1,2)=tab(2)/tempofr(l,1,1,1)
         call terp1(lrzmax,tr(1),presofr(1,1,1,1),presofr(1,1,1,3), &
              tr(l),1,tab,itab)
         presofr(l,1,1,2)=tab(2)/presofr(l,1,1,1)
      enddo

      endif

!.......................................................................
!     Thermal ion logrithmic derivatives
!     (Use density weighted ion temperature).
!     Find profiles.
!     Use coeff1/terp1 to obtain radial derivative.
!.......................................................................

      if (nionm.ne.0) then

      do k=ngen+1,ntotal
      do kk=1,nionm
         if (k.eq.kionm(kk)) then
            do l=1,lrzmax
               tempofr(l,2,1,1)=tempofr(l,2,1,1)+reden(k,l)
               presofr(l,2,1,1)=presofr(l,2,1,1)+reden(k,l)*temp(k,l)
            enddo
         endif
      enddo
      enddo

      do l=1,lrzmax
         tempofr(l,2,1,1)=presofr(l,2,1,1)/tempofr(l,2,1,1)
         tr(l)=rz(l)*.01
      enddo



      call coeff1(lrzmax,tr(1),tempofr(1,2,1,1),tempofr(1,2,1,3), &
           i1p,1,workk)
      call coeff1(lrzmax,tr(1),presofr(1,2,1,1),presofr(1,2,1,3), &
           i1p,1,workk)

      itab(1)=0
      itab(2)=1
      itab(3)=0
      do l=1,lrzmax

         call terp1(lrzmax,tr(1),tempofr(1,2,1,1),tempofr(1,2,1,3), &
              tr(l),1,tab,itab)
         tempofr(l,2,1,2)=tab(2)/tempofr(l,2,1,1)
         call terp1(lrzmax,tr(1),presofr(1,2,1,1),presofr(1,2,1,3), &
              tr(l),1,tab,itab)
         presofr(l,2,1,2)=tab(2)/presofr(l,2,1,1)
      enddo

      endif


!..................................................................
!     Nonthermal (general species) electron logrithmic derivatives
!     Find profiles.
!     Use coeff1/terp1 to obtain radial derivative.
!..................................................................

      if (kelecg.ne.0) then

      do k=1,ngen
      if (k.eq.kelecg) then
        do l=1,lrzmax
           tempofr(l,1,2,1)=(2./3.)*energy(k,l)
           presofr(l,1,2,1)=reden(k,l)*(2./3.)*energy(k,l)
           tr(l)=rz(l)*.01
        enddo
      endif
      enddo

      call coeff1(lrzmax,tr(1),tempofr(1,1,2,1),tempofr(1,1,2,3), &
           i1p,1,workk)
      call coeff1(lrzmax,tr(1),presofr(1,1,2,1),presofr(1,1,2,3), &
           i1p,1,workk)

      itab(1)=0
      itab(2)=1
      itab(3)=0
      do l=1,lrzmax
         call terp1(lrzmax,tr(1),tempofr(1,1,2,1),tempofr(1,1,2,3), &
              tr(l),1,tab,itab)
         tempofr(l,1,2,2)=tab(2)/tempofr(l,1,2,1)
         call terp1(lrzmax,tr(1),presofr(1,1,2,1),presofr(1,1,2,3), &
              tr(l),1,tab,itab)
         presofr(l,1,2,2)=tab(2)/presofr(l,1,2,1)
      enddo

      endif

!..................................................................
!     Nonthermal (general species) ion logrithmic derivatives
!     (Use density weighted ion temperature).
!     Find profiles.
!     Use coeff1/terp1 to obtain radial derivative.
!..................................................................

      if (niong.ne.0) then

      do k=1,ngen
      do kk=1,niong
         if (k.eq.kiong(kk)) then
            do l=1,lrzmax
               tempofr(l,2,2,1)=tempofr(l,2,2,1)+reden(k,l)
               presofr(l,2,2,1)=presofr(l,2,2,1) &
                            +reden(k,l)*(2./3.)*energy(k,l)
            enddo
         endif
      enddo
      enddo

      do l=1,lrzmax
         tempofr(l,2,2,1)=presofr(l,2,2,1)/tempofr(l,2,2,1)
         tr(l)=rz(l)*.01
      enddo

      call coeff1(lrzmax,tr(1),tempofr(1,2,2,1),tempofr(1,2,2,3), &
           i1p,1,workk)
      call coeff1(lrzmax,tr(1),presofr(1,2,2,1),presofr(1,2,2,3), &
           i1p,1,workk)

      itab(1)=0
      itab(2)=1
      itab(3)=0
      do l=1,lrzmax

         call terp1(lrzmax,tr(1),tempofr(1,2,2,1),tempofr(1,2,2,3), &
              tr(l),1,tab,itab)
         tempofr(l,2,2,2)=tab(2)/tempofr(l,2,2,1)
         call terp1(lrzmax,tr(1),presofr(1,2,2,1),presofr(1,2,2,3), &
              tr(l),1,tab,itab)
         presofr(l,2,2,2)=tab(2)/presofr(l,2,2,1)
      enddo

      endif

!.......................................................................
!     Loop over flux surfaces calculating bootstrap current
!     <j.B>/<|B|> (Amps/cm**2).
!.......................................................................


!     Minus in following should probably be plus.
!     (10**8 webbers/(guass*cm**2)/(10**2 cm/m))*10**6 cm**3/m**3)*
!        10**3 eV/keV * 1.6e-19 joules/ev * 10**-4 m**2/cm**2
!        /10**2 cm/m, for cgs fpsi/bmidplne). = 1.6022e-10
      cnst=-1.6022e-10


      do l=1,lrzmax

!        k refers to electrons or ions.  kk to thermal or nonthermal.

         do kk=1,2

         fj0oB=-fpsi(l)*presofr(l,1,kk,1)*cnst/ &
                           (dpsidrho(l)*psiavg(1,l)*bmidplne(l))

         if (jhirsh.eq.88) then

         zz=zeff(l)
         xx=trapfrac(l)/(1.-trapfrac(l))

         tmp1=a31(xx,zz)/ddf(xx,zz)
         tmp2=a32(zz)/ddf(xx,zz)

         k=1
         bscurm(l,k,kk)=fj0oB*xx*(a31(xx,zz)*presofr(l,k,kk,2) &
                       -a32(zz)*tempofr(l,k,kk,2))/ddf(xx,zz)

         k=2
         bscurm(l,k,kk)=fj0oB*xx*a31(xx,zz)/ddf(xx,zz)* &
                       tempofr(l,k,kk,1)/(tempofr(l,1,kk,1)*zz)* &
                       (presofr(l,k,kk,2)+alphai(xx)*tempofr(l,k,kk,2))

         elseif(jhirsh.eq.99) then

         zz=zeff(l)
         xx=trapfrac(l)

         tmp1=sa31(xx,zz)/xx
         tmp2=sa32(xx,zz)/xx

         k=1
         bscurm(l,k,kk)=fj0oB*(sa31(xx,zz)*presofr(l,k,kk,2) &
                       +sa32(xx,zz)*tempofr(l,k,kk,2))

         k=2
         bscurm(l,k,kk)=fj0oB*sa31(xx,zz)*(presofr(l,k,kk,1) &
                       /presofr(l,1,kk,1)*presofr(l,k,kk,2) &
                       +salpha(xx)*tempofr(l,k,kk,2))


         endif


         enddo

      enddo


      return
      end subroutine tdboothi
      
end module tdboothi_mod
