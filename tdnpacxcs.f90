module tdnpacxcs_mod

  !---BEGIN USE

  use iso_c_binding, only : c_double

  !---END USE

contains

       REAL FUNCTION cxcs(enin)
       IMPLICIT integer (i-n), real*8 (a-h,o-z)
       REAL*8,DIMENSION(9)::a
       REAL*8::enmin,enmax,x,t1,t2,t3,t4,t5,t6,t7,t8,en
!vt    use this function to get CX cross sections, V. Tang 9-25-05
!vt    Version 1.0, originally finalized 10/14/06, tested and compared
!vt    with plots from Redbook.

! gives cx cross sections for H-H like CX events, given a collisional energy;
! data is from oak ridge Atomic Data for Fusion, Vol1.
! Input: en is given in keV/amu of the moving species.
! Output: cs is in cm^2

! BH100423:  See http://www-amdis.iaea.org/,
! BH100423:  go to ALLADIN, Atomic & Molecular Data
! BH100423:  Substituted alcheb subroutine, below.
! BH100423:  Had problem here.

        en=enin*1000 !converts to eV
! uses chebyshev polynomials
! ORNL Cheb fitting parameters for cx
       a(1)=-72.6656
       a(2)=-5.49142
       a(3)=-3.42948
       a(4)=-1.98377
       a(5)=-.878009
       a(6)=-.198932
       a(7)=.0837431
       a(8)=.121252
       a(9)=.0827182
       enmin=1.2e-1
       enmax=6.3e5
!       print *,en, enmin,enmax
       dum1=log(en)-log(enmin)
       dum2=log(enmax)-log(en)
       dum3=log(enmax)-log(enmin)
       x=(dum1-dum2)/dum3
!       print *,x
! The chebyshev polynomials
       t1=x
       t2=2*x**2-1
       t3=4*x**3-3*x
       t4=8*x**4-8*x**2+1
       t5=16*x**5-20*x**3+5*x
       t6=32*x**6-48*x**4+18*x**2-1
       t7=64*x**7-112*x**5+56*x**3-7*x
       t8=128*x**8-256*x**6+160*x**4-32*x**2+1
       dum1=a(1)/2+a(2)*t1
       dum2=a(3)*t2+a(4)*t3+a(5)*t4
       dum3=a(6)*t5+a(7)*t6+a(8)*t7+a(9)*t8
! cxcs=exp(a(1)/2+a(2)*t1+a(3)*t2+a(4)*t3+a(5)*t4+a(6)*t5+a(7)*t6+a(8)*t7+a(9)*t8)/100/100
       cxcs=exp(dum1+dum2+dum3)
!       cxcs=1e-16 !testing
       END FUNCTION cxcs





! BH100423: From http://www-amdis.iaea.org/
! BH100423: Following subroutine calculates CS cross-section from
! BH100423: from the Chebyshev coefficients.  Replaces above calc.
!
!######################################################################
!

      subroutine alcheb(pet, pcf, kncf, pfit, kermsg)
!
!     this is an ornl:cfadc subroutine to calculate cross sections in
!     (cm[2]) versus energy in (ev/amu) or rate coefficients in
!     (cm[3]/s) versus maxwellian temperature in (ev) from chebyshev
!     polynomial fitting coefficients
!
!     these fits are valid only between the limits emin and emax,
!     which are coefficients pcf(10) and pcf(11) in the entry data field
!
!     pet = collision energy in ev/amu or maxwellian temperature in ev
!
!     kermsg = blank if no errors
!
!     pfit = cross section in cm[2] or rate coefficient in cm[3]/s
!
!     written by h. hunter, cfadc oak ridge national laboratory
!     (modified to aladdin calling structure 4/21/88 r.a. hulse)
!
!------------------------------------------------------------------------
!
      real(c_double) pet, pcf, pfit
      real(c_double) emin, emax, cheb, eminl, emaxl, enl, xnorm
      real(c_double) twox, prev, prev2
      dimension pcf(11)
      character*(*) kermsg
      emin = pcf(10)
      emax = pcf(11)
      if(pet .ge. emin .and. pet .le. emax) then
        kermsg = ' '
      else
        kermsg = 'outside range of fit in alcheb'
        return
      endif
!
!     calculate polynomial using recursion relation
!
      k = 9
      cheb = pcf(k)
      eminl = dlog(emin)
      emaxl = dlog(emax)
      enl= dlog(pet)
      k = k-1
      xnorm = (enl-eminl-(emaxl-enl)) / (emaxl-eminl)
      twox = 2.0d0 *  xnorm
      prev2 = 0.0d+00
   10 prev = cheb
      if(k .ne. 1) then
        cheb = pcf(k) + twox*prev - prev2
        prev2 = prev
        k = k-1
        go to 10
      endif
      cheb = 0.5d0*pcf(1) + xnorm*prev - prev2
      pfit = dexp(cheb)
  100 return
!
      end
end module tdnpacxcs_mod
