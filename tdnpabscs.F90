module tdnpabscs_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE


contains

       subroutine bscs(eng,ne,te,zeff,stopcs)
!........................................................
! This function gives the hydrogen beam stopping cross section based
! on Janev's 1989 NF Paper (pg. 2125); assumes a Carbon impurity
! V. Tang 9-27-05
! Inputs: eng=energy in keV/amu,
!         ne: electron density in 1/cm3,
!         te=electron temp in keV
!         zeff=zeff of plasma
! Outputs: stopcs: stopping cross section in cm**2
! Checks out with paper plots, 9-28-05
! Version 1.0, NPA mod
! BH100423:  Problem a low temperature.  NEED TO REPLACE.
!......................................................
       IMPLICIT NONE
       real(c_double),INTENT(IN)::eng,ne,te,zeff
       real(c_double),INTENT(OUT)::stopcs
       real(c_double),DIMENSION(3,3,3)::a,b
       real(c_double)::si,sz,i1,j1,k1,ne2,no
       INTEGER::i,j,k
!       print*,eng,ne,te,zeff,stopcs
!       ne2=ne/1e6 !converts electron density to 1/cm3
       no=1.0e13  !constant in fit
       ne2=ne/no !divide by no

!       print *,eng,ne2,te,zeff

! A coefficients
       a(1,1,1)=4.40
       a(1,1,2)=-2.49e-2
       a(1,2,1)=7.46e-2
       a(1,2,2)=2.27e-3
       a(1,3,1)=3.16e-3
       a(1,3,2)=-2.78e-5
       a(2,1,1)=2.30e-1
       a(2,1,2)=-1.15e-2
       a(2,2,1)=-2.55e-3
       a(2,2,2)=-6.2e-4
       a(2,3,1)=1.32e-3
       a(2,3,2)=3.38e-5

! B coefficients for Carbon
       b(1,1,1)=-1.49
       b(1,1,2)=-1.54e-2
       b(1,2,1)=-1.19e-1
       b(1,2,2)=-1.50e-2
       b(2,1,1)=5.18e-1
       b(2,1,2)=7.18e-3
       b(2,2,1)=2.92e-2
       b(2,2,2)=3.66e-3
       b(3,1,1)=-3.36e-2
       b(3,1,2)=3.41e-4
       b(3,2,1)=-1.79e-3
       b(3,2,2)=-2.04e-4

       si=0.0
       sz=0.0

       do 30 i=1,2
       do 20 j=1,3
       do 10 k=1,2
       i1=i-1.0
       j1=j-1.0
       k1=k-1.0
       si=si+a(i,j,k)*(log(eng))**i1*((log(ne2))**j1*(log(te))**k1)
 10        continue
 20        continue
 30        continue

       do 60 i=1,3
       do 50 j=1,2
       do 40 k=1,2
           i1=i-1.0
           j1=j-1.0
           k1=k-1.0
       sz=sz+b(i,j,k)*(log(eng))**i1*((log(ne2))**j1*(log(te))**k1)
 40         continue
 50         continue
 60         continue
! stopping cross section in cm^2
       stopcs=exp(si)*(1+(zeff-1)*sz)/eng*1e-16
! testing
!       stopcs=1e-30
       END SUBROUTINE bscs
end module tdnpabscs_mod
