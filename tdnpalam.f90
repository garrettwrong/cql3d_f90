module tdnpalam_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdnpabscs_mod, only : bscs

  !---END USE


contains

      subroutine lam(en_m,ene,stoplamda)
      use param_mod
      use comm_mod
!.............................................
! This subroutine determines the neutral mean free path as a function of
! its energy.
! It uses beam-stopping cross sections from tdnpabscs.f and assumes that
! the energtic neutral has energy greater than Ti and Te.
! The results are stored in array stoplamda(1:lrzmax,1:n_en_length) in cm
! Inputs: en_ is an array of energetic neutral energies with length nen
! en_ is the left edge of energy bin....
! Version 1.0
! Outputs checked and looks reasonable.
!..............V.Tang 9-25-05...............................................
       implicit integer (i-n), real(c_double) (a-h,o-z)
       real(c_double),dimension(lrzmax,nen_npa)::stoplamda
       real(c_double),dimension(nen_npa)::en_m
       real(c_double)::stopcs
       real(c_double),dimension(lrzmax)::ene

!       print*, "lam called"
       print*,lrzmax, nen_npa
       do 200 nn=1,lrzmax !flux surface
       do 100 mm=1,nen_npa !energy
!$$$       print*,'lam: loop enetered'
!$$$       print*,'lam: electron den is:',ene(nn)
!$$$       print*,'lam: temp(kelec,nn)is:',temp(kelec,nn)
!$$$       print*,'lam: zeff(nn)is:',zeff(nn)
!$$$       print*,'lam: fast neut engy:',en_m(mm)/bnumb(1)
       stopcs=0.
       call bscs(en_m(mm)/bnumb(1),ene(nn), &
        temp(kelec,nn),zeff(nn),stopcs) ! calls bscs to get cross section
!$$$       print*,'cross section:', stopcs
       stoplamda(nn,mm)=1/(ene(nn)*stopcs) ! mean-free path=1/(ne*sigma_stop) in cm
!$$$       print *,"lam: flux surface",nn,"Energy index",mm
!$$$       print *,"lam: Mean free path in cm",stoplamda(nn,mm)
100    continue
200    continue
       return
       END
end module tdnpalam_mod
