module micgetr_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine micgetr(jx,xmax,h,r,ksingul)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This subroutine computes a mesh ratio for the purpose of creating
!     a geometric mesh.
!     The solution is returned in r.
!     Mesh ratio spacing factor h must be bounded by 10.and 0.1
!     Accuracy is to within 1 part in 14 places
!     ksingul=1, if problem.
!..................................................................

      em90=1.d-90
      em14=1.d-14
      if (h .gt. em90 .and. xmax .gt. h .and. jx .ge. 3) go to 10
      ksingul=1
      return
 10   continue
      ehl=h/xmax
      jxm1=jx-1
      jxm2=jx-2
      toll=10.*(1.+em14)
      toli=1./toll
      r=1.
      delr=1.e-01
      sum=ehl
      do 20 k=1,jxm2
 20   sum=r*sum+ehl
      if (abs(sum-1.) .lt. em14) go to 60
      sign=1.
      if (sum .gt. 1.) sign=-1.
 30   r=r+sign*delr
      suml=sum
      if (r .gt. toll .or. r .lt. toli) go to 50
      sum=ehl
      do 40 k=1,jxm2
 40   sum=r*sum+ehl
      if (abs(sum-1.) .lt. em14) go to 60
      if ((sum-1.)*(suml-1.) .gt. 0.) go to 30
      delr=1.d-1*delr
      sign=-sign
      if (delr .lt. .9*em14) go to 60
      go to 30
 50   continue
      ksingul=1
      return
 60   continue
      return
      end subroutine micgetr

end module micgetr_mod
