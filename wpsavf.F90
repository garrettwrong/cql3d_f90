module wpsavf_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine wpsavf
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     This routine saves various versions of the distribution
!     function.
!..............................................................

!.......................................................................

      if (relaxtsp.ne."enabled" .or. n.lt.nonavgf .or. &
        n.gt.nofavgf) go to 200

!.......................................................................
!l    1. relaxtsp=enabled: mix f_n and f_n+1
!.......................................................................

      do 100 k=1,ngen
        do 110 j=0,jx+1
          do 120 l=0,setup0%ls+1
            do 130 i=0,iy_(l)+1
              fnp0(i,j,k,l)=0.5*(fnp0(i,j,k,l)+fnp1(i,j,k,l))
              fnp1(i,j,k,l)=fnp0(i,j,k,l)
 130        continue
 120      continue
 110    continue
 100  continue
      call dcopy(iyjx2*ngen*setup0%ls,fnp1(0:iy+1,0:jx+1,1:ngen,1:setup0%ls),1, &
                                  f(0:iy+1,0:jx+1,1:ngen,1:setup0%ls),1)
      call dcopy(iyjx2*ngen*setup0%ls,fnp1(0:iy+1,0:jx+1,1:ngen,1:setup0%ls),1, &
                                 f_(0:iy+1,0:jx+1,1:ngen,1:setup0%ls),1)

      return

!.......................................................................
!l    2.  f=f_n+1
!.......................................................................

 200  continue

      call dcopy(iyjx2*ngen*(setup0%ls+2),fnp1(0:iy+1,0:jx+1,1:ngen,0:setup0%ls+1),1, &
                                   fnp0(0:iy+1,0:jx+1,1:ngen,0:setup0%ls+1),1)
      call dcopy(iyjx2*ngen*setup0%ls,    fnp1(0:iy+1,0:jx+1,1:ngen,1:setup0%ls  ),1, &
                                      f(0:iy+1,0:jx+1,1:ngen,1:setup0%ls  ),1)
      call dcopy(iyjx2*ngen*setup0%ls,    fnp1(0:iy+1,0:jx+1,1:ngen,1:setup0%ls  ),1, &
                                     f_(0:iy+1,0:jx+1,1:ngen,1:setup0%ls  ),1)

      return
      end subroutine wpsavf


end module wpsavf_mod
