!
!
      subroutine frsplft(lrz,oldx,oldf,npts,ynewx,ynewf)
      use param_mod
      use zcunix_mod, only : coeff1
      use zcunix_mod, only : terp1
      implicit none

!..................................................................
!     Interpolates with splines between cql3d radial mesh and the
!     (finer) NFREYA radial mesh.
!..................................................................

      integer nwka ! local
      parameter (nwka=3*lrza+1)
      real(c_double) :: work(nwka) ! local
      integer l ! local
      integer lrz,npts ! arguments
      real(c_double) :: oldx(*),oldf(*),ynewx(*),ynewf(*) ! arguments
      integer i1p(2) ! local, and input for coeff1()
      real(c_double) :: secondd(lrza),tab(3) ! local, and arg for terp1()
      integer itab(3) ! local, and input for terp1()

      i1p(1)=4
      i1p(2)=4
      call coeff1(lrz,oldx,oldf,secondd,i1p,1,work)
      itab(1)=1
      itab(2)=0
      itab(3)=0
      do 10 l=1,npts
        call terp1(lrz,oldx,oldf,secondd,ynewx(l),1,tab,itab)
        ynewf(l)=tab(1)
 10   continue
      return
      end subroutine frsplft
