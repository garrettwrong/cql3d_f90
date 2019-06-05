module siggy_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      real(c_double) function siggy(ee)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!---------------------------------------------------------------------
!    this routine is used in conjunction with array svtab defined
!    in subroutine sigsetup. it is the table look-up.
!---------------------------------------------------------------------

      els=ee+em90
      inum=(els-elmin)/delegy+1.5
      inum=max0(inum,1)
      inum=min0(inum,mtab)
      val=svtab(inum)
      siggy=val
      return
      end function siggy
      
end module siggy_mod
