module siggy_mod

!
!

contains

      real*8 function siggy(ee)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
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
      end
end module siggy_mod
