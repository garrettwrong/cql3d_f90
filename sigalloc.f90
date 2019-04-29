module sigalloc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast

  !---END USE

!
!

contains

      subroutine sigalloc
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!dir$ nobounds

!..................................................................
!     A check on allocations is sucessful entering then exiting
!     the subroutine.
!..................................................................
      write(*,*)'sigalloc:  Entering sigalloc'



      msig=0
      do 10 i=1,6
        msig=msig+isigmas(i)
 10   continue

      jxis=(jx*(jx+1))/2
      mtab=mtaba

!.......................................................................
!     Allocate allocatable arrays for sigma-v modules
!.......................................................................

      lnln=jxis*(mmsv+1)*msig
      lntab=mtab

      lndumsg=lnln+lntab

      allocate(csv(jxis,0:mmsv,msig),STAT=istat)
      allocate(svtab(mtab),STAT=istat)
      allocate(tamm1(0:mmsv),STAT=istat)
      allocate(iind(jx),STAT=istat)
      call bcast(csv,zero,SIZE(csv))
      call bcast(svtab,zero,SIZE(svtab))
      call bcast(tamm1,zero,SIZE(tamm1))
      call ibcast(iind,0,SIZE(iind))

      write(*,*)'sigalloc:  Leaving sigalloc'

      return
      end
end module sigalloc_mod
