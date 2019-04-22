module eqalloc_mod

!
!

contains

      subroutine eqalloc
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!...................................................................
!     Allocate allocatable arrays
!...................................................................

!dir$ nobounds

!..................................................................
!     A check on allocations is sucessful entering then exiting
!     the subroutine.
!..................................................................
      write(*,*)'eqalloc:  Entering eqalloc'

      lnlfield=lfielda*lrzmax
      lndumeq=4*lnlfield
      allocate(drpmconz(lrzmax),STAT=istat)
      call bcast(drpmconz,zero,SIZE(drpmconz))
      allocate(eqdell(lfielda,lrzmax),STAT=istat)
      call bcast(eqdell,zero,SIZE(eqdell))
      allocate(eqbpol(lfielda,lrzmax),STAT=istat)
      call bcast(eqbpol,zero,SIZE(eqbpol))
      allocate(solr(lfielda,lrzmax),STAT=istat)
      call bcast(solr,zero,SIZE(solr))
      allocate(solz(lfielda,lrzmax),STAT=istat)
      call bcast(solz,zero,SIZE(solz))

      write(*,*)'eqalloc:  Leaving eqalloc'

      return
      end
end module eqalloc_mod