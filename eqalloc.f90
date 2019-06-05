module eqalloc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine eqalloc
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

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
      end subroutine eqalloc


end module eqalloc_mod
