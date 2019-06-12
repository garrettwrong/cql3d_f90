module vlfalloc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ccast
  use bcast_mod, only : ibcast

  !---END USE

!
!

contains

      subroutine vlfalloc
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      complex*16 czero
!dir$ nobounds

!..................................................................
!     A check on allocations is sucessful entering then exiting
!     the subroutine.
!..................................................................
      write(*,*)'vlfalloc:  Entering vlfalloc'

!.......................................................................
!     Called from subroutine vlfsetup
!     Allocate allocatable arrays for rf modules
!.......................................................................
      czero = (0.0,0.0)

      lniylz=iy*lz*setup0%lrzmax
      lnj=jx
      lni=iy
      lnjj=jjx

!...................................................................
!     Note three complex*16 arrays, cosz1, sinz1, sinz2 must be 2x longer
!     Same for cwexde,cweyde and cwezde.
!...................................................................

!BH031009      lnurfdum=lniylz+10*lnj+3*2*lni+2*lnjj
      lnurfdum=lniylz+8*lnj+(3*2+1)*lni+4*lnjj
      allocate(cosmz(iy,lz,setup0%lrzmax),STAT=istat)
      call bcast(cosmz,zero,SIZE(cosmz))
      allocate(alfag(jx),STAT=istat)
      call bcast(alfag,zero,SIZE(alfag))
      allocate(argmnt(jx),STAT=istat)
      call bcast(argmnt,zero,SIZE(argmnt))
      allocate(ilim1d(jx),STAT=istat)
      call ibcast(ilim1d,0,SIZE(ilim1d))
      allocate(ilim2d(jx),STAT=istat)
      call ibcast(ilim2d,0,SIZE(ilim2d))
      allocate(ilim1dd(jx),STAT=istat)
      call ibcast(ilim1dd,0,SIZE(ilim1dd))
      allocate(ilim2dd(jx),STAT=istat)
      call ibcast(ilim2dd,0,SIZE(ilim2dd))
      allocate(sx(jx),STAT=istat)
      call bcast(sx,zero,SIZE(sx))
      allocate(xmdx(jx),STAT=istat)
      call bcast(xmdx,zero,SIZE(xmdx))
      allocate(cosz1(iy),STAT=istat)
      call ccast(cosz1,czero,SIZE(cosz1))
      allocate(sinz1(iy),STAT=istat)
      call ccast(sinz1,czero,SIZE(sinz1))
      allocate(sinz2(iy),STAT=istat)
      call ccast(sinz2,czero,SIZE(sinz2))
      allocate(alfi(iy),STAT=istat)
      call bcast(alfi,zero,SIZE(alfi))
      allocate(ilim1(jjx),STAT=istat)
      call ibcast(ilim1,0,SIZE(ilim1))
      allocate(ilim2(jjx),STAT=istat)
      call ibcast(ilim2,0,SIZE(ilim2))
      allocate(ifct1(jjx),STAT=istat)
      call ibcast(ifct1,0,SIZE(ifct1))
      allocate(ifct2(jjx),STAT=istat)
      call ibcast(ifct2,0,SIZE(ifct2))

      write(*,*)'vlfalloc:  Leaving vlfalloc'

      return
      end subroutine vlfalloc

end module vlfalloc_mod
