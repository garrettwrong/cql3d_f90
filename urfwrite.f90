module urfwrite_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use urfwrite__mod, only : urfwrite_
  !---END USE

!
!

contains

      subroutine urfwrite
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine performs the converse operation on urfread:
!     i.e., a disk file "rayop" in written or updated with ray tracing
!     data.
!..................................................................

!MPIINSERT_INCLUDE

!MPIINSERT_IF_RANK_NE_0_RETURN

      krf=0
      if (lh.eq."enabled") then
        krf=krf+1
        open(unit=20,file='raylh',status='old')
        call urfwrite_(krf,20)
      endif
      if (ech.eq."enabled") then
        krf=krf+1
        open(unit=23,file='rayech',status='old')
        call urfwrite_(krf,23)
      endif
      if (fw.eq."enabled") then
        krf=krf+1
        open(unit=24,file='rayfw',status='old')
        call urfwrite_(krf,24)
      endif
      return
      end
end module urfwrite_mod
