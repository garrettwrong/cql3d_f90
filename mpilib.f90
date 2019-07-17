module mpilib_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use r8subs_mod, only : dcopy

  !---END USE

  !YuP[2019-06-26] All subroutines from here are moved into mpins_par.f (in /mpi/)
  !  to be inserted into a_cqlp.f90, tdchief.f90,  diagimpd.f90.
  !  Do not need mpilib.f90 anymore.

contains

!      subroutine init_mpi  !used in a_cqlp.f90 (only)  !content is moved to mpins_par.f 
!      include 'mpilib.h'
!      return
!      end subroutine init_mpi
!-------------------------------------------------------
!      subroutine close_mpi  !used in a_cqlp.f90 (only)  !content is moved to mpins_par.f
!      include 'mpilib.h'
!      return
!      end subroutine close_mpi
!-------------------------------------------------------

end module mpilib_mod
