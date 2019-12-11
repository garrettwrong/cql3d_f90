! Copyright 2019 Garrett Wright, Princeton Plasma Physics Laboratory,
!    contracted by the U.S. Department of Energy (DE-AC02-09CH11466).
!
! This file is part of cql3d_f90. See LICENSE.
!
! cql3d_f90 is free software: you can redistribute it and/or modify it
! under the terms of the GNU Affero General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! cql3d_f90 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with cql3d_f90.  If not, see <https://www.gnu.org/licenses/>.


!***********************************************************************
!
! For those interested in the historical development of cql3d_f90/CQL3D
!   please see compx_archive/compx_LICENSE and compx_archive/a_change
!
!***********************************************************************


program a_cql3d
  use abchief_mod, only : abchief
  use impavnc0_mod, only : de_alloc
  use impavnc0_mod, only : it3dalloc
  use impavnc0_mod, only : it3ddalloc
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_float
#ifdef __MPI
  use cqlmpilib_mod, only : cql_mpi_init
#endif


  real(c_float) tarray(2)
  character(len=*),parameter :: nml_file = 'cqlinput'

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

mpirank=0 ! When MPI is used, mpirank is set in init_mpi below

!     Initialize MPI:
#ifdef __MPI
      call cql_mpi_init()
!     It will insert:
!      call init_mpi ! Initialize MPI and set mpitime0=MPI_WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
#endif

call cpu_time(tarray(1))    !This is an f95 intrinsic subroutine
!------------!
call abchief(nml_file) !-> calls tdchief (only)
!------------!
call cpu_time(tarray(2))

#ifdef __MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
      if(mpirank.eq.0) then
#endif

if(setup0%verbose>0) WRITE(*,'(a,i5,f10.3)') ' a_cqlp: rank, Exec.time tarray(2)-tarray(1)', mpirank, tarray(2)-tarray(1)
!      WRITE(*,'(a)') ' a_cqlp: END of CQL3D, just before MPI_FINISH'

#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

call it3ddalloc ! Deallocate it3d related storage
call de_alloc   ! Deallocate other arrays

#ifdef __MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
!     close MPI (print 'MPI Full time =',MPI_WTIME()-mpitime0
!                then - MPI_FINALIZE )
      if(mpirank.eq.0) then
         if(setup0%verbose>0) WRITE(*,*) 'MPI Full time =',MPI_WTIME()-mpitime
      endif
      call MPI_FINALIZE(mpiierr)
      !PRINT *,'close_mpi:  mpirank===',mpirank
#endif


call exit(0)
stop
end program a_cql3d
