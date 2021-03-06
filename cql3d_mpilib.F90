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

module cqlmpilib_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use r8subs_mod, only : dcopy

  !---END USE

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

contains

  subroutine cql_mpi_init(skip_mpi_init)
    logical, intent(in), optional :: skip_mpi_init
    logical ::  do_mpi_init

    ! in standalone, we are repsonsible for MPI init
    do_mpi_init = .TRUE.
    if(present(skip_mpi_init)) then
       if(skip_mpi_init) do_mpi_init=.FALSE.
    end if

    if(do_mpi_init) call MPI_INIT(mpiierr)

    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpiierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpiierr)
    if(mpisize.le.1) stop '===   Run with number of cores >1   ==='
    if(mpirank.eq.0) then
       mpitime = MPI_WTIME()
    endif
  end subroutine cql_mpi_init

  subroutine cql3d_set_comm(comm)
    !     set communicator
    use cqlcomm_mod, only :cql3d_comm
#ifdef __MPI
    ! was for env_comm_add use mpi_proc_data         ! (from portlib library)
    implicit none
    include 'cql3d_mpilib.h'

    integer, intent(in) :: comm ! communicator

    !------------------------
    integer :: ierloc
    !------------------------

    cql3d_comm = comm

    !GBW, we should probably keep this in transp...
    ! call mpi_env_comm_add(comm,ierloc) ! ierloc not checked at present
#endif

  end subroutine cql3d_set_comm

  subroutine cql3d_set_myid(procid)
    ! set proc id
    implicit none
#ifdef __MPI
    include 'cql3d_mpilib.h'
#endif
    integer, intent(in) :: procid

#ifdef __MPI
    mpirank = procid
#endif

  end subroutine cql3d_set_myid

  subroutine cql3d_set_nprocs(nprocs)
    ! set proc id
    implicit none
#ifdef __MPI
    include 'cql3d_mpilib.h'
#endif
    integer, intent(in) :: nprocs

#ifdef __MPI
    mpisize = nprocs
#endif

  end subroutine cql3d_set_nprocs

end module cqlmpilib_mod
