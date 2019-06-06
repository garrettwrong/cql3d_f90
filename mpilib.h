      include 'mpif.h'

      integer mpirank, mpisize, mpiierr, mpitag, mpisrc
      integer mpicmd, mpimode, mpiworker, mpill, mpil_
      integer mpistatus(MPI_STATUS_SIZE)

      real(c_double) mpitime, mpitime1

      common /mpi/ mpirank, mpisrc, mpisize, mpiierr, mpistatus
      common /mpi/ mpicmd, mpimode, mpiworker, mpill
      common /mpi/ mpitime, mpitime1
