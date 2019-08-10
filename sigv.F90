module sigv_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use sigalloc_mod, only : sigalloc
  use sigsetup_mod, only : sigsetup
  use sigv5d_mod, only : sigv5d
  use tdnflxs_mod, only : tdnflxs

  !---END USE

!
!

contains

      subroutine sigv(icall)
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     sets up calls to sigma-v calculator
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

      character*8 icall
      real(c_double) buff(16) ! local buffer array for mpi

      if (icall.eq."first") then
        call sigalloc
        call sigsetup
!       return, if no nuclear fusion general species
        if (sigmamod.eq.'disabled') return
      endif

      call bcast(fuspwrvt,zero,4)
      call bcast(fuspwrmt,zero,4)
      call bcast(fuspwrv,zero,4*lrorsa)
      call bcast(fuspwrm,zero,4*lrorsa)
      call bcast(sigft,zero,4)
      call bcast(sigmt,zero,4)

      do 100  i=1,setup0%lrz           ! keep index 'i' - used by mpiworker_i below
#ifdef __MPI
         if(mpisize.gt.1) then
            mpiworker= MOD(i-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
#endif
#ifdef __MPI
      if(mpirank.eq.mpiworker) then

#endif
         call tdnflxs(i)        !-> determine flux surface index lr_
!        Compute reaction rate of general species with self and others.
         call sigv5d
!        Compute reaction rates for "equivalent" maxwellians.
! YuP[06-2016] Why do we need to call this sub again?
! All cases of isigsgv1 are treated inside of sigv5d*.
! Commenting out:
!YuP            if (isigsgv1.eq.1) call sigv5d
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif
#ifdef __MPI
      if(mpirank.eq.0) then !-------------------------------------------
        mpisz=4
        call MPI_RECV(buff, mpisz*4,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,mpiierr)
        mpitag=mpistatus(MPI_TAG)
        mpil_=mpitag ! determine which radial surface sent the data
        call dcopy(mpisz, buff(0*mpisz+1),1, fuspwrv(1:mpisz,mpil_),1)
        call dcopy(mpisz, buff(1*mpisz+1),1, fuspwrm(1:mpisz,mpil_),1)
        call dcopy(mpisz, buff(2*mpisz+1),1, sigf(1:mpisz,mpil_),   1)
        call dcopy(mpisz, buff(3*mpisz+1),1, sigm(1:mpisz,mpil_),   1)
      endif !-----------------------------------------------------------

#endif
#ifdef __MPI
      if(mpirank.eq.mpiworker) then !-----------------------------------
        mpisz=4
        call dcopy(mpisz, fuspwrv(1:mpisz,lr_), 1,buff(0*mpisz+1),1)
        call dcopy(mpisz, fuspwrm(1:mpisz,lr_), 1,buff(1*mpisz+1),1)
        call dcopy(mpisz, sigf(1:mpisz,lr_),    1,buff(2*mpisz+1),1)
        call dcopy(mpisz, sigm(1:mpisz,lr_),    1,buff(3*mpisz+1),1)
        mpitag= lr_ ! over flux surfaces =1,setup0%lrz
        call MPI_SEND(buff, mpisz*4,MPI_DOUBLE_PRECISION,0,mpitag,MPI_COMM_WORLD,mpiierr)
      endif !-----------------------------------------------------------
#endif
 100  continue                  ! on i=1,setup0%lrz

#ifdef __MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
#endif
#ifdef __MPI
      call MPI_BCAST(fuspwrv,4*lrorsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(fuspwrm,4*lrorsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sigf,4*lrorsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sigm,4*lrorsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
#endif

!
!         Four fusion reactions are treated (see cqlinput_help)
!         fuspwrv = fusion power density (from general species)
!         fuspwrm = fusion power density (Maxwellian with same density
!                                         energy as the gen. species)
      do ll=1,setup0%lrz
          call tdnflxs(ll)
          do 90 knumb=1,4
            fuspwrvt(knumb)=fuspwrvt(knumb) &
              +fuspwrv(knumb,lr_)*dvol(lr_)
            fuspwrmt(knumb)=fuspwrmt(knumb) &
              +fuspwrm(knumb,lr_)*dvol(lr_)
            sigft(knumb)=sigft(knumb) &
              +sigf(knumb,lr_)*dvol(lr_)
            sigmt(knumb)=sigmt(knumb) &
              +sigm(knumb,lr_)*dvol(lr_)
 90       continue
      enddo ! ll


!     Store time dependent data. nch has been set by ntdstore.
        if (n .eq. nchec*(n/nchec)) then
          do 80 knumb=1,4
            sigftt(nch(1),knumb)=sigft(knumb)
            sigmtt(nch(1),knumb)=sigmt(knumb)
 80       continue
        endif

!BH110401      endif  ! on icall

      return
      end subroutine sigv


end module sigv_mod
