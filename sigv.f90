module sigv_mod

  !---BEGIN USE

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
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     sets up calls to sigma-v calculator
!..................................................................

!MPIINSERT_INCLUDE

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

      do 100  i=1,lrz           ! keep index 'i' - used by mpiworker_i below
!MPIINSERT_MPIWORKER_I
!MPIINSERT_IF_RANK_EQ_MPIWORKER
         call tdnflxs(i)        !-> determine flux surface index lr_
!        Compute reaction rate of general species with self and others.
         call sigv5d
!        Compute reaction rates for "equivalent" maxwellians.
! YuP[06-2016] Why do we need to call this sub again?
! All cases of isigsgv1 are treated inside of sigv5d*.
! Commenting out:
!YuP            if (isigsgv1.eq.1) call sigv5d
!MPIINSERT_ENDIF_RANK
!MPIINSERT_RECV_FUS
!MPIINSERT_SEND_FUS
 100  continue                  ! on i=1,lrz

!MPIINSERT_BARRIER
!MPIINSERT_BCAST_FUS

!
!         Four fusion reactions are treated (see cqlinput_help)
!         fuspwrv = fusion power density (from general species)
!         fuspwrm = fusion power density (Maxwellian with same density
!                                         energy as the gen. species)
      do ll=1,lrz
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
      end
end module sigv_mod
