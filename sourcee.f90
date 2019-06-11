module sourcee_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use sounorm_mod, only : sounorm
  use sourc0_mod, only : sourc0
  use sourcef_mod, only : sourcef
  use sourceko_mod, only : sourceko
  use sourcpwr_mod, only : sourcpwr

  !---END USE

!
!

contains

      subroutine sourcee
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Subroutine sourcee is the controlling routine for the
!     analytic  SOURCE routines, all of which begin with "sou".
!     This is a facile model, it simply computes a source
!     profile (usually Gaussian in nature) with a specified
!     current. A more sophisticated model for ions utilizes NFREYA, the
!     Monte Carlo Beam deposition code. These are the "fr" routines with
!     frmod="enabled"
!     The "sou" routines are independent of the "fr" routines and vice-
!     versa.
!     Also, calls Knock On source modules, if specified.
!..................................................................

      save


!..................................................................
!     return if source is not to be recomputed this time step
!..................................................................

      if (nsou.eq.1 .or. n.eq.0) then
        continue
      elseif (mod(n,nsou).eq.1 .and. n.ne.1) then
        continue
      else
        return
      endif

!..................................................................
!     Return if the number of sources per species (nso) = 0
!..................................................................

      if (nso.eq.0) then
        return
      endif

!..................................................................
!     Initialization...
!..................................................................

      if (n.eq.0) then

!..................................................................
!     Compute constants used to determine Gaussian source profiles.
!     If soucoord="cart" specify cartesian species parameters,
!     and if soucoord="polar" specify polar parameters.
!     If soucoord="disabled", no gaussian sources.
!..................................................................

      if (soucoord .ne. "disabled") then
        do 20 k=1,ngen
!DIR$ NEXTSCALAR
          do 10 m=1,nso
            if (soucoord.eq."cart") then
!990131              sxllm1(k,m,lr_)=sign(1.,sellm1(k,m))*
              sxllm1(k,m,lr_)=sign(one,sellm1(k,m))* &
                sqrt(abs(sellm1(k,m))/fions(k))
              sxllm2(k,m,lr_)=sellm2(k,m)/fions(k)
              sxppm1(k,m,lr_)=sqrt(seppm1(k,m)/fions(k))
              sxppm2(k,m,lr_)=seppm2(k,m)/fions(k)
            else
              xem1(k,m,lr_)=sqrt(sem1(k,m)/fions(k))
              xem2(k,m,lr_)=sem2(k,m)/fions(k)
              cosm1(k,m,lr_)=cos(sthm1(k,m)*pi/180.)
              cosm2(k,m,lr_)=scm2(k,m)
            endif
            zm1(k,m,lr_)=zmax(lr_)*szm1(k,m)
            zm2(k,m,lr_)=(zmax(lr_)*szm2(k,m))**2
 10       continue
 20     continue
      endif  ! On soucoord

!..................................................................
!     sounor will contain normalization constants after the call
!     to sounorm
!..................................................................

        call bcast(sounor(1:ngen,1:nsoa,1:lz,lr_),one,ngen*nsoa*lz)
        if (soucoord .ne. "disabled") then
        call sounorm
        do 40 k=1,ngen
          do 30 m=1,nso
            if (nonso(k,m) .eq. 1) nonso(k,m)=0
 30       continue
 40     continue
        endif  ! On soucoord
      endif  ! On n.eq.0

!..................................................................
!     Initialize the source profile to zero.
!..................................................................

      call bcast(source(0:iy+1,0:jx+1,1:ngen,indxlr_),zero,iyjx2*ngen)

!..................................................................
!     xlncur will contain the source current (/cm**2/sec).
!     In general asor*zmaxpsi(lr_)=xlncur  (asor in units particles/cc)
!..................................................................

      call bcast(xlncur(1:ngen,lr_),zero,ngen)

!..................................................................
!     Determine Guassian+knock-on  source(i,j,k,indxlr_) array
!..................................................................

      call sourcef

      call sourceko

!..................................................................
!     define source profile uniquely at x=0.
!..................................................................

      call sourc0

!..................................................................
!     Compute the source power.
!..................................................................

      do 70 k=1,ngen
        call sourcpwr(k)
 70   continue
      return
      end subroutine sourcee


end module sourcee_mod
