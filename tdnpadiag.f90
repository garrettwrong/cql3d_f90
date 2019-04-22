module tdnpadiag_mod

!
!

contains

      subroutine tdnpadiag(icall)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      character*8 icall

!..................................................................
!     sets up call to NPA diagnostic, Version 1.0
!..................................................................

      character*8 iplotnbi

      write(*,*)
      write(*,*)'tdnpadiag, time step ',n


!     Call npa routines to calc and plot output....

      if (noplots.eq."enabled1") then
         iplotnbi='no'
      else
         iplotnbi='yes'
      endif

      do 1 l=1,lrzmax
        tr1(l)=reden(kelec,l)
 1    continue
      call tdnpa0(rrz,tr1(1),icall,iplotnbi)

      return
      end

end module tdnpadiag_mod