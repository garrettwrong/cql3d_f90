module bcast_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

! XXX every thing in this file and related to it should be deleted
! YuP: basically I agree. I guess the intent was to use a fast procedure with unrolled loops,
! similar to dcopy().
! Never accomplished.
! But now there are too many lines with "bcast" in the source.

!
!

contains

      subroutine bcast(a,val,n)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Temporary bcast routine until I can find UNICOS equivalent
!..................................................................

      dimension a(n)
      do 100 i=1,n
        a(i)=val
 100  continue
      return
      end subroutine bcast
!
!
      subroutine ibcast(ia,ival,n)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Temporary bcast routine until I can find UNICOS equivalent
!..................................................................

      dimension ia(n)
      do 100 i=1,n
        ia(i)=ival
 100  continue
      return
      end subroutine ibcast

! NME bcast routine for complex arrays
      subroutine ccast(c,cval,n)
      implicit integer (i-n), complex*16 (c)
      dimension c(n)
      do 100 i=1,n
         c(i)=cval
 100  continue
      return
      end subroutine ccast

end module bcast_mod
