module bcast_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  !---END USE

! XXX every thing in this file and related to it should be deleted

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
      end
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
      end

! NME bcast routine for complex arrays
      subroutine ccast(c,cval,n)
      implicit integer (i-n), complex*16 (c)
      dimension c(n)
      do 100 i=1,n
         c(i)=cval
 100  continue
      return
      end
end module bcast_mod
