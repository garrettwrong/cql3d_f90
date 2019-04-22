module hpalloc0_mod

!
!

contains

      subroutine hpalloc0(iptr,ln)
      use param_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!     NO_OP:  Not needed after 081216, after changing from
!             CRAY pointers to f90 pointers.

!..................................................................
!     Allocates memory AND zeroes out the array.
!..................................................................
!990131      integer err
!$$$CPTR>>>DELETE PTR-HPALLOC
!$$$      pointer(iptr,a(1))
!$$$
!$$$c_cray      call hpalloc(iptr,ln,err,1)
!$$$c_pc Unix libU77.a library.  malloc allocates bytes.
!$$$c_pc      iptr=malloc(8*ln)
!$$$c_pc      if (iptr.eq.0) stop 'Problem with malloc in hpalloc0'
!$$$      iptr=malloc(8*ln)
!$$$      if (iptr.eq.0) then
!$$$         write(*,*)'hpalloc0: ln=',ln
!$$$         stop 'Problem with malloc in hpalloc0'
!$$$      endif
!$$$
!$$$      zero=0.d0
!$$$      call bcast(a,zero,ln)
!$$$CPTR<<<END PTR-HPALLOC
      return
      end
end module hpalloc0_mod
