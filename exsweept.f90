module exsweept_mod

  !---BEGIN USE

  use diagwrng_mod, only : diagwrng

  !---END USE

!
!

contains

      subroutine exsweept(k)
      use param_mod
      use comm_mod
      use advnce_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!..................................................................
!     This routine checks the solution obtained after the velocity split
!     by differentiation of the solution in time and comparison with
!     the right hand side. temp4(i,j) will contain the l.h.s.and will b
!     multiplied by cint2(j)*cynt2(i,l_) so that it will contain the
!     the change in the local particle number density. temp4(i,j) will
!     be the r.h.s.  and the two arrays should (ideally) be identical.
!     Each array is summed for purposes of comparison.
!..................................................................


      sumleft=0.
      sumright=0.
      do 90 i=1,iy

!..................................................................
!     Pass/trapped boundary condition (itl and itu)
!..................................................................

        if (i .eq. itl .or. i .eq. itu) go to 90
        do 100 j=1,jx

!..................................................................
!     differentiate the flux hfu
!..................................................................

          temp6(i,j)=hfu(i,j)
          temp3(i,j)=(hfu(i,j)-hfu(i-1,j))*dx(j)*twopi &

!..................................................................
!     Particle source
!..................................................................

            +.5*so(i,j)*vptb(i,lr_)*cint2(j)*cynt2(i,l_)
          temp3(i,j)=temp3(i,j)*one_*dtr

!..................................................................
!     differentiate l.h.s.
!..................................................................

          temp4(i,j)=(temp2(i,j)-temp1(i,j)) &
            *vptb(i,lr_)*cint2(j)*cynt2(i,l_) &
            *one_
          sumleft=sumleft+temp4(i,j)
          sumright=sumright+temp3(i,j)
 100    continue
 90   continue

!..................................................................
!     Pass/trapped contribution...
!..................................................................

!dir$ ivdep
      do 200 j=1,jx
        temp6(itl,j)=hfu(itl,j)
        temp6(itu,j)=hfu(itu,j)
        temp3(itl,j)=-(hfu(itl-1,j)-2.*hfu(itl,j)-hfu(itu,j))*pi*dx(j) &
          +.5*cah(itl,j)*temp2(itl,j)*vptb(itl,lr_) &
          *cynt2(itl,l_)*cint2(j) &
          +.5*so(itl,j)*vptb(itl,lr_)*cint2(j)*cynt2(itl,l_)
        temp3(itl,j)=temp3(itl,j)*one_*dtr
        temp3(itu,j)=temp3(itl,j)
        temp4(itl,j)=(temp2(itl,j)-temp1(itl,j))*vptb(itl,lr_)*cint2(j) &
          *cynt2(itl,l_)*one_
        temp4(itu,j)=temp4(itl,j)
        sumleft=sumleft+2.*temp4(itl,j)
        sumright=sumright+2.*temp3(itl,j)
 200  continue
      error=(sumleft-sumright)/xlndn(k,lr_)
      if (iactst.eq."abort") then
        if (error.gt.1.e-8) call diagwrng(15)
      endif
      return
      end
end module exsweept_mod
