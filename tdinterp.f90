module tdinterp_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use zcunix_mod, only : coeff1
  use zcunix_mod, only : terp1

  !---END USE




contains

      subroutine tdinterp(ilow,iup,x,y,nold,xnew,ynew,nnew)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This subroutine uses the NCAR cubic spline routines to calculate
!     ynew as a function of xnew given y as a function of x.  The
!     boundary conditions are specified as follows:
!     ilow : 'zero', set derivative at lower x value to zero
!     'free', let derivative at lower x value be free
!     iup  : 'zero', set derivative at upper x value to zero
!     'free', let derivative at upper x value be free
!..................................................................
! YuP[2018-09-19] Noticed that sometimes Te is getting to
! negative values at plasma edge.
! It is related to interpolation by subroutine tdinterp().
! When tdinterp() is called by profiles.f,
! at the last (upper) point in ryain-coordinate,
! it uses iup="free" boundary condition, which means
! the first derivative is calculated by fitting a cubic to 4 last points
! in original arrays. This procedure may give the interpolated value
! of Te() at the last point in rya() a zero or negative value.
! To avoid such condition, an additional option is added
! for calling tdinterp() with iup="linear",
! which means that the derivative
! is set from the last two points in original array,
!   cd2(nold)= (y(nold)-y(nold-1))/(x(nold)-x(nold-1)),
! where x() corresponds to ryain() array, nold==njene.

      parameter(kjx=5000)
      parameter(kjx3p1=3*kjx+1)
      dimension xnew(nnew),ynew(nnew),x(nold),y(nold)
      dimension tab(3),itab(3),i0p(2),d2(kjx),work(kjx3p1)
      dimension xt(kjx),yt(kjx)
      character*(*) ilow,iup

!..................................................................
!     must assure that x(1).le.x(2)....
!..................................................................


      if(x(2).le.x(1)) go to 10
      iord=0
      if(nold .gt. kjx) go to 100
 60   continue

      if(ilow.eq."zero") then
        i0p(1)=2
        d2(1)=0.d0
      elseif(ilow.eq."free") then
        ! i0p(1)= 4  the first derivative at x(1) is
        ! calculated by fitting a cubic to points
        ! x(1) through x(4).
        i0p(1)=4
      else ! "linear" !YuP[2018-09-19] added
        i0p(1)=2
        ! i0p(1)= 2  first derivative given at x(1).
        ! Get it from first two points in (x,y) arrays:
        d2(1)=(y(2)-y(1))/(x(2)-x(1))
      endif

      if (iup.eq."zero") then
        i0p(2)=2
        d2(nold)=0.d0 ! zero deriv.
      elseif(iup.eq."free") then
        ! i0p(2)= 4  the first derivative at x(iup) is
        ! calculated by fitting a cubic to points
        ! x(iup-3) through x(iup).
        i0p(2)=4
      else ! "linear" !YuP[2018-09-19] added
        i0p(2)=2
        ! i0p(2)= 2  first derivative given at x(iup).
        ! Get it from last two points in (x,y) arrays:
        d2(nold)=(y(nold)-y(nold-1))/(x(nold)-x(nold-1))
      endif

      call coeff1(nold,x,y,d2,i0p,1,work)
      itab(1)=1
      itab(2)=0
      itab(3)=0
      do 500 i=1,nnew
        call terp1(nold,x,y,d2,xnew(i),1,tab,itab)
        ynew(i)=tab(1)
 500  continue
      if(iord.eq.0) return
      do 70 i=1,nold
        x(i)=xt(i)
        y(i)=yt(i)
 70   continue
      return
 10   continue
      iord=1
      do 30 i=1,nold
        xt(i)=x(i)
        yt(i)=y(i)
 30   continue
      do 50 i=1,nold
        x(i)=xt(nold+1-i)
        y(i)=yt(nold+1-i)
 50   continue
      go to 60
 100  write(*,1000)
 1000 format("failure in tdinterp: parameter kjx too small")
      call exit(1)
      end subroutine tdinterp
      
end module tdinterp_mod
