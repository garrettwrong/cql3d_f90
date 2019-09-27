! Copyright 2019 Garrett Wright, Princeton Plasma Physics Laboratory,
!    contracted by the U.S. Department of Energy (DE-AC02-09CH11466).
!
! This file is part of cql3d_f90. See LICENSE.
!
! cql3d_f90 is free software: you can redistribute it and/or modify it
! under the terms of the GNU Affero General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! cql3d_f90 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with cql3d_f90.  If not, see <https://www.gnu.org/licenses/>.

module zcunix_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      real(c_double) function taper(x,x0,x1,x2)
      implicit integer (i-n), real(c_double) (a-h,o-z)

!     A multiplication factor between 0.and 1., giving a tapered
!     (i.e., smoothed) box function.
!
!     For abs(x-x0) .lt. x1/2., gives 1.
!     For (x1+x2)/2 .gt. abs(x-x0) .gt. x1/2., gives
!                                  0.5*(1.+cos(abs(x-x1)/x2*pi))
!     For abs(x-x0) .ge. (x1+x2)/2., gives 0.
!
!     x1 and x2 are presumed .ge.0.

      save pi
      data pi/3.141592653589793/

      xx=abs(x-x0)

      if (xx.lt. x1/2.) then
        taper=1.
      elseif (xx.gt.x1/2. .and. xx.lt.(x1+x2)/2.) then
        taper=0.5*(1.+cos(xx/x2*pi))
      else
        taper=0.0
      endif

      return
      end function taper
!
!     subroutine aminmx(array,ifirst,ilast,istride,amin,amax,
!     *indmin,indmax)
!     real array(1)
!     length = 1+(ilast+1-ifirst)/istride
!     if(ilast.lt.ifirst+(length-1)*istride) length=length-1
!c    k=istride
!     m=ismin(length,array(ifirst),istride)
!     m=ifirst+(m-1)*istride
!c    m=ifirst
!c    do 1 i=ifirst+k,ilast,k
!c1   if(array(i).le.array(m)) m=i
!     amin=array(m)
!     m=ismax(length,array(ifirst),istride)
!     m=ifirst+(m-1)*istride
!c    m=ifirst
!c    do 2 i=ifirst+k,ilast,k
!c2   if(array(i).ge.array(m)) m=i
!     amax=array(m)
!     return
!     end

!*****************************************************************************
!     SPLINE PACKAGE
!*****************************************************************************
!     package cubspl         note--documentation for individual routines
!     follows the general package information
!
!     latest revision        january 1985
!
!     purpose                to perform one and two-dimensional cubic spline
!     interpolation with choice of boundary
!     conditions.  the function and selected
!     derivatives may be evaluated at any point where
!     interpolation is required.  in the
!     two-dimensional case, the given data points
!     must be on a rectangular grid, which need not
!     be equally spaced.  the package cubspl contains
!     seven routines.
!
!     subroutine coeff1
!     computes the coefficients for one-dimensional
!     cubic spline interpolation using one of
!     the following boundary conditions at
!     each end of the range.
!     . second derivative given at boundary.
!     . first derivative given at boundary.
!     . periodic boundary condition.
!     . first derivative determined by fitting a
!     cubic to the four points nearest to the
!     boundary.
!
!     subroutine terp1
!     using the coefficients computed by coeff1,
!     this routine evaluates the function and/or
!     first and second derivatives at any point
!     where interpolation is required.
!
!     subroutine coeff2
!     computes the coefficients for two-dimensional
!     bicubic spline interpolation with the same
!     choice of boundary conditions as for coeff1.
!
!     function terp2
!     using the coefficients produced by coeff2,
!     this subroutine evaluates the function or a
!     selected derivative at any point where
!     two-dimensional interpolation is required.
!
!     subroutine trip
!     a simple periodic, tridiagonal linear
!     equation solver used by coeff1.
!
!     subroutine search
!     performs a binary search in a one-dimensional
!     floating point table arranged in ascending
!     order.  this routine is called by terp1 and
!     terp2.
!
!     subroutine intrp
!     given coefficients provided by coeff1 and the
!     position of the interpolation point in the
!     independent variable table, this routine
!     performs one-dimensional interpolation for
!     the function value, first and second
!     derivative, as desired.  this routine is
!     called by terp1 and terp2.
!
!     usage                  for one-dimensional cubic spline interpolation,
!     the user first calls coeff1 by
!
!     call coeff1 (n,x,f,w,iop,int,wk)
!
!     this subroutine returns the coefficients
!     needed for the subsequent interpolation
!     in the array w.  the user then calls
!     subroutine terp1 by
!
!     call terp1 (n,x,f,w,y,int,tab,itab)
!
!     in order to compute the value of the
!     function and/or its derivatives.  the user
!     specifies y, the value of the independent
!     variable where the interpolation is to be
!     performed.  the interpolated values are
!     returned in tab.  the parameters
!     n,x,f,w, and int must not be changed
!     between the subroutine calls.
!
!     for two-dimensional cubic spline interpolation
!     the user first calls coeff2 by
!
!     call coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,
!     idm,ibd,wk)
!
!     this subroutine returns the coefficients
!     needed for the subsequent interpolation in
!     the arrays fxx, fyy, fxxyy.  the user then
!     calls the routine terp2 by
!
!     r = terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,
!     idm,ixd,iyd)
!
!     in order to perform the interpolation at the
!     point specified by the values of xb and yb.
!     depending on the values input in ixd and iyd,
!     the routine returns an interpolated value
!     for the function or one of its partial
!     derivatives.  the parameters nx,x,ny,y,f,fxx,
!     fyy,fxxyy, and idm must not be changed
!     between the calls to coeff2 and terp2.
!
!     special conditions     tables of independent variables must be
!     arranged in ascending order.  for
!     two-dimensional interpolation, the data points
!     must lie on a rectangular mesh.
!
!     i/o                    none
!
!     precision              single
!
!     files                  cray machines.
!
!     language               fortran
!
!     history                this package is based on the routines
!     la e 102a, spl1d1
!     la e 103a, spl1d2
!     la e 104a, spl2d1
!     la e 105a, spl2d2
!     of the los alamos cubic spline package written
!     by thomas j. jordan and bertha fagan, 1968.
!     the routines have been streamlined and
!     standardized.  the algorithm for
!     two-dimensional interpolation is considerably
!     modified.
!
!     algorithm              for one-dimensional interpolation, the cubic
!     spline is expressed in terms of the function
!     values and second derivatives at the data
!     points.  the second derivatives are
!     determined from a tridiagonal linear system
!     which describes the continuity of the first
!     derivative and incorporates the given
!     boundary conditions.  coeff1  sets up this
!     system and calls subroutine trip to solve it.
!
!     the cubic segment between two adjacent
!     tabular points is constructed from the
!     function values and second derivatives at
!     these points.  these provide the four
!     constants needed to define the cubic
!     uniquely.  from this cubic, values of the
!     function and its first and second
!     derivatives are readily determined at any
!     intermediate point.  one-dimensional
!     interpolation is performed by the routine
!     terp1.  for two-dimensional interpolation,
!     the bicubic spline is described in terms of
!     values of f,fxx,fyy, and fxxyy  at each
!     point on the given two-dimensional
!     rectangular grid of data points.  here f
!     is the function value,
!
!     fxx = (d/dx)**2*f
!
!     and so on.  the coefficients are determined
!     by coeff2, which uses successive
!     applications of coeff1.
!
!     1.  the array fxx is determined by applying
!     coeff1 to f along each line in the
!     x-direction.
!
!     2.  the array fyy is determined by applying
!     coeff1 to f along each line in the
!     y-direction.
!
!     3.  fxxyy is determined on the constant y
!     boundaries by applying coeff1 to fyy
!     along these boundaries.
!
!     4.  the remainder of the array fxxyy is
!     determined by applying coeff1 to fxx
!     along each line in the y-direction.
!
!     the bicubic within any rectangular element
!     of the grid is constructed from the values
!     of f,fxx,fyy,fxxyy at the four corners.
!     these provide the 16 constants necessary
!     to define the bicubic uniquely.  to find
!     the value of f corresponding to a point
!     (xb,yb) within the elementary rectangle,
!     (x(i),y(j)),(x(i+1),y(j)),(x(i),y(j+1)),
!     (x(i+1),y(j+1)), five one dimensional
!     interpolations are performed.
!
!     1.  f at (xb,y(j)) and (xb,y(j+1)) are
!     found by interpolating f in the
!     x-direction using fxx. (two interpolations)
!
!     2.  fyy at (xb,y(j)) and (xb,y(j+1)) are
!     found by interpolating fyy in the
!     x-direction using fxxyy. (two
!     interpolations.)
!
!     3.  finally f at (xb,yb) is found by
!     interpolating between f(xb,y(j)) and
!     f(xb,y(j+1)) in the y-direction using
!     values of fyy(xb,y(j)) and fyy(xb,y(j+1))
!     obtained above. (one interpolation).
!
!     two-dimensional interpolation is performed
!     in terp2.
!
!     references             for greater detail, see j.l.walsh,
!     j.h.ahlberg, e.n.nilsen, best approximation
!     properties of the spline fit, journal of
!     mathematics and mechanics, vol.ii(1962),
!     225-234.
!
!     t.l. jordan, smoothing and multivariable
!     interpolation with splines, los alamos
!     report, la-3137, 1965.
!
!     accuracy               near machine accuracy was obtained when
!     interpolating a cubic in one dimension
!     or a bicubic in two dimensions.
!
!     portability            fully portable with respect to fortran 66.
!***********************************************************************
!
!     subroutine coeff1 (n,x,f,w,iop,int,wk)
!
!     dimension of           x(n),f(int*(n-1)+1),w(int*(n-1)+1),iop(2),
!     arguments              wk(3*n+1)
!
!     purpose               :subroutine coeff1 computes the coefficients
!     for one-dimensional cubic spline
!     interpolation using one of the following
!     boundary conditions at each end of the
!     range
!
!     .  second derivatives given at boundary
!     .  first derivative given at boundary
!     .  periodic boundary conditions
!     .  first derivative calculated by fitting a
!     cubic to the four points nearest to the
!     boundary
!
!     note that terp1 must be called to perform
!     the interpolation.
!
!     usage                  call coeff1 (n,x,f,w,iop,int,wk)
!
!     arguments
!
!     on input               n
!     the number of data points.  n must be at
!     least 4.
!
!     x
!     table of n independent variable values in
!     ascending order.  dimension of x in calling
!     program must be at least n.
!
!     f
!     table of n corresponding dependent variable
!     values.  the values are separated by interval
!     int.  this is usually unity for
!     one-dimensional interpolation.  other values
!     are useful when coeff1 is incorporated in a
!     two-dimensional interpolation scheme (see
!     coeff2).  dimension of f in the calling
!     program must be at least (int*(n-1)+1).
!
!     iop
!     two element integer array defining boundary
!     conditions at x(1) and x(n) according to the
!     following code.
!
!     for iop(1)
!     = 1  second derivative given at x(1).  place
!     value of second derivative in w(1)
!     before call to coeff1.
!     = 2  first derivative given at x(1).  place
!     value of first derivative in w(1) before
!     call.
!     = 3  periodic boundary condition.  x(1) and
!     x(n) are equivalent points.  f(1) and
!     f(int*(n-1)+1) are equal.
!     = 4  the first derivative at x(1) is
!     calculated by fitting a cubic to points
!     x(1) through x(4).
!     similarly, iop(2) defines the boundary
!     condition at x(n).  when iop(2) = 1 (or 2),
!     the value of the second (or first) derivative
!     must be placed in w(int*(n-1)+1).  note that
!     if iop(1) = 3, consistency demands that
!     iop(2) = 3 also.
!
!     int
!     spacing in f and w tables.  for
!     one-dimensional interpolation this will
!     usually be unity.
!
!     wk
!     work area of dimension at least (3*n+1).
!
!     on output              w
!     table of second derivatives corresponding to
!     given x and f values.  the separation of
!     tabular entries is int (see above).
!     dimension of w in the calling program must be
!     at least (int*(n-1)+1).
!
!     the arrays x, f, w are used as input for the
!     routine terp1, which performs interpolation
!     at a given value of the independent variable.
!
!     timing                 the timing is linearly proportional to n, the
!     number of data points.
!***********************************************************************
!
!     subroutine terp1 (n,x,f,w,y,int,tab,itab)
!
!
!     dimension of           x(n),f(int*(n-1)+1),w(int*(n-1)+1),tab(3),
!     arguments              itab(3)
!
!     purpose                using the coefficients computed by coeff1,
!     this routine evaluates the function and/or
!     first and second derivatives at any point.
!
!     usage                  call terp1 (n,x,f,w,y,int,tab,itab)
!
!     arguments
!
!     on input               n
!     the number of data points.  n must be at
!     least 4.
!
!     x
!     table of n independent variable values in
!     ascending order.  dimension of x in the
!     calling program must be at least n.
!
!     f
!     table of n corresponding dependent variable
!     values separated by interval int, usually
!     unity for one-dimensional interpolation.
!     dimension of f in the calling program must be
!     at least (int*(n-1)+1).
!
!     w
!     table of second derivatives computed by
!     coeff1.  the separation of tabular entries is
!     int.  dimension of w in the calling program
!     must be at least (int*(n-1)+1).
!
!     y
!     value of the independent variable at which
!     interpolation is required.  if y lies outside
!     the range of the table, extrapolation takes
!     place.
!
!     int
!     spacing of tabular entries in f and w arrays.
!     this is usually unity for one-dimensional
!     interpolation.
!
!     itab
!     three element integer array defining
!     interpolation to be performed at y.
!     if itab(1) = 1, the function value is
!     returned in tab(1).
!     if itab(2) = 1, the first derivative is
!     returned in tab(2).
!     if itab(3) = 1, the second derivative is
!     returned in tab(3).
!     if itab(i) = 0 for i = 1, 2 or 3, the
!     corresponding function value or derivative is
!     not computed and tab(i) is not referenced.
!
!     on output              tab
!     three element array in which interpolated
!     function value, first and second derivatives
!     are returned as dictated by itab (see above).
!
!     timing                 this procedure is fast.  the maximum time for
!     the binary search is proportional to alog(n).
!     the time for function and derivative evaluation
!     is independent of n.
!***********************************************************************
!
!     subroutine coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,wk)
!
!
!     dimension of           x(nx),y(ny),f(idm,ny),fxx(idm,ny),fyy(idm,ny),
!     arguments              fxxyy(idm,ny),ibd(4),wk(3*max0(nx,ny)+1)
!     (idm must be .ge. nx)
!
!     purpose               :subroutine coeff2 computes the coefficients
!     for two-dimensional bicubic spline
!     interpolation with the same choice of
!     boundary conditions as for coeff1.  terp2
!     is called to perform interpolation.
!
!     usage                  call coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,
!     wk)
!
!     arguments
!
!     on input               nx
!     number of grid points in the x-direction.  nx
!     must be at least 4.
!
!     x
!     table of nx values of the first independent
!     variable arranged in ascending order.
!     dimension of x in the calling program must be
!     at least nx.
!
!     ny
!     number of grid points in the y-direction.  ny
!     must be at least 4.
!
!     y
!     table of ny values of the second independent
!     variable arranged in ascending order.
!     dimension of y in the calling program must be
!     at least ny.
!
!     f
!     two dimensional array of function values at
!     the grid points defined by the arrays x and
!     y.  dimension of f in the calling program is
!     (idm, nyy) where
!     idm .ge. nx
!     nyy .ge. ny
!
!     idm
!     first dimension in the calling program of
!     arrays f, fxx, fyy, fxxyy.  idm must be at
!     least nx.
!
!     ibd
!     four element integer array defining boundary
!     conditions according to the following code.
!     for ibd(1)
!     = 1  the second derivative of f with respect
!     to x is given at (x(1),y(j)) for
!     j = 1,ny,1.  values of this second
!     derivative must be placed in fxx(1,j)
!     for j = 1,ny,1, before calling coeff2.
!     = 2  the first derivative of f with respect
!     to x is given at (x(1),y(j)) for
!     j = 1,ny,1.  values of the derivative
!     must be placed in fxx(1,j) for
!     j = 1,ny,1 before calling coeff2.
!     = 3  periodic boundary condition in the
!     x-direction.  (x(1),y(j)) and
!     and (x(nx),y(j)) are equivalent points
!     for j = 1,ny,1.  f(1,j) and f(nx,j) are
!     equal.
!     = 4  the first derivative of f with respect
!     to x at (x(1),y(j)) is computed by
!     fitting a cubic to f(1,j) through f(4,j)
!     for j = 1,ny,1.
!
!     similarly, ibd(2) defines the boundary
!     condition at (x(nx),y(j)) for j = 1,ny,1.
!     when ibd(2) = 1 (or 2) the values of the
!     second (or first) derivative of f with
!     respect to x are placed in fxx(nx,j) for
!     j = 1,ny,1.
!     note that if ibd(1) = 3, consistency
!     requires that ibd(2) = 3 also.
!     for ibd(3)
!     = 1  the second derivative of f with respect
!     to y is given at (x(i),y(1)).  place
!     values of the derivative in fyy(i,1) for
!     i = 1,nx,1 before calling coeff2.
!     = 2  the first derivative of f with respect
!     to y is given at (x(i),y(1)).  values of
!     this derivative must be placed in
!     fyy(i,1) for i = 1,nx,1 before calling
!     coeff2.
!     = 3  periodic boundary condition in the
!     y-direction.  (x(i),y(1)) and
!     (x(i),y(ny)) are equivalent points.
!     f(i,1) and f(i,ny) are equal.
!     = 4  the first derivative of f with respect
!     to y at (x(i),y(1)) is computed by
!     fitting a cubic to f(i,1) through f(i,4)
!     for i = 1,nx,1.
!
!     similary, ibd(4) defines the boundary
!     condition at (x(i),y(ny)) for i = 1,nx,1 and
!     given derivative values are placed in
!     fyy(i,ny).
!     note that consistency demands that if
!     ibd(3) = 3, then ibd(4) = 3 also.
!
!     wk
!     work area of dimension at least
!     (3*max0(nx,ny)+1)
!
!     on output              fxx
!     array of second derivatives of f with respect
!     to x computed by coeff2.  fxx(i,j) is
!     derivative at (x(i),y(j)).  as for f,
!     dimension of fxx in the calling program is
!     (idm,nyy).
!
!     fyy
!     array of second derivatives of f with respect
!     to y computed by coeff2.  dimension of fyy in
!     the calling program is (idm,nyy).
!
!     fxxyy
!     array of fourth derivatives
!     (d/dx)**2*(d/dy)**2*f, computed by coeff2.
!     dimension of fxxyy in the calling program is
!     (idm,nyy).
!
!     the arrays x, y, f, fxx, fyy, fxxyy are used as
!     input for the routine terp2 which performs
!     interpolation at required values of the two
!     independent variables.
!
!     timing                 the timing is proportional to nx*ny.
!***********************************************************************
!
!     function terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd)
!
!
!     dimension of           x(nx),y(ny),f(idm,ny),fxx(idm,ny),fyy(idm,ny),
!     arguments              fxxyy(idm,ny))
!     (idm must be .ge. nx)
!
!     purpose                using the coefficients produced by coeff2,
!     this routine evaluates the function on a
!     selected derivative of any point where
!     two-dimensional interpolation is required.
!
!     usage                  r = terp2 (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,
!     ixd,iyd)
!
!     arguments
!
!     on input               xb, yb
!     values of the independent variables, x and y,
!     at which interpolation is required.
!
!     nx
!     number of grid points in the x-direction.  nx
!     must be at least 4.
!
!     x
!     table of nx values of the independent
!     variable, x, arranged in ascending order.
!     dimension of x in the calling program must be
!     at least nx.
!
!     ny
!     number of grid points in the y-direction.  ny
!     must be at least 4.
!
!     y
!     table of ny values of the independent
!     variable, y, arranged in ascending order.
!     dimension of y in the calling program must be
!     at least ny.
!
!     f
!     two-dimensional array of function values at
!     grid points defined by the arrays x and y.
!     dimension of f in the calling program is
!     (idm,nyy), where
!     idm .ge. nx
!     nyy .ge. ny
!
!     fxx
!     array of second derivatives of f with respect
!     to x computed by coeff2.  dimension of fxx in
!     the calling program is (idm,nyy).  see under
!     f above.
!
!     fyy
!     array of second derivatives of f with respect
!     to y computed by coeff2.  dimension of fyy in
!     the calling program is (idm,nyy).
!
!     fxxyy
!     array of fourth derivatives,
!     (d/dx)**2*(d/dy)**2*f, computed by coeff2.
!     dimension of fxxyy in the calling program is
!     (idm,nyy).
!
!     idm
!     first dimension in the calling program of
!     arrays f, fxx, fyy and fxxyy,
!     idm .ge. nx
!
!     ixd, iyd
!     define derivative to be returned by the
!     function terp2.  ixd, iyd may each take the
!     the values 0, 1, 2.  the derivative returned
!     is (d/dx)**ixd*(d/dy)**iyd*f.
!     note that if ixd = iyd = 0, the function
!     value itself is returned.
!
!     timing                 this procedure is fast.  the maximum
!     time for the binary search is proportional to
!     alog(nx*ny).  the time for function evaluation
!     is independent of n.
!***********************************************************************
      subroutine coeff1 (n,x,f,w,iop,int,wk)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      dimension       x(n)       ,f(n*int)       ,w(n*int)      ,iop(2), &
        wk(n,*)
!dir$ nobounds
      logical q8q4
      save q8q4
      data q8q4 /.true./
!
!     arithmetic statement function used to locate entries in f and w arrays
!
      ii(index)=(index-1)*int+1
!
!
!
!
!
!
!     the following call is for gathering statistics on library use at ncar
      if (q8q4) then
        q8q4 = .false.
      endif
!
!     start to set up tridiagonal system
!
      j0 = 1
      do 101 i=2,n
        jm = j0
        j0 = j0+int
        wk(i,1) = x(i)-x(i-1)
        wk(i,2) = (f(j0)-f(jm))/wk(i,1)
        wk(i,3) = wk(i,1)/6.
        wk(i,1) = wk(i,1)/3.
 101  continue
      nn = n
      mk = iop(1)
      ml = iop(2)
!
!     apply boundary conditions at boundary 1
!
      go to (102,103,104,105),mk
!
!     second derivative given at boundary 1
!
 102  continue
      wk(2,2) = wk(3,2)-wk(2,2)-wk(2,3)*w(1)
      wk(2,3) = 0.
      wk(2,1) = wk(2,1)+wk(3,1)
      i1 = 2
      nn = nn-1
      go to 106
!
!     first derivative given at boundary 1
!
 103  continue
      wk(1,2) = wk(2,2)-w(1)
      wk(2,2) = wk(3,2)-wk(2,2)
      wk(1,3) = 0.
      wk(1,1) = wk(2,1)
      wk(2,1) = wk(2,1)+wk(3,1)
      i1 = 1
      go to 106
!
!     periodic boundary condition
!
 104  continue
      y2 = wk(2,2)
      b2 = wk(2,1)
      wk(2,2) = wk(3,2)-wk(2,2)
      wk(2,1) = wk(3,1)+wk(2,1)
      i1 = 2
      nn = nn-1
      go to 106
!
!     first derivative at boundary 1 from 4 point interpolation.
!
 105  continue
      a12 = x(1)-x(2)
      a13 = x(1)-x(3)
      a14 = x(1)-x(4)
      a23 = x(2)-x(3)
      a24 = x(2)-x(4)
      a34 = x(3)-x(4)
      j1 = 1
      j2 = j1+int
      j3 = j2+int
      j4 = j3+int
      w(1)    = (1./a12+1./a13+1./a14)*f(j1)- &
        a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)- &
        a12*a13/(a14*a24*a34)*f(j4)
      go to 103
!     compute tridiagonal arrays
 106  continue
      i2 = n-2
      do 107 i=3,i2
        wk(i,2) = wk(i+1,2)-wk(i,2)
        wk(i,1) = wk(i+1,1)+wk(i,1)
 107  continue
!
!     apply boundary conditions at boundary 2.
!
      in = ii(n)
      go to (108,109,110,111),ml
!
!     second derivative given at boundary 2.
!
 108  continue
      wk(n-1,2) = wk(n,2)-wk(n-1,2)-wk(n,3)*w(in)
      wk(n,3) = 0.
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      nn = nn-1
      go to 112
!
!     first derivative given at boundary 2.
!
 109  continue
      wk(n-1,2) = wk(n,2)-wk(n-1,2)
      wk(n,2) = -wk(n,2)+w(in)
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      wk(1,4) = 0.
      go to 112
!
!     periodic boundary condition
!
 110  continue
      wk(n-1,2) = wk(n,2)-wk(n-1,2)
      wk(n,2) = y2-wk(n,2)
      wk(n-1,1) = wk(n-1,1)+wk(n,1)
      wk(n,1) = wk(n,1)+b2
      wk(1,4) = wk(2,3)
      go to 112
!
!     first derivative at boundary 2 from 4 point interpolation.
!
 111  continue
      a12 = x(n)-x(n-1)
      a13 = x(n)-x(n-2)
      a14 = x(n)-x(n-3)
      a23 = x(n-1)-x(n-2)
      a24 = x(n-1)-x(n-3)
      a34 = x(n-2)-x(n-3)
      j1 = in
      j2 = j1-int
      j3 = j2-int
      j4 = j3-int
      w(in)   = (1./a12+1./a13+1./a14)*f(j1)- &
        a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)- &
        a12*a13/(a14*a24*a34)*f(j4)
      go to 109
 112  continue
      iw1 = ii(i1)
      call trip (nn,wk(i1,3),wk(i1,1),wk(i1+1,3),wk(i1,2),w(iw1),int)
      go to (114,114,113,114),mk
 113  continue
      w(1) = w(in)
 114  continue
!dir$ bounds
      return
      end subroutine coeff1

      subroutine coeff2 (nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ibd,wk)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!BH070805      dimension       x(nx)       ,y(ny)       ,f(idm,ny)  ,fxx(idm,ny),
!BH070805     1  fyy(idm,ny) ,fxxyy(idm,ny)           ,ibd(4)     ,
!BH070805     2  iloc(2)    ,jloc(2)
      dimension       x(*)       ,y(*)       ,f(idm,*)  ,fxx(idm,*), &
        fyy(idm,*) ,fxxyy(idm,*)           ,ibd(4)     , &
        iloc(2)    ,jloc(2)      ,wk(*)
      logical q8q4
      save q8q4
      save iloc,jloc
      data q8q4 /.true./
      data iloc(1),iloc(2),jloc(1),jloc(2)/1,1,4,4/
!     the following call is for gathering statistics on library use at ncar
      if (q8q4) then
        q8q4 = .false.
      endif
!
!     compute fxx
!
      do 101 j=1,ny
        call coeff1 (nx,x,f(1,j),fxx(1,j),ibd(1),1,wk)
 101  continue
!
!     compute fyy
!
      do 102 i=1,nx
        call coeff1 (ny,y,f(i,1),fyy(i,1),ibd(3),idm,wk)
 102  continue
!
!     check for periodic boundary condition in both directions
!
      if (ibd(1) .eq. 3) go to 103
      if (ibd(3) .eq. 3) go to 105
!
!     calculate fxxyy along left and right boundaries
!
      call coeff1 (ny,y,fxx(1,1),fxxyy(1,1),jloc,idm,wk)
      call coeff1 (ny,y,fxx(nx,1),fxxyy(nx,1),jloc,idm,wk)
      go to 106
 103  continue
!
!     periodic in x direction . calculate fxxyy along lower and upper
!     boundaries.
!
      call coeff1 (nx,x,fyy(1,1),fxxyy(1,1),ibd(1),1,wk)
      call coeff1 (nx,x,fyy(1,ny),fxxyy(1,ny),ibd(1),1,wk)
!
!     calculate remaining fxxyy
!
      do 104 i=1,nx
        call coeff1 (ny,y,fxx(i,1),fxxyy(i,1),iloc,idm,wk)
 104  continue
      go to 108
 105  continue
!
!     periodic in y direction. calculate fxxyy along left and right
!     boundaries.
!
      call coeff1 (ny,y,fxx(1,1),fxxyy(1,1),ibd(3),idm,wk)
      call coeff1 (ny,y,fxx(nx,1),fxxyy(nx,1),ibd(3),idm,wk)
 106  continue
!
!     calculate remaining fxxyy
!
      do 107 j=1,ny
        call coeff1 (nx,x,fyy(1,j),fxxyy(1,j),iloc,1,wk)
 107  continue
 108  continue
      return
      end subroutine coeff2

      subroutine intrp (n,x,f,w,y,i,int,tab,itab)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      dimension       x(i+1)    ,f(i*int+1)    ,w(i*int+1)  ,tab(3) &
        ,itab(3)
!
!     arithmetic statement function used to locate entries in f and w arrays
!
      ii(index)=(index-1)*int+1
!
!     perform interpolation or extrapolation
!
      flk = x(i+1)-x(i)
      flp = x(i+1)-y
      fl0 = y-x(i)
      i0 = ii(i)
      ip = i0+int
      if (itab(1) .ne. 0) go to 101
      go to 102
 101  continue
!
!     calculate f(y)
!
      a = (w(i0)*flp**3+w(ip)*fl0**3)/(6.*flk)
      b = (f(ip)/flk-w(ip)*flk/6.)*fl0
      c = (f(i0)/flk-w(i0)*flk/6.)*flp
      !-YuP  Note: If w==0 (set all 2nd derivatives to zero),
      !-YuP  then a+b+c = f(i) + [y-x(i)]*[f(i+1)-f(i)]/[x(i+1)-x(i)]
      !-YuP  which is just a linear interpolation.
      tab(1) = a+b+c
 102  continue
      if (itab(2) .ne. 0) go to 103
      go to 104
 103  continue
!
!     calculate first derivative at y
!
      a = (w(ip)*fl0**2-w(i0)*flp**2)/(2.*flk)
      b = (f(ip)-f(i0))/flk
      c = (w(i0)-w(ip))*flk/6.
      tab(2) = a+b+c
 104  continue
      if (itab(3) .ne. 0) go to 105
      go to 106
 105  continue
!
!     calculate second derivative at y
!
      tab(3) = (w(i0)*flp+w(ip)*fl0)/flk
 106  continue
      return
      end subroutine intrp

      subroutine search (xbar,x,n,i)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      dimension       x(n)
      save b
      data b/.69314718/
!
!     if xbar is outside range of x table extrapolate
!
      if (xbar .gt. x(2)) go to 101
      i = 1
      return
 101  continue
      if (xbar .lt. x(n-1)) go to 102
      i = n-1
      return
 102  continue
!
!     find maximum power of two less than n
!
      m = int((alog(float(n)))/b)
      i = 2**m
      if (i .ge. n) i = i/2
      k = i
      nm1 = n-1
!
!     conduct binary search.
!
 103  continue
      k = k/2
      if (xbar .ge. x(i)) go to 104
      i = i-k
      go to 103
 104  continue
      if (xbar .le. x(i+1)) return
      i = min0(i+k,nm1)
      go to 103
      end subroutine search

      subroutine terp1 (n,x,f,w,y,int,tab,itab)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
      dimension       x(n)       ,f(n*int)       ,w(n*int)    ,tab(3), &
        itab(3)
!     the following call is for gathering statistics on library use at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
!
!     perform search
!
      call search (y,x,n,i)
!
!     carry out interpolation (or extrapolation)
!
      call intrp (n,x,f,w,y,i,int,tab,itab)
      return
      end subroutine terp1


      real(c_double) function terp2 &
           (xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
      dimension       x(nx)      ,y(ny)      ,f(idm,ny)  ,fxx(idm,ny), &
        fyy(idm,ny) ,fxxyy(idm,ny)           ,ff(2)      , &
        ww(2)      ,tab(3)     ,itab(3)
!     the following call is for gathering statistics
!       on library use at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
!
!     search in x and y arrays.
!
      call search (xb,x,nx,i)
      call search (yb,y,ny,j)
!
!     interpolate in x direction
!
      do 101 i1=1,3
        itab(i1) = 0
 101  continue
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
        jj = j+j1-1
        call intrp (n,x,f(1,jj),fxx(1,jj),xb,i,1,tab,itab)
        ff(j1) = tab(i1)
        call intrp (n,x,fyy(1,jj),fxxyy(1,jj),xb,i,1,tab,itab)
        ww(j1) = tab(i1)
 102  continue
!
!     interpolate in y direction
!
      do 103 j1=1,3
        itab(j1) = 0
 103  continue
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(j),ff,ww,yb,1,1,tab,itab)
      terp2 = tab(j1)
      return
!
!     revision history---
!
!     june 1977        replaced non-standard statement functions and
!     subscripts to enhance portability.
!
!     january 1978     deleted references to the  *cosy  cards, moved
!     the revision histories to appear before the
!     final end card, and moved the initial comment
!     cards to appear after the first subroutine card
!     and changed  itab  from logical to integer in
!     subroutine intrp and corrected problem with
!     version numbers in one statistics call
!-----------------------------------------------------------------------
      end function terp2

      subroutine searche (xbar,x,n,i,dx)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      dimension       x(n)
      save b
      data b/.69314718/
!
!     if xbar is outside range of x table extrapolate
!
      if (xbar .gt. x(2)) go to 101
      i = 1
      return
 101  continue
      if (xbar .lt. x(n-1)) go to 102
      i = n-1
      return
 102  continue
!..................................................................
!     This version knows data is evenly spaced with spacing dx
!..................................................................

      i=(xbar-x(1))/dx+1

      return
      end subroutine searche

      subroutine terp1e (n,x,f,w,y,int,tab,itab,dx)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
      dimension       x(n)       ,f(n*int)   ,w(n*int)    ,tab(3)     , &
        itab(3)
!     the following call is for gathering statistics on library use at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
!
!     perform search
!
      call searche (y,x,n,i,dx)
!
!     carry out interpolation (or extrapolation)
!
      call intrp (n,x,f,w,y,i,int,tab,itab)
      return
      end subroutine terp1e


!..................................................................
!     This version knows data is evenly spaced with spacing dx
!..................................................................
      real(c_double) function t2 &
           (dx,dy,xb,yb,nx,x,ny,y,f,fxx,fyy,fxxyy,idm,ixd,iyd)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
      dimension       x(nx)     ,y(ny)     ,f(idm,ny)  ,fxx(idm,ny) , &
        fyy(idm,ny)  ,fxxyy(idm,ny)        ,ff(2)      , &
        ww(2)      ,tab(3)     ,itab(3)
!     the following call is for gathering statistics
!       on library use at ncar
      logical q8q4
      save q8q4
      data q8q4 /.true./
      if (q8q4) then
        q8q4 = .false.
      endif
!
!     search in x and y arrays.
!
      call searche (xb,x,nx,i,dx)
      call searche (yb,y,ny,j,dy)
!
!     interpolate in x direction
!
      do 101 i1=1,3
        itab(i1) = 0
 101  continue
      i1 = ixd+1
      itab(i1) = 1
      do 102 j1=1,2
        jj = j+j1-1
        call intrp (n,x,f(1,jj),fxx(1,jj),xb,i,1,tab,itab)
        ff(j1) = tab(i1)
        call intrp (n,x,fyy(1,jj),fxxyy(1,jj),xb,i,1,tab,itab)
        ww(j1) = tab(i1)
 102  continue
!
!     interpolate in y direction
!
      do 103 j1=1,3
        itab(j1) = 0
 103  continue
      j1 = iyd+1
      itab(j1) = 1
      call intrp (2,y(j),ff,ww,yb,1,1,tab,itab)
      t2 = tab(j1)
      return
!
!     revision history---
!
!     june 1977        replaced non-standard statement functions and
!     subscripts to enhance portability.
!
!     january 1978     deleted references to the  *cosy  cards, moved
!     the revision histories to appear before the
!     final end card, and moved the initial comment
!     cards to appear after the first subroutine card
!     and changed  itab  from logical to integer in
!     subroutine intrp and corrected problem with
!     version numbers in one statistics call
!-----------------------------------------------------------------------
      end function t2

      subroutine trip (n,a,b,c,y,z,int)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      dimension       a(n)       ,b(n)       ,c(n)       ,y(n)       , &
        z(n*int)
!
!     arithmetic statement function used to locate entries in array z.
!
      ii(index)=(index-1)*int+1
!
!     gaussian elimination
!
      bn = b(n)
      yn = y(n)
      v = c(n)
      y(1) = y(1)/b(1)
      a(1) = a(1)/b(1)
      b(1) = c(1)/b(1)
      nm2 = n-2
      do 101 j=2,nm2
        den = b(j)-a(j)*b(j-1)
        b(j) = c(j)/den
        y(j) = (y(j)-a(j)*y(j-1))/den
        a(j) = -a(j)*a(j-1)/den
        bn = bn-v*a(j-1)
        yn = yn-v*y(j-1)
        v = -v*b(j-1)
 101  continue
      den = b(n-1)-a(n-1)*b(n-2)
      b(n-1) = (c(n-1)-a(n-1)*a(n-2))/den
      y(n-1) = (y(n-1)-a(n-1)*y(n-2))/den
      bn = bn-v*a(n-2)
      yn = yn-v*y(n-2)
      v = a(n)-v*b(n-2)
!     back substitution
      iin = ii(n)
      z(iin) = (yn-v*y(n-1))/(bn-v*b(n-1))
      iin2 = ii(n-1)
      z(iin2) = y(n-1)-b(n-1)*z(iin)
      nm1 = n-1
      in = ii(n)
      do 102 j=2,nm1
        k = n-j
        ik = ii(k)
        ikt = ik+int
      z(ik) = y(k)-b(k)*z(ikt)-a(k)*z(in)
      if(j.eq.(-n)) write(*,140) b(k), y(k), a(k), z(ikt),z(ik)
 102  continue
 140  format( '         cj_debug: trip', 5e16.8)
      return
      end subroutine trip

!************************************************************************
!     END OF SPLINES
!************************************************************************
!
!

!######date01jan1984     copyright ukaea, harwell.
!######aliasim01ad
      real(c_double) function im01ad(la,a,inc)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      real(c_double) a(la*inc)
      zero=0.d0
      kount = 0
      do 100 k=1,la
        ipos = (k-1)*inc + 1
        if (a(ipos).eq.zero) kount = kount + 1
 100  continue
      im01ad = kount
      return
      end function im01ad

!
!
!**********************************************************
      subroutine zzbeslri(x,nb,ize,b,ncalc)
      implicit integer (i-n), real(c_double) (a-h,o-z) ! DOUBLE NOW !
!mnt  this routine calculates Bessel functions I and J of real
!mnt  argument and integer order.


!mnt  explanation of variables in the calling sequence

!mnt  x (was single) DOUBLE precision real argument for which i,s or j,s
!mnt  are to be calculated.  if i,s are to be calculated,
!mnt  abs(x) must be less than exparg (which see below).

!mnt  nb      integer type.  1+highest order to be calculated.
!mnt  it must be positive.

!mnt  ize     integer type.  zero if j,s are to be calculated,
!mnt  1 if i,s are to be calculated.

!mnt  b (was single) DOUBLE precision vector of length nb, need not be
!mnt  initialized by user.  if the routine terminates
!mnt  normally (ncalc=nb), it returns j(or i) -sub-zero
!mnt  through j(or i) -sub-nb-minus-one of x in this
!mnt  vector.

!mnt  ncalc  integer type, need not be initialized by user.
!mnt  before using the results, the user should check that
!mnt  ncalc=nb.  i.e. all orders have been calculated to
!mnt  the desired accuracy.  see error returns below.


!mnt  explanation of machine-dependent constants

!mnt  nsig    decimal significance desired.  should be set to
!mnt  ifix(alog10(2)*nbit+1), where nbit is the number of
!mnt  bits in the mantissa of a double precision variable.
!mnt  setting nsig lower will result in decreased accuracy
!mnt  while setting nsig higher will increase cpu time
!mnt  without increasing accuracy.  the truncation error
!mnt  is limited to t=.5*10**-nsig for j,s of order less
!mnt  than argument, and to a relative error of t for
!mnt  i,s and the j,s.

!mnt  nten    largest integer k such that 10**k is machine-
!mnt  representable in single precision.

!mnt  largex  upper limit on the magnitude of x.  bear in mind
!mnt  that if abs(x)=n, then at least n iterations of the
!mnt  backward recursion will be executed.

!mnt  exparg  largest single precision argument that the library
!mnt  exp  routine can handle.


!mnt  error returns


!mnt  let g denote either i or j.
!mnt  in case of an error, ncalc.ne.nb, and not all g,s
!mnt  are calculated to the desired accuracy.
!mnt  if ncalc.lt.0, an argument is out of range.  nb.le.0
!mnt  or ize is neither 0 nor 1 or ize=1 and abs(x).ge.exparg.
!mnt  in this case, the b-vector is not calculated, and ncalc
!mnt  is set to min0(nb,0)-1 so ncalc.ne.nb.

!mnt  nb.gt.ncalc.gt.0 will occur if nb.gt.magx and abs(g-
!mnt  sub-nb-of-x/g-sub-magx-of-x).lt.10.**(nten/2), i.e.
!mnt  nb is much greater than magx.  in this case, b(n) is calcu
!mnt  lated to the desired accuracy for n.le.ncalc, but for
!mnt  ncalc.lt.n.le.nb, precision is lost.  if n.gt.ncalc and
!mnt  abs(b(ncalc)/b(n)).eq.10.**-k, then only the first nsig-k
!mnt  significant figures of b(n) may be trusted.  if the user
!mnt  wishes to calculate b(n) to higher accuracy, he should use
!mnt  an asymptotic formula for large order.



      dimension b(nb)
      save nsig,nten,largex,exparg
!BH150620      data nsig,nten,largex,exparg/14,307,100000,7.e2/
!BH150620:  Change error tolerance, for very high harmonic cases.
      data nsig,nten,largex,exparg/14,307,100000,7.e2/
      one=1.d0
      ten=10.d0
      tempa= abs(x)
!mnt  magx=ifix(sngl(tempa))
      magx=tempa
      if (nb.gt.0.and.magx.le.largex.and.(ize.eq.0.or. &
        (ize.eq.1.and.tempa.le.exparg)))goto 10
!mnt  error return -- x,nb,or ize is out of range.
      ncalc=min0(nb,0)-1
      return
 10   sign=1-2*ize
      ncalc=nb
!mnt  use 2-term ascending series for small x
!990131      if (tempa**4.lt..1e0**nsig) goto 230
      if (tempa**4.lt.(.1e0*one)**nsig) goto 230
!mnt  initialize the calculation of p,s
      nbmx=nb-magx
      n=magx+1
      plast=1.d0
      p=(n+n)/tempa
!mnt  calculate general significance test
!990131      test=2.e0*1.e1**nsig
      test=2.e0*ten**nsig
      if (ize.eq.1.and.2*magx.gt.5*nsig) test= sqrt(test*p)
      one585=1.585
!990131      if (ize.eq.1.and.2*magx.le.5*nsig) test=test/1.585**magx
      if (ize.eq.1.and.2*magx.le.5*nsig) test=test/one585**magx
      m=0
      if (nbmx.lt.3) goto 30
!mnt  calculate p,s until n=nb-1.  check for possible overflow.
!990131      tover=1.e1**(nten-nsig)
      tover=ten**(nten-nsig)
      nstart=magx+2
      nend=nb-1
      do 20 n=nstart,nend
        pold=plast
        plast=p
        p=(n+n)*plast/tempa-sign*pold
        if (p-tover) 20,20,40
 20   continue
      n=nend
!mnt  calculate special significance test for nbmx.gt.2.
!990131      test=amax1(test, sqrt(plast*1.e1**nsig)* sqrt(2.e0*p))
      test=max(test, sqrt(plast*ten**nsig)* sqrt(2.e0*p))
!mnt  calculate p,s until significance test passes.
 30   n=n+1
      pold=plast
      plast=p
      p=(n+n)*plast/tempa-sign*pold
      if (p.lt.test) goto 30
      if (ize.eq.1.or.m.eq.1) goto 90
!mnt  for j*s, a strong variant of the test is necessary.
!mnt  calculate it, and calculate p*s until this test is passed.
      m=1
      tempb=p/plast
      tempc=(n+1)/tempa
      if (tempb+1.e0/tempb.gt.2.e0*tempc) tempb=tempc+ sqrt &
        (tempc**2-1.e0)
      test=test/ sqrt(tempb-1.e0/tempb)
      if (p-test) 30,90,90
!mnt  to avoid overflow, divide p*s by tover.  calculate p*s
!mnt  until abs(p).gt.1).
 40   tover=ten**nten
      p=p/tover
      plast=plast/tover
      psave=p
      psavel=plast
      nstart=n+1
 50   n=n+1
      pold=plast
      plast=p
      p=(n+n)*plast/tempa-sign*pold
      if (p.le.1.e0) goto 50
      tempb=(n+n)/tempa
      if (ize.eq.1) goto 60
      tempc=.5e0*tempb
      tempb=plast/pold
      if (tempb+1.e0/tempb.gt.2.e0*tempc) tempb = tempc + sqrt(tempc**2-one)
!990131     1  (tempc**2-1.e0)
!mnt  calculate backward test, and find ncalc, the highest n
!mnt  such that the test is passed.
!990131 60   test=.5e0*pold*plast*(1.e0-1.e0/tempb**2)/1.e1**nsig
 60   test=.5e0*pold*plast*(1.e0-1.e0/tempb**2)/ten**nsig
      p=plast*tover
      n=n-1
      nend=min0(nb,n)
      do 70 ncalc=nstart,nend
        pold=psavel
        psavel=psave
        psave=(n+n)*psavel/tempa-sign*pold
        if (psave*psavel-test) 70,70,80
 70   continue
      ncalc=nend+1
 80   ncalc=ncalc-1
!mnt  the sum b(1)+2b(3)+2b(5)... is used to normalize.  m, the
!mnt  coefficient of b(n), is initialized to 2 or 0.
 90   n=n+1
      m=2*n-4*(n/2)
!mnt  initialize the backward recursion and the normalization
!mnt  sum.
      tempb=0.e0
      tempa=1.e0/p
      sum=m*tempa
      nend=n-nb
      if (nend) 140,120,100
!mnt  recur backward via difference equation.  calculating (but
!mnt  not storing) b(n), until n=nb.
 100  do 110 l=1,nend
        n=n-1
        tempc=tempb
        tempb=tempa
        tempa=(n+n)*tempb/x-sign*tempc
        m=2-m
 110  sum=sum+m*tempa
!mnt  store b(nb)
 120  b(n)=tempa
      if (nb.gt.1) goto 130
!mnt  nb=1.  since 2*tempa was added to the sum, tempa must be
!mnt  subtracted.
      sum=sum-tempa
      goto 200
!mnt  calculate and store b(nb-1)
 130  n=n-1
      b(n)=(n+n)*tempa/x-sign*tempb
      if (n.eq.1) goto 190
      m=2-m
      sum=sum+m*b(n)
      goto 160
!mnt  n.lt.nb, so store b(n) and set higher orders to zero
 140  b(n)=tempa
      nend=-nend
      do 150 l=1,nend
 150  b(n+l)=0.e0
 160  nend=n-2
      if (nend.eq.0) goto 180
!mnt  calculate via difference equation and store b(n),
!mnt  until n=2.
      do 170 l=1,nend
        n=n-1
        b(n)=(n+n)*b(n+1)/x-sign*b(n+2)
        m=2-m
 170  sum=sum+m*b(n)
!mnt  calculate b(1)
 180  b(1)=2.e0*b(2)/x-sign*b(3)
 190  sum=sum+b(1)
!mnt  normalize--if ize=1, divide sum by cosh(x).  divide all
!mnt  b(n) by sum.
 200  if (ize.eq.0) goto 210
      tempa= exp( abs(x))
      sum=2.e0*sum/(tempa+1.e0/tempa)
 210  do 220 n=1,nb
 220  b(n)=b(n)/sum
      return
!mnt  two-term ascending series for small x
 230  tempa=1.e0
      tempb=-.25e0*x*x*sign
      b(1)=1.e0+tempb
      if (nb.eq.1) goto 250
      do 240 n=2,nb
        tempa=tempa*x/(n+n-2)
 240  b(n)=tempa*(1.e0+tempb/n)
 250  return
      end subroutine zzbeslri

!***********************************************************************
!$$$      subroutine zzechk(nchars,narray)
!$$$      implicit integer (i-n), real(c_double) (a-h,o-z)
!$$$      save
!$$$
!$$$cmnt  abstract
!$$$cmnt  zzechk is a companion routine of zzrchk.  it is called
!$$$cmnt  just like zzrchk, and messages from it may be suppressed
!$$$cmnt  by an appropriate call to zzxset.  it differs from zzrchk
!$$$cmnt  in that each call to zzechk will produce no more than one
!$$$cmnt  printed message, regardless of how many times that call is
!$$$cmnt  executed, and zzechk never terminates execution.
!$$$cmnt  its purpose is to provide one-time-only informative
!$$$cmnt  diagnostics.
!$$$
!$$$cmnt  description of arguments
!$$$cmnt  nchars - number of characters in the message.
!$$$cmnt  if negated, the message will be printed (once) even
!$$$cmnt  if nfatal has been set to 0 (see zzxset).
!$$$cmnt  narray - same as in zzrchk
!$$$
!$$$
!$$$
!$$$cmnt  zzechk uses subroutines zzrget, zzrprt, zzxset, zzstgt
!$$$cmnt  compile decks zzrchk
!$$$
!$$$      dimension narray(14)
!$$$      data nflag/4h.$,*/
!$$$      if (narray(1).eq.nflag) return
!$$$      call zzrget(nf,nt)
!$$$      if ((nf.eq.0).and.(nchars.gt.0)) return
!$$$      call zzrprt (59,59hthe following informative diagnostic is
!$$$     1 printed  only once.)
!$$$      call zzrprt(iabs(nchars),narray)
!$$$      if (nf.gt.0) nf = nf-1
!$$$
!$$$cmnt  calculates psi!!, the second derivative of psi with
!$$$cmnt  respect to arc length, as a functionn of psi:  d(dpsi/dz)/dz
!$$$
!$$$      call zzxset(nf,nt)
!$$$      narray(1) = nflag
!$$$      end

!***********************************************************************
      subroutine zzrchk(nchars,narray)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!mnt  sandia mathematical program library
!mnt  applied mathematics division 2642
!mnt  sandia laboratories
!mnt  albuquerque, new mexico 87115

!mnt  simplified version for stand-alone use.     april 1977

!mnt  abstract
!mnt  the routines zzrchk, zzxset, and zzrget together provide
!mnt  a uniform method with several options for the processing
!mnt  of diagnostics and warning messages which originate
!mnt  in the mathematical program library routines.
!mnt  zzrchk is the central routine, which actually processes
!mnt  messages.

!mnt  description of arguments
!mnt  nchars - number of characters in hollerith message.
!mnt  if nchars is negated, zzrchk will unconditionally
!mnt  print the message and stop execution.  otherwise,
!mnt  the behavior of zzrchk may be controlled by
!mnt  an appropriate call to zzxset.
!mnt  narray - name of array or variable containing the message,
!mnt  or else a literal hollerith constant containing
!mnt  the message.  by convention, all messages should
!mnt  begin with *in subnam, ...*, where subnam is the
!mnt  name of the routine calling zzrchk.

!mnt  examples
!mnt  1. to allow control by calling zzxset, use
!mnt  call zzrchk(30,30hin quad, invalid value of err.)
!mnt  2. to unconditionally print a message and stop execution, use
!mnt  call zzrchk(-30,30hin quad, invalid value of err.)



!mnt  zzrchk uses subroutines zzrget, zzrprt, zzxset, zzstgt
!mnt  compile decks zzrchk

      dimension narray(14)

      call zzrget(nf,nt)
!mnt  if zzrchk was called with negative character count, set fatal flag
      if (nchars.lt.0) nf = -1
!mnt  if messages are to be suppressed, return
      if (nf.eq.0) return
!mnt  if character count is invalid, stop
      if (nchars.eq.0) print 10010
      if (nchars.eq.0) stop 'zcunix: zzrchk'
!mnt  print message
      call zzrprt(iabs(nchars),narray)
!mnt  if last message, say so
      if (nf.eq.1) print 10020
!mnt  print trace-back if asked to
!mnt  if ((nt.gt.0).or.(nf.lt.0)) call system routine for traceback
!mnt  decrement message count
      if (nf.gt.0) nf = nf-1
      call zzxset(nf,nt)
!mnt  if all is well, return
      if (nf.ge.0) return
!mnt  if this message is suppressable by an zzxset call,
!mnt  then explain zzxset usage.
      if (nchars.gt.0) print 10030
      print 10040
      stop
10010 format(/31h zzrchk was called incorrectly.)
10020 format (30h zzrchk message limit reached.)
10030 format (/13h *** note ***, &
        'to make the error message printed above be nonfatal, '&
        'or to suppress the message completely, '&
        'insert an appropriate zzxset, '&
        'at the start of your program.' &
        'for example, to print up to 10 nonfatal warning messages, '&
        '  use     call zzxset(10,0)    ')
10040 format (/28h program abort due to error.)
      end subroutine zzrchk

!***********************************************************************
      subroutine zzrget(nfatal,ntrace)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!mnt  abstract
!mnt  zzrget is a companion routine to subroutine zzrchk.
!mnt  zzrget assigns to nfatal and ntrace respectively the values
!mnt  of nf and nt in common block mlblk0 thereby ascertaining the
!mnt  state of the options which control the execution of zzrchk.

!mnt  description of arguments
!mnt  both arguments are output arguments of data type integer.
!mnt  nfatal - current value of nf (see description of zzxset.)
!mnt  ntrace - current value of nt (see description of zzxset.)

      call zzstgt(1,nfatal,ntrace)
      return
      end subroutine zzrget

!***********************************************************************
      subroutine zzrprt(nchars,narray)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!mnt  utility routine to simply print the hollerith message in narray,
!mnt  whose length is nchars characters.

      dimension narray(14)

!mnt  note - nch must be the number of hollerith characters stored
!mnt  per word.  if nch is changed, format 1 must also be
!mnt  changed correspondingly.

      nch = 10
!mnt  for line printers, use
!mnt  for data terminals, use
!mnt  1 format (1x,7a10)
      nwords = (nchars+nch-1)/nch
      print 10010,(narray(i),i=1,nwords)
      return
10010 format (1x,13a10)
      end subroutine zzrprt

!***********************************************************************
      subroutine zzstgt(k,nfatal,ntrace)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!mnt  this routine is a slave to zzrget and errset which keeps
!mnt  the flags as local variables.

!mnt  *** if local variables are not normally retained between
!mnt  calls on this system, the variables lnf and lnt can be
!mnt  placed in a common block and preset to the following
!mnt  values in the main program.

      data lnf/-1/,lnt/0/
      if (k.le.0) lnf = nfatal
      if (k.le.0) lnt = ntrace
      if (k.gt.0) nfatal = lnf
      if (k.gt.0) ntrace = lnt
      return
      end subroutine zzstgt

!***********************************************************************
      subroutine zzxset(nfatal,ntrace)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!mnt  abstract
!mnt  zzxset is a companion routine to subroutine zzrchk.
!mnt  zzxset assigns the values of nfatal and ntrace respectively
!mnt  to nf and nt in common block mlblk0 thereby specifying the
!mnt  state of the options which control the execution of zzrchk.

!mnt  description of arguments
!mnt  both arguments are input arguments of data type integer.
!mnt  nfatal - is a fatal-error / message-limit flag. a negative
!mnt  value denotes that detected difficulties are to be
!mnt  treated as fatal errors.  nonnegative means nonfatal.
!mnt  a nonnegative value is the maximum number of nonfatal
!mnt  warning messages which will be printed by zzrchk,
!mnt  after which nonfatal messages will not be printed.
!mnt  (default value is -1.)
!mnt  ntrace - .ge.1 will cause a trace-back to be given,
!mnt  if this feature is implemented on this system.
!mnt  .le.0 will suppress any trace-back, except for
!mnt  cases when execution is terminated.
!mnt  (default value is 0.)

!mnt  *note* -- some calls to zzrchk will cause unconditional
!mnt  termination of execution.  zzxset has no effect on such calls.

!mnt  examples
!mnt  1. to print up to 100 messages as nonfatal warnings use
!mnt  call zzxset(100,0)
!mnt  2. to suppress all mathlib warning messages use
!mnt  call zzxset(0,0)



!mnt  zzxset uses subroutines zzstgt
!mnt  compile decks zzrchk

      call zzstgt(0,nfatal,ntrace)
      return
      end subroutine zzxset



!     Following two subroutines from ONETWO
      subroutine allocate_error(var,myid,istat)
!------------------------------------------------------------------
      character *(*) var
      integer istat,myid
      write(*,1)var,istat,myid
 1    format(2x,"Memory Allocation error encountered",/, &
             2x,"Unable to allocate ",a,/, &
             2x,"status =",i5)
      istat =0 !reset for next case
      return
      end subroutine allocate_error



      subroutine deallocate_error(var,myid,istat)
!---------------------------------------------------------------------
      character *(*) var
      integer istat, myid
      write(*,1)var,istat,myid
 1    format(2x,"Memory DE-Allocation error encountered",/, &
             2x,"Unable to deallocate ",a,/, &
             2x,"status =",i5, " process rank =",i5)
      istat =0 !reset for next case
      return
      end subroutine deallocate_error

end module zcunix_mod
