!     advnce.f90
!..................................................................
!     advnce contains the statement functions utilized by splitting
!     or implicit time advancement routines and by their
!     diagnostic routines.
!..................................................................

module advnce_mod

  !---BEGIN USE

  use bsu_mod, only : bsu
  use bsl_mod, only : bsl
  use comm_mod
  use r8subs_mod, only : cvmgt

  !---END USE

  implicit none
  save

  integer, public, target ::  advnce_k

contains

  !..................................................................
  !     Define some integration coefficients.
  !..................................................................

  real(c_double) function qz(j)
    integer :: j
    qz = 1./cint2(j)
  end function qz

  real(c_double) function ry(i,j)
    integer :: i,j
    ry = dx(j)*twopi/(cynt2(i,l_)*cint2(j))
  end function ry

  real(c_double)function cl(i,j)
    integer :: i,j
    cl = .25*vptb(itl,lr_)/vptb(i,lr_)*dc(i,j)
  end function cl

  real(c_double) function r2y(j)
    integer :: j
    r2y = ry(itl,j)*.5
  end function r2y

  real(c_double) function dithta(i,j,l_)
    integer :: i,j, l_
    dithta = 0.5-0.5*(1/(i+1))+0.5*(1/(iy_(l_)+1-i))
  end function dithta

  !###########################################################
  !     Statement functions for explicit time advancement (splitting)
  !     follow.
  !###########################################################

  !..................................................................
  !     Using Chang-Cooper weights dj(i,j,k,l_) and di(i,j,k,l_) define
  !     the weighted averages of the distribution function:
  !     f1* - before split;    f2* - after split
  !..................................................................

  !     f1j(i,j)=temp1(i,j+1)*(1.-dj(i,j,k,l_))+temp1(i,j)*dj(i,j,k,l_)
  !     f2j(i,j)=temp2(i,j+1)*(1.-dj(i,j,k,l_))+temp2(i,j)*dj(i,j,k,l_)
  !     f1i(i,j)=temp1(i+1,j)*(1.-di(i,j,k,l_))+temp1(i,j)*di(i,j,k,l_)
  !     f2i(i,j)=temp2(i+1,j)*(1.-di(i,j,k,l_))+temp2(i,j)*di(i,j,k,l_)
  ! YuP-101228: same as above, but re-arranged to have one '*'
  real(c_double) function f1j(i,j)
    integer :: i,j
    f1j = temp1(i,j+1) + (temp1(i,j)-temp1(i,j+1))*dj(i,j,advnce_k,l_)
  end function f1j

  real(c_double) function f2j(i,j)
    integer :: i,j
    f2j = temp2(i,j+1) + (temp2(i,j)-temp2(i,j+1))*dj(i,j,advnce_k,l_)
  end function f2j

  real(c_double) function f1i(i,j)
    integer :: i,j
    f1i = temp1(i+1,j) + (temp1(i,j)-temp1(i+1,j))*di(i,j,advnce_k,l_)
  end function f1i

  real(c_double) function f2i(i,j)
    integer :: i,j
    f2i = temp2(i+1,j) + (temp2(i,j)-temp2(i+1,j))*di(i,j,advnce_k,l_)
  end function f2i

  !..................................................................
  !     The first half of the splitting scheme consists of forward and
  !     backward sweeps in the "x" or velocity direction. The forward
  !     sweep requires the coefficients alpx betx gamx and the r.h.s.
  !     delx of the equation:
  !
  !     -alpx(i,j)*f(i,j+1,l_) +betx(i,j)*f(i,j,l_) -gamx(i,j)*f(i,j-1,l_)
  !
  !     =delx(i,j)
  !
  !     Boundary conditions at x=0  (x=xmax)  that is j=1 (j=jx)
  !     automatically force gamx(i,1) (alpx(i,jx)) to be zero.
  !..................................................................

  real(c_double) function alpx(i,j)
    integer :: i,j
    alpx =  (da(i,j)*(1.-dj(i,j,advnce_k,l_))+db(i,j)*exp5(j)) * qz(j)
  end function alpx

  real(c_double) function betx(i,j)
    integer :: i,j
    betx =  (db(i,j)*exp5(j)+db(i,j-1)*exm5(j)-da(i,j) &
         *dj(i,j,advnce_k,l_)+da(i,j-1) &
         *(1.-dj(i,j-1,advnce_k,l_))) &
         * qz(j)-vptb(i,lr_) &
         *cah(i,j)+vptb(i,lr_)*rbgn
  end function betx

  real(c_double) function gamx(i,j)
    integer :: i,j
    gamx =  (-da(i,j-1)*dj(i,j-1,advnce_k,l_)+db(i,j-1)*exm5(j)) * qz(j)
  end function gamx

  !..................................................................
  !     The quantity delx is complicated by the averaging done at the
  !     pass/trapped boundary (i=itl).
  !
  !     Define the appropriate average at the p/t boundary
  !..................................................................

  real(c_double) function cdf(j)
    integer :: j
    cdf =  (cl(itl-1,j)*(f1j(itl,j)-f1j(itl-1,j)) &
         *eyp5(itl-1,l_)+2. &
         *cl(itl+1,j)*(f1j(itl+1,j)-f1j(itl,j))*eyp5(itl,l_) &
         +cl(itu+1,j)*(f1j(itu+1,j)-f1j(itu,j))*eyp5(itu,l_))
  end function cdf

  real(c_double) function delx(i,j)
    use iso_c_binding, only : c_double
    integer :: i,j
    !XXXX
    delx = cvmgt( DBLE(&
         (dc(i,j)*(f1j(i+1,j)-f1j(i-1,j))*0.5*dyi(i,l_) &
         -dc(i,j-1)*(f1j(i+1,j-1)-f1j(i-1,j-1))*0.5*dyi(i,l_)) *qz(j)), &
         DBLE((cdf(j)-cdf(j-1)) * qz(j)*ident(i)), &
         iota(i).ne.itl .and. iota(i).ne.itu) &
         +vptb(i,lr_)*.5*so(i,j) &
         +vptb(i,lr_)*temp1(i,j)*rbgn
  end function delx

  !..................................................................
  !     The sweep (y) requires the coefficients alpy, bety, gamy, and the
  !     r.h.s. dely to the equation:
  !
  !     -alpy(i,j)*f(i+1,j,l_) +bety(i,j)*f(i,j,l_) -gamy(i,j)*f(i-1,j,l_)
  !
  !     = dely(i,j)
  !
  !     Boundary conditions at y=0 and y=pi automatically force gamy(1,j)
  !     and alpy(iy,j) to be equal to zero.
  !..................................................................

  real(c_double) function alpy(i,j)
    integer :: i,j
    alpy =  (dd(i,j)*(1.-di(i,j,advnce_k,l_))+df(i,j)*eyp5(i,l_)) *ry(i,j)
  end function alpy

  real(c_double) function bety(i,j)
    integer :: i,j
    bety =  (-dd(i,j)*di(i,j,advnce_k,l_) &
         +df(i,j)*eyp5(i,l_) &
         +dd(i-1,j)*(1.-di(i-1,j,advnce_k,l_)) &
         +df(i-1,j)*eym5(i,l_)) *ry(i,j) &
         +vptb(i,lr_)*rbgn
  end function bety

  real(c_double) function gamy(i,j)
    integer :: i,j
    gamy =  -ry(i,j)*(dd(i-1,j)*di(i-1,j,advnce_k,l_) - df(i-1,j)*eym5(i,l_))
  end function gamy

  real(c_double) function dely(i,j)
    integer :: i,j
    dely =  ry(i,j)*0.5*dxi(j)*(de(i,j)*(f1i(i,ifp(j))-f1i(i,j-1)) &
         -de(i-1,j)*(f1i(i-1,ifp(j))-f1i(i-1,j-1))) &
         +vptb(i,lr_)*rbgn*temp1(i,j) &
          +vptb(i,lr_)*.5*so(i,j)
  end function dely

  !..................................................................
  !     Define the flux related quantities G and H the are used in
  !     the r.h.s. of the Fokker-Planck equation.
  !..................................................................
  !     -YuP: Moved gfu to diagentr.f, to avoid cvmgt() construct
  !     gfu(i,j)=da(i,j)*f2j(i,j)
  !     1  +db(i,j)*(temp2(i,j+1)-temp2(i,j))*exp5(j)
  !     1  +cvmgt(  dc(i,j)*(f1j(i+1,j)-f1j(i-1,j))*0.5*dyi(i,l_),
  !     1  cdf(j),
  !     1  i .ne. itl  .and. i .ne. itu)

  real(c_double) function hfu(i,j)
    integer :: i,j
    hfu = dd(i,j)*f2i(i,j) &
         +de(i,j)*(f1i(i,ifp(j))-f1i(i,j-1))*0.5*dxi(j) &
         +df(i,j)*(temp2(i+1,j)-temp2(i,j))*eyp5(i,l_)
  end function hfu

  !..................................................................
  !     End splitting scheme time advancement statement functions...
  !..................................................................


  !######################################################
  !     Begin implicit time advancement statement functions..
  !######################################################

  !..................................................................
  !     advnce contains all of the parameters, arrays, and function
  !     statements necessary to create the sparse matrix which represents
  !     the implicit set of Fokker-Planck equations. The matrix is inverted
  !     via Gaussian elimination (White routines ZSGBFA,ZSGBZL)

  !
  !     Below are statement functions which are used to determine
  !     the matrix to be inverted. The i and j are the indices of the
  !     equation considered - i being the theta index, j the velocity
  !     index. On the left hand side the "m" stands for minus, the "p"
  !     for plus and the "u" for upper pass/trapped mesh point itu and
  !     "0" for neutral. The "x" stands for coefficients that hold
  !     everywhere but the pass/trapped boundary and the "t" for
  !     the equation at the pass/trapped boundary. For example,
  !     the equation which represents the mesh point (i,j) involves
  !     a sum of products of the distribution function on the left
  !     hand side. The quantity which multiplies f(i-1,j,l_) is
  !     xm0(i,j).
  !
  !     In terms of coefficients in Killeen et al.(1986) book
  !     within a multiplicative constant:
  !     da(i,j)=A_i,j+1/2      db(i,j)=B_i,j+1/2      dc(i,j)=C_i,j+1/2
  !     dd(i,j)=D_i+1/2,j      de(i,j)=D+i+1/2,j      df(i,j)=F_i+1/2,j
  !     di(i,j)=delta_i+1/2,j  dj(i,j)=delta_i,j+1/2
  !
  !.......................................................................

  real(c_double) function xmm(i,j)
    integer :: i,j
    xmm = (-qz(j)*dc(i,j-1)*dj(i-1,j-1,advnce_k,l_)*dyi(i,l_) &
         -ry(i,j)*de(i-1,j)*di(i-1,j-1,advnce_k,l_)*dxi(j))*.5
  end function xmm

  real(c_double) function x0m(i,j)
    integer :: i,j
    x0m = qz(j)*da(i,j-1)*dj(i,j-1,advnce_k,l_)+ry(i,j)*(de(i,j)* &
          di(i,j-1,advnce_k,l_)-de(i-1,j)*(1.-di(i-1,j-1,advnce_k,l_)))*0.5*dxi(j) &
          -qz(j)*db(i,j-1)*exm5(j)
  end function x0m

  real(c_double) function xpm(i,j)
    integer :: i,j
    xpm = (qz(j)*dc(i,j-1)*dj(i+1,j-1,advnce_k,l_)*dyi(i,l_)+ &
         ry(i,j)*de(i,j)*(1.-di(i,j-1,advnce_k,l_))*dxi(j))*.5
  end function xpm

  real(c_double) function xm0(i,j)
    integer :: i,j
  !  integer, parameter :: l_=1 !XXXXXX
    xm0 = qz(j)*(dc(i,j)*dj(i-1,j,advnce_k,l_) &
         -dc(i,j-1)*(1.-dj(i-1,j-1,advnce_k,l_))) &
         *0.5*dyi(i,l_)+ry(i,j)*(dd(i-1,j)*di(i-1,j,advnce_k,l_) &
         -df(i-1,j)*eym5(i,l_)) &
         +cthta(i,j)*dithta(i-1,j,l_) !Added since 1992
  end function xm0

  real(c_double) function x00(i,j)
    integer :: i,j
    x00 = qz(j)* &
         (-da(i,j)*dj(i,j,advnce_k,l_)+da(i,j-1)*(1.-dj(i,j-1,advnce_k,l_)) &
         +db(i,j)*exp5(j)+db(i,j-1)*exm5(j)) &
         +ry(i,j)*(-dd(i,j)*di(i,j,advnce_k,l_) &
         +dd(i-1,j)*(1.-di(i-1,j,advnce_k,l_))+df(i,j)*eyp5(i,l_)+df(i-1,j) &
         *eym5(i,l_)) &
         -vptb(i,lr_)*(cah(i,j)-1./dtreff) &
         +cthta(i,j)*(1.-dithta(i-1,j,l_)-dithta(i,j,l_)) !Added since 1992
  end function x00

  real(c_double) function xp0(i,j)
    integer :: i,j
    xp0 = qz(j)*(-dc(i,j)*dj(i+1,j,advnce_k,l_) &
         +dc(i,j-1)*(1.-dj(i+1,j-1,advnce_k,l_)))*0.5*dyi(i,l_) &
         -ry(i,j)*(dd(i,j) &
         *(1.-di(i,j,advnce_k,l_))+df(i,j)*eyp5(i,l_)) &
         -cthta(i,j)*(1.-dithta(i,j,l_)) !Added since 1992
  end function xp0

  real(c_double) function xmp(i,j)
    integer :: i,j
    xmp = qz(j)*dc(i,j)*(1.-dj(i-1,j,advnce_k,l_))*.5*dyi(i,l_)+ &
         ry(i,j)*de(i-1,j)*.5*dxi(j)*di(i-1,j+1,advnce_k,l_)
  end function xmp

  real(c_double) function x0p(i,j)
    integer :: i,j
    x0p = qz(j)*(-da(i,j)*(1.-dj(i,j,advnce_k,l_))-db(i,j)*exp5(j)) &
         +ry(i,j)*(-de(i,j)*di(i,j+1,advnce_k,l_) &
         +de(i-1,j)*(1.-di(i-1,j+1,advnce_k,l_)))*0.5*dxi(j)
  end function x0p

  real(c_double) function xpp(i,j)
    integer :: i,j
    xpp = -qz(j)*dc(i,j)*(1.-dj(i+1,j,advnce_k,l_))*0.5*dyi(i,l_) &
         -ry(i,j)*de(i,j)*(1.-di(i,j+1,advnce_k,l_))*0.5*dxi(j)
  end function xpp

  !.......................................................................
  !     z00 is the right hand side of the equation, and holds the
  !     explicit-in-time rhs of the FP difference equations.
  !
  !     The terms involving the factors bsl, bsu , x**_ and t0**_
  !     are related to calculation of the bootstrap effect.
  !     We assume virtually that the distribution is skewed asymetrically
  !     in the trapped region...that is we assume (virtually) that
  !     f(itl).ne.f(itu) and that the difference is driven by
  !     a df/dr term through bsl and bsu. Since this term involves f at
  !     different radial positions, it cannot figure into the solution
  !     implicitly, that is, it is differenced explicitly. The resulting
  !     contributions appear below. There will be contributions from
  !     i=itl-1, itu+1, itl and itu only.
  !     All contributions are zero elsewhere, and are zero everywhere
  !     if bootcalc= "disabled".   (Refer to Harvey et al, 1993 Sherwood
  !     Theory Mtg; E. Westerhof and A.G. Peters, Computer Physics Comm.,
  !     Vol. 95, p. 131-138 (1996).)
  !.......................................................................


  !cc   z00f(i,j)=vptb(i,lr_)*(f_(i,j,k,l_)/dtreff+so(i,j)) +
  !cc   +  spasou(i,j,k,l_)

  !     itl-1 case:

  !cc   z00itl1(i,j)=z00f(i,j)
  !cc   1           -xpm(i,j)*bsl(j-1,k,l_)-xp0(i,j)*bsl(j,k,l_)
  !cc   2           -xpp(i,j)*bsl(j+1,k,l_)

  !     itu+1 case:

  !cc   z00itu1(i,j)=z00f(i,j)
  !cc   1           -xmm(i,j)*bsu(j-1,k,l_)-xm0(i,j)*bsu(j,k,l_)
  !cc   2           -xmp(i,j)*bsu(j+1,k,l_)

  !     itl or itu case:

  !cc   t0ml_(j)=qz(j)*(
  !cc   1 cl(itl-1,j-1)*dj(itl,j-1,k,l_)*eym5(itl,l_))
  !cc   1 +r2y(j)*(-de(itl-1,j)*(1.-di(itl-1,j-1,k,l_)))
  !cc   1 *0.5*dxi(j)

  !cc   t00l_(j)=
  !cc   1 +qz(j)*(
  !cc   1 -cl(itl-1,j)*dj(itl,j,k,l_)*eym5(itl,l_)
  !cc   1 +cl(itl-1,j-1)*(1.-dj(itl,j-1,k,l_))*eym5(itl,l_))
  !cc   1 +r2y(j)*(dd(itl-1,j)*(1.-di(itl-1,j,k,l_))
  !cc   1 +df(itl-1,j)*eym5(itl,l_))

  !cc   t0pl_(j)=qz(j)*(
  !cc   1 -cl(itl-1,j)*eym5(itl,l_)*(1.-dj(itl,j,k,l_)))
  !cc   1 +r2y(j)*de(itl-1,j)*0.5*dxi(j)*(1.-di(itl-1,j+1,k,l_))


  !cc   t0mu_(j)=qz(j)*(
  !cc   1 -cl(itu+1,j-1)*dj(itu,j-1,k,l_)*eyp5(itu,l_))
  !cc   1 +r2y(j)*(
  !cc   1 +de(itu,j)*di(itu,j-1,k,l_))*0.5*dxi(j)

  !cc   t00u_(j)=
  !cc   1 +qz(j)*(
  !cc   1 +cl(itu+1,j)*dj(itu,j,k,l_)*eyp5(itu,l_)
  !cc   1 -cl(itu+1,j-1)*(1.-dj(itu,j-1,k,l_))*eyp5(itu,l_))
  !cc   1 +r2y(j)*(
  !cc   1 -dd(itu,j)
  !cc   1 *di(itu,j,k,l_)
  !cc   1 +df(itu,j)*eyp5(itu,l_))

  !cc   t0pu_(j)=qz(j)*(
  !cc   1 +cl(itu+1,j)*(1.-dj(itu,j,k,l_))*eyp5(itu,l_))
  !cc   1 +r2y(j)*(-de(itu,j)*di(itu,j+1,k,l_)*0.5*dxi(j))

  !cc   z00itl(i,j)=z00f(i,j)
  !cc   1           -(t0ml_(j)*bsl(j-1,k,l_)+t00l_(j)*bsl(j,k,l_)
  !cc   2           +t0pl_(j)*bsl(j+1,k,l_)+t0mu_(j)*bsu(j-1,k,l_)
  !cc   3           +t00u_(j)*bsu(j,k,l_)+t0pu_(j)*bsu(j+1,k,l_))

  !cc   z00ff(i,j)=cvmgt(z00itu1(i,j),z00itl(i,j),i.eq.(itu+1))

  !cc   z00t(i,j)=cvmgt(z00itl1(i,j),z00ff(i,j),i.eq.(itl-1))

  !cc   z00(i,j)=cvmgt(z00t(i,j),z00f(i,j),bootcalc.ne."disabled".and.
  !cc   1                                  (i.eq.(itl-1).or.i.eq.itl.or.
  !cc   2                                   i.eq.itu.or.i.eq.(itu+1)))

  !.......................................................................
  !     Pass/Trapped boundary statement functions follow..
  !.......................................................................


  real(c_double) function tmm(j)
    integer :: j
    tmm = -qz(j)*cl(itl-1,j-1)*dj(itl-1,j-1,advnce_k,l_)*eym5(itl,l_) &
         -r2y(j)*di(itl-1,j-1,advnce_k,l_)*de(itl-1,j)*0.5*dxi(j)
  end function tmm

  real(c_double) function tm0(j)
    integer :: i,j
    tm0 = qz(j)*cl(itl-1,j)*dj(itl-1,j,advnce_k,l_)*eym5(itl,l_) &
         -qz(j)*cl(itl-1,j-1)*(1.-dj(itl-1,j-1,advnce_k,l_))*eym5(itl,l_) &
         +r2y(j)*(dd(itl-1,j)*di(itl-1,j,advnce_k,l_) &
         -df(itl-1,j)*eym5(itl,l_))
  end function tm0

  real(c_double) function tmp(j)
    integer :: j
    tmp = qz(j)*cl(itl-1,j)*(1.-dj(itl-1,j,advnce_k,l_))*eym5(itl,l_) &
         +r2y(j)*di(itl-1,j+1,advnce_k,l_)*de(itl-1,j)*0.5*dxi(j)
  end function tmp

  real(c_double) function t0m(j)
    integer :: j
    t0m = qz(j)*(da(itl,j-1)*dj(itl,j-1,advnce_k,l_)-db(itl,j-1)*exm5(j)+ &
         cl(itl-1,j-1)*dj(itl,j-1,advnce_k,l_)*eym5(itl,l_)-2.*cl(itl+1,j-1)* &
         dj(itl,j-1,advnce_k,l_)*eyp5(itl,l_) &
         -cl(itu+1,j-1)*dj(itu,j-1,advnce_k,l_)*eyp5(itu,l_)) &
         +r2y(j)*(-de(itl-1,j)*(1.-di(itl-1,j-1,advnce_k,l_))+2.*de(itl,j) &
         *di(itl,j-1,advnce_k,l_)+de(itu,j)*di(itu,j-1,advnce_k,l_))*0.5*dxi(j)
  end function t0m

  real(c_double) function t00(j)
    integer :: j
    t00 = vptb(itl,lr_)/dtreff &
         +qz(j)*(-da(itl,j)*dj(itl,j,advnce_k,l_)+db(itl,j)* &
         exp5(j)-cl(itl-1,j)*dj(itl,j,advnce_k,l_)*eym5(itl,l_) &
         +2.*cl(itl+1,j)*dj(itl,j,advnce_k,l_) &
         *eyp5(itl,l_) &
         +cl(itu+1,j)*dj(itu,j,advnce_k,l_)*eyp5(itu,l_) &
         +da(itl,j-1)*(1.-dj(itl,j-1,advnce_k,l_)) &
         +db(itl,j-1)*exm5(j) &
         +cl(itl-1,j-1)*(1.-dj(itl,j-1,advnce_k,l_))*eym5(itl,l_) &
         -2.*cl(itl+1,j-1)*(1.-dj(itl,j-1,advnce_k,l_))*eyp5(itl,l_) &
         -cl(itu+1,j-1)*(1.-dj(itu,j-1,advnce_k,l_))*eyp5(itu,l_)) &
         +r2y(j)*(dd(itl-1,j)*(1.-di(itl-1,j,advnce_k,l_)) &
         +df(itl-1,j)*eym5(itl,l_) &
         -2.*dd(itl,j)*di(itl,j,advnce_k,l_) &
         +2.*df(itl,j)*eyp5(itl,l_)-dd(itu,j) &
         *di(itu,j,advnce_k,l_) &
         !XXXXXXXXXXXXXXXXXXXXXXx bug? +df(itu,j)*eyp5(itu,l_))-vptb(i,lr_)*cah(itl,j)
         +df(itu,j)*eyp5(itu,l_))-vptb(itl,lr_)*cah(itl,j)
  end function t00

  real(c_double) function t0p(j)
    integer :: j
    t0p = qz(j)*(-da(itl,j)*(1.-dj(itl,j,advnce_k,l_))-db(itl,j)*exp5(j) &
         -cl(itl-1,j)*eym5(itl,l_)*(1.-dj(itl,j,advnce_k,l_))+2.*cl(itl+1,j) &
         *eyp5(itl,l_)*(1.-dj(itl,j,advnce_k,l_)) &
         +cl(itu+1,j)*(1.-dj(itu,j,advnce_k,l_))*eyp5(itu,l_)) &
         +r2y(j)*(de(itl-1,j)*0.5*dxi(j)*(1.-di(itl-1,j+1,advnce_k,l_))- &
         2.*de(itl,j)*di(itl,j+1,advnce_k,l_)*0.5*dxi(j) &
         -de(itu,j)*di(itu,j+1,advnce_k,l_)*0.5*dxi(j) )
  end function t0p

  real(c_double) function tpm(j)
    integer :: j
    tpm = 2.*qz(j)*cl(itl+1,j-1)*eyp5(itl,l_)*dj(itl+1,j-1,advnce_k,l_) &
          +2.*r2y(j)*de(itl,j)*(1.-di(itl,j-1,advnce_k,l_))*0.5*dxi(j)
  end function tpm

  real(c_double) function tp0(j)
    integer :: j
    tp0 = -2.*qz(j)*(cl(itl+1,j)*dj(itl+1,j,advnce_k,l_)*eyp5(itl,l_)- &
         cl(itl+1,j-1)*(1.-dj(itl+1,j-1,advnce_k,l_))*eyp5(itl,l_)) &
         -2.*r2y(j)*df(itl,j)*eyp5(itl,l_) &
         -2.*r2y(j)*dd(itl,j)*(1.-di(itl,j,advnce_k,l_))
  end function tp0

  real(c_double) function tpp(j)
    integer :: j
    tpp = -2*qz(j)*cl(itl+1,j)*eyp5(itl,l_)*(1.-dj(itl+1,j,advnce_k,l_)) &
         -2.*r2y(j)*de(itl,j)*0.5*dxi(j)*(1.-di(itl,j+1,advnce_k,l_))
  end function tpp

  real(c_double) function tum(j)
    integer :: j
    tum = qz(j)*cl(itu+1,j-1)*eyp5(itu,l_)*dj(itu+1,j-1,advnce_k,l_) &
         +r2y(j)*de(itu,j)*0.5*dxi(j)*(1.-di(itu,j-1,advnce_k,l_))
  end function tum

  real(c_double) function tu0(j)
    integer :: j
    tu0 = -qz(j)*cl(itu+1,j)*dj(itu+1,j,advnce_k,l_)*eyp5(itu,l_) &
         +qz(j)*cl(itu+1,j-1)*eyp5(itu,l_)*(1.-dj(itu+1,j-1,advnce_k,l_)) &
         -r2y(j)*(dd(itu,j)*(1.-di(itu,j,advnce_k,l_)) &
         +df(itu,j)*eyp5(itu,l_))
  end function tu0

  real(c_double) function tup(j)
    integer :: j
    tup = -qz(j)*cl(itu+1,j)*eyp5(itu,l_)*(1.-dj(itu+1,j,advnce_k,l_)) &
         -r2y(j)*de(itu,j)*(1.-di(itu,j+1,advnce_k,l_))*0.5*dxi(j)
  end function tup

  !..................................................................
  !     Express the Chang-Cooper weighted average f(i,j+1/2,l_): fpj
  !..................................................................

  real(c_double) function fpj(i,j,k)
    integer :: i,j, k
    !k bug
    fpj = f(i,j+1,k,l_)+ (f(i,j,k,l_)-f(i,j+1,k,l_))*dj(i,j,advnce_k,l_)
  end function fpj

  real(c_double) function fpjp(i,j, k)
    integer :: i,j, k
    ! k bug
    fpjp = fpj(i+1,j,k) + cvmgt(bsl(j,k,l_),zero,(i+1).eq.itl)
  end function fpjp

  real(c_double) function fpj0(i,j,k)
    integer :: i,j, k
    ! k bug
    fpj0 = fpj(i,j,k) + cvmgt(bsu(j,k,l_),zero,i.eq.itu)
  end function fpj0

  !YuP-110106: error corrected:  l_
  !     Note: no need to check bootcalc="disabled" or not,
  !     because when bootcalc="disabled",  bsl==0 and bsu==0.
  !..................................................................
  !     Express the velocity flux at (i,j+1/2)
  !..................................................................
  !     -YuP: Moved gft,gfi to diagentr.f, to avoid cvmgt() construct
  !     gft(j)=
  !     1  +cl(itl-1,j)*eyp5(itl-1,l_)*(fpjp(itl-1,j)-fpj(itl-1,j))
  !     1  +2.*cl(itl+1,j)*eyp5(itl,l_)*(fpj(itl+1,j)-fpj(itl,j))
  !     1  +cl(itu+1,j)*eyp5(itu,l_)*(fpj(itu+1,j)-fpj0(itu,j))
  !
  !     gfi(i,j)=da(i,j)*fpj(i,j)
  !     1  +db(i,j)*exp5(j)*(f(i,j+1,k,l_)-f(i,j,k,l_))
  !     1  +cvmgt(dc(i,j)*0.5*dyi(i,l_)*(fpjp(i,j)-fpj0(i-1,j)),
  !     1  gft(j),
  !     1  (i.ne.itl .and. i.ne.itu) .or. symtrap.ne."enabled")

  !..................................................................
  !     Express the Chang-Cooper weighted average f(i+1/2,j,l_): fpi
  !..................................................................

  real(c_double) function fpip(i,j)
    !XXX integer :: k=1, l_=1, we should specifiy.....
    integer :: i,j
    fpip = f(i+1,j,advnce_k,l_) + cvmgt(bsl(j,advnce_k,l_),zero,(i+1).eq.itl)
  end function fpip

  real(c_double) function fpi0(i,j)
    integer :: i,j
    fpi0 = f(i,j,advnce_k,l_) + cvmgt(bsu(j,advnce_k,l_),zero,i.eq.itu)
  end function fpi0

  !     Note: no need to check bootcalc="disabled" or not,
  !     because when bootcalc="disabled",  bsl==0 and bsu==0.

  real(c_double) function fpi(i,j)
    integer :: i,j
    fpi = fpip(i,j)*(1.-di(i,j,advnce_k,l_)) + fpi0(i,j)*di(i,j,advnce_k,l_)
  end function fpi

  !..................................................................
  !     Express the theta flux at (i+1/2,j)
  !..................................................................

  real(c_double) function hfi(i,j)
    integer :: i,j
    hfi = dd(i,j)*fpi(i,j) &
         +de(i,j)*0.5*dxi(j)*(fpi(i,ifp(j))-fpi(i,j-1)) &
         +df(i,j)*eyp5(i,l_)*(fpip(i,j)-fpi0(i,j))
  end function hfi

!..................................................................
!     End of statement functions used for implicit time advancement.
!..................................................................

end module advnce_mod
