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

!     trans.h

!..............................................................
!     If the transport model is utilized, define the alpha, beta
!     gamma and delta for the tridiagonal sweep option next.
!     h_r defined =H*rho=del_V/del_rho/(4*pi**2*R0)
!     For new soln_method="it3drv" option:
!     BH080307:
!     In alprp,betrp,gamrp, a zmaxpsi(l) factor is added
!     BH080410:
!     vptb/zmaxpsi factor at appropriate l absorbed into each
!     alprp,betrp,gamrp.
!     BH171226:
!     For difus_io(k)= "drrin" and "drrdrin":
!     Added scale factor drrt(k) for d_rr and drt(k) for d_r.  These
!     scale factors are defaulted to 1d0, and set with function
!     difus_io_scale(,).  The d_rr (and d_r for difus_io(k).eq."drrdrin"
!     are not changed by the scale factor, rather the scale factors
!     are applied as appropriate.
!..............................................................

      ztr(i,l)=cosovb(idx(i,l),l)/dvol(setup0%lrindx(l))*4.*pi**2*radmaj
      ztrp(i,l)=cosovb(idx(i,l),l)/dvol(setup0%lrindx(l))*4.*pi**2*radmaj* &
                zmaxpsi(setup0%lrindx(l))
      ytr(i,l)=h_r(setup0%lrindx(l))*drrt(k)*d_rr(idx(i,l),j,k,indxlr(l))/ &
        drp5(setup0%lrindx(l))*bovcos(idx(i,l),l)
! BH080410:  Changing sign of xtr (compatible with usual positive
! BH080410:  radial drift in the positive radial direction.
!            Need to check have changed sign for
!            soln_method=direct, transp=enabled case.
!     The drt(k) scale factor is only applied for difus_io(k).eq.
!     "drrdrin"; else d_r is calculated from the scaled d_rr.
      xtr(i,l)=cvmgt(-h_r(setup0%lrindx(l))*drt(k)*d_r(idx(i,l),j,k,indxlr(l))* &
        bovcos(idx(i,l),l),-h_r(setup0%lrindx(l))*d_r(idx(i,l),j,k,indxlr(l))* &
        bovcos(idx(i,l),l),difus_io(k).eq."drrdrin")

      alpr(i,l)=ztr(i,l)*(ytr(i,l)+xtr(i,l)*(1.-dl(idx(i,l),j,k,l)))
      alprp(i,l)=ztrp(i,l)*(ytr(i,l)+xtr(i,l)*(1.-dl(idx(i,l),j,k,l)))* &
        vptb_(idx(i,l+1),setup0%lrindx(l+1))*zmaxpsii(setup0%lrindx(l+1))
      betr(i,l)=ztr(i,l)*(ytr(i,l)-xtr(i,l)*dl(idx(i,l),j,k,l)+ &
        ytr(i,l-1)+xtr(i,l-1)*(1.-dl(idx(i,l-1),j,k,l-1)))+1./dttr
! BH071109:  Removing 1./dttr for soln_method="it3drv"
      betrp(i,l)=ztrp(i,l)*(ytr(i,l)-xtr(i,l)*dl(idx(i,l),j,k,l)+ &
        ytr(i,l-1)+xtr(i,l-1)*(1.-dl(idx(i,l-1),j,k,l-1)))* &
        vptb_(idx(i,l),setup0%lrindx(l))*zmaxpsii(setup0%lrindx(l))
      gamr(i,l)=ztr(i,l)*(ytr(i,l-1)-xtr(i,l-1)*dl(idx(i,l-1),j,k,l-1))
      gamrp(i,l)=ztrp(i,l)*(ytr(i,l-1)-xtr(i,l-1)* &
        dl(idx(i,l-1),j,k,l-1))* &
        vptb_(idx(i,l-1),setup0%lrindx(l-1))*zmaxpsii(setup0%lrindx(l-1))
      delr(i,l)=frn(idx(i,l),j,k,l)/dttr*vptb_(idx(i,l),setup0%lrindx(l))* &
        zmaxpsii(setup0%lrindx(l)) + velsou(idx(i,l),j,k,l)* &
        zmaxpsii(setup0%lrindx(l))

      f1l(i,j,k,l)=frn(idx(i,l),j,k,l)*vptb_(idx(i,l),setup0%lrindx(l))* &
        zmaxpsii(setup0%lrindx(l))*dl(idx(i,l),j,k,l)+ &
        frn(idx(i,l+1),j,k,l+1)*vptb_(idx(i,l+1),setup0%lrindx(l)+1) &
        *zmaxpsii(setup0%lrindx(l)+1)*(1.-dl(idx(i,l),j,k,l))

      sfu(i,j,k,l)=drrt(k)*d_rr(idx(i,l),j,k,indxlr(l)) &
        *(frn(idx(i,l+1),j,k,l+1) &
        *vptb_(idx(i,l+1),setup0%lrindx(l)+1)*zmaxpsii(setup0%lrindx(l)+1)- &
        frn(idx(i,l),j,k,l)* &
        vptb_(idx(i,l),setup0%lrindx(l))*zmaxpsii(setup0%lrindx(l)))/drp5(setup0%lrindx(l))+ &
        cvmgt(drt(k)*d_r(idx(i,l),j,k,indxlr(l))*f1l(i,j,k,l), &
        d_r(idx(i,l),j,k,indxlr(l))*f1l(i,j,k,l), &
        difus_io(k).eq."drrdrin")

!BH080425:  Changed sign on d_r, consistent above BH080410
      sfup(i,j,k,l)=drrt(k)*d_rr(idx(i,l),j,k,indxlr(l))* &
        (frn(idx(i,l+1),j,k,l+1)*vptb_(idx(i,l+1),setup0%lrindx(l)+1)* &
        zmaxpsii(setup0%lrindx(l)+1)-frn(idx(i,l),j,k,l)* &
        vptb_(idx(i,l),setup0%lrindx(l))*zmaxpsii(setup0%lrindx(l)))/drp5(setup0%lrindx(l))- &
        cvmgt(drt(k)*d_r(idx(i,l),j,k,indxlr(l))*f1l(i,j,k,l), &
        d_r(idx(i,l),j,k,indxlr(l))*f1l(i,j,k,l), &
        difus_io(k).eq."drrdrin")
!*************************************************************
