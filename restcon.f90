module restcon_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine restcon
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Connor formula for toroidal resitivity,
!     J.W. Connor, R.C. Grimm, R.J. Hastie, P.M. Keeping,
!     Nucl. Fus. 13, 211 (1973).
!     The following calculation of the transiting particle fraction,
!     xconn (=I in Connor et al.),
!     is not entirely clear to me, but for a multi-flux-surface
!     test case, 1.-xconn agrees to almost
!     3 significant figures with trapfrac(lr_) from eqfndpsi.f.
!     (BobH, 010803).
!..................................................................

      s=0.
      ssg=0.
      u=0.
      uu=0.
      ilzhfs=lz
      if (numclas .eq. 1) ilzhfs=lz/2+1
      do 10 l=1,ilzhfs
        uu=uu+dz(l,lr_)*bbpsi(l,lr_)**2
 10   u=u+dz(l,lr_)*bbpsi(l,lr_)
      qz=.75*u/pi
      qqg=.75*uu/pi
      do 11 i=1,itl-2
        qx=coss(i,l_)*sinn(i,l_)**2*cynt2(i,l_)
        s=s+qx/(tau(i,lr_)*(psiiv(i,lr_)-sinn(i,l_)**2))
        ssg=ssg+qx/(tau(i,lr_)*(1.-psiba(i,lr_)*sinn(i,l_)**2))
 11   continue
      qx=cynt2(itl-1,l_)+pi*sinn(itl,l_)*(y(itl,l_)-y(itl-1,l_))
      qx=coss(itl-1,l_)*sinn(itl-1,l_)**2*qx
      s=s+qx/(tau(itl-1,lr_)*(psiiv(itl-1,lr_)-sinn(itl-1,l_)**2))
      ssg=ssg+qx/(tau(itl-1,lr_)*(1.-psiba(itl-1,lr_)* &
        sinn(itl-1,l_)**2))
      xconn=s*qz
      rovsc(lr_)=(1.+.471*(1.-xconn))/xconn/(1.+.039*(1.-xconn))
      fc=ssg*qqg
      cr=.56/bnumb(kionn)*(3.-bnumb(kionn))/(3.+bnumb(kionn))
      eovsz=1./((1.-cr*(1.-fc))*fc)
      return
      end subroutine restcon
      
end module restcon_mod
