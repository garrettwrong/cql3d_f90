module cfpleg_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use r8subs_mod, only : dscal

  !---END USE

!
!

contains

      subroutine cfpleg(m,l,jzval)
        use cqlconf_mod, only : setup0
        use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dscal
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!.......................................................................
!     This routine evaluates the Legendre coefficients given that the
!     the midplane or the local distribution is stored in temp3.
!
!     Note:dcofleg is weight for interval [i-1/2,i+1/2],
!     therefore temp3 should have the value of f on the nodes and
!     not at i+1/2 as previously.
!
!     This will return in tam1(j) the coefficient evalutated at z(l,lr_).
!     The first non-zero coefficient will be indexed by jzval - the
!     others presumably not being required.
!     The method for integrating over theta to obtain the Leg. coeff.
!     depends on the value of ngauss (and subsequently of analegco)
!.......................................................................


      call bcast(tam1,zero,jx) !-> Reset all Legendre coefficients to zero

!     Modification introduced by YuP, for relativ='fully':  There are
!     problems with subtraction of two large terms at low v for
!     nonrelativistic Te.  The following jzmin is introduced for a case
!     with jx=300, Te=5 keV.  A general solution would be to combine
!     the relativ="enabled", quasi-relativistic expressions at lower v
!     with relativ="fully" treatment at higher v  (YuP: Oct'09).
!
      jzmin = jzval  !-YuP-> first non-zero coefficient
      !-YuP-> Keep high-m Legendre coefficients for small j (v~0) equal zero:

      if(relativ.eq."fully".and. m.ge.1) jzmin=MAX(3,jzval) ! For temp=5keV:
                                                             ! 33 ok (10% of jx=330);
                                                             ! 25 still bad

      !-YuP-> For high-m, first non-zero coeff. is j=33 or jzval, if jzval>33

!.......................................................................
!     dcofleg(i) is the weight between i-1/2 and i+1/2 such that the
!     Legendre coefficients tam1(j) = sum(i) [f(i)*dcofleg(i)], i=1,imax
!     and iy+1-imax,iy
!.......................................................................

      iend=itl
      if (setup0%cqlpmod .eq. "enabled") iend=iyh

!     passing particles
      do 100 i=1,iend
        ia=iy-i+1
        do 110 j=jzmin,jx
          tam1(j)=tam1(j)+temp3(i,j)*dcofleg(i,l,m,indxlr_)+ &
            temp3(ia,j)*dcofleg(ia,l,m,indxlr_)
 110    continue
 100  continue

      if (setup0%cqlpmod .eq. "enabled") go to 200
      if (l.eq.lz .or. imax(l,lr_).eq.itl) go to 200

!     particles in trapped region but passing poloidal position l
!     Note: dcofleg(imax) includes contribution up to turning point

      do 120 i=iend+1,imax(l,lr_)
        ia=iy-i+1
        do 121 j=jzmin,jx
          tam1(j)=tam1(j)+temp3(i ,j)*dcofleg(i,l,m,indxlr_)+ &
            temp3(ia,j)*dcofleg(ia,l,m,indxlr_)
 121    continue
 120  continue

 200  continue
      sa=DBLE(2*m+1)*0.5  !Per Eq. 3.1.63 in Killeen book for V_m
      call dscal(jx,sa,tam1,1)

      return
      end subroutine cfpleg


end module cfpleg_mod
