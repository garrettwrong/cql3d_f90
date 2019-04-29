module eqvolpsi_mod

  !---BEGIN USE

  use eqorbit_mod, only : eqorbit

  !---END USE

!
!

contains

      subroutine eqvolpsi(epsicon_,volum,areac)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)


!..................................................................
!     This routine returns the volume and area of a flux surface given
!     that the psi label of the flux surface is epsi. If we
!     can assume that eqorbit has been previously called for
!     this flux surface, set eqorb="disabled" before the call, otherwise
!     set eqorb="enabled". Note if eqorb is "disabled" then it is
!     assumed that the field line data exists in arrays such as
!     solr_, solz_, etc.
!..................................................................

      if (eqorb.eq."enabled") then
        call eqorbit(epsicon_)
      endif


!......................................................................
!     The volume (2*pi*Z*R*dR) is computed by a simple mid-point rule
!     For the non-updown-symm case, measure volume and area to the
!     the contour from the diagonal from (rmcon_,zmcon) to (rpcon_,zpcon_).
!......................................................................

      volum=0.
      areac=0.
         do  l=2,lorbit_
            ravg=0.5*(solr_(l)+solr_(l-1))
            zdiag=zmcon_+(ravg-rmcon_)*(zpcon_-zmcon_)/(rpcon_-rmcon_)
            volum=volum+abs((solr_(l)-solr_(l-1))* &
                 (solz_(l)+solz_(l-1)-2.*zdiag)*ravg)
            areac=areac+abs((solr_(l)-solr_(l-1))* &
                 0.5*(solz_(l)+solz_(l-1)-2.*zdiag))
         enddo
!      else
!         do l=2,lorbit_
!            volum=volum+(solr_(l)-solr_(l-1))*(solz_(l)+solz_(l-1))*
!     1           (solr_(l)+solr_(l-1))
!            areac=areac+(solr_(l)-solr_(l-1))*(solz_(l)+solz_(l-1))
!         enddo
!      endif
      if (eqsym.ne."none") then ! up-down symmetry
          !-YuP: same as non-sym case,
          ! but only half of flux surface was traced
          volum=volum*2.
          areac=areac*2.
      endif

      volum=abs(volum*pi)
      areac=abs(areac)

      !YuP[2015/05/02] Added a more accurate def. of volume and area
      ! near magnetic axis.
      !This definition is not dependent/sensitive to tracing the surface.
      ! Gives better convergence in eqfndpsi in cases for
      ! radcoord.eq."sqtorflx", when 'rhonew' is proportional to
      ! sqrt of 'volum' (near magn.axis).
      if(abs(epsicon_-psimag).le.psimag*1.d-3)then
         ! Very close to magn.axis: assume circular flux surfaces
         areac= pi*(solr_(1)-rmag)**2 ! solr_(1)-rmag == minor radius
         volum= 2*pi*rmag*areac
      endif

      return
      end
end module eqvolpsi_mod
