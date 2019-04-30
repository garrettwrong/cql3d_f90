module tdfinterp_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use lookup_mod, only : lookup

  !---END USE

!
!

contains

      subroutine tdfinterp(u,pitch,rho,polang,fdist)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!
!     Interpolate f for value at u,pitch,rho,polang ==> fdist
!          u=momentum-per-rest-mass (cm/sec)
!          pitch=pitch angle at polang [0,pi]
!          rho=normalized radial coord
!          polang=poloidal angle, in [0,2*pi)
!                                    [or [0,pi] if eqsym.ne.'none']
!
!     f will be in code units (i.e., f_cgs*vnorm**3). where
!     f_cgs is in cgs units 1/(volume*velocity**3).
!     Normalization is thus such that integral{f*d**3x}=ne(#/cm**3).
!
!     Zero-orbit-width (readily generalized).
!
!     Bob Harvey, June, 2000/April, 2010
!
!..................................................................

      real(c_double), dimension(:), allocatable ::  rovera_
      data ifirst/1/, icount/0/


!.......................................................................
!
!    CQL3D has a facility to computer FP solutions on a subset of
!      the full radial mesh on which plasma parameters are specified.
!
!    In the following, we are only interested in the given data
!    which is at flux surfaces numbered ll=1:lrz
!    (thus we can consider lrindx(ll)=ll).
!
!
!     distribution function normalized so that
!         integral( (dx)**3 f) = density at minumum B point.
!
!.......................................................................
!
!.......................................................................
!     Fill rovera_, containing FP'd flux surfaces (lrz.le.lrzmax).
!     f is given on this restricted set of surfaces f(,,,1:lrz).
!.......................................................................

      if (ifirst.eq.1 ) then
         allocate (rovera_(lrz), STAT=istat)
         do ll=1,lrz
            rovera_(ll)=rovera(lrindx(ll))
         enddo
         ifirst=0
      endif
!
!.......................................................................
!     For given poloidal angle and pitch, determine equatorial
!     plane pitch angle.
!     ll= radial index in rovera_(1:lrz) grid
!     lr= radial index in rovera(1:lrzmax) grid
!     lookup() adjusts weights for endpoints.
!.......................................................................

      call lookup(rho,rovera_,lrz,rweightu_,rweightl_,ll)
      call lookup(rho,rovera,lrzmax,rweightu,rweightl,lr)

!     Lookup polang on surrounding pol[1:lz,1:lrzmax] meshes,
!     interpolate for B/B0, and obtain equatorial plane pitch angle.

      if (eqsym.ne."none" .and. polang.gt.pi) &
           polang=2*pi-polang   !assumes up-down symm

      call lookup(polang,pol(1:lz,lr-1),lz,pweightu1,pweightl1,lt1)
      call lookup(polang,pol(1:lz,lr),lz,pweightu2,pweightl2,lt2)
      bdb01=pweightl1*bbpsi(lt1,lr-1)+pweightu1*bbpsi(lt1,lr-1)
      bdb02=pweightl2*bbpsi(lt2,lr)+pweightu2*bbpsi(lt2,lr)
      bdb0=rweightl*bdb01+rweightu*bdb02

!BH110817      sin_pitch02=sin(pitch)**2/bdb0
!BH110817      pitch0=asin(sqrt(sin_pitch02))  ! Corresponding midplane pitch
!BH110817: Fix for pitch.gt.pio2.  Since tdfinterp only used so far
!BH110817: for perp viewing NPA in C-Mod, this hasn't caused a problem.
      sin_pitch0=sin(pitch)/sqrt(bdb0)
      pitch0=asin(sin_pitch0)  ! Corresponding midplane pitch
      if (pitch.gt.pio2) pitch0=pi-pitch0


!.......................................................................
!     BEGIN: Interpolate f for value at u,pitch0,rho
!     Single species, although easy to generalize
!.......................................................................

      k=1
      xvel=u/vnorm

!     First, determine if the particle at (u,pitch0,radius)
!       is trapped.
!     Use lookup to find weights for distribution at neighboring
!       tabulated points.  Use closest pitch angle mesh to
!       determine whether trapped.

!$$$      if (rweightl.ge.rweightu) llp=lr-1
!$$$      lpitch=luf(pitch,y(1,llp),iy_(llp))
!$$$      if (lpitch.gt.(iy_(llp)/2)) lpitch=lpitch-1 ! Gives symmetry.
!$$$      itrap=0                                   ! Trapping flag.
!$$$      if (lpitch.gt.itl_(llp) .and. lpitch.lt.itu_(llp)) itrap=1
!$$$
!$$$c     Shift radial location of trapped particles by a half banana width
!$$$c     This will be an outward/inward shift for pos/neg bthr*bnumb
!$$$c     (bnumb is charge number, including sign).
!$$$c     Transiting particles are unshifted.
!$$$
!$$$      rr=rho
!$$$      charge=4.8032d-10
!$$$      if (itrap.eq.1) then
!$$$         qb_mc=bnumb(k)*charge*bthr(lr)/(fmass(k)*clight)
!$$$         rban=xvel*cos(pitch)*vnorm/qb_mc
!$$$c     Shift radial location of trapped particles by a half banana width
!$$$c     This will be an outward/inward shift for pos/neg bthr*bnumb
!$$$c     (bnumb is charge number, including sign).
!$$$         rrr=rho+rban/radmin
!$$$         if (rr.lt. 0.d0) rr=0.d0
!$$$         if (rr.gt. 1.d0)  rr=1.d0
!$$$      endif

!..................................................................
!     Use lookup to find weights for distribution at neighboring
!     tabulated points.  Interpolate to velocity and pitch angle
!     to the specified (xvel,pitch) values.
!..................................................................

!     lookup pitch0 on y-grid corresponding to neighboring f()
      lrm=lrindx(ll-1)
      call lookup(pitch0,y(1:iy_(lr),lr),iy_(lr),pweightu,pweightl,iyy)
      call lookup(pitch0,y(1:iy_(lrm),lrm),iy_(lrm),pweightum,pweightlm,iyym)
      call lookup(xvel,x,jx,xweightu,xweightl,jxx)

!     Tri-linear interpolate (Eliminate f<0 values).

!$$$      fdist=
!$$$     +     rweightl_*(pweightlm*xweightl*max(0d0,f(iyym-1,jxx-1,k,ll-1))
!$$$     +               +pweightum*xweightl*max(0d0,f(iyym,jxx-1,k,ll-1))
!$$$     +               +pweightlm*xweightu*max(0d0,f(iyym-1,jxx,k,ll-1))
!$$$     +               +pweightum*xweightu*max(0d0,f(iyym,jxx,k,ll-1)))
!$$$     +     +rweightu_*(pweightl*xweightl*max(0d0,f(iyy-1,jxx-1,k,ll))
!$$$     +               +pweightu*xweightl*max(0d0,f(iyy,jxx-1,k,ll))
!$$$     +               +pweightl*xweightu*max(0d0,f(iyy-1,jxx,k,ll))
!$$$     +               +pweightu*xweightu*max(0d0,f(iyy,jxx,k,ll)))

!     Tri-linear interpolate log(f) (and eliminate log(f)<-200 values).
!     (The log of constant temperature distributions are linear vs vel,
!      making linear interpolation appropriate.)
      fmin= 1.383896526737250d-87  ! log(1.383896526737250d-87)= -200.
!      fmin= 3.720075976020916d-44  ! log(3.720075976020916d-44)= -100.


!     Adapt to calc of fdist from either favg or f:
      if (f4d_out.eq."tavg" .and. n.eq.nstop) then
         fmmm=favg(iyym-1,jxx-1,k,ll-1)
         f0mm=favg(iyym,  jxx-1,k,ll-1)
         fm0m=favg(iyym-1,jxx,  k,ll-1)
         f00m=favg(iyym,  jxx,  k,ll-1)
         fmm0=favg(iyym-1,jxx-1,k,ll)
         f0m0=favg(iyym,  jxx-1,k,ll)
         fm00=favg(iyym-1,jxx,  k,ll)
         f000=favg(iyym,  jxx,  k,ll)
      else
         fmmm=f(iyym-1,jxx-1,k,ll-1)
         f0mm=f(iyym,  jxx-1,k,ll-1)
         fm0m=f(iyym-1,jxx,  k,ll-1)
         f00m=f(iyym,  jxx,  k,ll-1)
         fmm0=f(iyym-1,jxx-1,k,ll)
         f0m0=f(iyym,  jxx-1,k,ll)
         fm00=f(iyym-1,jxx,  k,ll)
         f000=f(iyym,  jxx,  k,ll)
      endif

      fdist= rweightl_* &
          (  pweightlm*xweightl*log(max(fmmm,fmin)) &
            +pweightum*xweightl*log(max(f0mm,fmin)) &
            +pweightlm*xweightu*log(max(fm0m,fmin)) &
            +pweightum*xweightu*log(max(f00m,fmin)) ) &
            +rweightu_* &
          (  pweightl*xweightl*log(max(fmm0,fmin)) &
            +pweightu*xweightl*log(max(f0m0,fmin)) &
            +pweightl*xweightu*log(max(fm00,fmin)) &
            +pweightu*xweightu*log(max(f000,fmin)) )
      fdist=exp(fdist)

      return
      end

!
end module tdfinterp_mod
