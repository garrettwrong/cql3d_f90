module wpchgdy_mod

!
!

contains

      subroutine wpchgdy
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!.......................................................................
!     Adapt dyp5(iyh) and dy(iyh) and related quantities, such that
!
!     sum(i=1,iyh) [cynt2(i,l_)] = sum(i=1,iyh) [cynt2(i,ll)]
!     where ll=1 for CQLP and lrz for CQL3D. This may be useful for
!     meshy=fixed_mu option
!.......................................................................

!.......................................................................
!l    1. Re-define dyp5(iyh,l_)
!.......................................................................

      ilref=1
      if (cqlpmod .ne. "enabled") ilref=lrz
      if (l_ .eq. ilref) go to 999

      zsum0=0.0
      do 100 i=1,iyh_(ilref)
        zsum0=zsum0+cynt2(i,ilref)/twopi
 100  continue
      zsuml=0.0
      do 110 i=1,iyh-1
        zsuml=zsuml+cynt2(i,l_)/twopi
 110  continue
      dyp5(iyh,l_)=2.*(zsum0-zsuml)/sinn(iyh,l_) - dym5(iyh,l_)

!.......................................................................
!l    2. Correct related quantities
!.......................................................................

      dy(iyh,l_)=0.5*(dym5(iyh,l_)+dyp5(iyh,l_))
      dy(iyh+1,l_)=dy(iyh,l_)
      eyp5(iyh,l_)=1./dyp5(iyh,l_)
      dym5(iyh+1,l_)=dyp5(iyh,l_)
      eym5(iyh+1,l_)=eyp5(iyh,l_)
      cynt2(iyh,l_)=dy(iyh,l_)*twopi*sinn(iyh,l_)
      cynt2(iyh+1,l_)=cynt2(iyh,l_)
      twoint(l_)=0.
      do 200 i=1,iy
        twoint(l_)=twoint(l_)+cynt2(i,l_)
 200  continue

 999  return
      end
end module wpchgdy_mod
