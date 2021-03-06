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

module wparsou_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

  subroutine wparsou
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This routine computes the spatial source term, due to the parallel
!     transport operator applied on f_n, needed for the velocity split
!     with the ADI (Alternating Direction Implicit) method.
!.......................................................................

      dimension zdns(lrorsa)
      dimension zdnspa(lrorsa),zdnvel(lrorsa),zdnp1(lrorsa)
      dimension zdnspa2(iy),zdnvel2(iy),zspa2i(iy/2),zvel2i(iy/2)
      dimension zdnuspa(jx),zdnuvel(jx)

      include 'wpadvnc.h'
!.......................................................................
!     with new try if sbdry=periodic .and. laddbnd=1 (see wpbdry)
!.......................................................................

!     restarted runs read spasou from restart file
      if (n.eq.0 .and. (setup0%nlrestrt.ne."disabled")) return

!%OS
      if (n .eq. 0) return
!%OS

      !call bcast(spasou(0:iy+1,0:jx+1,1:ngen,1:lrors),zero,iyjx2*ngen*lrors)
      spasou = zero
      call bcast(zdns,zero,lrors)
      call bcast(zdnspa,zero,lrors)
      call bcast(zdnvel,zero,lrors)
      call bcast(zdnp1,zero,lrors)
      call bcast(zdnspa2,zero,iymax)
      call bcast(zdnvel2,zero,iymax)
      call bcast(zspa2i,zero,iymax/2)
      call bcast(zvel2i,zero,iymax/2)
      call bcast(zdnuspa,zero,jx)
      call bcast(zdnuvel,zero,jx)

!     standard case: numindx=0 or 1
      ipshft2=0
      ipshft1=-1
      imshft2=0
      imshft1=-1
!     l loop between 1+i.strt and l_upper(i)+i.end, with .=p or m for i< or >iyh
      ipstrt=+1
      ipend=0
      imstrt=+1
      imend=0
      if (numindx .eq. 2) then
        ipshft2=(numixts+1)/2
        ipshft1=(numixts-1)/2
        imshft2=-(numixts-1)/2
        imshft1=-(numixts+1)/2
        ipstrt=-(numixts-1)/2
        ipend=-(numixts+1)/2
        imstrt=(numixts+1)/2
        imend=(numixts-1)/2
      else if (numindx .eq. 3) then
        ipshft2=+1
        ipshft1=0
        imshft2=+1
        imshft1=0
        ipstrt=0
        ipend=-1
        imstrt=0
        imend=-1
      else if (numindx .eq. 4) then
        ipshft2=+1
        ipshft1=-1
        imshft2=+1
        imshft1=-1
        ipstrt=+1
        ipend=-1
        imstrt=+1
        imend=-1
      endif
      if (sbdry .eq. "periodic") then
        ipstrt=0
        ipend=0
        imstrt=0
        imend=0
      endif
      ip1t1m2=ipshft1*(1-ipshft2)
      im1t1m2=imshft1*(1-imshft2)
      if (lmidvel .eq. 0) ip1t1m2=0
      if (lmidvel .eq. 0) im1t1m2=0

!.......................................................................
!l    1. Loop over general species
!.......................................................................

      do 100 k=1,ngen

!.......................................................................
!l    2. Loop over theta in [0,pi/2]
!.......................................................................

        do 200 i=1,iymax/2
          if (l_upper(i) .eq. 1) go to 200

!.......................................................................
!l    3. Loop over all possible l for given i and numerical scheme
!.......................................................................

          do 300 l=1+ipstrt*imstrt,l_upper(i)-ipend*imend

!.......................................................................
!l    3.0 new try
!.......................................................................

            if (sbdry.eq."periodic" .and. laddbnd.eq.1 .and. &
              l.eq.l_upper(i) .and. l_upper(i).ne.setup0%ls) then
              do 401 j=2,jx
!     points A_i and A_ii (notes p. Num(13))
                spasou(i,j,k,l)= &
                  (fnp1(i,j,k,l)-fnhalf(i,j,k,l))/dtreff-velsou(i,j,k,l)
                zdns(l)=zdns(l)+spasou(i,j,k,l)*cynt2(i,l)*cint2(j)
                iief=iy_(l)+1-i
                spasou(iief,j,k,l)= &
                  (fnp1(iief,j,k,l)-fnhalf(iief,j,k,l))/dtreff &
                  -velsou(iief,j,k,l)
                zdns(l)=zdns(l)+spasou(iief,j,k,l)*cynt2(iief,l) &
                  *cint2(j)

!     points B_i and B_ii (notes p. Num(13))
                ll=setup0%ls+2-l_upper(i)
                spasou(i,j,k,ll)= &
                  (fnp1(i,j,k,ll)-fnhalf(i,j,k,ll))/dtreff &
                  -velsou(i,j,k,ll)
                zdns(ll)=zdns(ll)+spasou(i,j,k,ll)*cynt2(i,ll)*cint2(j)
                iief=iy_(ll)+1-i
                spasou(iief,j,k,ll)= &
                  (fnp1(iief,j,k,ll)-fnhalf(iief,j,k,ll))/dtreff &
                  -velsou(iief,j,k,ll)
                zdns(ll)=zdns(ll)+spasou(iief,j,k,ll)*cynt2(iief,ll)* &
                  cint2(j)
 401          continue
              go to 300
            endif

!.......................................................................
!l    3.1 l in [1,l_upper(i)]
!.......................................................................

            if (l.lt.1+ipstrt .or. l.gt.l_upper(i)+ipend .or. &
              ilpm1ef(i,l,ipshft1).eq.-999 .or. &
              ilpm1ef(i,l,ipshft2).eq.-999) go to 312

!.......................................................................
!l    3.1.1 cos(theta)>0
!.......................................................................

            l1=lpm1eff(l,ipshft1)
            l2=lpm1eff(l,ipshft2)
            l12=lpm1eff(l,ipshft1+ipshft2)
            ztrders=vnorm*coss(i,l)/(ipshft2*dszp5(l)-ipshft1*dszm5(l))
            do 411 j=2,jx
              spasou(i,j,k,l+ip1t1m2)=-x(j)*ztrders* &
                (fnp1(i,j,k,l2)-fnp1(i,j,k,l1))
              zdns(l)=zdns(l)+spasou(i,j,k,l+ip1t1m2)*cynt2(i,l) &
                *cint2(j)
 411        continue

!.......................................................................
!l    3.1.2 cos(theta) < 0
!.......................................................................

 312        continue

            if (l.lt.1+imstrt .or. l.gt.l_upper(i)+imend .or. &
              ilpm1ef(i,l,imshft1).eq.-999 .or. &
              ilpm1ef(i,l,imshft2).eq.-999) go to 320

            iieff=iy_(l)+1-i
            l1=lpm1eff(l,imshft1)
            il1=ilpm1ef(iieff,l,imshft1)
            l2=lpm1eff(l,imshft2)
            il2=ilpm1ef(iieff,l,imshft2)
            l12=lpm1eff(l,imshft1+imshft2)
            il12=ilpm1ef(iieff,l,imshft1+imshft2)
            ztrders=vnorm*coss(iieff,l)/ &
              (imshft2*dszp5(l)-imshft1*dszm5(l))
            do 412 j=2,jx
              spasou(iieff,j,k,l+im1t1m2)=-x(j)*ztrders* &
                (fnp1(il2,j,k,l2)-fnp1(il1,j,k,l1))
              zdns(l)=zdns(l)+spasou(iieff,j,k,l+im1t1m2)* &
                cynt2(iieff,l)*cint2(j)
 412        continue

 320        if (sbdry.ne."periodic" .or. l_upper(i).eq.setup0%ls .or. l.eq.1) &
              go to 300

!.......................................................................
!l    3.2 point ll=setup0%ls+2-l, i.e. bottom half cross-section. Not yet treated
!     in 3.1 if particle is trapped.
!     Thus: l in [setup0%ls+2-l_upper(i),setup0%ls]
!.......................................................................

            ll=setup0%ls+2-l

            if (ilpm1ef(i,ll,ipshft1).eq.-999 .or. &
              ilpm1ef(i,ll,ipshft2).eq.-999) go to 322

!.......................................................................
!l    3.2.1 cos(theta)>0, ll>setup0%ls/2+1
!.......................................................................

            l1=lpm1eff(ll,ipshft1)
            l2=lpm1eff(ll,ipshft2)
            l12=lpm1eff(ll,ipshft1+ipshft2)
            ztrders=vnorm*coss(i,ll) &
              /(ipshft2*dszp5(ll)-ipshft1*dszm5(ll))
            do 421 j=2,jx
              spasou(i,j,k,ll+ip1t1m2)=-x(j)*ztrders* &
                (fnp1(i,j,k,l2)-fnp1(i,j,k,l1))
              zdns(ll)=zdns(ll)+spasou(i,j,k,ll+ip1t1m2)*cynt2(i,ll)* &
                cint2(j)
 421        continue

!.......................................................................
!l    3.2.2 cos(theta) < 0, ll>setup0%ls/2+1
!.......................................................................

 322        continue

            if (ilpm1ef(i,ll,imshft1).eq.-999 .or. &
              ilpm1ef(i,ll,imshft2).eq.-999) go to 300

            iieff=iy_(ll)+1-i
            l1=lpm1eff(ll,imshft1)
            il1=ilpm1ef(iieff,ll,imshft1)
            l2=lpm1eff(ll,imshft2)
            il2=ilpm1ef(iieff,ll,imshft2)
            l12=lpm1eff(ll,imshft1+imshft2)
            il12=ilpm1ef(iieff,ll,imshft1+imshft2)
            ztrders=vnorm*coss(iieff,ll)/ &
              (imshft2*dszp5(ll)-imshft1*dszm5(ll))
            do 422 j=2,jx
              spasou(iieff,j,k,ll+im1t1m2)=-x(j)*ztrders* &
                (fnp1(il2,j,k,l2)-fnp1(il1,j,k,l1))
              zdns(ll)=zdns(ll)+spasou(iieff,j,k,ll+im1t1m2)* &
                cynt2(iieff,ll)*cint2(j)
 422        continue

 300      continue
 200    continue
 100  continue

!%OSc interpolate on node points
!%OS  if (mod(nummods,10) .le. 4) then
!%OS  do 210 i=1,iy_(setup0%ls)
!%OS  spasou(i,j,k,setup0%ls)=spasou(ilpm1ef(i,setup0%ls,-1),j,k,setup0%ls-1)+
!%OS  /             dszm5(setup0%ls)/(dszm5(setup0%ls)+dszm5(setup0%ls-1))*
!%OS  1           (spasou(ilpm1ef(i,setup0%ls,-1),j,k,setup0%ls-1)-
!%OS  1               spasou(ilpm1ef(ilpm1ef(i,setup0%ls,-1),setup0%ls-1,-1),j,k,setup0%ls-2))
!%OS  210        continue
!%OS  do 211 l=setup0%ls-1,2,-1
!%OS  do 211 i=1,iy_(l)
!%OS  spasou(i,j,k,l)=(1.-dls(i,j,k,l))*spasou(i,j,k,l)+
!%OS  +                      dls(i,j,k,l)*spasou(ilpm1ef(i,l,-1),j,k,l-1)
!%OS  211        continue
!%OS  do 212 i=1,iy_(1)
!%OS  spasou(i,j,k,1)=2.*spasou(i,j,k,1)-
!%OS  1                                     spasou(ilpm1ef(i,1,+1),j,k,2)
!%OS  212        continue
!%OS  endif
!%OS  200    continue
!%OS  100  continue

!%OS
!     check spasou and velsou total densities
      do 500 k=1,ngen
        do 501 j=1,jx
          do 502 l=1,lrors
            do 503 i=1,iy_(l)
              zdnp1(l)=zdnp1(l)+fnp1(i,j,k,l)*cynt2(i,l)*cint2(j)
              zdnspa(l)=zdnspa(l)+spasou(i,j,k,l)*cynt2(i,l)*cint2(j)
              zdnvel(l)=zdnvel(l)+velsou(i,j,k,l)*cynt2(i,l)*cint2(j)
 503        continue
            do 504 i=1,iyh_(l)
              ii=iymax+1-i
              ief=iy_(l)+1-i
              zdnspa2(i)=zdnspa2(i)+spasou(i,j,k,l)*cynt2(i,l)*cint2(j) &
                *dsz(l)/dsz(1)/lrors
              zdnspa2(ii)=zdnspa2(ii)+spasou(ief,j,k,l)*cynt2(ief,l)* &
                cint2(j)*dsz(l)/dsz(1)/lrors
              zdnvel2(i)=zdnvel2(i)+velsou(i,j,k,l)*cynt2(i,l)*cint2(j) &
                *dsz(l)/dsz(1)/lrors
              zdnvel2(ii)=zdnvel2(ii)+velsou(ief,j,k,l)*cynt2(ief,l)* &
                cint2(j)*dsz(l)/dsz(1)/lrors
              zdnuspa(j)=zdnuspa(j)+(spasou(i,j,k,l)*cynt2(i,l)+ &
                spasou(ief,j,k,l)*cynt2(ief,l))* &
                cint2(j)*dsz(l)/dsz(1)/lrors
              zdnuvel(j)=zdnuvel(j)+(velsou(i,j,k,l)*cynt2(i,l)+ &
                velsou(ief,j,k,l)*cynt2(ief,l))* &
                cint2(j)*dsz(l)/dsz(1)/lrors
 504        continue
 502      continue
 501    continue
 500  continue
!
      zp1to=0.0
      zspato=0.0
      zvelto=0.0
      do 510 l=1,lrors
        zp1to=zp1to+dsz(l)/dsz(1)*zdnp1(l)/lrors
        zspato=zspato+dsz(l)/dsz(1)*zdnspa(l)/lrors
        zvelto=zvelto+dsz(l)/dsz(1)*zdnvel(l)/lrors
 510  continue

      zspa2to=0.0
      zvel2to=0.0
      do 511 i=1,iymax
        zspa2to=zspa2to+zdnspa2(i)
        zvel2to=zvel2to+zdnvel2(i)
 511  continue
      do 512 i=1,iymax/2
        zspa2i(i)=zdnspa2(i)+zdnspa2(iymax+1-i)
        zvel2i(i)=zdnvel2(i)+zdnvel2(iymax+1-i)
 512  continue
!%OS

      return
      end subroutine wparsou

end module wparsou_mod
