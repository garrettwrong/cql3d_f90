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

module tdtransp_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use r8subs_mod, only : dcopy
  use tdtravct_mod, only : tdtravct
  use tdtrchk_mod, only : tdtrchk
  use tdtrrtov2_mod, only : tdtrrtov2
  use tdtrrtov_mod, only : tdtrrtov
  use tdtrsym_mod, only : tdtrsym
  use tdtrvtor2_mod, only : tdtrvtor2
  use tdtrvtor_mod, only : tdtrvtor
  use tdtrwtl_mod, only : tdtrwtl

  !---END USE


contains

      subroutine tdtransp
      use param_mod
      use cqlcomm_mod
      use cqlconf_mod, only : setup0
      use r8subs_mod, only : cvmgt, dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
!..............................................................
!     Time advancement splitting step for radial transport.
!..............................................................
      include 'trans.h'
      nobind="disabled"
!.......................................................................
     !YuP: frn,fvn,fvn_1,frn_1,frn_2 are dimensioned as 0:iyp1,0:jxp1,1:ngen,0:lrors
     !                          and f is dimensioned as 0:iy+1,0:jx+1,1:ngen,1:lrors
      call dcopy(iyjx2*ngen*lrors,f(0:iy+1,0:jx+1,1:ngen,1:lrors),1, &
                                fvn(0:iy+1,0:jx+1,1:ngen,1:lrors),1)
!..............................................................
!     The velocity split for this time step is done and we assume
!     that the results of this split resides in fvn(i,j,k,l) and
!     that this is defined on the velocity or complete mesh. Since
!     the tranport equation is not advanced in the vicinity of the
!     pass/trapped boundary we first interpolate onto a new mesh
!     that that does not include these mesh points, and this is
!     accomplished in a density conserving manner.
!..............................................................
      if (nonadi .eq. 6) then
         !XXX YuP:This subr. uses internal loops in 0:iyp1,0:jxp1,1:ngen,1:lrors(or setup0%lrz)
        call tdtrvtor2(fvn(0:iyp1,0:jxp1,1:ngen,1), &
                       frn(0:iyp1,0:jxp1,1:ngen,1), vpint,vpint_,1)
      else
        !     YuP:This subr. uses internal loops in 0:iyp1,0:jxp1,1:ngen,0:lrors(or setup0%lrz)
        call tdtrvtor(fvn,frn)
      endif
!..............................................................
!     Generate advection coefficients - then determine the
!     Chang-Cooper weights (radial direction).
!     Presently only set up for 1 general species,
!       except for pinch.eq."disabled" (adv()=0.).
!..............................................................
      kprofile=kelecm
      if (niong.ne.0) kprofile=kionm(1)
      if (colmodl.eq.0) kprofile=1
      call tdtravct(frn,kprofile,1)  !No-op if pinch="disabled"
      call tdtrwtl
!..............................................................
!     Keep a copy for accuracy check later...
!..............................................................
     !YuP: frn,fvn,fvn_1,frn_1,frn_2 are dimensioned as 0:iyp1,0:jxp1,1:ngen,0:lrors
      call dcopy(iyjx2*ngen*(lrors+1),frn(0:iy+1,0:jx+1,1:ngen,0:lrors),1, &
                                    frn_2(0:iy+1,0:jx+1,1:ngen,0:lrors),1)
!..............................................................
!     Loop over species index.
!..............................................................
      do 100 k=1,ngen
!..............................................................
!     Loop over momentum
!..............................................................
        do 200 ic=1,iytr(lrors)/2
          if (l_lower(ic).eq.lrors) go to 200
          i_=iytr(lrors)+1-ic
          do 201 ii=1,2
            i=ic
            if (ii.eq.2) i=i_
            lu=l_lower(i)
!..............................................................
!     Initialize the calculation (zero flux) at the innermost
!     radial mesh point seen by a particle with a given adiabatic
!     invariant if meshy="fixed_mu" or at l=1 if meshy="fixed_y".
!     In the event that there is a binding equation ((lpt(i).ne.
!     lrors) meaning a particle with a given "i" is passing at the center
!     (l=1) and at some point (lpt(i_)) becomes trapped) we will
!     also initialize at l=lrors. There will then be two initial
!     sweeps, instead of just one. One sweep will be up in "l" the
!     other down, and both terminate at lpt(i).
!..............................................................
            do 222 j=1,jx
              fg_(j,ii,lu)=delr(i,lu)/betr(i,lu)
              eg_(j,ii,lu)=alpr(i,lu)/betr(i,lu)
              frn(idx(i,lrors),j,k,lrors)=vptb_(idx(i,lrors), &
                setup0%lrindx(lrors))/zmaxpsi(setup0%lrindx(lrors)) &
                *frn(idx(i,lrors),j,k,lrors)
              eg_(j,ii,lrors)=0.
              fg_(j,ii,lrors)=frn(idx(i,lrors),j,k,lrors)
 222        continue
!..................................................................
!     Now complete the initial sweep over 1 (or 2) regions.
!..................................................................
            do 240 l=l_lower(i)+1,lpt(i)-1
              do 230 j=1,jx
                eg_(j,ii,l)=alpr(i,l)/(betr(i,l)-gamr(i,l) &
                  *eg_(j,ii,l-1))
                fg_(j,ii,l)=(delr(i,l)+gamr(i,l)*fg_(j,ii,l-1)) &
                  /(betr(i,l)-gamr(i,l)*eg_(j,ii,l-1))
 230          continue
 240        continue
            if (lpt(i).eq.lrors) go to 256
!..................................................................
!     second region, if necessary...
!..................................................................
            do 250 l=lrors-1,lpt(i)+1,-1
              do 255 j=1,jx
                eg_(j,ii,l)=gamr(i,l)/(betr(i,l)-alpr(i,l) &
                  *eg_(j,ii,l+1))
                fg_(j,ii,l)=(delr(i,l)+alpr(i,l) &
                  *fg_(j,ii,l+1)) &
                  /(betr(i,l)-alpr(i,l)*eg_(j,ii,l+1))
 255          continue
 250        continue
 256        continue
 201      continue
!..................................................................
!     The binding equation (at lpt(i)) is next..
!..................................................................
          i=ic
          if (lpt(i).eq.lrors) go to 260
          lv=lpt(i)
          if (nobind.eq."disabled") then

            do 265 j=1,jx
              tam1(j)=-(gamr(i,lv)*eg_(j,1,lv-1)+gamr(i_,lv) &
                *eg_(j,2,lv-1))*.5 &
                +(betr(i,lv)+betr(i_,lv))*.5 - alpr(i,lv)*eg_(j,1,lv+1)
              tam2(j)=(gamr(i,lv)*fg_(j,1,lv-1)+gamr(i_,lv) &
                *fg_(j,2,lv-1))*.5 &
                +alpr(i,lv)*fg_(j,1,lv+1)+delr(i,lv)
              frn(idx(i,lpt(i)),j,k,lpt(i))=tam2(j)/tam1(j)
              frn(idx(i_,lpt(i_)),j,k,lpt(i_))=tam2(j)/tam1(j)
 265        continue
          else
            do 266 j=1,jx
              tam1(j)=-(gamr(i,lv)*eg_(j,1,lv-1)) &
                +(betr(i,lv)) - alpr(i,lv)*eg_(j,1,lv+1)
              tam2(j)=(gamr(i,lv)*fg_(j,1,lv-1)) &
                +alpr(i,lv)*fg_(j,1,lv+1)+delr(i,lv)
              frn(idx(i,lpt(i)),j,k,lpt(i))=tam2(j)/tam1(j)
              tam3(j)=-(gamr(i_,lv)*eg_(j,2,lv-1)) &
                +(betr(i_,lv)) - alpr(i_,lv)*eg_(j,2,lv+1)
              tam4(j)=(gamr(i_,lv)*fg_(j,2,lv-1)) &
                +alpr(i_,lv)*fg_(j,2,lv+1)+delr(i_,lv)
              frn(idx(i_,lpt(i_)),j,k,lpt(i_))=tam4(j)/tam3(j)
 266        continue
          endif
 260      continue
!..................................................................
!     solve for the new distribution...
!..................................................................
          do 285 ij=1,2
            i=ic
            if (ij.eq.2) i=i_
            do 280 l=lpt(i)-1,l_lower(i),-1
              do 270 j=1,jx
                frn(idx(i,l),j,k,l)=eg_(j,ij,l)*frn(idx(i,l+1),j,k,l+1) &
                  +fg_(j,ij,l)
 270          continue
 280        continue
            do 295 l=lpt(i)+1,lrors
              do 290 j=1,jx
                frn(idx(i,l),j,k,l)=eg_(j,ij,l)*frn(idx(i,l-1),j,k,l-1) &
                  +fg_(j,ij,l)

 290          continue
 295        continue

            do 310 l=l_lower(i),lrors
              ii=idx(i,l)
              do 312 j=1,jx
                frn(ii,j,k,l)=frn(ii,j,k,l)/vptb_(ii,setup0%lrindx(l))* &
                  zmaxpsi(setup0%lrindx(l))
 312          continue
 310        continue

 285      continue

 200    continue
 100  continue

      call tdtrchk

      if (nobind.eq."enabled") then
        call tdtrsym
      endif
!.................................................................
!     Redefine frn at v=0 so it is unique.
!.................................................................
      if (nonadi .ne. 3) then
        do 400 k=1,ngen
          do 410 l=1,lrors
            zs=0.
            zt=0.
            do 420 i=1,iytr(l)
              zs=zs+vpint_(idx(i,l),setup0%lrindx(l))
              zt=zt+vpint_(idx(i,l),setup0%lrindx(l))*frn(idx(i,l),1,k,l)
 420        continue
            do 430 i=1,iytr(l)
              frn(idx(i,l),1,k,l)=zt/zs
 430        continue
 410      continue
 400    continue
      endif
!......................................................................
!     The distribution frn is currently defined on the transport
!     mesh - interpolate onto the velocity mesh, returning it in frn.
!......................................................................
      if (nonadi .eq. 6) then
        !     YuP:This subr. uses internal loops in 0:iyp1,0:jxp1,1:ngen,1:lrors(or setup0%lrz)
        call tdtrrtov2(frn(0:ipy1,0:jxp1,1:ngen,1), &
                       frn(0:ipy1,0:jxp1,1:ngen,1), vpint,vpint_,1)
      else
        call tdtrrtov(frn)
      endif
      return
      end subroutine tdtransp


end module tdtransp_mod
