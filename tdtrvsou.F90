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

module tdtrvsou_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use diagentr_mod, only : gfi
  use tdtrvtor3_mod, only : tdtrvtor3

  !---END USE

!
!

contains

      subroutine tdtrvsou(k)
      use param_mod
      use cqlcomm_mod
      use advnce_mod !here: in tdtrvsou.  To get gfi(),hfi(),gfu(),hfu(),etc.
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     compute source term, due to velocity operator evaluated on
!     f_n+1/2, needed for the transport equation solved with ADI.
!
!     For up/down symmetry cases, also compute the value of f
!     at l_=1 and setup0%ls (if sbdry.ne."periodic")
!..............................................................

      dimension zdns(lrorsa),zdns1(lrorsa)

      fpithta(i,j)=f(i+1,j,k,l_)*(1.-dithta(i,j,l_)) + &
                   f(i  ,j,k,l_)*dithta(i,j,l_)
!.......................................................................

!     include cthta in vsou at n=nontran-1/2
!%OS  if (n .eq. nontran) call wpcthta

      zdns(l_)=0.0
      zdns1(l_)=0.0

!.......................................................................
!l    1. Compute velsou for main mesh points
!.......................................................................

      do 100 i=1,iy
        if ((i.eq.itl.or.i.eq.itu) .and. setup0%cqlpmod.ne."enabled") go to 100
!%OS  should be j=2,jx-1? may depend on lbdry
!%OS  do 110 j=1,jx
        do 110 j=2,jx-1
          velsou(i,j,k,l_)= qz(j)*(gfi(i,j,k)-gfi(i,j-1,k))+ &
                          ry(i,j,l_)*(hfi(i,j,k,l_)-hfi(i-1,j,k,l_))+ &
            vptb(i,lr_)*(cah(i,j)*f(i,j,k,l_)+so(i,j))+ &
            cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
          velsou2(i,j,k,l_)=cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
          zdns(l_)=zdns(l_)+velsou(i,j,k,l_)*cynt2(i,l_)*cint2(j)
          zdns1(l_)=zdns1(l_)+velsou2(i,j,k,l_)*cynt2(i,l_)*cint2(j)
 110    continue
        if (setup0%cqlpmod.eq."enabled" .and. updown.eq."symmetry" .and. &
          sbdry.ne."periodic") then
          if (l_ .le. 2) then
            do 111 j=2,jx-1
!%OS  fedge(i,j,k,l_)=velsou(i,j,k,l_)*dtreff-
!%OS  -                      (dtreff+0.5*dszp5(1)/vnorm/x(j)/coss(i,l_))*
!%OS  *                          cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
              fedge(i,j,k,l_)=-dszp5(1)/vnorm/x(j)/coss(i,l_)* &
                cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
 111        continue
          else if (l_ .ge. lrors-1) then
            do 112 j=2,jx-1
!%OS  fedge(i,j,k,l_-lrors+4)=velsou(i,j,k,l_)*dtreff-
!%OS  -                     (dtreff+0.5*dszm5(setup0%ls)/vnorm/x(j)/coss(i,l_))*
!%OS  *                          cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
              fedge(i,j,k,l_-lrors+4)=-dszm5(setup0%ls)/vnorm/x(j)/ &
                coss(i,l_)*cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
 112        continue
          endif
        endif
 100  continue

!..................................................................
!l    2. Compute velsou at pass/trapped boundary
!..................................................................

      if (setup0%cqlpmod .ne. "enabled") then
        do 200 j=1,jx
          velsou(itl,j,k,l_)=qz(j)*(gfi(itl,j,k)-gfi(itl,j-1,k))+ &
            r2y(j,l_)*(-hfi(itl-1,j,k,l_)+2.*hfi(itl,j,k,l_)+hfi(itu,j,k,l_))+ &
            vptb(itl,lr_)*(cah(itl,j)*f(itl,j,k,l_)+so(itl,j))
 200    continue
      endif

!..................................................................
!l    2.1 Symmetry about pi/2. in trapped region.
!..................................................................

      if (symtrap .eq. "enabled") then
        do 210 i=iyh+1,itu
          ii=iy+1-i
          do 215 j=1,jx
            velsou(i,j,k,l_)=velsou(ii,j,k,l_)
 215      continue
 210    continue
      endif

!.......................................................................
!l    3. Interpolate on radial velocity mesh, keeping integral of velsou
!l    over velocity conserved (for radial transport case)
!.......................................................................

      if (setup0%cqlpmod .ne. "enabled") &
        call tdtrvtor3(velsou(0:iyp1,0:jxp1,1:ngen,1), &
                       velsou(0:iyp1,0:jxp1,1:ngen,1), cynt2,cynt2_,2,k)
      !YuP:This subr. uses internal loops in 0:iyp1,0:jxp1,1:ngen, [for given l_]

      return
      end subroutine tdtrvsou


end module tdtrvsou_mod
