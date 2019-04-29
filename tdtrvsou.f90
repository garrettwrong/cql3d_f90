module tdtrvsou_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use diagentr_mod, only : gfi
  use tdtrvtor3_mod, only : tdtrvtor3

  !---END USE

!
!

contains

      subroutine tdtrvsou(k)
      use param_mod
      use comm_mod
      use advnce_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     compute source term, due to velocity operator evaluated on
!     f_n+1/2, needed for the transport equation solved with ADI.
!
!     For up/down symmetry cases, also compute the value of f
!     at l_=1 and ls (if sbdry.ne."periodic")
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
        if ((i.eq.itl.or.i.eq.itu) .and. cqlpmod.ne."enabled") go to 100
!%OS  should be j=2,jx-1? may depend on lbdry
!%OS  do 110 j=1,jx
        do 110 j=2,jx-1
          velsou(i,j,k,l_)=qz(j)*(gfi(i,j,k)-gfi(i,j-1,k))+ &
            ry(i,j)*(hfi(i,j)-hfi(i-1,j))+ &
            vptb(i,lr_)*(cah(i,j)*f(i,j,k,l_)+so(i,j))+ &
            cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
          velsou2(i,j,k,l_)=cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
          zdns(l_)=zdns(l_)+velsou(i,j,k,l_)*cynt2(i,l_)*cint2(j)
          zdns1(l_)=zdns1(l_)+velsou2(i,j,k,l_)*cynt2(i,l_)*cint2(j)
 110    continue
        if (cqlpmod.eq."enabled" .and. updown.eq."symmetry" .and. &
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
!%OS  -                     (dtreff+0.5*dszm5(ls)/vnorm/x(j)/coss(i,l_))*
!%OS  *                          cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
              fedge(i,j,k,l_-lrors+4)=-dszm5(ls)/vnorm/x(j)/ &
                coss(i,l_)*cthta(i,j)*(fpithta(i,j)-fpithta(i-1,j))
 112        continue
          endif
        endif
 100  continue

!..................................................................
!l    2. Compute velsou at pass/trapped boundary
!..................................................................

      if (cqlpmod .ne. "enabled") then
        do 200 j=1,jx
          velsou(itl,j,k,l_)=qz(j)*(gfi(itl,j,k)-gfi(itl,j-1,k))+ &
            r2y(j)*(-hfi(itl-1,j)+2.*hfi(itl,j)+hfi(itu,j))+ &
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

      if (cqlpmod .ne. "enabled") &
        call tdtrvtor3(velsou(0:iyp1,0:jxp1,1,1),velsou(0:iyp1,0:jxp1,1,1),cynt2,cynt2_,2,k)

      return
      end
end module tdtrvsou_mod
