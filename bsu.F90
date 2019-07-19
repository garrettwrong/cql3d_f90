! Copyright 2019 Garrett Wright, Princeton Plasma Physics Laboratory,
!    contracted by the U.S. Department of Energy (putnumberhere).
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

module  bsu_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use lookup_mod, only : lookup

  !---END USE
  use iso_c_binding, only : c_double

contains

  real(c_double) function bsu(jj,kk,ll)
    use param_mod
    use cqlcomm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save
    !      dimension f_(0:iy+1,0:jx+1,ngen,setup0%lrz)
    !      f_ is passed in common
    allocatable :: bsu_s(:,:)

    !.........................................................................
    !     This routine is used to compute the skewing effect at the p/t bndry
    !     due to the bootstrap effect. It is not used if bootcalc="disabled",
    !     or if advnce="explicit" or if setup0%lrz=1
    !     Subroutine bsu is for the upper theta tp-bndry at i=itu.
    !     There is also a separate subroutine bsl for the lower itl bndry.
    !.........................................................................

    bsu=0.
    if (bootcalc.eq."disabled") return
    if (implct.ne."enabled" .or. setup0%lrz.eq.1 .or. n.lt.nonboot) return

    if (.NOT. ALLOCATED(bsu_s)) then
       allocate( bsu_s(0:jx+1,setup0%lrz) )
       call bcast(bsu_s,zero,(jx+2)*setup0%lrz)
    endif

    qb_mc=bnumb(kk)*charge*bthr(ll)/(fmass(kk)*clight)

    jjj=max(jj,1)   ! to limit from below, for x(jjj)
    jjj=min(jjj,jx) ! to limit from above, for x(jjj)
    rban=cursign*x(jjj)*coss(itu_(ll),ll)*vnorm/qb_mc

    if (bootcalc.eq."method1") then

       if (n.eq.nonboot.or.bootupdt.eq."enabled") then
          if (ll.eq.1) then
             !            dfdr=(f_(itu_(ll+1),jj,kk,ll+1)-f_(itu_(ll),jj,kk,ll))/
             !     1      (rz(ll+1)-rz(ll))
             p1=rz(ll+1)-rz(ll)
             p2=rz(ll+2)-rz(ll+1)
             p3=rz(ll+2)-rz(ll)
             dfdr=-(p1+p3)/(p1*p3)*f_(itu_(ll),jj,kk,ll) &
                  +p3/(p1*p2)*f_(itu_(ll+1),jj,kk,ll+1) &
                  -p1/(p2*p3)*f_(itu_(ll+2),jj,kk,ll+2)
          elseif (ll.eq.setup0%lrz) then
             !            dfdr=(f_(itu_(ll),jj,kk,ll)-f_(itu_(ll-1),jj,kk,ll-1))/
             !     1      (rz(ll)-rz(ll-1))
             p1=rz(ll-1)-rz(ll-2)
             p2=rz(ll)-rz(ll-1)
             p3=rz(ll)-rz(ll-2)
             dfdr=+p2/(p1*p3)*f_(itu_(ll-2),jj,kk,ll-2) &
                  -p3/(p1*p2)*f_(itu_(ll-1),jj,kk,ll-1) &
                  +(p2+p3)/(p2*p3)*f_(itu_(ll),jj,kk,ll)
          else
             !            dfdr=(f_(itu_(ll+1),jj,kk,ll+1)-f_(itu_(ll-1),jj,kk,ll-1))/
             !     1      (rz(ll+1)-rz(ll-1))
             p1=rz(ll)-rz(ll-1)
             p2=rz(ll+1)-rz(ll)
             p3=rz(ll+1)-rz(ll-1)
             dfdr=-p2/(p1*p3)*f_(itu_(ll-1),jj,kk,ll-1) &
                  -(p1-p2)/(p1*p2)*f_(itu_(ll),jj,kk,ll) &
                  +p1/(p2*p3)*f_(itu_(ll+1),jj,kk,ll+1)
          endif
          bsu_=-bootsign*dfdr*rban
          bsu_s(jj,ll)=bsu_
          !          if(jj.eq.30 .and. (n.eq.nonboot .or. n.eq.4) )
          !     +      write(*,'(a,2i3,e11.3,4e12.4)')
          !     +     'kk,ll, dfdr, f_(ll=2,4,6,8)=',
          !     +      kk,ll, dfdr, f_(itu_(2),jj,kk,2),f_(itu_(4),jj,kk,4),
          !     +                   f_(itu_(6),jj,kk,6),f_(itu_(8),jj,kk,8)
       else
          bsu_=bsu_s(jj,ll)
       endif

       !     Limit the jump at trapped-passing boundary to 0.2*(distn functn).
       !     If getting larger values, should probably at least consider
       !     method2.  (BH).
       bsum=abs(0.2*f_(itu_(ll),jj,kk,ll))
       bsu=sign(one,bsu_)*min(bsum,abs(bsu_))

    elseif (bootcalc.eq."method2") then

       !       Assuming positive current here (should be generalized).
       rrr=rz(ll)-rban
       if (rrr.gt.rz(setup0%lrzmax)) rrr=rz(setup0%lrzmax)
       !BH080714:  rz dimensioned 0:lrza
       call lookup(rrr,rz(1),setup0%lrzmax,weightu,weightl,irrr)
       !BH051101        if (irrr.eq.lzrmax+1) then
       !$$$        if (irrr.ge.setup0%lrzmax+1) then
       !$$$          irrr=setup0%lrzmax+1
       !$$$          bsu=f_(itu_(irrr-1),jj,kk,irrr-1)-f_(itu_(ll),jj,kk,ll)
       !$$$        else
       bsu=weightl*f_(itu_(irrr-1),jj,kk,irrr-1) &
            +weightu*f_(itu_(irrr),jj,kk,irrr) &
            -f_(itu_(ll),jj,kk,ll)
       !$$$        endif

       bsu_s(jj,ll)=bsu

    endif

    return
  end function bsu

end module bsu_mod
