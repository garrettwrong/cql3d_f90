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

module tdtry_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdwrng_mod, only : tdwrng

  !---END USE

!
!

contains

  subroutine tdtry
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)


!.......................................................................
!     Called by micxinit, which is called repeatedly from tdinit with
!     ll=lrors,1,-1;  tdnflxs(ll).
!
!     This routine computes the y (pitch angle) mesh
!     from the normalized mesh determined in tdtrmuy or wptrmuy.
!     In general the y mesh will vary from flux surface to flux surface.
!.......................................................................

      if (eqsym.eq."none") then  !Only if eqmod.eq."enabled",
                                 !        eqsource="eqdsk"
         ilzhfs=lz_bmax(lr_)
      else
         ilzhfs=lz
      endif
      if (numclas .eq. 1) ilzhfs=lz/2+1

      write(*,*)'tdtry:  ipacktp =',ipacktp

      if (setup0%cqlpmod .ne."enabled") then
        zsntrp2=1./bbpsi(ilzhfs,lr_)
!BH070419        ipacktp=3
!BH070419 NOTE:  ipacktp now set in tdtrmuy, to coordinate with iytr.
!BH070419        Now can be 0 or 3 here, depending on tfac sign.
      else
        zsntrp2=psis(l_)/bbpsi(ilzhfs,lr_)
!BH070419        ipacktp=0
!BH070419 NOTE:  ipacktp now set in tdtrmuy, to coordinate with iytr.
      endif
      thb(l_)=asin(sqrt(zsntrp2))
      iadd=0
      iu=0
      mark=0

!..............................................................
!     constant mu mesh first.
!..............................................................

      if (meshy.eq."fixed_mu") then
        do 10 i=1,iymax/2
          if (setup0%cqlpmod .ne."enabled") then
            sinsq=mun(i)*bmidplne(lr_)/bmidplne(setup0%lrzmax)
          else
            sinsq=mun(i)*psis(l_)
          endif
          if (sinsq.ge.1.) then
            if(mark.eq.0 .and. setup0%cqlpmod.ne."enabled") then
              call tdwrng(8)
            else
              it=i-1
              iyh_(l_)=it+ipacktp
              iyh=it+ipacktp
              go to 11
            endif
          elseif (sinsq.ge.zsntrp2 .and. mark.eq.0) then
            if (setup0%cqlpmod .ne. "enabled") then
              if (sinsq.eq.zsntrp2) call tdwrng(9)
              iu=i-1
            else
              iu=i-1
            endif
            mark=1
            iadd=ipacktp
          endif
          ix=i+iadd
          idx(i,l_)=ix
          y(ix,l_)=asin(sqrt(sinsq))
 10     continue
        if (mark.eq.0) call tdwrng(8)
        it=iymax/2-ipacktp
        iyh_(l_)=iymax/2
        iyh=iyh_(l_)
 11     continue
        iyy=2*iyh   !-YuP-101215: Don't use iy=; it's in common /params/
                    ! Don't let overwrite the cqlinput value!
        iy_(l_)=iyy
!BH070419        if (setup0%cqlpmod .ne. "enabled") then
!BH070419          itl=iu+2
!BH070419        else
!BH070419          if (iu .eq. 0) iu=iyh
!BH070419          itl=iu
!BH070419        endif
        if (ipacktp.eq.3) then
           itl=iu+2
        elseif (ipacktp.eq.0) then
           if (iu .eq. 0) iu=iyh
           itl=iu
        else
           WRITE(*,*)'tdtry:  ipacktp =',ipacktp
           WRITE(*,*)'STOP in tdtry:  Check ipacktp'
           stop
        endif
        itl_(l_)=itl
        itu=iyy+1-itl
        itu_(l_)=itu

!..............................................................
!     Now for the fixed theta mesh case.
!..............................................................

      elseif (meshy.eq."fixed_y") then
        iyh=iyh_(lrors)
        iyh_(l_)=iyh
        iyy=iymax  !-YuP-101215: Don't use iy=; it's in common /params/
                       ! Don't let overwrite the cqlinput value!
        iy_(l_)=iyy
        do 13 i=1,iytr(lrors)/2
          if (mun(i).gt.thb(l_) .and. mark.eq.0) then
            if (setup0%cqlpmod .ne. "enabled") iu=i-1
            if (setup0%cqlpmod .eq. "enabled") iu=i-1
            mark=1
            iadd=ipacktp
          endif
          ix=i+iadd
          y(ix,l_)=mun(i)
          idx(i,l_)=ix
 13     continue
!BH070419        if (setup0%cqlpmod .ne. "enabled") itl=iu+2
!BH070419        if (setup0%cqlpmod .eq. "enabled") itl=iu
        if (ipacktp.eq.3) then
           itl=iu+2
        elseif (ipacktp.eq.0) then
           if (iu .eq. 0) iu=iyh
           itl=iu
        else
           WRITE(*,*)'tdtry:  ipacktp =',ipacktp
           WRITE(*,*)'STOP in tdtry:  Check ipacktp'
           stop
        endif

        itl_(l_)=itl
        itu=iyy+1-itl
        itu_(l_)=itu
        if (mark.eq.0) call tdwrng(10)

      else if (meshy.eq."free") then
        return
      endif

      iytr(l_)=iyy-2*ipacktp
      iyjx=iyy*jx
      iyjx_(l_)=iyjx

!..............................................................
!     Make sure the passed/trapping interface is between two meshpoints.
!..............................................................

!BH070419      if (setup0%cqlpmod .ne. "enabled") then
      if (ipacktp.eq.3) then
        diff1=thb(l_)-y(iu,l_)
        if (diff1.lt.tbnd(l_)) tbnd(l_)=diff1/3.
        ytop=y(iu+4,l_)
        if (iyh.eq.itl+1) ytop=pi/2.
        diff2=ytop-thb(l_)
        if (diff2.lt.tbnd(l_)) tbnd(l_)=diff2/3.
        y(iu+1,l_)=thb(l_)-tbnd(l_)
        y(iu+2,l_)=thb(l_)
        y(iu+3,l_)=thb(l_)+tbnd(l_)
      endif

!BH070419 for tfac.lt.0., have set ipacktp=0, i.e., pack zero additional
!BH070419 y-mesh points around the t-p boundary.   Will simply adjust
!BH070419 a y-mesh point to be at t-p boundary (not sure this is
!BH070419 essential) and space surrounding mesh points equidistant
!BH070419 away.   This is to give same velocity and radial transport
!BH070419 related y-meshes.
!BH070419 The y-meshes will not be identical at t-p bndry, but
!BH070419 will ignore this for now [Probably have to have y(itl)
!BH070419 exactly at thb(l_), but not certain. Check this out later....]

      if (setup0%cqlpmod .ne. "enabled" .and. ipacktp.eq.0) then

         deltay1=abs(thb(l_)-y(iu,l_))
         deltay2=abs(y(iu+1,l_)-thb(l_))
         if (deltay1.le.deltay2) then
            y(itl,l_)=thb(l_)
            y(itl-1,l_)=y(itl,l_)-deltay2  ! i.e., y(itl+/-1) equidistant
         else
            y(itl+1,l_)=thb(l_)
            y(itl+2,l_)=y(itl+1,l_)+deltay1
            itl=itl+1                !reset itl,itu
            itl_(l_)=itl
            itu=iyy+1-itl
            itu_(l_)=itu
         endif
      endif

!..............................................................
!     Define the mesh on the other side of pi/2
!..............................................................

      do 12 i=1,iyh
        ii=iyy+1-i
        y(ii,l_)=pi-y(i,l_)
 12   continue
      do 14 i=1,iytr(l_)/2
        ii=iytr(lrors)+1-i
        idx(ii,l_)=iy_(l_)+1-idx(i,l_)
 14   continue
      if (l_.eq.1) then
        do 15 i=1,iymax
          idx(i,0)=1
 15     continue
      endif

!cc      write(*,*)'tdtry: l_,iytr(l_),idx()=',
!cc     +                  l_,iytr(l_),idx(1:iytr(l_),l_)


      return
      end subroutine tdtry


end module tdtry_mod
