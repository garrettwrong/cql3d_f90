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

module pack21_mod

  !---BEGIN USE

  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_float

  !---END USE

!
!

!contains
!YuP: if empty, gives a compilation error in IntelFortran

end module pack21_mod


!XXXX this is a mess (relates to netcdf files)
      subroutine pack21(a,ibot,itop,jbot,jtop,b,iy,jx)
        use iso_c_binding, only : c_double
        use iso_c_binding, only : c_float
        use bcast_mod, only : bcast
        use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     It sometimes becomes necessary to take a
!     2-D array dimensioned ibot:itop by jbot:jtop
!     and repack it as though it were
!     dimensioned 1:iy by 1:jx, starting at a(1,1).
!     This routine does this, transfering relevant data
!     from array a to b.
!.......................................................................

      save
      dimension a(ibot:itop,jbot:jtop)
      dimension b(iy*jx)
      b=0.d0
      do 1 j=1,jx
        i1=(j-1)*iy+1
        call dcopy(iy,a(1:iy,j),1,b(i1:i1+iy-1),1) ! a-->b
 1    continue
      return
      end
!
!
      subroutine unpack21(a,ibot,itop,jbot,jtop,b,iy,jx)
        use iso_c_binding, only : c_double
        use iso_c_binding, only : c_float
        use bcast_mod, only : bcast
        use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This is the inverse of subroutine pack21.
!     It sometimes becomes necessary to take a
!     2-D array b dimensioned 1:iy,1:jx starting at b(1,1)
!     and repack it as though it were dimensioned
!     ibot:itop by jbot:jtop starting at a(ibot,jbot).
!     This routine does this, transfering relevant data
!     from array b to a.
!.......................................................................

      save
      dimension a(ibot:itop,jbot:jtop)
      dimension b(iy*jx)
      zero=0.
      a=zero !YuP[2019-06-08]was call bcast(a,zero,(itop-ibot+1)*(jtop-jbot+1))
      do 1 j=1,jx
        i1=(j-1)*iy+1
        call dcopy(iy,b(i1:i1+iy-1),1,a(1:iy,j),1) ! b-->a
 1    continue
      return
      end
!
!
      subroutine ipack21(ia,ibot,itop,jbot,jtop,ib,iy,jx)
        use iso_c_binding, only : c_double
        use iso_c_binding, only : c_float
        implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     It sometimes becomes necessary to take a
!     2-D array dimensioned ibot:itop by jbot:jtop
!     and repack it as though it were
!     dimensioned 1:iy by 1:jx, starting at a(1,1).
!     This routine does this, transfering relevant data
!     from array a to b.
!     This routine does this, transfering relevant data from
!     array ia to ib.
!.......................................................................

      save
      dimension ia(ibot:itop,jbot:jtop)
      dimension ib(iy*jx)
      do j=1,jx
        i1=(j-1)*iy+1
        do i=0,iy-1
           ib(i1+i)=ia(1+i,j)
        enddo
      enddo
      return
      end
