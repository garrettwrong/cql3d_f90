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

module sounorm_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use soup_mod, only : soup

  !---END USE


!*****************************************************************

!
!

contains

      subroutine sounorm
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine establishes an array of normalization constants
!     which allows the code to force a particular poloidal source
!     profile. It is called only at initialization.
!..................................................................




!..................................................................
!     Set a flag.
!..................................................................

      isounor=1
      do 200 l=1,lz
        do 100 k=1,ngen
          call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)
          do 50 m=1,nso

!..................................................................
!     Call a routine which determines the non-normalized source
!     distribution for a given i and all j.
!..................................................................

            do 40 i=1,iy
              call soup(coss(i,l_),l,k,m)
              do 30 j=1,jx
                temp1(i,j)=soupp(j,lr_)
 30           continue
 40         continue

!..................................................................
!     Set the array sounor such that when the computed source
!     profile is multiplied by sounor the resulting current
!     density is unity.
!..................................................................

            s=0.
            do 10 i=1,iy
              do 20 j=1,jx
                s=s+temp1(i,j)*cynt2(i,l_)*cint2(j)
 20           continue
 10         continue
            if (s.ne.zero) sounor(k,m,l,lr_)=1./(s*one_)
 50       continue
 100    continue
 200  continue

!..................................................................
!     reset the flag (subroutine soup)
!..................................................................

      isounor=0
      return
      end subroutine sounorm


end module sounorm_mod
