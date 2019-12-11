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

module sourcpwr_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine sourcpwr(k)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif


!..................................................................
!     This routine computes the Neutral Beam power in Watts/cc
!..................................................................

      call bcast(tam1,zero,jx)
      do 501 i=1,iy
        do 502 j=1,jx
         if(gone(i,j,k,lr_).eq.zero)then !YuP[2017-11-21] added.
          !YuP: added gone() check: ONLY COUNT not-lost particles.
          !YuP: This is similar to the CQL3D-FOW version (added on 03/23/2015)
          tam1(j)=tam1(j)+source(i,j,k,indxlr_)*cynt2(i,l_)*vptb(i,lr_)
          !Note: tam1 is needed for sorpw_nbi below, which is only used
          !for plots and NetCDF file.
          !It has no effect on the source operator or solution.
          !The consequence of adding if(gone..) condition is
          !the drop of NBI power profile at plasma edge
          !(if lossmode is not disabled) because of banana+Larm.rad losses;
          !see tdpltjop, plot 'FSA SOURCE POWER DEN: ...';
          !Also, 'Power from NBI', and 'Power integrated over radius'
          !become smaller (values are printed in *.ps file).
         endif
 502    continue
 501  continue
      s=0.
      do 503 j=1,jx
        s=s+tam1(j)*tcsgm1(j)*cint2(j)
 503  continue

!..................................................................
!     Source power of species k in Watts/cc averaged over the flux surfa
!..................................................................

      sorpw_nbi(k,lr_)=s*fions(k)*one_*1.6022e-16*zmaxpsii(lr_) !-YuP->NBI source

#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      if(setup0%verbose>0) WRITE(*,'(a,2i5,e12.4)') &
            'sourcpwr: k,lr, sorpw_nbi', &
                       k,lr_,sorpw_nbi(k,lr_)
      if(setup0%verbose>0) WRITE(*,*)'----------------------------------------------------'
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

      return
      end subroutine sourcpwr


end module sourcpwr_mod
