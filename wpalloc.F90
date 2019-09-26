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

module wpalloc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast

  !---END USE

!
!

contains

  subroutine wpalloc
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     Allocates arrays used by parallel transport routines.
!..............................................................

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif
!.......................................................................

!..................................................................
!     A check on allocations is sucessful entering then exiting
!     the subroutine.
!..................................................................
#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      WRITE(*,*)'wpalloc:  Entering wpalloc'
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

      ! YuP-101220: allocation of wcqlb-wcqlf is moved to vlf.f

      lnyxgs2=(iyp1+1)*(jxp1+1)*ngen*(setup0%ls+2)
      lny2gx=iy*jx*ngen*4
      lnsbn2y=(lsa+nbanda+2)*max(iy,jx)*2
      lnys2bn=max(iy,jx)*(lsa+2)*nbanda*2

      lndums=5*lnyxgs2+lny2gx+lnsbn2y+lnys2bn

      if (vlfmod.eq."enabled") then
         lnyxms=iy*jx*nmodsa*setup0%ls
         lndums=lndums+4*lnyxms
      endif

      allocate(l_upper(1:iy),STAT=istat)
      allocate(ilpm1ef(0:iy+1,0:lsa1,-1:+1),STAT=istat)
      call ibcast(l_upper,0,SIZE(l_upper))
      call ibcast(ilpm1ef,0,SIZE(ilpm1ef))

      allocate(fnhalf(0:iy+1,0:jx+1,ngen,0:setup0%ls+1),STAT=istat)
      call bcast(fnhalf,zero,SIZE(fnhalf))
      allocate(fnp0(0:iy+1,0:jx+1,ngen,0:setup0%ls+1),STAT=istat)
      call bcast(fnp0,zero,SIZE(fnp0))
      allocate(fnp1(0:iy+1,0:jx+1,ngen,0:setup0%ls+1),STAT=istat)
      call bcast(fnp1,zero,SIZE(fnp1))
      allocate(dls(0:iy+1,0:jx+1,ngen,0:ls+1),STAT=istat)
      call bcast(dls,zero,SIZE(dls))
      allocate(fg(0:iy+1,0:jx+1,ngen,0:ls+1),STAT=istat)
      call bcast(fg,zero,SIZE(fg))
      allocate(fh(0:iy+1,0:jx+1,ngen,0:ls+1),STAT=istat)
      call bcast(fh,zero,SIZE(fh))

      allocate(fedge(iy,jx,ngen,4),STAT=istat)
      call bcast(fedge,zero,SIZE(fedge))

      allocate(rhspar(0:lsa+nbanda+1,jx,2),STAT=istat)
      call bcast(rhspar,zero,SIZE(rhspar))
      allocate(bndmats(jx,0:lsa+1,nbanda,2),STAT=istat)
      call bcast(bndmats,zero,SIZE(bndmats))


#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      WRITE(*,*)'wpalloc:  Leaving wpalloc'
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

      return
      end subroutine wpalloc

end module wpalloc_mod
