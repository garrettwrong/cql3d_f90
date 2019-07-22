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

module rdc_bplt_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use pltdf_mod, only : pltcont

  !---END USE


contains

      subroutine rdc_bplt(krf)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     plots rf rdcb coefficient as a contour plot on cql3d grid.
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

      character*8 pltvlhb
      character*8 pltovlp

      data pltvlhb /'enabled'/
      data pltovlp /'enabled'/


#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      if (setup0%noplots.eq."enabled1") return
      if (pltvlhb.ne."enabled") return
!$$$      if (pltovlp.eq."enabled".and. mrfn.gt.1) then
!cc      call bcast(temp1(1,0),zero,iy*(jx+1)) ! YuP-101215: error?
      call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)

      write(*,*)'rdc_bplt(krf): mrfn =',mrfn,' krf=',krf

!$$   do 560 k=1,mrfn
      do 561 j=1,jx
         do 562 i=1,iy
            temp1(i,j)=rdcb(i,j,lr_,krf)
 562     continue
 561  continue
#ifndef NOPGPLOT
      CALL PGPAGE
#endif
      itype=7 ! means: plots are made for rdcb
      call pltcont(nrdcspecies(krf),1, &
           'Contours of RdcB vs. v_parallel,v_perp',7)
      write(t_,552) lr_
 552  format(" Flux surface number",i3,";   all modes, krf=",i2)
#ifndef NOPGPLOT
      CALL PGMTXT('B',10.,0.,0.,t_)
#endif

!$$$  560    continue

!$$$  endif


      return
      end subroutine rdc_bplt


end module rdc_bplt_mod
