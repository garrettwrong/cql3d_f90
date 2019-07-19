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

module vlhbplt_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use pltdf_mod, only : pltcont

  !---END USE


contains

  subroutine vlhbplt
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     plots rf cqlb coefficient as a contour plot.
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

      character*8 pltvlhb
      character*8 pltovlp

      save pltvlhb,pltovlp
      data pltvlhb /'enabled'/
      data pltovlp /'enabled'/


#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      if (setup0%noplots.eq."enabled1") return
      if (pltvlhb.ne."enabled") return
      if (pltovlp.eq."enabled".and. mrfn.gt.1) then
!cc        call bcast(temp1(1,0),zero,iy*(jx+1))  ! YuP-101215: error?
        call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)


        do 560 k=1,mrfn
          do 561 j=1,jx
            do 562 i=1,iy
              temp1(i,j)=cqlb(i,j,indxlr_,k)
 562        continue
 561      continue
#ifndef NOPGPLOT
          CALL PGPAGE
#endif
          itype=6 ! means: plots are made for vlhb
          call pltcont(1,1,'Contours of CqlB vs. v_parallel,v_perp', &
                       itype)
!$$$          call gstxno(80.)
!$$$          call gscpvs(.15,.35)
          write(t_,552) lr_
 552      format(" Flux surface number",i3,";   all modes")
#ifndef NOPGPLOT
          CALL PGMTXT('B',10.,0.,0.,t_)
#endif

 560    continue

      endif

      do 680 k=1,mrfn
!..................................................................
!     Compute vpar21/vte and vpar11/vte normalized velocities
!     Following vlh.f
!..................................................................
        if (vprprop .eq. "enabled") then
          xvpr=1.
        else
          xvpr=0.
        endif

!
!     determine the vparallel interval at R=R0
!
        vpmax=vparmax(k)*clight
        vpmin=vparmin(k)*clight
!
!     determine the vparallel range of nonzero D at outside
!     of flux surface.
!
        rovr0=(radmaj+xvpr*radmin*rovera(lr_))/radmaj
        vmin=vpmin*rovr0
        vmax=vpmax*rovr0

        vpar21dv=vmin/vth(1,lr_)
        vpar11dv=vmax/vth(1,lr_)


        do  j=1,jx
           do  i=1,iy
              temp1(i,j)=cqlb(i,j,indxlr_,k)
           enddo
        enddo

#ifndef NOPGPLOT
        CALL PGPAGE
#endif
        itype=6 ! means: plots are made for vlhb
        call pltcont(1,1,'Contours of CqlB vs. v_parallel,v_perp',itype)
!$$$        call gstxno(80.)
!$$$        call gscpvs(.15,.35)
        write(t_,660) lr_,k
#ifndef NOPGPLOT
        CALL PGMTXT('B',10.,0.,0.,t_)
#endif
        write(t_,661) vpar21dv,vpar11dv
#ifndef NOPGPLOT
        CALL PGMTXT('B',11.,0.,0.,t_)
#endif

 680  continue
 660  format("Flux surface number",i3," mode=",i1)
 661  format("vpar21/vth=",1pe15.7,"   vpar11/vth=",1pe15.7)

      return
      end subroutine vlhbplt


end module vlhbplt_mod
