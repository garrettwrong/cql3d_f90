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

module urfbplt_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use pltdf_mod, only : pltcont

  !---END USE


contains

  subroutine urfbplt
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     plots rf urfb coefficient as a contour plot.
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

      character*8 pltovlp

      data pltovlp /'enabled'/


#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      if (setup0%noplots.eq."enabled1") return
!BH081105      iiplt3d=0
      iiplt3d=-1
      do i=1,nplota
         if (n.eq.nplt3d(i)) iiplt3d=n
      enddo

!      if (n/nplt3d*nplt3d.ne.n .and. n.ne.1) return
!BH081105      if (iiplt3d.eq.0 .and. n.ne.1) return
      if (iiplt3d.eq.-1) return
      if (mplot(l_).eq."disabled") return
      if (urfmod.eq."disabled") return
      if (plturfb.eq."disabled") return

      if (pltovlp.eq."enabled".and. mrfn.gt.1) then

!       This do 560 loop was overlapping the modes on one plot.
        !YuP[2016] instead of overlapping contour plots,
        !sum them and then plot contours for the total
        !UrfB(all modes) at a given surface (and for each k-species).
        do k=1,ngen ! YuP[10-2016] scan general species: can be more than one
           !initialize for each species: one plot for each k
           call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)
           do 560 krf=1,mrfn
             if (nrfspecies(krfn(krf)) .eq.k) then
               !sum-up modes for a given species only
               do j=1,jx
               do i=1,iy
                 !temp1(i,j)=urfb(i,j,indxlr_,krf) ! YuP: original
                 temp1(i,j)= temp1(i,j)+urfb(i,j,indxlr_,krf) ! YuP: new version: sum-up
               enddo
               enddo
               !call pltcont(1,1,'Contours of UrfB vs. v_parallel,v_perp')  !YuP:original
             endif
 560       continue ! krf mode (usually = harmonic number)
           if( MAXVAL(temp1)-MINVAL(temp1) .gt. 0.d0 ) then
#ifndef NOPGPLOT
             CALL PGPAGE ! new page for each k
#endif
               itype=4 ! means: plots are made for urfb
             call pltcont(k,1,'Contours of UrfB vs. v_parallel,v_perp', &
               itype) !YuP:summed-up
             write(t_,552)
#ifndef NOPGPLOT
             CALL PGMTXT('B',10.,0.,0.,t_)
#endif
             write(t_,553) lr_
#ifndef NOPGPLOT
             CALL PGMTXT('B',11.,0.,0.,t_)
#endif
             write(t_,692) MAXVAL(temp1) !YuP[10-2016] max value for this krf
#ifndef NOPGPLOT
             CALL PGMTXT('B',12.,0.,0.,t_)
#endif
             write(t_,693) k
#ifndef NOPGPLOT
             CALL PGMTXT('B',13.,0.,0.,t_)
#endif
           endif
        enddo ! k species

      endif

!     This do 680 loop plots the individual mode contributions:
      do k=1,ngen ! YuP[10-2016] scan general species
      do 680 krf=1,mrfn
         do j=1,jx
            do i=1,iy
               temp1(i,j)=urfb(i,j,indxlr_,krf)
            enddo
         enddo
         if (nrfspecies(krfn(krf)) .eq. k) then
         if( MAXVAL(temp1)-MINVAL(temp1) .gt. 0.d0 ) then
#ifndef NOPGPLOT
          CALL PGPAGE ! opens new page for each krf-mode
#endif
          itype=4 ! means: plots are made for urfb
          call pltcont(k,1,'Contours of UrfB vs. v_parallel,v_perp', &
            itype)
          write(t_,690)
#ifndef NOPGPLOT
          CALL PGMTXT('B',10.,0.,0.,t_)
#endif
          ! write flux surface number and mode number;
          ! also harmonic number and species number (added YuP[10-2016])
          write(t_,691) lr_ ,krf,nharm(krf),k
#ifndef NOPGPLOT
          CALL PGMTXT('B',11.,0.,0.,t_)
#endif
          write(t_,692) MAXVAL(temp1) !YuP[10-2016] max value for this krf
#ifndef NOPGPLOT
          CALL PGMTXT('B',12.,0.,0.,t_)
#endif
         endif
         endif
 680  continue ! krf
      enddo ! k species

 552  format("Contours of the rf (v,v) diffusion coefficient, urfb")
 553  format(" Flux surface number",1x,i3,"; all modes")
 690  format("Contours of the rf (v,v) diffusion coefficient, urfb")
 691  format("Flux surf.N",i3,";  mode,nharm=",2i5,";  Species k=",i1)
 692  format("Max value for this surface/mode:",e13.3)
 693  format("Species k=",i1)
      return
      end subroutine urfbplt


end module urfbplt_mod
