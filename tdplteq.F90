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

module tdplteq_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use frplteq_mod, only : micfrplt

  !---END USE

!
!

contains

  subroutine tdplteq(krf)
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      use frplteq_mod, only : micfrplt, textt
      implicit integer (i-n), real(c_double) (a-h,o-z)
#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

      character*8 pltrays
      character*8 plteq

      REAL rbot,zbot,rtop,ztop
      REAL RTAB1(LFIELDA),RTAB2(LFIELDA)
      REAL PGER1,PGERNNR,PGZERO,PGZMIN,PGZMAX,PGEZNNZ
      REAL wk_ray_r(nrayelts), wk_ray_z(nrayelts)

!..................................................................
!     This routine plots out the contours (flux surfaces) on
!     which the calculations are done.
!..................................................................

#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      pltrays="enabled" ! for plotting rays over cross-sectional view.
                        ! Can be moved to namelist later.

      if (setup0%noplots.eq."enabled1") return
      plteq="enabled"
      if (plteq.eq."disabled") return

!     Set up textt to contain number 1 to 200:
      call micfrplt
!     Would be nice to label the flux surfaces, but
!     it isn't done at the moment....

      if (eqsym.ne."none") then ! Up-dn symm
         ztop=.75*ez(nnz)/(er(nnr)-er(1))+.05
      else ! Not up-dn symm
         ztop=0.95
      endif

      ! YuP[03-2016][07-2017] Added plotting rays in cross-sectional view
      if (urfmod.ne."disabled" .and. pltrays.eq.'enabled') then
         ! over-write setting for page size - plot whole cross-section.
         ! For up-dn symmetrical case only half of surfaces are plotted
         ! but rays could be in the other hemisphere.
         ztop=.95
      endif


      ! YuP: min and max value of Z-coord over all flux surfaces:
      solz_min=MINVAL(solz)
      solz_max=MAXVAL(solz)
      PGZMIN=real(solz_min)
      PGZMAX=real(solz_max)

#ifndef NOPGPLOT
      CALL PGPAGE
#endif
#ifndef NOPGPLOT
      CALL PGSVP(.15,.85,.15,ztop)
#endif

      PGER1=er(1)
      PGERNNR=er(nnr)
      PGZERO=0.
      PGEZNNZ=ez(nnz)
      if( (eqsym.ne."none") .and. &
          (urfmod.eq."disabled" .or. pltrays.eq.'disabled') ) then
         !plot half only
#ifndef NOPGPLOT
!-YuP         CALL PGSWIN(PGER1,PGERNNR,PGZERO,PGEZNNZ)
#endif
#ifndef NOPGPLOT
!-YuP         CALL PGWNAD(PGER1,PGERNNR,PGZERO,PGEZNNZ)
#endif
#ifndef NOPGPLOT
         CALL PGSWIN(PGER1,PGERNNR,PGZMIN,PGZMAX)
#endif
#ifndef NOPGPLOT
         CALL PGWNAD(PGER1,PGERNNR,PGZMIN,PGZMAX)
#endif
      else !eqsym=none; and/or urfmod='enabled',
         ! plot upper and lower halves:
#ifndef NOPGPLOT
         CALL PGSWIN(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
#endif
#ifndef NOPGPLOT
         CALL PGWNAD(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
#endif
      endif

#ifndef NOPGPLOT
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
#endif
      if ( (urfmod.ne."disabled") .and. (pltrays.eq.'enabled') &
                                  .and. (krf.gt.0)            ) then
#ifndef NOPGPLOT
         CALL PGLAB('Major radius (cms)','Vert height (cms)', &
              'Fokker-Planck Flux Surfaces + Rays')
#endif
      else ! urfmod='disabled', or krf=0
#ifndef NOPGPLOT
         CALL PGLAB('Major radius (cms)','Vert height (cms)', &
              'Fokker-Planck Flux Surfaces')
#endif
      endif

      IF (setup0%lrzmax.GT.200) STOP 'TDPLTEQ: CHECK DIM OF TEXTT'

      do 10 l=1,setup0%lrzmax
         IF (LORBIT(L).GT.LFIELDA) STOP 'TDPLTEQ: CHECK DIM OF RTAB1/2'
        do 20 j=1,lorbit(l)
           RTAB1(j)=solr(lorbit(l)+1-j,l)
           RTAB2(j)=solz(lorbit(l)+1-j,l)
 20     continue
        text(1)=textt(l)
#ifndef NOPGPLOT
        CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
#endif
        ! YuP[03-2016] Added plotting rays in cross-sectional view
        if ( (urfmod.ne."disabled") .and. (pltrays.eq.'enabled') &
             .and. (eqsym.ne.'none')) then
         ! Add surfaces in whole cross-section.
         ! For up-dn symmetrical case only half of surfaces are plotted
         ! but rays could be in the other hemisphere,
         ! so plot the other half:
#ifndef NOPGPLOT
         CALL PGLINE(LORBIT(L),RTAB1,-RTAB2)
#endif
        endif
 10   continue

      if((eqsym.eq.'none').or. &
       ((urfmod.ne."disabled") .and. (pltrays.eq.'enabled')) ) then

        if(ncontr.gt.1) then
          ! YuP[2015/05/03] Add LCFS, if available
          ncontr_= min(ncontr,LFIELDA)
          do ilim=1,ncontr_
             RTAB1(ilim)=rcontr(ilim)
             RTAB2(ilim)=zcontr(ilim)
          enddo
#ifndef NOPGPLOT
          CALL PGLINE(ncontr_,RTAB1,RTAB2)
#endif
        endif

        if(nlimiter.gt.1) then
          ! YuP[2019-02-22] Add "last surface" (plasma border), if available
          nline= min(nlimiter,LFIELDA)
          do ilim=1,nline
             RTAB1(ilim)=rlimiter(ilim)
             RTAB2(ilim)=zlimiter(ilim)
          enddo
#ifndef NOPGPLOT
          CALL PGSLW(setup0%lnwidth*2) ! bold
#endif
#ifndef NOPGPLOT
          CALL PGLINE(nline,RTAB1,RTAB2)
#endif
          if(machine.eq."mirror") then
#ifndef NOPGPLOT
          CALL PGLINE(nline,-RTAB1,RTAB2) !mirror area to the left of Z-axis
#endif
          endif
#ifndef NOPGPLOT
          CALL PGSLW(setup0%lnwidth) ! restore
#endif
        endif

      endif !eqsym.eq.'none' or (urfmod.ne."disabled")&(pltrays.eq.'enabled')



!..................................................................
! YuP[03-2016] Added plotting rays in cross-sectional view
      if ((urfmod.ne."disabled") .and. (pltrays.eq.'enabled') .and. &
          (krf.gt.0)  ) then
        !do krf=1,mrfn ! YuP[2019-02-22] separate page for each krf mode now
        do iray=1,nray(krf)  !Loop over rays
           nrayelt00=nrayelt(iray,krf)
!           write(*,*)'tdplteq: krf, iray, nrayelt00=',
!     +                         krf, iray, nrayelt00
!           write(*,'(a,3i6)')
!     +      'tdplteq: iray,lloc(nrayelt00),llray(nrayelt00)=',
!     +      iray,lloc(nrayelt00,iray,krf),llray(nrayelt00,iray,krf) !local to ray element point.
           wk_ray_r=0.0 ! reset
           wk_ray_z=0.0 ! reset
           do is=1,nrayelt00
             wk_ray_r(is)=wr(is,iray,krf)
             wk_ray_z(is)=wz(is,iray,krf)
           enddo
!           write(*,*)'tdplteq: krf, minval(wk_ray_r),maxval(wk_ray_r)=',
!     +                         krf, minval(wk_ray_r),maxval(wk_ray_r)
#ifndef NOPGPLOT
           CALL PGLINE(nrayelt00,wk_ray_r,wk_ray_z)
#endif
        enddo  ! iray
        !Add some info on rays:
        write(t_,'(a,i3, a,2i3, a,e12.5)') 'krf=',krf, &
                              '  nharm,nharms=', nharm(krf),nharms(krf), &
                              '    f[Hz]=', freqcy(krf)
#ifndef NOPGPLOT
        CALL PGMTXT('T',0.9,0.0,0.,t_) ! 0.9=just outside of viewport
#endif
        !enddo
      endif
!..................................................................

      return
      end subroutine tdplteq


end module tdplteq_mod
