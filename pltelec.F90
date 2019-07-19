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

module pltelec_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx

  !---END USE

!
!
!

contains

  subroutine pltelec
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!mnt  this routine plots electron density as a function of poloidal angl
!

      REAL RTAM1(LZA),RTAM2(LZA)
      REAL RPGMIN,RPGMAX
      REAL RILIN

      if (setup0%noplots.eq."enabled1") return
!$$$      call gxglfr(0)
      call aminmx(densz(1:lz,ngen+1,negyrg,lr_), &
       1,lz,1,fmin,fmax,kmin,kmax)
      if (fmin .eq. fmax) fmin=.9*fmax-1.e-20
!$$$      call gswd2d("linlin$",pol(1,lr_),pol(lz,lr_),fmin,fmax)
!$$$      call gsvp2d(.2,.8,.25,.95)
!$$$      call gpgr80("linlin$")

#ifndef NOPGPLOT
      CALL PGPAGE
#endif
#ifndef NOPGPLOT
      CALL PGSVP(.2,.8,.45,.95)
#endif

      DO L=1,LZ
         RTAM1(L)=pol(L,lr_)
      ENDDO

      RPGMIN=fmin
      RPGMAX=fmax

      IF ( RPGMAX-RPGMIN .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPGMAX= RPGMIN+1.e-16
      ENDIF
#ifndef NOPGPLOT
      CALL PGSWIN(RTAM1(1),RTAM1(LZ),RPGMIN,RPGMAX)
#endif

!$$$      do 3001 l=1,lz
!$$$        tz1(l)=densz(l,ngen+1,negyrg,lr_)
!$$$ 3001 continue
!$$$      call gpcv2d(pol(1,lr_),tz1,lz)

      do 3001 l=1,lz
         RTAM2(L)=densz(l,ngen+1,negyrg,lr_)
 3001 continue
#ifndef NOPGPLOT
      CALL PGLINE(LZ,RTAM1,RTAM2)
#endif

 3002 continue
!$$$      call gscvlb(0)
!$$$      call gstxno(100.)
!$$$      call gscpvs(.5,.2)

!$$$      write(t_,610) kelec,n,timet,xlndnz(ngen+1,negyrg)
!$$$      call gptx2d(t_)
      RILIN=1.
#ifndef NOPGPLOT
      CALL PGMTXT(B,RILIN,0.,0.,T_)
#endif
      write(t_,611) kelec,n,timet
      RILIN=RILIN+1.
#ifndef NOPGPLOT
      CALL PGMTXT(B,RILIN,0.,0.,T_)
#endif
      write(t_,612) xlndnz(ngen+1,negyrg)
      RILIN=RILIN+1.
#ifndef NOPGPLOT
      CALL PGMTXT(B,RILIN,0.,0.,T_)
#endif
 610  format("Density as a function of poloidal angle(=pi*z/zmax)")
 611  format("species ",i3, " (electrons)   n= ",i5,"  time= ",1pe14.4)
 612  format("density (line-integration) =",1pe16.5)
      return
      end subroutine pltelec


end module pltelec_mod
