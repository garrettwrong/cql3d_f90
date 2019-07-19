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

module pltdnz_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use pltelec_mod, only : pltelec

  !---END USE

!
!

contains

  subroutine pltdnz
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!     plot the density as a function of poloidal angle for a given
!     set of energy ranges.
!

      REAL RTAM1(LZA),RTAM2(LZA)
      REAL RPGMIN,RPGMAX
      REAL RILIN,PGCOORD

#ifndef NOPGPLOT
      CALL PGSAVE
#endif

!
      if (setup0%noplots.eq."enabled1") return
      if (pltdn .eq. "disabled") return
      do 100 k=1,ngen
        fu=.99999
        do 3006 ny=1,negyrg
          call aminmx(densz(1:lz,k,ny,lr_),1,lz,1,fmin,fmax,kmin,kmax)
!cc          write(*,*) 'pltdnz: fmin,fmax',fmin,fmax
          fmax=fmax+em90
          tam1(ny)=fmax
          tam2(ny)=fmin
          if (fmin/fmax .lt. fu) fu=fmin/fmax
 3006   continue

#ifndef NOPGPLOT
        CALL PGPAGE
#endif
#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.45,.95)
#endif

        DO L=1,LZ
           RTAM1(L)=pol(l,lr_)
        ENDDO
        RPGMIN=fu
        RPGMAX=1.


#ifndef NOPGPLOT
        CALL PGSWIN(RTAM1(1),RTAM1(LZ),RPGMIN,RPGMAX)
#endif
#ifndef NOPGPLOT
        CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
#endif


        xu=negyrg
        do 3002 ny=1,negyrg
          if (jegy(ny,1,k,lr_) .eq. 0 .or. eegy(ny,2,k,lr_) .lt. 1.e-15) &
            go to 3002
          xv=ny-1
          do 3001 l=1,lz
            tz1(l)=densz(l,k,ny,lr_)/tam1(ny)
 3001     continue

          DO L=1,LZ
             RTAM2(L)=tz1(l)
          ENDDO
#ifndef NOPGPLOT
          CALL PGSLS(MOD(NY,5))
#endif
#ifndef NOPGPLOT
          CALL PGLINE(LZ,RTAM1,RTAM2)
#endif

 3002   continue
#ifndef NOPGPLOT
        CALL PGLAB('Poloidal angle (radians)','Normalized density', &
             'Density as a function of poloidal angle(=pi*z/zmax) ')
#endif

        PGCOORD=-.15
!        RILIN=5.
#ifndef NOPGPLOT
!        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
#endif
        write(T_,611)
 611    FORMAT("(curves are normalized to a maximum of 1.)")
        RILIN=5.
#ifndef NOPGPLOT
        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
#endif
        write(T_,6111)
 6111   FORMAT( &
        "(Line order: full,dashed,dot-dash,dotted,dash-dot-dot-dot)")
        RILIN=RILIN+1.
#ifndef NOPGPLOT
        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
#endif
        write(t_,612) k,rovera(lr_),n,timet
 612    FORMAT( &
        "species",i2,"  r/a=",1pe8.2,"  n=",i5,"  time= ",1pe11.4)
        RILIN=RILIN+2.
#ifndef NOPGPLOT
        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
#endif
!        write(t_,613) xlndnz(k,negyrg)
! 613    FORMAT("line density (line-integration) =",1pe16.5)
!        RILIN=RILIN+1.
#ifndef NOPGPLOT
!        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
#endif


        RILIN=RILIN+1.
        do 3005 ny=1,negyrg
           RILIN=RILIN+1.
           write(t_,3003) eegy(ny,1,k,lr_),eegy(ny,2,k,lr_)
#ifndef NOPGPLOT
           CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
#endif
           RILIN=RILIN+1.
           write(t_,3004) tam1(ny)
#ifndef NOPGPLOT
           CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
#endif
 3005      continue

 3003   format("lower egy =",1pe10.2," kev;  upper egy =",1pe10.2,"kev")
 3004   format("maximum density on this curve  =",1pe16.5,"/cm**3")

        if (k.eq.ngen.and.kelecg.eq.0.and.locquas.eq."enabled$") &
          call pltelec

 100  continue

#ifndef NOPGPLOT
      CALL PGUNSA
#endif

      return
      end subroutine pltdnz


end module pltdnz_mod
