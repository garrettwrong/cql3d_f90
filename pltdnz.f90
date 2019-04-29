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
      use param_mod
      use comm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!     plot the density as a function of poloidal angle for a given
!     set of energy ranges.
!

      REAL RTAM1(LZA),RTAM2(LZA)
      REAL RPGMIN,RPGMAX
      REAL RILIN,PGCOORD

      CALL PGSAVE

!
      if (noplots.eq."enabled1") return
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

        CALL PGPAGE
        CALL PGSVP(.2,.8,.45,.95)

        DO L=1,LZ
           RTAM1(L)=pol(l,lr_)
        ENDDO
        RPGMIN=fu
        RPGMAX=1.


        CALL PGSWIN(RTAM1(1),RTAM1(LZ),RPGMIN,RPGMAX)
        CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)


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
          CALL PGSLS(MOD(NY,5))
          CALL PGLINE(LZ,RTAM1,RTAM2)

 3002   continue
        CALL PGLAB('Poloidal angle (radians)','Normalized density', &
             'Density as a function of poloidal angle(=pi*z/zmax) ')

        PGCOORD=-.15
!        RILIN=5.
!        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
        write(T_,611)
 611    FORMAT("(curves are normalized to a maximum of 1.)")
        RILIN=5.
        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
        write(T_,6111)
 6111   FORMAT( &
        "(Line order: full,dashed,dot-dash,dotted,dash-dot-dot-dot)")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
        write(t_,612) k,rovera(lr_),n,timet
 612    FORMAT( &
        "species",i2,"  r/a=",1pe8.2,"  n=",i5,"  time= ",1pe11.4)
        RILIN=RILIN+2.
        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
!        write(t_,613) xlndnz(k,negyrg)
! 613    FORMAT("line density (line-integration) =",1pe16.5)
!        RILIN=RILIN+1.
!        CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)


        RILIN=RILIN+1.
        do 3005 ny=1,negyrg
           RILIN=RILIN+1.
           write(t_,3003) eegy(ny,1,k,lr_),eegy(ny,2,k,lr_)
           CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
           RILIN=RILIN+1.
           write(t_,3004) tam1(ny)
           CALL PGMTXT('B',RILIN,PGCOORD,0.,T_)
 3005      continue

 3003   format("lower egy =",1pe10.2," kev;  upper egy =",1pe10.2,"kev")
 3004   format("maximum density on this curve  =",1pe16.5,"/cm**3")

        if (k.eq.ngen.and.kelecg.eq.0.and.locquas.eq."enabled$") &
          call pltelec

 100  continue

      CALL PGUNSA

      return
      end
end module pltdnz_mod
