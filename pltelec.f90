module pltelec_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx

  !---END USE

!
!
!

contains

      subroutine pltelec
      use param_mod
      use comm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!mnt  this routine plots electron density as a function of poloidal angl
!

      REAL RTAM1(LZA),RTAM2(LZA)
      REAL RPGMIN,RPGMAX
      REAL RILIN

      if (noplots.eq."enabled1") return
!$$$      call gxglfr(0)
      call aminmx(densz(1:lz,ngen+1,negyrg,lr_), &
       1,lz,1,fmin,fmax,kmin,kmax)
      if (fmin .eq. fmax) fmin=.9*fmax-1.e-20
!$$$      call gswd2d("linlin$",pol(1,lr_),pol(lz,lr_),fmin,fmax)
!$$$      call gsvp2d(.2,.8,.25,.95)
!$$$      call gpgr80("linlin$")

      CALL PGPAGE
      CALL PGSVP(.2,.8,.45,.95)

      DO L=1,LZ
         RTAM1(L)=pol(L,lr_)
      ENDDO

      RPGMIN=fmin
      RPGMAX=fmax

      IF ( RPGMAX-RPGMIN .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPGMAX= RPGMIN+1.e-16
      ENDIF
      CALL PGSWIN(RTAM1(1),RTAM1(LZ),RPGMIN,RPGMAX)

!$$$      do 3001 l=1,lz
!$$$        tz1(l)=densz(l,ngen+1,negyrg,lr_)
!$$$ 3001 continue
!$$$      call gpcv2d(pol(1,lr_),tz1,lz)

      do 3001 l=1,lz
         RTAM2(L)=densz(l,ngen+1,negyrg,lr_)
 3001 continue
      CALL PGLINE(LZ,RTAM1,RTAM2)

 3002 continue
!$$$      call gscvlb(0)
!$$$      call gstxno(100.)
!$$$      call gscpvs(.5,.2)

!$$$      write(t_,610) kelec,n,timet,xlndnz(ngen+1,negyrg)
!$$$      call gptx2d(t_)
      RILIN=1.
      CALL PGMTXT(B,RILIN,0.,0.,T_)
      write(t_,611) kelec,n,timet
      RILIN=RILIN+1.
      CALL PGMTXT(B,RILIN,0.,0.,T_)
      write(t_,612) xlndnz(ngen+1,negyrg)
      RILIN=RILIN+1.
      CALL PGMTXT(B,RILIN,0.,0.,T_)
 610  format("Density as a function of poloidal angle(=pi*z/zmax)")
 611  format("species ",i3, " (electrons)   n= ",i5,"  time= ",1pe14.4)
 612  format("density (line-integration) =",1pe16.5)
      return
      end
end module pltelec_mod
