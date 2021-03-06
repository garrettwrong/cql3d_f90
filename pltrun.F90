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

module pltrun_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use tdnflxs_mod, only : tdnflxs

  !---END USE

!
!

contains

  subroutine pltrun
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : rbound
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      dimension yg(nonch),xg(nonch)

!     PGPLOT REAL Variables:
      real(c_float) RPG1,RPG2
      real(c_float) RNONCHA1(nonch),RNONCHA2(nonch)
      !XXX real(c_float) RBOUND
!-----------------
!   This subroutine writes and plots runaway population and current
!------------------

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      if (pltra.eq."disabled") return

      iounit=35
      open(unit=iounit,file='runaway.out',status='unknown')

!BH070405      if (nch(l_).gt.noncha .or. nch(l_).gt.500)
      if (nch(1).gt.nonch) stop 'check dimensions in pltrun'

      write (iounit,20001)
20001 format(3x,'time',12x,'runaway x'//)
      do ll=1,setup0%lrz
        write(iounit,20010) setup0%lrindx(ll)
20010   format(//2x,'flux surface ',i3)
        write(iounit,20015)
20015   format(//19x,'denra:'/)
!BH070405        do ntm=1,nch(l_)
        do ntm=1,nch(l_)
          write(iounit,20020) ptime(ntm,ll),pdenra(ntm,ll)
20020     format(2x,1pe12.4,2x,1pe18.10)
        enddo
        write(iounit,20016)
20016   format(//19x,'curra:'/)
!BH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),pcurra(ntm,ll)
        enddo
        write(iounit,20017)
20017   format(//19x,'frac. ra-dens.:'/)
!BH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),pfdenra(ntm,ll)
        enddo
        write(iounit,20018)
20018   format(//19x,'frac. ra-cur.:'/)
!BH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),pfcurra(ntm,ll)
        enddo
        write(iounit,20013)
20013   format(//19x,'ucrit/c.:'/)
!BH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),pucrit(ntm,ll)/cnorm
        enddo
        write(iounit,20014)
20014   format(//19x,'e/ecrit0.:'/)
!BH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),peoe0(ntm,ll)
        enddo
      enddo

      write(iounit,20050)
20050 format(//2x,'Pitch Angle Averaged Distribution At Time Slices:'//)

      nframep=min(iplot,5)
      do ns=1,nframep
        write(iounit,20060) tplot(ns)
20060   format(///2x,'Time=', f9.4,'secs.'/)
        do ll=1,setup0%lrz
          write(iounit,20070) ll
20070     format(2x,i3,'th flux surface'/, &
              5x,'GAMMA',15x,'F(GAMMA)'/)
          do j=1,jx
            write(iounit,20080) gamma(j),f_aveth(j,1,ll,ns)
20080       format(1x,e16.8,2x,e20.10)
          enddo
        enddo
      enddo

      close(unit=iounit)

!-----------------------------------------------------------------------
!     Plotting
!-----------------------------------------------------------------------

      if (setup0%noplots.eq."enabled1") return

      dgts=1.e-8

      do ll=1,setup0%lrz
        call tdnflxs(ll)
        rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
!BH070405        do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pdenra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))

        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
!        if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
!        call pltcycl(iymin,iymax,ymin,ymax)

        ymax=1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)



!      Following plot is with pdena,... with 1st dimension up to nch(l_).

#ifndef NOPGPLOT
        CALL PGPAGE
        CALL PGSVP(.2,.8,.6,.9)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
#endif
        IF (RPG1.ne.RPG2) THEN
#ifndef NOPGPLOT
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
#endif
        ENDIF
#ifndef NOPGPLOT
        CALL PGLAB('time (secs)','RA density (cm\u-3\d)', &
             'Runaway Density and Current vs. Time')
#endif

!BH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pcurra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))

        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)

!        if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)

        ymax=1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.2,.5)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
#endif
        IF (RPG1.ne.RPG2) THEN
#ifndef NOPGPLOT
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
#endif
        ENDIF
#ifndef NOPGPLOT
        CALL PGLAB('time (secs)','RA curr den (A/cm\u2\d)',' ')
#endif


        write(t_,10010) n,timet
#ifndef NOPGPLOT
        CALL PGMTXT('B',6.,-.1,0.,t_)
#endif
        write(t_,10011) rovera(lr_),ll,rr
#ifndef NOPGPLOT
        CALL PGMTXT('B',7.,-.1,0.,t_)
#endif

!-----------------------------------------------------------------------


!BH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pfdenra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))

        enddo
!        write(*,*) 'ptime(ijk,1),pfdenra(ijk,ll),ijk=1,5',
!     +              (ptime(ijk,1),pfdenra(ijk,ll),ijk=1,5)
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
!        if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
!        call pltcycl(iymin,iymax,ymin,ymax)

        ymax=1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


#ifndef NOPGPLOT
        CALL PGPAGE
        CALL PGSVP(.2,.8,.6,.9)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
#endif
        IF (RPG1.ne.RPG2) THEN
#ifndef NOPGPLOT
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
#endif
        ENDIF
#ifndef NOPGPLOT
        CALL PGLAB('time (secs)','Fraction RA density', &
             'Runaway Fraction of Density and Current vs. Time')
#endif

        do nt=1,nch(l_)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pfcurra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))

         enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
!        if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
!        call pltcycl(iymin,iymax,ymin,ymax)

        ymax=1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)



#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.2,.5)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
#endif
        IF (RPG1.ne.RPG2) THEN
#ifndef NOPGPLOT
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
#endif
#ifndef NOPGPLOT
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
#endif
        ENDIF
#ifndef NOPGPLOT
        CALL PGLAB('time (secs)','Fraction RA curr den',' ')
#endif

        write(t_,10010) n,timet
#ifndef NOPGPLOT
        CALL PGMTXT('B',6.,-.1,0.,t_)
#endif
        write(t_,10011) rovera(lr_),ll,rr
#ifndef NOPGPLOT
        CALL PGMTXT('B',7.,-.1,0.,t_)
#endif

!-----------------------------------------------------------------------




!     Plots of ucrit and e/e_dreicer, versus time:

!BH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
           xg(nt)=ptime(nt,1)
           yg(nt)=abs(pucrit(nt,ll))
           RNONCHA1(nt)=RBOUND(XG(nt))
           RNONCHA2(nt)=RBOUND(YG(nt))
           RNONCHA2(nt)=LOG10(RNONCHA2(nt))
        enddo


      if (pltlim.eq.'u/c') then
         do nt=1,nch(l_)
!YuP            yg(nt)=yg(nt)*cnormi !cnormi=0 when relativ.eq."disabled"
          yg(nt)=yg(nt)/cnorm ! YuP[07-2016]
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
         enddo
      elseif (pltlim.eq.'energy') then
!BH070405         do nt=1,nch(l_)
         do nt=1,nch(ll)
            yg2=yg(nt)*yg(nt)
            if(yg2*cnorm2i.lt.1.e-8 .or. relativ.eq."disabled") then
               g1=.5*yg2/cnorm2
            else
               g1=sqrt(1.+yg2*cnorm2i)-1.
            endif
            yg(nt)=g1*restmkev
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
         enddo
      endif

      call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
!     if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
!     call pltcycl(iymin,iymax,ymin,ymax)


      ymin=0.97*ymin
      ymax=1.e3*ymin

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


      if (pltlim.eq.'u/c') then
         write(t_,1020)
 1020    format("Critical runaway vel/c")
      elseif (pltlim.eq.'energy') then
         write(t_,1021)
 1021    format("Critical runaway energy (keV)")
      else
         write(t_,1022)
 1022    format("Critical runaway vel/vnorm")
      endif

#ifndef NOPGPLOT
        CALL PGPAGE
#endif
#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.6,.9)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
#endif
        IF (RPG1.ne.RPG2) THEN
#ifndef NOPGPLOT
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
#endif
#ifndef NOPGPLOT
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
#endif
        ENDIF
#ifndef NOPGPLOT
        CALL PGLAB('time (secs)',t_, &
             'Critical vel(energy) and E/E\dDreicer\u vs. Time')
#endif

!BH070405      do nt=1,nch(l_)
      do nt=1,nch(ll)
         xg(nt)=ptime(nt,1)
         yg(nt)=abs(peoed(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
      enddo
      call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
!     if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
!     call pltcycl(iymin,iymax,ymin,ymax)

      ymax=1.03*ymax
      ymin=ymax/1.e3

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)

      write(t_,1023)
 1023 format("E-field/E_Dreicer")
      write(t_,1024)
 1024 format("time(sec)")


#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.2,.5)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
#endif
        IF (RPG1.ne.RPG2) THEN
#ifndef NOPGPLOT
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
#endif
#ifndef NOPGPLOT
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
#endif
        ENDIF
#ifndef NOPGPLOT
        CALL PGLAB('time (secs)','E-field/E\dDreicer\u',' ')
#endif

        write(t_,10010) n,timet
#ifndef NOPGPLOT
        CALL PGMTXT('B',6.,-.1,0.,t_)
#endif
        write(t_,10011) rovera(lr_),ll,rr
#ifndef NOPGPLOT
        CALL PGMTXT('B',7.,-.1,0.,t_)
#endif

!-----------------------------------------------------------------------


!     Plots of e/e0 and KO source, versus time:

!BH070405      do nt=1,nch(l_)
      do nt=1,nch(ll)
         xg(nt)=ptime(nt,1)
         yg(nt)=abs(peoe0(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
      enddo

      call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
!     if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
!     call pltcycl(iymin,iymax,ymin,ymax)

!$$$      call gxglfr(0)

      ymax=1.03*ymax
      ymin=ymax/1.e3

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


#ifndef NOPGPLOT
        CALL PGPAGE
#endif
#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.6,.9)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
#endif
        IF (RPG1.ne.RPG2) THEN
#ifndef NOPGPLOT
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
#endif
#ifndef NOPGPLOT
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
#endif
        ENDIF
#ifndef NOPGPLOT
        CALL PGLAB('time (secs)','E-field/E0', &
             'E/(Critical E0) and KO Source vs. Time')
#endif


!BH070405      do nt=1,nch(l_)
      do nt=1,nch(ll)
         xg(nt)=ptime(nt,1)
         yg(nt)=abs(psrc(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
      enddo
      call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
!     if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
!     call pltcycl(iymin,iymax,ymin,ymax)

      ymax=1.03*ymax
!BH180714      ymin=ymax/1.e3
      ymin=ymax/1.e8

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.2,.5)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
#endif
        IF (RPG1.ne.RPG2) THEN
#ifndef NOPGPLOT
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
#endif
#ifndef NOPGPLOT
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
#endif
        ENDIF
#ifndef NOPGPLOT
        CALL PGLAB('time (secs)','KO Src(electrons/cm\u3\d/sec)',' ')
#endif

        write(t_,10010) n,timet
#ifndef NOPGPLOT
        CALL PGMTXT('B',6.,-.1,0.,t_)
#endif
        write(t_,10011) rovera(lr_),ll,rr
#ifndef NOPGPLOT
        CALL PGMTXT('B',7.,-.1,0.,t_)
#endif

!-----------------------------------------------------------------------


      enddo   ! on ll:1,setup0%lrz

10010 format("After time step(n)=",i4,5x,"time=",1pe14.6," secs")
10011 format("r/a=",f7.4,",  radial bin=",i3, &
             ",  radial position(R)=",1pe10.3," cm")


      return
      end subroutine pltrun
end module pltrun_mod
