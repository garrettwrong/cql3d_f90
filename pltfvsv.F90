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

module pltfvsv_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine pltfvsv
        use cqlconf_mod, only : setup0
        use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : luf
      use pltdf_mod, only : JXQ
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

      REAL RTAM1(jx),RTAM2(jx)
      REAL REMAX,REMIN
      CHARACTER*64 TX_
      real(c_double) wkd(jx)

!...
!mnt  This routine plots f as a function of velocity or energy
!mnt  for various angles.
!...
      if (setup0%noplots.eq."enabled1") return


      do 120 k=1,ngen ! ---------------------------------------------

         if (pltlim.eq."disabled") then
            jxq=jx
            xmaxq=x(jxq)
            do j=1,jxq
               tam1(j)=x(j)
            enddo
            TX_='u/vnorm'
         endif
         if (tandem.eq."enabled" .and. fmass(k).gt.1.e-27) then
            jxq=jlwr
            xmaxq=xlwr
            do j=1,jxq
               tam1(j)=x(j)
            enddo
            TX_='u/vnorm'
         ! If pltlim .ne. "disabled", plot versus
         ! 'x', 'u/c', or 'energy', up to maximum pltlimm.
         elseif (pltlim.eq.'x') then
            jxq=min(luf(pltlimm,x,jx),jx)
            xmaxq=x(jxq)
            do j=1,jxq
               tam1(j)=x(j)
            enddo
            TX_='u/vnorm'
         elseif (pltlim.eq.'u/c') then
            jxq=min(luf(pltlimm,uoc,jx),jx)
            xmaxq=uoc(jxq)
            do j=1,jxq
               tam1(j)=uoc(j)
            enddo
            TX_='u/c\dlight\u'
         elseif (pltlim.eq.'energy') then
            pltlimmm=pltlimm
            wkd(1:jx)=enerkev(1:jx,k)
            jxq=min(luf(pltlimmm,wkd,jx),jx)
            xmaxq=enerkev(jxq,k) !YuP[2018-01-08] added 2nd index (k)
            do j=1,jxq
               tam1(j)=enerkev(j,k) !YuP[2018-01-08] added 2nd index (k)
            enddo
            TX_='Energy (keV)'
         endif

         imsh(1)=1
         imsh(2)=itl
         imsh(3)=iyh
         imsh(4)=iy
         imsh(5)=0

         call bcast(tam4,zero,jx)
         call bcast(tam5,zero,jx)
!     tam4 contains pitch angle averaged equatorial distribution
!          function.
!     tam5 to contain pitch angle integrated distribution <f> such
!          that integral{<f> d(gamma)} = fsa density.
!          tam5=2*pi*x*cnorm2*gamma*int_i{sinn*dy*vptb*f_code}/zmaxpsi.
         do 40 i=1,iy
            do 30 j=1,jx
               tam4(j)=tam4(j)+f(i,j,k,l_)*cynt2(i,l_)
               tam5(j)=tam5(j)+f(i,j,k,l_)*cynt2(i,l_)*vptb(i,l_)
 30         continue
 40      continue
         do 41 j=1,jx
            tam4(j)=tam4(j)/twoint(l_)
            tam5(j)=tam5(j)/zmaxpsi(lr_)*x(j)*gamma(j)*cnorm2
 41      continue
!-----------------------------
!  Save pitch angle averaged fsa distribution
!  (for seperate print and output). (SC Sept./96)
!bh970403        itstp=n_(l_)/nplot
!bh971027         itstp=n_(l_)/nplot+1
        if (iplot.gt.0 .and. iplot.le.5) then
           do 35 j=1,jx
              f_aveth(j,k,l_,iplot)=tam5(j)
 35        continue
        endif
!-----------------------------

        if (pltfvs.eq."enabled") then

        call aminmx(tam4,1,jxq,1,emin,emax,kmin,kmax)
        do 60 iu=1,5
           if (iu.ne.5) then
              i=imsh(iu)
              do j=1,jxq
                 tam2(j)=f(i,j,k,l_)
              enddo
           else
              do j=1,jxq
                 tam2(j)=tam4(j)
              enddo
           endif
           call aminmx(tam2,1,jxq,1,fmin,fmax,kmin,kmax)
           if (fmin .lt. emin) emin=fmin
           if (fmax .gt. emax) emax=fmax
 60     continue
        emax=emax*1.03
        emin=emax/1.e+12


        DO J=1,JXQ
           RTAM1(J)=TAM1(J)
           TAM2(J)=ABS(TAM2(J))
           TAM2(J)=MAX(em300,TAM2(j))
           RTAM2(J)=LOG10(TAM2(J))
        ENDDO
        REMIN=LOG10(EMIN)
        REMAX=LOG10(EMAX)

#ifndef NOPGPLOT
        CALL PGPAGE
#endif

#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.35,.9)
#endif
        IF ( ReMAX-ReMIN .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           ReMAX= ReMIN+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(Rtam1(1),Rtam1(jxq),Remin,Remax)
#endif
#ifndef NOPGPLOT
        CALL PGSAVE
#endif
#ifndef NOPGPLOT
        CALL PGSCH(0.8)
#endif
#ifndef NOPGPLOT
        CALL PGBOX('BCNST',0.,0,'BCNSTL',0.,0)
#endif
#ifndef NOPGPLOT
        CALL PGUNSA
#endif
#ifndef NOPGPLOT
        CALL PGSAVE
#endif
#ifndef NOPGPLOT
        CALL PGSCH(1.44)
#endif
#ifndef NOPGPLOT
        CALL PGLAB(TX_, 'f', 'Cuts of f vs. v, at cnst pitch angle')
#endif
#ifndef NOPGPLOT
        CALL PGUNSA
#endif

        do 110 iu=1,5
           if (iu.ne.5) then
              i=imsh(iu)
              do j=1,jxq
                 tam2(j)=f(i,j,k,l_)
              enddo
           else
              do j=1,jxq
                 tam2(j)=tam4(j)
              enddo
           endif
           if (iu .eq. 1) text(1)="pll$"
           if (iu .eq. 2) text(1)="trp/ps$"
           if (iu .eq. 3) text(1)="perp$"
           if (iu .eq. 4) text(1)="pll-pi$"
           if (iu .eq. 5) text(1)="avg$"
 90        continue
           do 100 j=1,jxq
              if(tam2(j).le.emin) tam2(j)=emin
              RTAM2(J)=LOG10(TAM2(J))
 100       continue
#ifndef NOPGPLOT
           CALL PGLINE(JXQ,RTAM1,RTAM2)
#endif
          xu=float(iu)
 110   continue





       rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
       if (pltlim.eq.'u/c') then
          write(t_,10021)
       elseif (pltlim.eq.'energy') then
          write(t_,10022)
       else
          write(t_,10020)
       endif

#ifndef NOPGPLOT
        CALL PGMTXT('B',7.,0.,0.,t_)
#endif
        write(T_,10023) K,ENORM
#ifndef NOPGPLOT
        CALL PGMTXT('B',8.,0.,0.,t_)
#endif



       write(t_,10010) n,timet
#ifndef NOPGPLOT
        CALL PGMTXT('B',9.,0.,0.,t_)
#endif
       write(t_,10011) rovera(lr_),rr
#ifndef NOPGPLOT
        CALL PGMTXT('B',10.,0.,0.,t_)
#endif

#ifndef NOPGPLOT
!        CALL PGEND
#endif
!        STOP
#ifndef NOPGPLOT
!        CALL PGPAGE
#endif

       endif


!     Plot pitch angle integrated distribution function <f>, normalized
!     so integral{<f> d(gamma)} = flux surface averaged density.
 120  continue
10010 format("time step (n)=",i5,5x,"time=",e14.6," secs")
10011 format("r/a=",1pe10.2,5x,"radial position(R)=",1pe10.3," cm")
10020 format("Distribution function vs. velocity for some angles")

10021 format("Distribution function vs. u/c for some angles")

10022 format("Distribution function vs. energy(kev) for some angles")

10023 format("Species number=",i1,",   enorm=",1pd10.2)
10030 format("pll ---- theta=0 radians",";", &
           "pll-pi ---- theta=pi radians",";", &
           "trp/ps ---- theta = ",e13.5," radians",";", &
           "perp ---- theta = ",e13.5," radians",";", &
           "avg ---- theta averaged over pi radians","$")

20010 format("time step (n) is",i5,5x,"time=",e14.6," secs",";", &
           "r/a=",e14.6,5x,"radial position (R) =",e14.6," cm","$")
20020 format("Pitch angle integrated electron distribution vs v/vnorm" &
           ,";","(normed so that int{<f> d(gamma)}=fsa n_e", &
           ";","$")
20021 format("Pitch angle integrated electron distribution vs u/c " &
           ,";","(normed so that int{<f> d(gamma)}=fsa n_e", &
           ";","$")
20022 format("Pitch angle integrated electron distribution vs energy " &
           ,";","(normed so that int{<f> d(gamma)}=fsa n_e", &
           ";","$")

      return
      end subroutine pltfvsv


end module pltfvsv_mod
