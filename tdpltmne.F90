module tdpltmne_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use tdpltjop_mod, only : tdpltjop

  !---END USE

!
!

contains

  subroutine tdpltmne
      use param_mod
      use cqlcomm_mod
      use cqlconf_mod, only : setup0
      use r8subs_mod, only : rbound
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!     This routine plots out data as a function of r/a.
!     (Also a time-dependent plot of fusion rate.)
!
#ifdef __MPI
      include 'mpilib.h'
#endif
      dimension trilr(lrza),tr1s(lrorsa),tr2s(lrorsa)
!     real(c_float) variables (and function rbound) for pgplot:
      real(c_float) RPG1,RPG2, RPGmin, RPGmax
      real(c_float) RLRZAP1(0:LRZA),RLRZAP11(0:LRZA),RLRZAP12(0:LRZA), &
           RLRZAP13(0:LRZA)
      real(c_float) RLRZAP(0:LRZA)
      real(c_float) RLRORSA(LRORSA),RLRORSA1(LRORSA),RLRORSA2(LRORSA)
      real(c_float) RNONCHA1(nonch),RNONCHA2(nonch)
      !XXX real(c_float) RBOUND
      character*16 t_horiz

      data em33/1.d-33/

!.......................................................................

#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      if (setup0%noplots.eq."enabled1") return
      if (plt3d .ne. "enabled") return
      if (setup0%lrzmax .le. 2) return
      lrzevn=(setup0%lrzmax/2)*2
      dgts=1.e-8
!..................................................................
!     Determine the x-axis for plots (psi or rho - both normalized).
!..................................................................

      if (pltvs.eq."psi") then
        do 20 l=1,setup0%lrzmax
          tr(l)=(equilpsi(0)-equilpsi(l))/equilpsi(0)
 20     continue
        write(t_horiz,'(a3)') 'psi'
      else
        do 30 l=1,setup0%lrzmax
          tr(l)=rz(l)/radmin
 30     continue
        write(t_horiz,'(a6,a8,a1)') 'rho (=', radcoord, ')'
      endif

      do 21 l=1,setup0%lrz
        trilr(l)=tr(setup0%lrindx(l))
 21   continue



      do 300 k=1,ntotal
        if (k.gt.ngen .and. k.ne.kelec .and. n.gt.1) go to 300
#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)
#endif

#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5,'DENSITIES (/CC) OF SPECIES')

        CALL PGUNSA
#endif

        write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
#ifndef NOPGPLOT
        CALL PGMTXT('T',3.,0.,0.,t_)
#endif
 3002   format("species no. ",i2,4x,a8,2x,a8,2x," time step n=",i4)
        do 19 l=1,setup0%lrzmax
          tr1(l)=reden(k,l)
          tr2(l)=0.0
          tr3(l)=0.0
          if (zmaxpsi(l).ne.0.0 .and. k.le.ngen) then
            tr2(l)=xlndnv(k,l)/zmaxpsi(l)
            tr3(l)=xlndnr(k,l)/zmaxpsi(l)
          endif
 19     continue
        fmin=0.
        fmax=0.
        call aminmx(tr1(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
        if (fmin .ge. .95*fmax) fmin=.95*fmax
        call aminmx(tr2(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin1,fmax1,kmin,kmax)
        if (fmin1.lt.fmin) fmin=fmin1
        if(fmax1.gt.fmax) fmax=fmax1
        call aminmx(tr3(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin1,fmax1,kmin,kmax)
        if (fmin1.lt.fmin) fmin=fmin1
        if(fmax1.gt.fmax) fmax=fmax1
        fmax=fmax*1.05 ! extend range, in case the profile is flat
        fmin=0.0 ! Set lower limit to 0.

        DO I=1,setup0%lrzmax
           RLRZAP1(I)=tr(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
#ifndef NOPGPLOT
        CALL PGSWIN(RLRZAP1(1),RLRZAP1(setup0%lrzmax),RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP11(I)=tr1(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))

        CALL PGSLS(2)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP12(I)=tr2(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))

        CALL PGSLS(3)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP13(I)=tr3(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax))

        CALL PGSLS(1)
#endif
        write(t_,4050) pltvs
 4050   format(a8)
#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGLAB(t_,'density (/cm\u3\d)',' ')

        CALL PGUNSA
#endif


 300  continue

      do 400 k=1,ntotal
        if (k.gt.ngen .and. n.ge.1) go to 400

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)
#endif

#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5,'ENERGIES OF SPECIES IN KEV')

        if(k.le.ngen) CALL PGMTXT('T',5.2,0.5,0.5,'(Solid: <..>_FSA)')

        CALL PGUNSA
#endif
        write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
#ifndef NOPGPLOT
        CALL PGMTXT('T',3.,0.,0.,t_)
#endif

        if (n .eq. 0) then
          do 18 l=1,setup0%lrzmax
            tr1(l)=energy(k,l) ! FSA, solid line
            tr2(l)=0
!     tr3 set to zero so no need to change following plot statements
            tr3(l)=0
 18       continue
          if (k .le. ngen) then
            do 181 ll=1,setup0%lrz
              tr2(setup0%lrindx(ll))=energym(k,ll) ! midplane, dashed line
 181        continue
          endif
        else if (k .le. ngen) then ! and n>0
          do 182 l=1,setup0%lrzmax
            tr1(l)=energy(k,l) ! FSA, solid line
            tr2(l)=0.0
            tr3(l)=0.0
 182      continue
          do 183 l=1,setup0%lrz
            ilr=setup0%lrindx(l)
            tr2(ilr)=energyv(k,ilr) !?? FSA, transport related,  dashed line
            tr3(ilr)=energyr(k,ilr) !?? FSA, transport related,  dot-dashed line
 183      continue
        endif

        fmin=0.
        fmax=0.
        call aminmx(tr1(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
        if (fmin .ge. .95*fmax) fmin=.95*fmax
        call aminmx(tr2(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin1,fmax1,kmin,kmax)
        if (fmin1.lt.fmin) fmin=fmin1
        if(fmax1.gt.fmax) fmax=fmax1
        call aminmx(tr3(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin1,fmax1,kmin,kmax)
        if (fmin1.lt.fmin) fmin=fmin1
        if(fmax1.gt.fmax) fmax=fmax1
        fmax=fmax*1.05 ! extend range, in case the profile is flat

        DO I=1,setup0%lrzmax
           RLRZAP1(I)=tr(i)
           RLRZAP11(I)=tr1(i)
           RLRZAP12(I)=tr2(i)
           RLRZAP13(I)=tr3(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
        RPG1=min(fmin,0.)
#ifndef NOPGPLOT
        CALL PGSWIN(RLRZAP1(1),RLRZAP1(setup0%lrzmax),RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif

#ifndef NOPGPLOT
        CALL PGSLS(1)  !solid

        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))
#endif

#ifndef NOPGPLOT
        CALL PGSLS(2)  !Set dashed line style

        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif

#ifndef NOPGPLOT
        CALL PGSLS(3)  !Set dot-dash line style

        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax))
#endif

#ifndef NOPGPLOT
        CALL PGSLS(1)  !Re-set solid line style for annotation
#endif

        write(t_,4050) pltvs
#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGLAB(t_,'energy (kev)',' ')

        CALL PGUNSA
#endif

 400  continue

      if (setup0%cqlpmod .eq. "enabled") then

        do 350 k=1,ntotal
          if (k.gt.ngen .and. k.ne.kelec .and. n.gt.1) go to 350

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)
#endif

#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5,'DENSITIES (/CC) OF SPECIES')

        CALL PGUNSA
#endif
          write(t_,3050) k,kspeci(1,k),kspeci(2,k)
#ifndef NOPGPLOT
        CALL PGMTXT('T',3.,0.,0.,t_)
#endif
          write(t_,3051) setup0%lrindx(1),rz(setup0%lrindx(1))/radmin
#ifndef NOPGPLOT
        CALL PGMTXT('T',2.,0.,0.,t_)
#endif
 3050     format("species no. ",i2,2x,a8,2x,a8)
 3051     format("r(",i2,")/a=",1pe11.4)


          do 352 ll=1,setup0%lsmax
            tr1s(ll)=denpar(k,ll)
 352      continue
          fmin=0.
          fmax=0.
          call aminmx(tr1s(1:setup0%lsmax),1,setup0%lsmax,1,fmin,fmax,kmin,kmax)
          if (fmin .ge. .95*fmax) fmin=.95*fmax
          fmin=0.d0 ! YuP: make lower limit =0.0
          fmax=fmax*1.05 ! extend range, in case the profile is flat

        DO I=1,setup0%lsmax
           RLRORSA(I)=z(i,setup0%lrindx(1))
           RLRORSA1(I)=tr1s(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRORSA(1)
      RPGmax=RLRORSA(setup0%lsmax)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.

        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

        CALL PGLINE(setup0%lsmax,RLRORSA(1:setup0%lsmax),RLRORSA1(1:setup0%lsmax))
#endif

          write(t_,4051)
 4051     format("s (cms)")

#ifndef NOPGPLOT
          CALL PGLAB(t_,'Density (/cm\u3\d)',' ')
#endif

 350    continue

        do 450 k=1,ntotal
          if (k.gt.ngen .and. n.gt.1) go to 450

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)
#endif

#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5,'ENERGIES OF SPECIES IN KEV')

        CALL PGMTXT('T',5.,0.,0., &
               'Solid: midplane;  Dashed: <..>_FSA')
        CALL PGUNSA
#endif
          write(t_,3050) k,kspeci(1,k),kspeci(2,k)
#ifndef NOPGPLOT
        CALL PGMTXT('T',3.,0.,0.,t_)
#endif
          write(t_,3051) setup0%lrindx(1),rz(setup0%lrindx(1))/radmin
#ifndef NOPGPLOT
        CALL PGMTXT('T',2.,0.,0.,t_)
#endif

          do 452 ll=1,setup0%lsmax
            tr1s(ll)=enrgypa(k,ll) ! Midplane
 452      continue
          fmin=0.
          fmax=0.
          call aminmx(tr1s(1:setup0%lsmax),1,setup0%lsmax,1,fmin,fmax,kmin,kmax)
          if (fmin .ge. .95*fmax) fmin=.95*fmax
          fmin=0.d0 ! YuP: make lower limit =0.0
          fmax=fmax*1.05 ! extend range, in case the profile is flat

        DO I=1,setup0%lsmax
           RLRORSA(I)=z(i,setup0%lrindx(1))
           RLRORSA1(I)=tr1s(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRORSA(1)
      RPGmax=RLRORSA(setup0%lsmax)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.

        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

        CALL PGLINE(setup0%lsmax,RLRORSA(1:setup0%lsmax),RLRORSA1(1:setup0%lsmax))
#endif

          write(t_,4051)

#ifndef NOPGPLOT
          CALL PGLAB(t_,'Energy (keV)',' ')
#endif


 450    continue

      endif

!.......................................................................
!     n > 0
!.......................................................................

      if (n .gt. 0) then

        do 500 k=1,ngen

 5001     format(2("curr(",e12.5,")=",e15.5," Amps per cm**2",3x),"$")
 5011     format(("curr(",e12.5,")=",e15.5," Amps per cm**2",3x),"$")

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)
#endif

#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5, &
                    'FLUX SURF. AV. CURNT. (AMPS/CM\u2\d)')
        CALL PGUNSA
#endif
          write(t_,5002) k, currza(k)
 5002     format("Species:",i2,"  Total current =",e16.6," Amps")
#ifndef NOPGPLOT
        CALL PGMTXT('T',3.,0.,0.,t_)
#endif



          do 69 l=1,setup0%lrzmax
            tr2(l)=currr(k,l)
            tr3(l)=currv_(k,l)
            tr1(l)=currz(k,l)
 69       continue
          fmin=0.
          fmax=0.
          call aminmx(tr1(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
          if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
          call aminmx(tr2(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin1,fmax1,kmin,kmax)
          if (fmin1.lt.fmin) fmin=fmin1
          if(fmax1.gt.fmax) fmax=fmax1
          call aminmx(tr3(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin1,fmax1,kmin,kmax)
          if (fmin1.lt.fmin) fmin=fmin1
          if(fmax1.gt.fmax) fmax=fmax1
          if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range


        DO I=1,setup0%lrzmax
           RLRZAP1(I)=tr(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
#ifndef NOPGPLOT
        CALL PGSWIN(RLRZAP1(1),RLRZAP1(setup0%lrzmax),RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP11(I)=tr1(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))

        CALL PGSLS(2)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP12(I)=tr2(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))

        CALL PGSLS(3)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP13(I)=tr3(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax))

        CALL PGSLS(1) ! restore solid line
#endif
        write(t_,4050) pltvs
#ifndef NOPGPLOT
        CALL PGLAB(t_,'Curr Den (A/cm\u2\d)',' ')
#endif







!
          if (setup0%cqlpmod.eq."enabled") then

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)
#endif

#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5, &
             'LOCAL IN s CURNT. DENS. (AMPS per CM**2)')
        CALL PGUNSA
#endif
          write(t_,3050) k,kspeci(1,k),kspeci(2,k)
#ifndef NOPGPLOT
        CALL PGMTXT('T',3.,0.,0.,t_)
          write(t_,3051) setup0%lrindx(1),rz(setup0%lrindx(1))/radmin
        CALL PGMTXT('T',2.,0.,0.,t_)
#endif





            do 551 ll=1,setup0%lsmax
              tr1s(ll)=pcurr(nch(ll),k,ll)
              tr2s(ll)=psis(ll)*pcurr(nch(1),k,1)
 551        continue
            ilsmx=setup0%lsmax
            if (sbdry .eq. "periodic") then
              ilsmx=setup0%lsmax+1
              tr1s(ilsmx)=tr1s(1)
              tr2s(ilsmx)=tr2s(1)
            endif
            fmin=0.
            fmax=0.
            call aminmx(tr1s(1:ilsmx),1,ilsmx,1,fmin,fmax,kmin,kmax)
            if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
            call aminmx(tr2s(1:ilsmx),1,ilsmx,1,fmin1,fmax1,kmin,kmax)
            if (fmin1.lt.fmin) fmin=fmin1
            if(fmax1.gt.fmax) fmax=fmax1
            if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range


        DO I=1,setup0%lsmax
           RLRORSA(I)=z(i,setup0%lrindx(1))
           RLRORSA1(I)=tr1s(i)
           RLRORSA2(I)=tr2s(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RLRORSA(1),RLRORSA(setup0%lsmax),RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

        CALL PGLINE(setup0%lsmax,RLRORSA(1:setup0%lsmax),RLRORSA1(1:setup0%lsmax))

        CALL PGSLS(2)

        CALL PGLINE(setup0%lsmax,RLRORSA(1:setup0%lsmax),RLRORSA2(1:setup0%lsmax))

        CALL PGSLS(1)
#endif
        write(t_,4051)
#ifndef NOPGPLOT
        CALL PGLAB(t_,'Curr Den (A/cm\u2\d',' ')
#endif



            go to 500
          endif
!
          if (nrf .ne. 0) then

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)
#endif

#ifndef NOPGPLOT
        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5, &
             'LOCAL RF POWER (WATTS per CC)')
        CALL PGUNSA
#endif
          write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
#ifndef NOPGPLOT
        CALL PGMTXT('T',3.,0.,0.,t_)
#endif

            write(t_,6002) rfpwrt(k)
 6002       format("total rf power =",1pe10.2," Watts")
#ifndef NOPGPLOT
        CALL PGMTXT('T',2.,0.,0.,t_)
#endif

            do 99 ll=1,setup0%lrzmax
 99         tr1(ll)=rfpwrz(k,ll)
            fmin=0.
            fmax=0.
            call aminmx(tr1(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
!BH090220            if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
        if (abs(fmin-fmax).lt.fmax*dgts) fmax=fmin+.001*abs(fmin)
        if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range

        DO I=1,setup0%lrzmax
!BH090220           RLRZAP1(I)=tr(i)
           RLRZAP1(I)=RBOUND(tr(i))
        ENDDO
!BH090220        RPG1=fmin
!BH090220        RPG2=fmax
        RPG1=RBOUND(fmin)
        RPG2=RBOUND(fmax)
        if (abs(RPG1-RPG2).le.(2.1*em33)) RPG2=RPG1+10.*em33
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RLRZAP1(1),RLRZAP1(setup0%lrzmax),RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP11(I)=tr1(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))

        CALL PGLAB('Radius','RF Power (W/cm\u3\d)',' ')
#endif


          endif
!
          if(syncrad.ne."disabled" .and. k.eq.kelecg) then

#ifndef NOPGPLOT
            CALL PGPAGE

            CALL PGSVP(.2,.8,.2,.6)

            CALL PGSAVE

            CALL PGSCH(1.44)

            CALL PGMTXT('T',6.,0.5,0.5, &
                 'LOCAL RF POWER (WATTS per CC)')
            CALL PGUNSA
#endif
            write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
#ifndef NOPGPLOT
            CALL PGMTXT('T',3.,0.,0.,t_)
            write(t_,6004) psynct
 6004       format(";","synchrotron radiated power =",e16.6," Watts")
            CALL  PGMTXT('T',2.,0.,0.,t_)
#endif

            do 79 ll=1,setup0%lrzmax
 79            tr1(ll)=psyncz(ll)
            fmin=0.
            fmax=0.
            call aminmx(tr1(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
            if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
            if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range

            DO I=1,setup0%lrzmax
               RLRZAP1(I)=tr(i)
            ENDDO
            RPG1=fmin
            RPG2=fmax
            IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
               RPG2= RPG1+1.e-16
            ENDIF
#ifndef NOPGPLOT
            CALL PGSWIN(RLRZAP1(1),RLRZAP1(setup0%lrzmax),RPG1,RPG2)

            CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
            DO I=1,setup0%lrzmax
               RLRZAP11(I)=tr1(i)
            ENDDO
#ifndef NOPGPLOT
            CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))

            CALL PGLAB('Radius','Synch Power (W/cm\u3\d)',' ')
#endif

          endif ! syncrad.ne."disabled" .and. k.eq.kelecg
!
          if(sigmamod.eq."enabled" .and. pltsig.eq."enabled" &
                                   .and. k.eq.kiong(1)) then
            do 800 lsig=1,4
               if(isigmas(lsig).eq.0) go to 800

#ifndef NOPGPLOT
               CALL PGPAGE

               CALL PGSVP(.2,.8,.2,.5)

               CALL PGSAVE

               CALL PGSCH(1.44)

               CALL PGMTXT('T',8.,0.5,0.5, &
                    'Fusion Power (Watts/cm\u3\d)')
               CALL PGUNSA
#endif
               write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
#ifndef NOPGPLOT
               CALL PGMTXT('T',7.,0.,0.,t_)
#endif
               write(t_,8001) lsig
 8001          format("Fusion reaction number ",i3)
#ifndef NOPGPLOT
               CALL PGMTXT('T',6.,0.,0.,t_)
#endif
               if (lsig.eq.1) then
                  write(t_,*) "d+t=>alpha(3.5MeV)+n(14.1MeV)"
               elseif (lsig.eq.2) then
                  write(t_,*) "d+he3=>alpha(3.6MeV)+p(14.7MeV)"
               elseif (lsig.eq.3) then
                  write(t_,*) "d+d=>t(1.01Mev)+p(3.02MeV)"
               elseif (lsig.eq.4) then
                  write(t_,*) "d+d=>he3(.82MeV)+n(2.45MeV)"
               endif
#ifndef NOPGPLOT
               CALL PGMTXT('T',5.,0.,0.,t_)
#endif

               write(t_,8002) fuspwrvt(lsig)
 8002          format("fusion power =",1pe16.6," Watts")
#ifndef NOPGPLOT
               CALL PGMTXT('T',3.,0.,0.,t_)
#endif

!BH120314:  Have removed the isigsgv2 option and references to it, since
!BH120314:  this functionality appears to have no physical use.
!BH120314               if (isigsgv2.eq.1) then
!BH120314                  write(t_,8003) fuspwrmt(lsig)
!BH120314 8003           format("fusion power(equiv. Maxwln) =",1pe16.6," Watts")
#ifndef NOPGPLOT
!BH120314                  CALL PGMTXT('T',2.,0.,0.,t_)
#endif
!BH120314               endif

               do 801 ll=1,setup0%lrzmax
                  tr1(ll)=fuspwrv(lsig,ll)
                  tr2(ll)=fuspwrm(lsig,ll)
 801           continue

               fmin=0.
               fmax=0.
               call aminmx(tr1(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
!BH120314               if (isigsgv2.eq.1) then
!BH120314                 call aminmx(tr2(1),1,setup0%lrzmax,1,fmin1,fmax1,kmin,kmax)
!BH120314                 fmin=min(fmin,fmin1)
!BH120314                 fmax=max(fmax,fmax1)
!BH120314               endif
               if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
               if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range

               DO I=1,setup0%lrzmax
                  RLRZAP1(I)=tr(i)
               ENDDO
               RPG1=fmin
               RPG2=fmax
               IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
                  RPG2= RPG1+1.e-16
               ENDIF
#ifndef NOPGPLOT
               CALL PGSWIN(RLRZAP1(1),RLRZAP1(setup0%lrzmax),RPG1,RPG2)

               CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
               DO I=1,setup0%lrzmax
                  RLRZAP11(I)=tr1(i)
               ENDDO
#ifndef NOPGPLOT
               CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))

               CALL PGSLS(2)
#endif
               DO I=1,setup0%lrzmax
                  RLRZAP12(I)=tr2(i)
               ENDDO
#ifndef NOPGPLOT
               CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))

               CALL PGSLS(1)

               CALL PGSAVE

               CALL PGSCH(1.44)

               CALL PGLAB(t_,'Fusion Power (W/cm/u3/d',' ')

               CALL PGUNSA
!...           mnt  Generate time-dependent plot of total fusion rates
               CALL PGPAGE
               CALL PGSVP(.2,.8,.2,.5)
               CALL PGSAVE
               CALL PGSCH(1.44)
               CALL PGMTXT('T',8.,0.5,0.5, &
                  'Fusion Rx Rate (/sec) Vs Time')
               CALL PGUNSA
               write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
               CALL PGMTXT('T',7.,0.,0.,t_)
               write(t_,8005) lsig
 8005          format('Fusion reaction number ',i3)
               CALL PGMTXT('T',6.,0.,0.,t_)
               write(t_,8006) sigftt(nch(1),lsig)
 8006          format("Reaction rate =",1pe16.6," /sec")
               CALL PGMTXT('T',4.,0.,0.,t_)
#endif

!BH120314               if (isigsgv2.eq.1) then
!BH120314                 write(t_,8007) sigmtt(nch(1),lsig)
!BH120314 8007           format("Reaction rate(equiv. Maxwln) =",1pe14.6," /sec")
#ifndef NOPGPLOT
!BH120314                 CALL PGMTXT('T',2.,0.,0.,t_)
#endif
!BH120314               endif

               call aminmx(sigftt(1:nch(1),lsig),1,nch(1),1, &
               emin,emax,kmin,kmax)
!BH120314               if (isigsgv2.eq.1) then
!BH120314                  call aminmx(sigmtt(1,lsig),1,nch(1),1,
!BH120314     ~            emin1,emax1,kmin,kmax)
!BH120314                  emin=min(emin,emin1)
!BH120314                  emax=max(emax,emax1)
!BH120314               endif
               dgts=1.e-8
               if (abs(emin-emax).lt.emax*dgts) emax=emin+.001*abs(emin)
               if(emax.gt.0.) emax=emax*1.05 ! extend the upper range
               DO I=1,NCH(1)
                  RNONCHA1(I)=ptime(i,1)
               ENDDO
               RPG1=emin
               RPG2=emax
               IF ( RPG2-RPG1 .le. 1.e-16 ) THEN
                   RPG2= RPG1+1.e-16
               ENDIF
#ifndef NOPGPLOT
               CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(1)),RPG1,RPG2)

               CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

               CALL PGLAB('Time (sec)','Fusion Rx Rate (/sec)',' ')
#endif
               DO I=1,NCH(1)
                  RNONCHA2(I)=sigftt(I,lsig)
               ENDDO
#ifndef NOPGPLOT
               CALL PGLINE(nch(1),RNONCHA1,RNONCHA2)
#endif
!BH120314               if (isigsgv2.eq.1) then
!BH120314                  DO I=1,NCH(1)
!BH120314                     RNONCHA2(I)=sigmtt(I,lsig)
!BH120314                  ENDDO
#ifndef NOPGPLOT
!BH120314                  CALL PGLINE(nch(1),RNONCHA1,RNONCHA2)
#endif
!BH120314               endif

 800        continue ! lsig=1,4
          endif !if(sigmamod="enabled" .and. pltsig="enabled" .and. k=kiong(1)

!
!$$$  if (tauegy(k,0).gt.0.) then
!$$$  call gxglfr(0)
!$$$            call gscpvs(.5,.95)
!$$$            call gstxjf("center","top")
!$$$            call gstxno(40.)
!$$$            call gitxft(ift)
!$$$            if (ift .ne. ift07a) call gstxft(ift07a)
!$$$            call gptx2d("LOCAL LOSS POWER (WATTS per CC);$")
!$$$            call gitxft(ift)
!$$$            if (ift.ne.ift21a) call gstxft(ift21a)
!$$$            call gstxno(100.)
!$$$            call gstxjf("left","top")
!$$$            call gscpvs(.01,.85)
!$$$            write(t_,3002) k,kspeci(1,k),kspeci(2,k)
!$$$            call gptx2d(t_)
!$$$            do 6030 l=1,setup0%lrzmax-1,2
!$$$              write(t_,6005) tr(l),pegyz(k,l),tr(l+1),pegyz(k,l+1)
!$$$              call gptx2d(t_)
!$$$ 6030       continue
!$$$            if (setup0%lrzmax.ne.lrzevn) then
!$$$              write(t_,6031) tr(setup0%lrzmax),pegyz(k,setup0%lrzmax)
!$$$              call gptx2d(t_)
!$$$            endif
!$$$ 6005       format(2("pegyz(",e12.5,")=",e15.5," Watts per cc",3x),"$")
!$$$ 6031       format(("pegyz(",e12.5,")=",e15.5," Watts per cc",3x),"$")
!$$$            write(t_,6006) pegyt(k)
!$$$            call gptx2d(t_)
!$$$ 6006       format(";","total phenomenological energy loss =",e16.6,
!$$$     1        " Watts","$")
!$$$
!$$$            do 89  l=1,setup0%lrzmax
!$$$ 89         tr1(l)=pegyz(k,l)
!$$$            fmin=0.
!$$$            fmax=0.
!$$$            call aminmx(tr1(1),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
!$$$            if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
!$$$            call gswd2d("linlin$",tr(1),tr(setup0%lrzmax),fmin,fmax)
!$$$            call gsvp2d(.2,.8,.2,.5)
!$$$            call gpgr80("linlin$")
!$$$            call gpcv2d(tr(1),tr1(1),setup0%lrzmax)
!$$$C%OS  call gxglfr(0)
!$$$          endif
!$$$
!$$$          if (torloss(k).ne."disabled") then
!$$$            call gxglfr(0)
!$$$            call gscpvs(.5,.95)
!$$$            call gstxjf("center","top")
!$$$            call gstxno(40.)
!$$$            call gitxft(ift)
!$$$            if (ift .ne. ift07a) call gstxft(ift07a)
!$$$            call gptx2d("Particle Loss Power (Watts per cc);$")
!$$$            call gitxft(ift)
!$$$            if (ift.ne.ift21a) call gstxft(ift21a)
!$$$            call gstxno(100.)
!$$$            call gstxjf("left","top")
!$$$            call gscpvs(.01,.85)
!$$$            write(t_,3002) k,kspeci(1,k),kspeci(2,k)
!$$$            call gptx2d(t_)
!$$$            do 6040 l=1,setup0%lrzmax-1,2
!$$$              write(t_,6007) tr(l),pplossz(k,l),tr(l+1),pplossz(k,l+1)
!$$$              call gptx2d(t_)
!$$$ 6040       continue
!$$$            if (setup0%lrzmax.ne.lrzevn) then
!$$$              write(t_,6041) tr(setup0%lrzmax),pplossz(k,setup0%lrzmax)
!$$$              call gptx2d(t_)
!$$$            endif
!$$$ 6007       format(2("pplossz(",e12.5,")=",e15.5," Watts per cc",3x),
!$$$     +        "$")
!$$$ 6041       format(("pplossz(",e12.5,")=",e15.5," Watts per cc",3x),"$")
!$$$            write(t_,6008) pplosst(k)
!$$$            call gptx2d(t_)
!$$$ 6008       format(";","total particle loss power =",e16.6," Watts","$")
!$$$
!$$$            do 109  l=1,setup0%lrzmax
!$$$ 109        tr1(l)=pegyz(k,l)
!$$$            fmin=0.
!$$$            fmax=0.
!$$$            call aminmx(tr1(1),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
!$$$            if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
!$$$            call gswd2d("linlin$",tr(1),tr(setup0%lrzmax),fmin,fmax)
!$$$            call gsvp2d(.2,.8,.2,.5)
!$$$            call gpgr80("linlin$")
!$$$            call gpcv2d(tr(1),tr1(1),setup0%lrzmax)
!$$$            call gxglfr(0)
!$$$          endif
!
 500    continue ! k




        call tdpltjop



!     endif n>0, before 500 loop
      endif

      return
      end subroutine tdpltmne
end module tdpltmne_mod
