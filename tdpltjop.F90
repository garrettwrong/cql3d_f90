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

module tdpltjop_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx

  !---END USE


!
!

contains

  subroutine tdpltjop
      use param_mod
      use cqlcomm_mod
      use cqlconf_mod, only : setup0
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
#ifdef __MPI
      include 'mpilib.h'
#endif

      REAL RILIN
      REAL RPG1,RPG2, RPGmin, RPGmax
      REAL RLRZAP1(0:LRZA),RLRZAP11(0:LRZA),RLRZAP12(0:LRZA), &
           RLRZAP13(0:LRZA),RLRZAP14(0:LRZA)
      real(c_float) RLRZAP(0:LRZA)
      character*16 t_horiz

!..................................................................
!     This routine plots a number of current drive diagnostics
!     for NB, LH and ECH excitation. The x axis of the plots are
!     rho/rhomax or psi/psimax depending on the value of bcd
!     word psival. Below "fi" means fast ions (usually in conjunction
!     with NBI. "fi+e" would include the electron contribution to
!     the fast ions, this would be either from Coulomb drag on the
!     electrons by the fast ions or RF excitation.
!..................................................................

#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      if (setup0%noplots.eq."enabled1") return

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
          tr(l)=rya(l) !YuP: was rz(l)/radmin
          !write(*,*)'lr,radmin,rz(lr),rya(lr)=',l,radmin,rz(l),rya(l)
 30     continue
        write(t_horiz,'(a6,a8,a1)') 'rho (=', radcoord, ')'
      endif

      do l=1,setup0%lrzmax ! for plots of profiles vs R
          RLRZAP(l)=rpcon(l) ! R_outermost of flux surf.
      enddo


!..................................................................
!     plots of current..
!..................................................................
      ! Choose which j_bs component to plot:
      ! thermal(=maxwellian) or non-thermal(=general species)
      if(kelecg.ne.0)then
         kke=2 ! I_bs for e_general (non-thermal)
      else
         kke=1 ! I_bs for e_maxw
      endif
      if(niong.ne.0)then
         kki=2 ! I_bs for i_general (non-thermal)
      else
         kki=1 ! I_bs for i_maxw
      endif
      if (jhirsh.eq.0) then
         ! In this case, j_bs is only calculated for maxwellian part.
         ! (In fact, it is done for electrons only)
         kke=1
         kki=1
      endif

      fmin=0.
      fmax=0.
      call aminmx(currtz(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
      call aminmx(currtpz(1:setup0%lrzmax),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
      if (gmin.lt.fmin) fmin=gmin
      if (fmax.lt.gmax) fmax=gmax
      call aminmx(bscurm(1:setup0%lrzmax,1,kke),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
      if (gmin.lt.fmin) fmin=gmin
      if (fmax.lt.gmax) fmax=gmax
      call aminmx(bscurm(1:setup0%lrzmax,2,kki),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
      if (gmin.lt.fmin) fmin=gmin
      if (fmax.lt.gmax) fmax=gmax
      if (fmax-fmin.lt.1.e-8) fmin=fmax-.1*abs(fmax)-1.e-5
      RPG1=min(fmin,0.) !Make lower limit 0 when fmin>0
      RPG2=max(fmax,0.)

      DO I=1,setup0%lrzmax
         RLRZAP1(I)=tr(i) ! rho or psi
         RLRZAP11(I)=currtz(i)
         RLRZAP12(I)=currtpz(i)
         RLRZAP13(I)=bscurm(i,1,kke) !bscurm(*,1,*) is for electrons
         RLRZAP14(I)=bscurm(i,2,kki) !bscurm(*,2,*) is for ions
      ENDDO

      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRZAP1(1)
      RPGmax=RLRZAP1(setup0%lrzmax)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.

!yup      call aminmx(totcurz(1),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
!yup      if (gmin.lt.fmin) fmin=gmin
!yup      if (fmax.lt.gmax) fmax=gmax
!yup [May-2014] Do not plot the total current anymore,
!yup because in present setup, it may include
!yup a bootstrap analytical current (jhirsh88/99) for electrons.
!yup But if bootcalc='method1' or bootcalc='method2',
!yup totcurz() will also include a numerical bootstrap current
!yup for a given general species. In such a case, totcurz() will
!yup include the bootstrap current from both analytical model
!yup and numerical model; makes no sense.
!yup For now, make plots of current density profiles from:
!yup general ions (if any), general electrons (it will be in "fi+e").
!yup Also, plotted profiles of bootstrap current (jhirsh88/99)
!yup separately for electrons and ions, based on maxwellian T,n profiles
!yup or, if available, for general electrons or ions;
!yup they are designated as "bs_e" and "bs_i" in the plots.

      write(t_,4040)
 4040 format("Current summed over all species",";", &
        "fi - fast ion current",3x, &
        "fi+e - fi + electrons" &
        ,";","bs - Bootstrap current",3x,"tot - total current","$")


#ifndef NOPGPLOT
        CALL PGPAGE
#endif

        ! FIRST PANEL: profiles vs rho
#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.6,.9)
#endif
#ifndef NOPGPLOT
        CALL PGSAVE
#endif
#ifndef NOPGPLOT
        CALL PGSCH(1.44)
#endif
#ifndef NOPGPLOT
        CALL PGMTXT('T',2.5,0.5,0.5,'CURRENT (AMPS/CM\u2\d)')
#endif
#ifndef NOPGPLOT
        CALL PGUNSA
#endif
!..................................................................
!     currents, printed and plotted next
!..................................................................
 4012 format("fi [solid]=  ",1pe10.3,5x,"fi+e[--]= ",  1pe10.3)
 4013 format("bs_e[-.-]= ",1pe10.3,5x,"bs_i[.....]= ", 1pe10.3, " Amps")
        write(t_,4012) currtza,currtpza
#ifndef NOPGPLOT
        CALL PGMTXT('T',2.,0.,0.,t_)
#endif
        write(t_,4013) bscurma(1,kke),bscurma(2,kki) ! (e,*), (i,*)
#ifndef NOPGPLOT
        CALL PGMTXT('T',1.,0.,0.,t_)
#endif
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)
#endif
#ifndef NOPGPLOT
        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
        ! fi (general ions):
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))
#endif
        ! fi+e [can be general ions (if any) + general electrons,
        ! or general ions + screening current from e; see eleccomp='enabled' option]:
#ifndef NOPGPLOT
        CALL PGSLS(2) ! ---
#endif
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif
        ! Bootstrap for e, based on jhirsh88/99 models:
#ifndef NOPGPLOT
        CALL PGSLS(3) ! -.-
#endif
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax))
#endif
        ! Bootstrap for ions, based on jhirsh88/99 models:
#ifndef NOPGPLOT
        CALL PGSLS(4) ! ...
#endif
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP14(1:setup0%lrzmax))
#endif
#ifndef NOPGPLOT
        CALL PGSLS(1) ! Restore solid line
#endif
!yup        ! Total: [YuP: do not plot anymore]
!yup        DO I=1,setup0%lrzmax
!yup           RLRZAP13(I)=totcurz(i)
!yup        ENDDO
#ifndef NOPGPLOT
!yup        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax))
#endif
#ifndef NOPGPLOT
        CALL PGSAVE
#endif
#ifndef NOPGPLOT
        CALL PGSCH(1.44)
#endif
#ifndef NOPGPLOT
        CALL PGLAB(' ','curr density (A/cm\u2\d)',' ')
#endif
#ifndef NOPGPLOT
        CALL PGMTXT('B',1.8,0.5,0.5,t_horiz)
#endif
#ifndef NOPGPLOT
        CALL PGUNSA
#endif

        ! SECOND PANEL: same profiles vs R
#ifndef NOPGPLOT
        CALL PGSVP(.2,.8,.2,.5)
#endif
        RPGmin=min(RLRZAP(1),rmag)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RPGmin,RLRZAP(setup0%lrzmax),RPG1,RPG2)
#endif
#ifndef NOPGPLOT
        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
        ! fi (general ions):
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))
#endif
        ! fi+e [can be general ions (if any) + general electrons,
        ! or general ions + screening current from e; see eleccomp='enabled' option]:
#ifndef NOPGPLOT
        CALL PGSLS(2) ! ---
#endif
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif
        ! Bootstrap for e, based on jhirsh88/99 models:
#ifndef NOPGPLOT
        CALL PGSLS(3) ! -.-
#endif
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax))
#endif
        ! Bootstrap for ions, based on jhirsh88/99 models:
#ifndef NOPGPLOT
        CALL PGSLS(4) ! ...
#endif
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP(1:setup0%lrzmax),RLRZAP14(1:setup0%lrzmax))
#endif
#ifndef NOPGPLOT
        CALL PGSLS(1) ! Restore solid line
#endif
#ifndef NOPGPLOT
        CALL PGSAVE
#endif
#ifndef NOPGPLOT
        CALL PGSCH(1.44)
#endif
#ifndef NOPGPLOT
        CALL PGLAB(' ','curr density (A/cm\u2\d)',' ')
#endif
#ifndef NOPGPLOT
        CALL PGMTXT('B',1.8,0.5,0.5,'R (=rpcon)  (cm)')
#endif
#ifndef NOPGPLOT
        CALL PGUNSA
#endif



!..................................................................
!     Partially integrated in rho (or psi)
!..................................................................

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)

        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5, &
             'CURRENT (AMPS)')

        CALL PGMTXT('T',5.,0.5,0.5, &
                    '(INTEGRATED UP TO RHO or PSI)')
        CALL PGUNSA
#endif

!..................................................................
!     plots of partial integration of currents
!..................................................................


!%OS  call gscpvs(.5,.95)
!%OS  call gstxjf("center","top")
!%OS  call gstxno(43.)
!%OS  call gitxft(ift)
!%OS  if (ift .ne. ift07a) call gstxft(ift07a)
!%OS  call gptx2d("CURRENT (AMPS/CM**2);$")
!%OS  call gscpvs(.5,.9)
!%OS  call gptx2d("(INTEGRATED UP TO RHO OR PSI);$")
!%OS  call gitxft(ift)
!%OS  if (ift.ne.ift21a) call gstxft(ift21a)
!%OS  call gstxno(100.)
!%OS  call gscvlb(1)
      fmin=0.
      fmax=0.

      call aminmx(currtzi(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
      call aminmx(currtpzi(1:setup0%lrzmax),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
      if (gmin.lt.fmin) fmin=gmin
      if (fmax.lt.gmax) fmax=gmax
!yup      call aminmx(totcurzi(1),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
!yup      if (gmin.lt.fmin) fmin=gmin
!yup      if (fmax.lt.gmax) fmax=gmax
      if (fmax-fmin.lt.1.e-8) fmin=fmax-.1*abs(fmax)-1.e-5
      if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range
      RPG1=fmin
      RPG2=fmax
      RPG1=min(fmin,0.) !Make lower limit 0 when fmin>0

      DO I=1,setup0%lrzmax
         RLRZAP1(I)=tr(i)
      ENDDO
      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRZAP1(1)
      RPGmax=RLRZAP1(setup0%lrzmax)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.

        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP11(I)=currtzi(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))

        CALL PGSLS(2)
#endif
        DO I=1,setup0%lrzmax
           RLRZAP12(I)=currtpzi(i)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))

!yup        CALL PGSLS(3)
#endif
!yup        DO I=1,setup0%lrzmax
!yup           RLRZAP13(I)=totcurzi(i)
!yup        ENDDO
#ifndef NOPGPLOT
!yup        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax))

        CALL PGSLS(1)

        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGLAB(' ','current  (Amps)',' ')

        CALL PGMTXT('B',1.8,0.5,0.5,t_horiz)

        CALL PGUNSA
#endif


!..................................................................
!     Source Power...
!..................................................................

!-----------------------------------------------------------------------
!     IF NO POWER NO PLOT
!-----------------------------------------------------------------------
!%OS
!     Adjust for kfrsou=0 (can occur when no NBI):
      if (kfrsou.ne.0) then
         kfrsou1=kfrsou
      else
         kfrsou1=1
      endif

      IF (sorpwtza .le. 1.0e-25) go to 809
!%OS
      ! PRINT OUT OF SOURCE POWER (WATTS/CC):
      ! sorpwt(l),sorpw_nbi(kfrsou,l),sorpw_rf(1,l),sorpw_rf(2,l),sorpw_rf(3,l)
#ifndef NOPGPLOT
      CALL PGPAGE

      CALL PGSVP(.05,.95,.05,.95)
#endif
      RILIN=0.
#ifndef NOPGPLOT
      CALL PGSAVE

      CALL PGSCH(1.44)

      CALL PGMTXT('T',-RILIN,.5,.5, &
           'SOURCE POWER: (WATTS/CC)')
#endif
      RILIN=RILIN+3.
#ifndef NOPGPLOT
      CALL PGSCH(1.)
#endif
      write(t_,6013) pltvs
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+2.
      write(t_,6014) pltvs
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+2.

 6013 format(4x,a8,"   NBI+RF      NBI", &
        "         RF(1)      RF(2)       RF(3)")
 6014 format(4x,a8,"   (sorpwt)   (sorpw_nbi)", &
        "      (sorpw_rf for gen.species 1,2,3)")

!     Start printing results on first page
      do 10  l=1,min(40,setup0%lrzmax)
        if (ngen.eq.1) then
          write(t_,6015) tr(l),sorpwt(l),sorpw_nbi(kfrsou1,l), &
            sorpw_rf(1,l)
        elseif (ngen.eq.2) then
          write(t_,6016) tr(l),sorpwt(l),sorpw_nbi(kfrsou1,l), &
            sorpw_rf(1,l),sorpw_rf(2,l)
        else
          write(t_,6017) tr(l),sorpwt(l),sorpw_nbi(kfrsou1,l), &
            sorpw_rf(1,l),sorpw_rf(2,l),sorpw_rf(3,l)
        endif
#ifndef NOPGPLOT
        CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
        RILIN=RILIN+1.
 10   continue

!     Continue printing results on second page, if setup0%lrzmax.gt.40
      if (setup0%lrzmax.gt.40) then
#ifndef NOPGPLOT
         CALL PGPAGE

         CALL PGSVP(.05,.95,.05,.95)
#endif
         RILIN=0.+3.
         do l=41,setup0%lrzmax
            if (ngen.eq.1) then
               write(t_,6015) tr(l),sorpwt(l),sorpw_nbi(kfrsou1,l), &
                    sorpw_rf(1,l)
            elseif (ngen.eq.2) then
               write(t_,6016) tr(l),sorpwt(l),sorpw_nbi(kfrsou1,l), &
                    sorpw_rf(1,l),sorpw_rf(2,l)
            else
               write(t_,6017) tr(l),sorpwt(l),sorpw_nbi(kfrsou1,l), &
                    sorpw_rf(1,l),sorpw_rf(2,l),sorpw_rf(3,l)
            endif
#ifndef NOPGPLOT
            CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
            RILIN=RILIN+1.
         enddo
      endif

 6015 format(1x, 1pe9.3, 3(1x,e10.3) )
 6016 format(1x, 1pe9.3, 4(1x,e10.3) )
 6017 format(1x, 1pe9.3, 5(1x,e10.3) )


      RILIN=RILIN+2.
      write(t_,6022) sorpwtza !sorpwtza=sorpwti(setup0%lrzmax) !NBI+RF, all gen.species
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+1.
      write(t_,6023) sorpw_nbii(kfrsou1,setup0%lrzmax)     ! NBI only
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+1.
      if (ngen.ge.1) then
         write(t_,6024) sorpw_rfi(1,setup0%lrzmax) ! RF(1st gen.species) only
#ifndef NOPGPLOT
         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
         RILIN=RILIN+1.
      endif
      if (ngen.ge.2) then
         write(t_,6025) sorpw_rfi(2,setup0%lrzmax) ! RF(2nd gen.species) only
#ifndef NOPGPLOT
         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
         RILIN=RILIN+1.
      endif
      if (ngen.ge.3) then
         write(t_,6026) sorpw_rfi(3,setup0%lrzmax) ! RF(3rd gen.species) only
#ifndef NOPGPLOT
         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
         RILIN=RILIN+1.
      endif
#ifndef NOPGPLOT
      CALL PGUNSA
#endif
 6022 format("Power integr. over radius (RF+NBI, all gen.species)=", &
             1pe12.4,"Watts")
 6023 format("Power from NBI (sorpw_nbii)=",1pe12.4,"Watts")
 6024 format("Power from RF  (sorpw_rfi) Gen.species no.1 =", &
             1pe12.4,"Watts")
 6025 format("Power from RF  (sorpw_rfi) Gen.species no.2 =", &
             1pe12.4,"Watts")
 6026 format("Power from RF  (sorpw_rfi) Gen.species no.3 =", &
             1pe12.4,"Watts")



      ! PRINT OUT OF DEPOSITED POWER (WATTS/CC):
      ! powrft(l), powrf(l,1),...,powrf(l,5)
#ifndef NOPGPLOT
      CALL PGPAGE

      CALL PGSVP(.05,.95,.05,.95)
#endif
      RILIN=0.
#ifndef NOPGPLOT
      CALL PGSAVE

      CALL PGSCH(1.3) !PGSCH(1.44) ! set character size; default is 1.

      CALL PGMTXT('T',-RILIN,.5,.5, &
           'DEPOSITED POWER: (WATTS/CC)')
#endif
      RILIN=RILIN+3.
#ifndef NOPGPLOT
      CALL PGSCH(1.) ! set character size; default is 1.
#endif
      write(t_,6113) pltvs
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+2.
      write(t_,6114) pltvs
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+2.

 6113 format(1x,a8,"TOTAL", &
        "       RF1         RF2        RF3         RF4        RF5")
 6114 format(1x,a8,"(powrft)", &
        "           (powrf(*,harmonic) for harmonics = 1-5)")

!     Start printing results on first page
      do l=1,min(40,setup0%lrzmax)
        if (mrfn.eq.1) then
          write(t_,6115) tr(l),powrft(l),powrf(l,1)
        elseif (mrfn.eq.2) then
          write(t_,6116) tr(l),powrft(l),powrf(l,1),powrf(l,2)
        elseif (mrfn.eq.3) then
          write(t_,6117) tr(l),powrft(l),powrf(l,1),powrf(l,2), &
                               powrf(l,3)
        elseif (mrfn.eq.4) then
          write(t_,6118) tr(l),powrft(l),powrf(l,1),powrf(l,2), &
                               powrf(l,3),powrf(l,4)
        else
          write(t_,6119) tr(l),powrft(l),powrf(l,1),powrf(l,2), &
                               powrf(l,3),powrf(l,4),powrf(l,5)
        endif
#ifndef NOPGPLOT
        CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
        RILIN=RILIN+1.
      enddo

!     Continue printing results on second page, if setup0%lrzmax.gt.40
      if (setup0%lrzmax.gt.40) then
#ifndef NOPGPLOT
         CALL PGPAGE

         CALL PGSVP(.05,.95,.05,.95)
#endif
         RILIN=0.+3.
         do l=41,setup0%lrzmax
           if (mrfn.eq.1) then
             write(t_,6115) tr(l),powrft(l),powrf(l,1)
           elseif (mrfn.eq.2) then
             write(t_,6116) tr(l),powrft(l),powrf(l,1),powrf(l,2)
           elseif (mrfn.eq.3) then
             write(t_,6117) tr(l),powrft(l),powrf(l,1),powrf(l,2), &
                                  powrf(l,3)
           elseif (mrfn.eq.4) then
             write(t_,6118) tr(l),powrft(l),powrf(l,1),powrf(l,2), &
                                  powrf(l,3),powrf(l,4)
           else
             write(t_,6119) tr(l),powrft(l),powrf(l,1),powrf(l,2), &
                                  powrf(l,3),powrf(l,4),powrf(l,5)
           endif
#ifndef NOPGPLOT
           CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
           RILIN=RILIN+1.
         enddo
      endif

 6115 format(f5.3, 2(1x,e9.2) )
 6116 format(f5.3, 3(1x,e9.2) )
 6117 format(f5.3, 4(1x,e9.2) )
 6118 format(f5.3, 5(1x,e9.2) )
 6119 format(f5.3, 6(1x,e9.2) )

      RILIN=RILIN+2.
      write(t_,6122) sorpwtza
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+1.
      write(t_,6123) powurf(0)
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+1.

      do krf=1,mrfn ! Print-out for all harmonics now (was 5 only)
         write(t_,6131) krf,nharm(krf),powurf(krf)
#ifndef NOPGPLOT
         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
         RILIN=RILIN+1.
      enddo
 6131 format("      mode/harmonic krf, nharm(krf), powurf(krf)=", &
                    2i4,1pe12.4)

!      if (mrfn.ge.1) then
!         write(t_,6124) powurf(1)
#ifndef NOPGPLOT
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
!         RILIN=RILIN+1.
!      endif
!      if (mrfn.ge.2) then
!         write(t_,6125) powurf(2)
#ifndef NOPGPLOT
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
!         RILIN=RILIN+1.
!      endif
!      if (mrfn.ge.3) then
!         write(t_,6126) powurf(3)
#ifndef NOPGPLOT
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
!         RILIN=RILIN+1.
!      endif
!      if (mrfn.ge.4) then
!         write(t_,6127) powurf(4)
#ifndef NOPGPLOT
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
!         RILIN=RILIN+1.
!      endif
!      if (mrfn.ge.5) then
!         write(t_,6128) powurf(5)
#ifndef NOPGPLOT
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
!         RILIN=RILIN+1.
!      endif

      write(t_,6129) powurfc(0)
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif
      RILIN=RILIN+1.
      write(t_,6130) powurfl(0)
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif

#ifndef NOPGPLOT
      CALL PGUNSA
#endif

 6122 format &
       ("Power sources integr. over rad. (RF+NBI, all gen.species)=", &
             1pe12.4,"W")
 6123 format("Power from intern ray diagnostic[powurf(0)]=",1pe12.4,"W")
 6124 format("                mode/harmonic 1 [powurf(1)]=",1pe12.4)
 6125 format("                mode/harmonic 2 [powurf(2)]=",1pe12.4)
 6126 format("                mode/harmonic 3 [powurf(3)]=",1pe12.4)
 6127 format("                mode/harmonic 4 [powurf(4)]=",1pe12.4)
 6128 format("                mode/harmonic 5 [powurf(5)]=",1pe12.4)
 6129 format("Power by collisions (from ray data)    =",1pe12.4,"W")
 6130 format("Power by linear damping (from ray data)=",1pe12.4,"W")


!..................................................................
!     plots of sources and deposited power density [W/cm^3]..
!..................................................................

      fmin=0.
      fmax=0.
      gmin=0.
      gmax=0.
      call aminmx(powrft(1:setup0%lrzmax),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
      call aminmx(sorpwt(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
      if (gmin.lt.fmin) fmin=gmin
      if (gmax.gt.fmax) fmax=gmax
      !if (fmax-fmin.lt.1.e-16) fmin=fmax-.1*abs(fmax)-1.e-5
      if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range
      RPG1=fmin ! could be negative because of numerical errors
      RPG2=fmax
      RPG1=min(fmin,0.) !Make lower limit 0 when fmin>0

      do l=1,setup0%lrzmax
        tr1(l)=powrft(l)
      enddo

#ifndef NOPGPLOT
      CALL PGPAGE ! START FSA SOURCE POWER DEN

      CALL PGSVP(.2,.8,.2,.6)

      CALL PGSAVE

      CALL PGSCH(1.44) ! set character size; default is 1.

      CALL PGSLW(setup0%lnwidth) ! line thickness/width

      CALL PGMTXT('T',6.,0.5,0.5, &
                    'FSA SOURCE POWER DEN: (WATTS/CM\u3\d)')
      CALL PGMTXT('T',5.,0.,0., &
         "Solid line: NBI+RF for all gen.species [sorpwt]")
      CALL PGMTXT('T',4.,0.,0., &
         "Dashed: NBI (beam ions) [sorpw_nbi]")
      CALL PGMTXT('T',3.,0.,0., &
         "Solid-bold: total absorbed RF power [powrft]")
      CALL PGMTXT('T',2.,0.,0., &
         "Other: RF general species (each) [sorpw_rf]")
      CALL PGUNSA
#endif

        DO I=1,setup0%lrzmax
           RLRZAP1(I)=tr(i)
        ENDDO
      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRZAP1(1)
      RPGmax=RLRZAP1(setup0%lrzmax)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.
      IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
         RPG2= RPG1+1.e-16
      ENDIF
#ifndef NOPGPLOT
      CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

        CALL PGSLS(1) ! 1-> solid
#endif
        DO I=1,setup0%lrzmax
           RLRZAP11(I)=sorpwt(i) ! solid: NBI+RF(all gen.species)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))

        CALL PGSLS(2) ! 2-> dashed
#endif
        DO I=1,setup0%lrzmax
           RLRZAP12(I)=sorpw_nbi(kfrsou1,I) ! dashed: NBI only
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif
      do k=1,ngen ! rf sources for general species
         DO I=1,setup0%lrzmax
            RLRZAP12(I)=sorpw_rf(k,I)
         ENDDO
#ifndef NOPGPLOT
         CALL PGSLS(k+2) ! 3-> -.-.- ;   4-> .....

         CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif
      enddo ! k=1,ngen
      !
#ifndef NOPGPLOT
      CALL PGSLS(1) ! solid

      CALL PGSLW(setup0%lnwidth+1) ! bold
#endif
      DO I=1,setup0%lrzmax
         RLRZAP13(I)=powrft(i)
      ENDDO
#ifndef NOPGPLOT
      CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax)) !solid bold: total rf

      CALL PGSLS(1)

      CALL PGSLW(setup0%lnwidth) !

      CALL PGSAVE

      CALL PGSCH(1.44)

      CALL PGLAB(' ','power density (W/cm\u3\d)',' ')

      CALL PGMTXT('B',1.8,0.5,0.5,t_horiz)

      CALL PGUNSA
#endif
      ! DONE FSA SOURCE POWER


      !YuP[08-2017] Make a separate page for only RF-power density
      !(because, if the power level is too small, the curve can be too low)
#ifndef NOPGPLOT
      CALL PGPAGE ! START RF-only POWER DEN

      CALL PGSVP(.2,.8,.2,.6)

      CALL PGSAVE

      CALL PGSCH(1.44) ! set character size; default is 1.

      CALL PGSLW(setup0%lnwidth) ! line thickness/width

      CALL PGMTXT('T',6.,0.5,0.5, &
                    'FSA RF POWER DEN: (WATTS/CM\u3\d)')
      CALL PGMTXT('T',3.,0.,0., &
         "Solid-bold: total absorbed RF power [powrft]")
      CALL PGMTXT('T',2.,0.,0., &
         "Other: RF general species (each) [sorpw_rf]")
      CALL PGUNSA
#endif
      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRZAP1(1)
      RPGmax=RLRZAP1(setup0%lrzmax)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.
      ! Vertical axis limits:
      call aminmx(powrft(1:setup0%lrzmax),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
      RPG1=gmin ! could be negative because of numerical errors
      RPG2=gmax*1.2 ! give 20% extra
      RPG1=min(gmin,0.) !Make lower limit 0 when gmin>0
      IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
         RPG2= RPG1+1.e-16
      ENDIF
#ifndef NOPGPLOT
      CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)

      CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
#endif
      do k=1,ngen ! rf sources for general species
         DO I=1,setup0%lrzmax
            RLRZAP12(I)=sorpw_rf(k,I)
         ENDDO
#ifndef NOPGPLOT
         CALL PGSLS(k+2) ! 3-> -.-.- ;   4-> .....

         CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif
      enddo ! k=1,ngen
      !
#ifndef NOPGPLOT
      CALL PGSLS(1) ! solid

      CALL PGSLW(setup0%lnwidth+1) ! bold
#endif
      DO I=1,setup0%lrzmax
         RLRZAP13(I)=powrft(i)
      ENDDO
#ifndef NOPGPLOT
      CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax)) !solid bold: total rf

      CALL PGSLS(1)

      CALL PGSLW(setup0%lnwidth) !

      CALL PGSAVE

      CALL PGSCH(1.44)

      CALL PGLAB(' ','power density (W/cm\u3\d)',' ')

      CALL PGMTXT('B',1.8,0.5,0.5,t_horiz)

      CALL PGUNSA
#endif
      ! DONE RF-only POWER


!..................................................................
!     plots of partial integration of powers [Watts]
!..................................................................

      fmin=0.
      fmax=0.
      call aminmx(sorpwti(1:setup0%lrzmax),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
      if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range

      do 50 kk=1,mrfn
        call aminmx(powurfi(1:setup0%lrzmax,kk),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
        if (gmin.lt.fmin) fmin=gmin
        if (fmax.lt.gmax) fmax=gmax
 50   continue
      if (fmax-fmin.lt.1.e-8) fmin=fmax-.1*abs(fmax)-1.e-5
      RPG1=fmin
      RPG2=fmax
      RPG1=min(fmin,0.) !Make lower limit 0 when fmin>0

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)

        CALL PGSAVE

        CALL PGSCH(1.44)

        CALL PGMTXT('T',6.,0.5,0.5, &
                    'SOURCE POWER (integr. up to rho or psi) (WATTS)')
        CALL PGMTXT('T',5.,0.,0., &
         "Solid line: NBI+RF for all gen.species [sorpwti]")
        CALL PGMTXT('T',4.,0.,0., &
         "Dashed: NBI (beam ions) [sorpw_nbii]")
        CALL PGMTXT('T',3.,0.,0., &
         "Solid-bold: total absorbed RF [powurfi(*,0)]")
        CALL PGMTXT('T',2.,0.,0., &
         "Other: RF general ions (each) [sorpw_rfi]")
        CALL PGUNSA
#endif

        DO I=1,setup0%lrzmax
           RLRZAP1(I)=tr(i)
        ENDDO

      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRZAP1(1)
      RPGmax=RLRZAP1(setup0%lrzmax)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.

        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
#ifndef NOPGPLOT
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)

        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

        CALL PGSLS(1) ! 1-> Solid line
#endif
        DO I=1,setup0%lrzmax
           RLRZAP11(I)=sorpwti(i) ! solid: NBI+RF(all gen.species)
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP11(1:setup0%lrzmax))

        CALL PGSLS(2) ! 2-> dashed
#endif
        DO I=1,setup0%lrzmax
           RLRZAP12(I)=sorpw_nbii(kfrsou1,I) ! dashed: NBI only
        ENDDO
#ifndef NOPGPLOT
        CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif
      do k=1,ngen ! rf sources for general species
         DO I=1,setup0%lrzmax
            RLRZAP12(I)=sorpw_rfi(k,I)
         ENDDO
#ifndef NOPGPLOT
         CALL PGSLS(k+2) ! 3-> -.-.- ;   4-> .....

         CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif
      enddo ! k=1,ngen
#ifndef NOPGPLOT
      CALL PGSLS(1) ! solid

      CALL PGSLW(setup0%lnwidth+1) ! bold
#endif
      DO I=1,setup0%lrzmax
         RLRZAP13(I)=powurfi(i,0) ! = SUM_harmonics{powurfi(l,harmonics)}
      ENDDO
#ifndef NOPGPLOT
      CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax)) !solid bold: total rf

      CALL PGSLS(1) ! restore: solid line

      CALL PGSLW(setup0%lnwidth) ! restore

      CALL PGSAVE

      CALL PGSCH(1.44)

      CALL PGLAB(' ','Int power: 0 to psi (or rho)',' ')

      CALL PGMTXT('B',1.8,0.5,0.5,t_horiz)

      CALL PGUNSA ! restore
#endif


!..................................................................
!     plots of partial integration of RF-only powers [Watts]
      !YuP[08-2017] Make a separate page for RF-power only
      !(because, if the power level is too small, the curve can be too low)
!..................................................................

      call aminmx(powurfi(1:setup0%lrzmax,0),1,setup0%lrzmax,1,gmin,gmax,kmin,kmax)
      if (gmax-gmin.lt.1.e-8) gmin=gmax-.1*abs(gmax)-1.e-5
      RPG1=gmin
      RPG2=gmax*1.2 ! give 20% more
      RPG1=min(gmin,0.) !Make lower limit 0 when gmin>0

#ifndef NOPGPLOT
        CALL PGPAGE

        CALL PGSVP(.2,.8,.2,.6)

        CALL PGSAVE

        CALL PGSCH(1.3)

        CALL PGMTXT('T',6.,0.5,0.5, &
                    'RF POWER (integr. up to rho or psi) (WATTS)')
        CALL PGMTXT('T',3.,0.,0., &
         "Solid-bold: total absorbed RF [powurfi(*,0)]")
        CALL PGMTXT('T',2.,0.,0., &
             "Other: RF general species (each) [sorpw_rfi]")
        CALL PGUNSA
#endif

      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRZAP1(1)
      RPGmax=RLRZAP1(setup0%lrzmax)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.

      IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
          RPG2= RPG1+1.e-16
      ENDIF
#ifndef NOPGPLOT
      CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)

      CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)

      CALL PGSLS(1) ! 1-> Solid line
#endif
      !
      do k=1,ngen ! rf sources for general species
         DO I=1,setup0%lrzmax
            RLRZAP12(I)=sorpw_rfi(k,I)
         ENDDO
#ifndef NOPGPLOT
         CALL PGSLS(k+2) ! 3-> -.-.- ;   4-> .....

         CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP12(1:setup0%lrzmax))
#endif
      enddo ! k=1,ngen
#ifndef NOPGPLOT
      CALL PGSLS(1) ! solid

      CALL PGSLW(setup0%lnwidth+1) ! bold
#endif
      DO I=1,setup0%lrzmax
         RLRZAP13(I)=powurfi(i,0) ! = SUM_harmonics{powurfi(l,harmonics)}
      ENDDO
#ifndef NOPGPLOT
      CALL PGLINE(setup0%lrzmax,RLRZAP1(1:setup0%lrzmax),RLRZAP13(1:setup0%lrzmax)) !solid bold: total rf

      CALL PGSLS(1) ! restore: solid line

      CALL PGSLW(setup0%lnwidth) ! restore

      CALL PGSAVE

      CALL PGSCH(1.44)

      CALL PGLAB(' ','Int power: 0 to psi (or rho)',' ')

      CALL PGMTXT('B',1.8,0.5,0.5,t_horiz)

      CALL PGUNSA ! restore
#endif

!..................................................................
!     J/P print out
!..................................................................

!$$$      call gxglfr(0)
!$$$      call gscpvs(.5,.99)
!$$$      call gstxjf("center","top")
!$$$      call gstxno(40.)
!$$$      call gitxft(ift)
!$$$      if (ift .ne. ift07a) call gstxft(ift07a)
!$$$      call gptx2d("J/P-(AMPS/WATT);$")
!$$$      call gitxft(ift)
!$$$      if (ift.ne.ift21a) call gstxft(ift21a)
!$$$      call gstxno(100.)
!$$$      call gstxjf("center","top")
!$$$      call gscpvs(.5,.93)
!$$$      write(t_,8010) pltvs
!$$$      call gptx2d(t_)
!$$$ 8010 format(4x,a8,2x," J/P (fi)   J/P (fi+e)","$")
!$$$      do 801 l=1,setup0%lrzmax
!$$$        write(t_,8011) tr(l),bdre(l),bdrep(l)
!$$$        call gptx2d(t_)
!$$$ 801  continue
!$$$ 8011 format((1x,1pe10.3,2x,e10.3,2x,e10.3),"$")
!%OS  call gxglfr(0)
!%OS  call gscpvs(.5,.95)
!%OS  call gstxjf("center","top")
!%OS  call gstxno(40.)
!%OS  call gitxft(ift)
!%OS  if (ift .ne. ift07a) call gstxft(ift07a)
!%OS  call gptx2d("J/P - (AMPS/WATT);$")
!%OS  call gitxft(ift)
!%OS  if (ift.ne.ift21a) call gstxft(ift21a)
!%OS  call gstxno(100.)
!%OS  call gscvlb(1)
!$$$      fmin=0.
!$$$      fmax=0.
!$$$      call aminmx(bdrep(1),1,setup0%lrzmax,1,fmin,fmax,kmin,kmax)
!$$$      call aminmx(bdre(1),1,setup0%lrzmax,1,fminn,fmaxx,kmin,kmax)
!$$$      if (fminn.lt.fmin) fmin=fminn
!$$$      if (fmaxx.gt.fmax) fmax=fmaxx
!$$$      if (fmax-fmin.lt.1.e-8) fmin=fmax-.1*abs(fmax)-1.e-5
!$$$      call gswd2d("linlin$",tr(1),tr(setup0%lrzmax),fmin,fmax)
!$$$      call gsvp2d(.2,.9,.15,.6)
!$$$      call gscvft(0.)
!$$$      text(1)="fi+e$"
!$$$      call gscvtx(loc(text))
!$$$      call gpgr80("linlin$")
!$$$      call gpcv2d(tr(1),bdrep(1),setup0%lrzmax)
!$$$      call gscvft(.3)
!$$$      text(1)="fi$"
!$$$      call gscvtx(loc(text))
!$$$      call gpcv2d(tr(1),bdre,setup0%lrzmax)
!$$$      call gstxan(90.)
!$$$      call gscpvs(.01,.375)
!$$$      call gscvlb(0)
!$$$      write(t_,8080)
!$$$      call gptx2d(t_)
!$$$ 8080 format("fi - J/P for fast ions - NBI",";",
!$$$     1  "fi+e - J/P (fast ions+electrons - NBI+RF)","$")
!$$$      call gstxan(0.)
!$$$      call gscpvs(.57,.1)
!$$$      call gptx2d(t_)

!..................................................................
!     Print out figures af merit and other average quantities.
!..................................................................

 809  CONTINUE

      return
      end subroutine tdpltjop

end module tdpltjop_mod
