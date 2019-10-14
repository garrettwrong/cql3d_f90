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

!
!
module eqrhopsi_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use eqflxavg_mod, only : eqflxavg
  use eqfpsi_mod, only : eqfpsi
  use eqonovrp_mod, only : eqonovrp
  use eqorbit_mod, only : eqorbit
  use eqvolpsi_mod, only : eqvolpsi
  use eqwrng_mod, only : eqwrng
  use exlin_mod, only : exlin
  use zcunix_mod, only : coeff1
  use zcunix_mod, only : coeff2
  use zcunix_mod, only : terp1
  use zcunix_mod, only : terp2

  !---END USE

  use param_mod
  use cqlcomm_mod
  use iso_c_binding, only : c_double
#ifdef __MPI
  include 'cql3d_mpilib.h'
#endif
  real(c_double), private :: btor00,bthr00,bmod00
  save

contains

  subroutine eqrhopsi(generate)
    use aminmx_mod, only : aminmx
    implicit integer (i-n), real(c_double) (a-h,o-z)
    character*8 generate

    parameter(nworka=3*nconteqa+1)
    dimension workk(nworka)

    !..................................................................
    !     This routine determines an array eqrho(eqpsi(k)) where
    !     k=1,...,nconteqn, eqpsi is a value of psi (equilibrium
    !     psi_poloidal) and eqrho is the associated  radial coordinate
    !     (sqrt(toroidal flux),sqrt(area), sqrt(volume),
    !     radial extent of flux surface, normalized poloidal flux,
    !     or normalized sqrt(pol flux), according to the value
    !     of radcoord.  See cqlinput_help.)
    !     The eqpsi array can be provided externally
    !     (generate="disabled") or it can be generated within
    !     this routine (generate="enabled"). [It is called with
    !     generate="enabled" from eqcoord when lr_.eq.setup0%lrzmax.]
    !     Also, determines some approximate geometric quantities
    !     of the eqdsk.
    !..................................................................

    !......................................................................
    !     General, non-circular flux-surface psi values, epsi, have been
    !     set up from eqdsk, topeol, or elliptical plasma models, depending
    !     on eqmod.ne.'disabled' and eqsource.  epsi, at this point is an
    !     decreasing function from the magnetic axis to the plasam edge.
    !     For eqmod.ne.'disabled', namelist eqsource gives the source of
    !     equilibrium flux surface data, giving flux (epsi) on a regular
    !     cylindrical coordinate R,Z-grid.  Also provided is f=R*B_phi.
    !     If eqsource="ellipse", then epsi, R, Z etc are determined
    !     from a subroutine eqelpse, and f is determined arbitrarily from
    !     an input parameter fpsimodl.
    !     In order to solve the orbit equations it will be necessary to
    !     know all the derivatives of epsi up to second order. A 2-D
    !     spline package will be used. Set up the NCAR spline package.
    !     Ordinarily the pol. flux psi value associated with the flux
    !     surface of interest will be labeled by common variable epsicon_.
    !
    !     091016: Setup of 2D splines and the rmaz,zmag re-calculation
    !             below has been moved up from subroutine eqorbit,
    !             to facilitate de-updown-symmetrization.
    !
    !     Calls to the current subroutine eqrhospi are from eqcoord, in
    !     descending order lr_=setup0%lrzmax,1,-1.
    !......................................................................


    !      if (eqcall.eq."enabled") then   !Called on first run through.
    if (lr_.eq.setup0%lrzmax) then   !Called on first run through.
       ibd(1)=4
       ibd(2)=4
       ibd(3)=4
       ibd(4)=4

       !       2D bicubic spline coeffs of poloidal flux. epsi has sign change
       !       after reading eqdsk, so psimag is max psi, psilim is less.
       call coeff2(nnr,er,nnz,ez,epsi,epsirr,epsizz,epsirz,nnra,ibd, &
            wkepsi)
       !..................................................................
       !     Now determine the exact location of the magnetic axis. From
       !     The point of view of the spline package it is somewhere
       !     between er(imag-1) and er(imag+1). Use Newton iteration.
       !     imag was determined in subroutine eqrhopsi.
       !..................................................................

       !        rmag=er(imag)
       !        zmag=0.
       !        do 5 iter=1,8
       !          dpsidr=terp2(rmag,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz,
       !     1      epsirz,nnra,1,0)
       !          d2psidrr=terp2(rmag,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz,
       !     1      epsirz,nnra,2,0)
       !          rmag=rmag-dpsidr/d2psidrr
       ! 5      continue
       !        psimag=terp2(rmag,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz,
       !     1    epsirz,nnra,0,0)
       !        eqpsi(1)=psimag
       !      endif


       !     Simple search for more accurate rmag, based in bi-cubic spline
       !     of the equilibrium data.  Above simple Newton-Raphson failed
       !     for high beta (shifted flux surface) equilibria (bobh, 960309).
       !     [BH091017: Maybe could be cured by better starting imag derived
       !     from the input equilibrium raxis,zaxis, as below?]
       !     It is assumed here that the magnetic axis is at a maximum
       !     of psi.   Up-down symmetry (zmag=0.) assumed.
       !
       !     BH091018:  Should touchup following calc of the magnetic
       !                axis with a Newton-Raphson iteration.  Inaccurate
       !                magnetic axis location probably limits minimum
       !                rya() which can be successfully used.
       rmag_old=rmag ! From eqdsk
       zmag_old=zmag
       psimag_old=psimag

       if (eqsym.ne."none") then  !i.e., assuming mag axis at 0.
          !BH091017           imag=nnr/2  !nnr set, e.g., by read of eqdsk
          drr=er(2)-er(1)
          dzz=ez(2)-ez(1)
          imag=nint((rmag-er(1))/drr)+1  !nint() is nearest integer
          jmag=nint((zmag-ez(1))/dzz)+1
          if(imag.gt.1)then
             er1=er(imag-1)
             er2=er(imag+1)
             rmag=er(imag-1)
          else ! imag=1 (could be for eqsource="mirror1"
             er1=er(1)
             er2=er(2)
             rmag=er(1)
          endif
          zmag=0.
          psi1=terp2(rmag,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
               epsirz,nnra,0,0)
          npoints=101
          dr=(er2-er1)/(npoints-1)
          do i=1,npoints
             erlocal=er1+(i-1)*dr
             psilocal=terp2(erlocal,zmag,nnr,er,nnz,ez,epsi,epsirr, &
                  epsizz,epsirz,nnra,0,0)
             if (psilocal.gt.psi1) then
                psi1=psilocal
                rmag=erlocal
             endif
          enddo
          psimag=terp2(rmag,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
               epsirz,nnra,0,0)
          eqpsi(1)=psimag

       else  !i.e.  No symmetrization: Uses entire equilibrium
          !      Perform 2D search for magnetic axis
          drr=er(2)-er(1)
          dzz=ez(2)-ez(1)
          imag=nint((rmag-er(1))/drr)+1  !nint() is nearest integer
          jmag=nint((zmag-ez(1))/dzz)+1
          !
          !        write(*,*)'eqrhopsi: er(imag+-1)',er(imag-1),er(imag),er(imag+1)
          !        write(*,*)'eqrhopsi: ez(jmag+-1)',ez(jmag-1),ez(jmag),ez(jmag+1)
          !        write(*,*)'epsi',
          !     +    epsi(imag-1,jmag-1),epsi(imag,jmag-1),epsi(imag+1,jmag-1)
          !        write(*,*)'epsi',
          !     +    epsi(imag-1,jmag),epsi(imag,jmag),epsi(imag+1,jmag)
          !        write(*,*)'epsi',
          !     +    epsi(imag-1,jmag+1),epsi(imag,jmag+1),epsi(imag+1,jmag+1)
          !
          if(imag.gt.1)then
             er1=er(imag-1)
             er2=er(imag+1)
             rmag=er(imag-1)
          else ! imag=1 (could be for eqsource="mirror1"
             er1=er(1)
             er2=er(2)
             rmag=er(1)
          endif
          zmag=ez(jmag-1)
          ez1=ez(jmag-1)
          ez2=ez(jmag+1)  ! YuP[2015/05/03] Bug fix: was imag+1
          psi1=terp2(rmag,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
               epsirz,nnra,0,0)
          npoints=101
          drrr=(er2-er1)/(npoints-1)
          dzzz=(ez2-ez1)/(npoints-1)
          do i=1,npoints
             do j=1,npoints
                erlocal=er1+(i-1)*drrr
                ezlocal=ez1+(j-1)*dzzz
                psilocal=terp2(erlocal,ezlocal,nnr,er,nnz,ez,epsi,epsirr, &
                     epsizz,epsirz,nnra,0,0)
                if (psilocal.gt.psi1) then
                   psi1=psilocal
                   rmag=erlocal
                   zmag=ezlocal
                endif
             enddo
          enddo
          psimag=terp2(rmag,zmag,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
               epsirz,nnra,0,0)
          eqpsi(1)=psimag
       endif  !On eqsym

       write(*,*)'eqrhopsi: psilim=',psilim
       !write(*,*)'eqrhopsi: ez(jmag+-1)',ez(jmag-1),ez(jmag),ez(jmag+1)

       write(*,'(a,i6,2e13.4)')'eqrhopsi: imag, rmag_old/new =', &
            imag, rmag_old, rmag
       write(*,'(a,i6,2e13.4)')'eqrhopsi: jmag, zmag_old/new =', &
            jmag, zmag_old, zmag
       write(*,'(a,2e13.4)')'eqrhopsi: psimag_old/new', &
            psimag_old, psimag
       write(*,'(a,e13.4)')'eqrhopsi: epsi(imag,jmag)=',epsi(imag,jmag)
       write(*,'(a,e13.4)')'eqrhopsi: er(imag)=',er(imag)
       write(*,'(a,e13.4)')'eqrhopsi: ez(jmag)=',ez(jmag)

       !      endif  !On eqcall
    endif   !On lr_.eq..rzmax



    !..................................................................
    !     This routine  is called with generate="enabled" from eqcoord
    !     when lr_.eq.setup0%lrzmax.  radcoord is namelist input indicating
    !     the specific definition of the radial coordinate.
    !..................................................................

      write(*,*)'eqrhopsi: generate,radcoord: ',generate,radcoord

    if (generate.eq."enabled") then     ! generate eqpsi
       !        nzc=(nnz-1)/2+1                   ! up-down symmetry assumed
       nzc=jmag
       if(eqsource.eq."ellipse") psilim=0.
       !        nrc=nnr/2
       nrc=imag
       amax=epsi(nrc,nzc) ! =psimag
       imag=nrc ! =1 in mirror machine
       !       Find jpsimnr such that psi~psilim
       do 1 j=nrc,nnr ! going from imag to larger R, along Z=zmag
          if (epsi(j,nzc).lt.psilim .and. amax.gt.psilim) then
             iminval=j-1
             eqpmn1=epsi(iminval,nzc)
             go to 2 ! done, found: epsi(j,nzc) is outside of psilim
          endif     !      while epsi(j-1,nzc) was inside of psilim
          if (epsi(j,nzc).gt.amax) then
             amax=epsi(j,nzc)
             imag=j
          endif
1     continue
      iminval=nnr
      eqpmn1=epsi(nnr,nzc)
2     continue
      jpsimnr=iminval + 1 - 1/(nnr+1-iminval)
      write(*,*)'eqrhopsi: iminval,jpsimnr=',iminval,jpsimnr
      !
      !       Find jpsimnl such that psi~psilim (at the inboard)
      jpsimnl= jpsimnr ! will be over-written in  a tokamak machine
      !(but in a mirror machine, the loop is skipped because nrc=1)
      do 3 j=nrc-1,1,-1 !going from er(imag-1) to inner smaller er
         jpsimnl=j
         if(epsi(j,nzc).lt.psilim) go to 4
         if (epsi(j,nzc).gt.amax) then
            amax=epsi(j,nzc)
            imag=j
         endif
3     continue
4     continue
      eqpsimin=eqpmn1
      eqpsimax=amax

      !..................................................................
      !     The eqpsi array will be chosen if possible so that it is spaced
      !     between the value of psi at the magnetic axis and the value of
      !     psi for a flux surface that passes through:
      !     the point (R=rmag,Z=radmin) if eqsource="ellipse"
      !     psilim if eqsource = "topeol" or "eqdsk"
      !..................................................................

      if (eqsource.eq."ellipse") then
         if (ez(nnz).lt.radmin) call eqwrng(9)
         do 5 i=nzc,nnz
            if (ez(i).ge.radmin) go to 6
5        continue
6        continue
         ival=i
         z1=ez(ival-1)
         z2=ez(ival)
         psi1=epsi(imag,ival-1)
         psi2=epsi(imag,ival)
         eqpmn2=(psi1-psi2)/(z1-z2)*(radmin-z1)+psi1
         eqpsimin=eqpmn2
         if (eqpmn1.gt.eqpmn2) eqpsimin=eqpmn1
         !*bh*931222elseif(eqsource.eq."eqdsk".or.eqsource.eq."tsc") then
      else if (eqsource.eq."eqdsk" .or. eqsource.eq."tsc") then
         !..................................................................
         !     If an eqdsk equilibrium or a ONETWO equilibrium is being used
         !     we already know psi at the limiter.
         !..................................................................
         eqpsimin=psilim
         eqpsimax=psimag
      else if (eqsource.eq."topeol") then
         !         Prior to moving calc of 2D spline of epsi to beginning of
         !         this subroutine, following call appears to be a bug (BH091016).
         !         [There must have been some prior rearrangement.]
         write(*,*)'eqrhopsi: rmaxcon,zmag,rmag =',rmaxcon,zmag,rmag
         psilim=terp2  (rmaxcon,zmag,nnr,er,nnz,ez,epsi,epsirr, &
              epsizz,epsirz,nnra,0,0)
         eqpsimin=psilim
         psimag=terp2(rmag   ,zmag,nnr,er,nnz,ez,epsi,epsirr, &
              epsizz,epsirz,nnra,0,0)
         eqpsimax=psimag
      endif

      !..................................................................
      !     Determine the array of psi values..
      !     (Positioned on the epsi grid if nconteq.eq."psigrid",
      !     equally spaced in psi if nconteq.ne."psigrid" and
      !     nconteqn is a positive number)
      !..................................................................
      eqpsi=eqpsimin !initialize (not for all indexes eqpsi is set below)
      if (nconteq .ne. "psigrid") then
         !         Equispaced from Mag axis to psilim
         !         (Could be problem in eqorbit if psilim surface has separatrix.
         !          Could add a psifactr here, or improve eqorbit LCFS calc. BH).
         !          Alternatively, we pull back eqpsi(nconteqn) slightly, as in
         !          the nconteq.eq."psigrid" case.
         delpsi=(eqpsimax-eqpsimin)/(nconteqn-1)
         do 10 j=1,nconteqn
            eqpsi(j)=-(j-1)*delpsi+eqpsimax !from eqpsimax down to eqpsimin
10       continue
            ! Adjust the last point:
            !          eqpsi(nconteqn)=eqpsi(nconteqn-1)+
            !     +                    0.9*(eqpsi(nconteqn)-eqpsi(nconteqn-1))
            !YuP: which is same as
         eqpsi(nconteqn)=eqpsi(nconteqn-1)-0.9*delpsi
         !So, instead of reducing by the whole delpsi,
         ! we reduce by 0.9*delpsi, to stay away from eqpsimin
      else ! nconteq = "psigrid"
         !         This calc of eqpsi values chooses eqpsi(nconteqn) slightly
         !         inside the psilim flux surface.
         nconteqn=0
         !for a mirror machine imag=1, so the whole range 1:nnr can be used
         do 11 j=imag,nnr
            ! as soon as epsi(j,nzc) got outside of psilim, then:
            if (epsi(j,nzc).le. eqpsimin) then
               !YuP[03-2016]  changed to .le. (was .lt.) see notes below.
               nconteqn=nconteqn+1
               ! !!eqpsi(nconteqn)=eqpsimin+(eqpsi(nconteqn-1)-eqpsimin)*.1
               !YuP[03-2016] The above line appears to have a flaw:
               !(Note: eqpsimin=psilim, and epsi(j,nzc) and eqpsi()
               ! are descending with j).
               !It can happen that at previous step,
               ! eqpsi(nconteqn-1) was exactly equal to eqpsimin.
               ! In this case, at present step, eqpsi(nconteqn)=eqpsimin
               ! and so we got two points with same value.
               !Correction made: changed .lt. to .le.,
               ! and the definition for the last point:
               if(     (eqpsi(nconteqn-1)-eqpsimin)  .gt. &
                    0.5*(eqpsi(nconteqn-2)-eqpsi(nconteqn-1)) )then
                  !The previous point was far enough from psilim;
                  !Form the new (last) point, slightly inside psilim :
                  del_psi= eqpsi(nconteqn-2)-eqpsi(nconteqn-1)
                  !Note: del_psi is positive (because eqpsi is descending)
                  ![this was corrected on 2017-12-05;
                  !the wrong sign got from 03-2016/cql3d-mirror-version]
                  eqpsi(nconteqn)= eqpsimin+del_psi*0.1
               else
                  !The previous point too close from psilim;
                  !Do not form the new point;
                  !consider the previous points as the last :
                  nconteqn=nconteqn-1 ! un-do the increment done above.
               endif
               !              write(*,*)'eqrhopsi: nconteqn,eqpsi(nconteqn-1),eqpsimin',
               !     +                             nconteqn,eqpsi(nconteqn-1),eqpsimin
               go to 12 ! done, quit the loop
            endif
            eqpsi(j-imag+1)=epsi(j,nzc)
            !So, epsi(imag,nzc)   --> eqpsi(1)
            !    epsi(imag+1,nzc) --> eqpsi(2)
            !    epsi(imag+2,nzc) --> eqpsi(3)
            !    ..... until epsi(j,nzc) .le. eqpsimin
            !         (got onto or outside of LCFS)
            !          Then, eqpsi(nconteqn)= eqpsimin+del_psi*0.1
            nconteqn=j-imag+1
 11      continue
 12      continue
         endif ! nconteq .ne. "psigrid"

         ! Check that eqpsi is strictly descending:
         do j=2,nconteqn
#ifdef __MPI
         ! for MPI, only print from master
         if(mpirank.eq.0) then
#endif
            write(*,'(a,i4,e16.8)')' j,eqpsi(j)-psilim=',j,eqpsi(j)-psilim
#ifdef __MPI
         endif !mpirank 0
#endif
            if(eqpsi(j)-eqpsi(j-1).ge.0.d0)then
               WRITE(*,*)'eqrhopsi: eqpsi is not descending for j,j-1=',j,j-1
               stop
            endif
         enddo

    endif  ! On generate

!..................................................................
!     Determine the volume,... of the flux surface.
!..................................................................

    eqvol(1)=0.
    eqarea(1)=0.
    do 20 j=2,nconteqn
       if (j.eq.2) then
          if (eqsource.ne."ellipse") then
             eqcall="enabled"
          else if (eqsource.eq."ellipse" .and. n.eq.0) then
             eqcall="enabled"
          else
             eqcall="disabled"
          endif
       else
          eqcall="disabled"
       endif
       epsicon_=eqpsi(j)

       call eqorbit(epsicon_)

       eqfopsi(j)=fpsi_
       eqrpcon(j)=rpcon_
       eqrmcon(j)=rmcon_
       eqzpcon(j)=zpcon_
       eqzmcon(j)=zmcon_
       eqorb="disabled"
       call eqvolpsi(epsicon_,volum,areac)
       eqvol(j)=volum
       eqarea(j)=areac
       !..................................................................
       !     Compute <1./R**2>, <1./R>
       !..................................................................

       call eqonovrp(epsicon_,onovrp1,onovrp2)
       eqovrp(j,1)=onovrp1
       eqovrp(j,2)=onovrp2
       !
       zmaxpsi_=0.
       do 50 l=2,lorbit_
          zmaxpsi_=zmaxpsi_+(es_(l)-es_(l-1)) &
               /(.5*bpsi_(l)+.5*bpsi_(l-1))
50     continue
       if (j.eq.nconteqn) then
          do 60 l=1,lorbit_
             tlorb1(l)=eqbpol_(l)**2
60        continue
          call eqflxavg(epsicon_,tlorb1,bpolsqa_,flxavgd_)
          bpolsqlm=bpolsqa_/bmidplne_*zmaxpsi_/eqdells_
          bpolsqlm=bpolsqlm**2
       endif
       q_(j)=2.*fpsi_*onovrp2*zmaxpsi_/(twopi*bmidplne_)
       !        write(*,*)'eqrhopsi: j,eqvol(j),eqrmcon(j),eqrpcon(j),q_(j) =',
       !     +                       j,eqvol(j),eqrmcon(j),eqrpcon(j),q_(j)
20     continue
       eqrpcon(1)=rmag
       eqrmcon(1)=rmag
       if(rmag.gt.em8)then ! rmag=0 in a mirror machine
          eqovrp(1,1)=1./rmag
          eqovrp(1,2)=1./rmag**2
       endif
       call eqfpsi(eqpsimax,fpsi_,fppsi_)
       eqfopsi(1)=fpsi_

       !..................................................................
       !     rmaxcon=point of intersection of "largest" flx surface with
       !     midplane at low B field side; rmincon analogous for high field.
       !     BH140122:  This calculation could be improved upon, since there
       !     BH140122:  there are sometimes problems as the LCFS is
       !     BH140122:  approached.  For example, consider ONETWO type calc
       !     BH140122:  of the LCFS.
       !..................................................................

       rmaxcon=rpcon_
       rmincon=rmcon_
       !      write(*,*)'eqrhospi: rmincon,rmaxcon',rmincon,rmaxcon

       !..................................................................
       !     Determine the toroidal flux array, eqrho.
       !     Also determine the area associated with the flux surface.
       !..................................................................

       eqrho(1)=0.0
       if (radcoord.eq."sqtorflx") then
          do 30 j=2,nconteqn
             dvolum=(eqvol(j)-eqvol(j-1))
             eqrph=(eqovrp(j,2)+eqovrp(j-1,2))*.5
             eqfh=(eqfopsi(j)+eqfopsi(j-1))*.5
             eqrho(j)=eqrho(j-1)+eqrph*dvolum*eqfh/pi*0.5
30        continue
          do 40 j=2,nconteqn
             eqrho(j)=sqrt(eqrho(j)/pi/btor)
40        continue

       elseif (radcoord.eq."sqarea") then
          do 31 j=2,nconteqn
             eqrho(j)=sqrt(eqarea(j)/pi)
31           continue

       elseif (radcoord.eq."sqvol") then
          do 32 j=2,nconteqn
             eqrho(j)=sqrt(eqvol(j)/(2.*pi**2*rmag))
32        continue

       elseif (radcoord.eq."rminmax") then
          do 33 j=2,nconteqn
             eqrho(j)=0.5*(eqrpcon(j)-eqrmcon(j))
33        continue

       elseif (radcoord.eq."polflx") then
          do j=2,nconteqn
             eqrho(j)=(eqpsi(j)-psimag)/(psilim-psimag)
          enddo

       elseif (radcoord.eq."sqpolflx") then
          do j=2,nconteqn
             eqrho(j)=(eqpsi(j)-psimag)/(psilim-psimag)
             if (eqrho(j).ne.zero) eqrho(j)=sqrt(eqrho(j))
          enddo

       endif

       write(*,*)'eqrhopsi: psilim,psimag',psilim,psimag
       write(*,*)'eqrhopsi: rmag,zmag',rmag,zmag
       write(*,*)'eqrhopsi: eqpsi(j),j=1,nconteqn', &
            (eqpsi(j),j=1,nconteqn)
       write(*,*)'eqrhopsi: eqrho(j),j=1,nconteqn', &
            (eqrho(j),j=1,nconteqn)

       rhomax=exlin(eqrho(nconteqn-1),eqrho(nconteqn), &
            eqpsi(nconteqn-1),eqpsi(nconteqn),eqpsimin)
       radmin=rhomax
       btor00=fpsi_/rpcon_
       bthr00=eqbpol_(1)
       bmod00=sqrt(btor00**2+bthr00**2)
       write(*,*)'eqrhospi: rhomax = ',rhomax

       !.................................................................
       !     Setup 1D spline coeffs for eqrho (determined according to
       !     radcoord) versus eqpsi(1:nconteqn).
       !..................................................................

       i1p(1)=4
       i1p(2)=4
       !     Change sign of eqpsi, to get asceding order necessary
       !     for coeff1.  (Remember this when using the spline coeffs.)
       do j=1,nconteqn
          eqpsi(j)=-eqpsi(j)
       enddo
       call coeff1(nconteqn,eqpsi,eqrho,d2eqrho,i1p,1,workk)
       !     Change back eqpsi sign.
       do j=1,nconteqn
          eqpsi(j)=-eqpsi(j)
       enddo

       !.................................................................
       !     Determine the maximum axial position for the plasma.
       !..................................................................
       !      write(*,*)'eqrhopsi:solz_(j),j=1,lorbit_', (solz_(j),j=1,lorbit_)
       call aminmx(solz_,1,lorbit_,1,zmincon,zmaxcon,kmin,kmax)
       if (eqsym.ne."none") then
          zmaxcon=max(abs(zmincon),abs(zmaxcon))
          zmincon=-zmaxcon
       endif
       write(*,*)'eqrhospi: zmincon,zmaxcon',zmincon,zmaxcon

       !-----------------------------------------------------------------------
       !     Determine some approx geometrical factors needed to calculate
       !     the aspect ratio and the elongation [Assumes up-down symm.]
       !     Could check for non-up-down symmetric case, BH091023.
       !-----------------------------------------------------------------------
       !     left and right radius of psi=psilim surface at midplane
       !     assume magnetic axis is at the midplane of the eqdsk.
       !      izmag=nnz/2 + 1
       izmag=jmag
       rgeom1=exlin(er(jpsimnl),er(jpsimnl+1), &
         epsi(jpsimnl,izmag),epsi(jpsimnl+1,izmag),eqpsimin)
       rgeom2=exlin(er(jpsimnr),er(jpsimnr-1), &
         epsi(jpsimnr,izmag),epsi(jpsimnr-1,izmag),eqpsimin)
       !     minor radius of plasma surface
       rgeomp=0.5 * (rgeom2 - rgeom1)
       !     geometrical center of plasma surface
       r0geomp=0.5 * (rgeom1 + rgeom2)
       !
       irc=nnr/2 + 1
       do 100 jz=izmag,nnz
          jzpsimx=jz
          if (epsi(irc,jz+1) .lt. eqpsimin) go to 101
100    continue
101    continue
!
       itestl=0
       itestr=0
       jrpsmxl=irc
       jrpsmxr=irc
       do 102 ir=1,nnr-irc
          if (epsi(irc-ir,jzpsimx).lt.eqpsimin .and. itestl.eq.0) then
             jrpsmxl=irc - ir + 1
             itestl=1
          endif
          if (epsi(irc+ir,jzpsimx).lt.eqpsimin .and. itestr.eq.0) then
             jrpsmxr=irc + ir - 1
             itestr=1
          endif
          if (itestl.eq.1 .and. itestr.eq.1) go to 103
 102  continue
 103  continue
!
      jrmaxp=irc
      jzmaxp=jzpsimx
      do 104 jr=jrpsmxl,jrpsmxr
         do 105 jz=jzpsimx+1,nnz
            if (epsi(jr,jz) .lt. eqpsimin) then
               if (jz .le. jzmaxp) go to 106
               jrmaxp=jr
               jzmaxp=jz
               go to 106
            endif
105      continue
106      continue
104  continue

     !  not very precise, but used only for diagnostic of plasma elongation
     zgeomp=exlin(ez(jzmaxp-1),ez(jzmaxp),&
          epsi(jrmaxp,jzmaxp-1),epsi(jrmaxp,jzmaxp),eqpsimin)

      return

    end subroutine

  end module eqrhopsi_mod
