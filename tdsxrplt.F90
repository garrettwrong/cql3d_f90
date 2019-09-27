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

module tdsxrplt_mod

  !---BEGIN USE

  use aminmx_mod, only : aminmx
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine tdsxrplt(en,eflux,nen,nenaa, &
                          efluxt,nv,inegsxr,softxry,lnwidth)
      use param_mod
      use r8subs_mod, only : rbound
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save


#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

      real(c_float) RTAM1(nena),RTAM2(nena)
      real(c_float) REMAX,REMIN
!BH092022: XXX
      !real(c_float) RBOUND

!..................................................................
!mnt  this routine plots SXR/NPA energy/particle flux spectra
!mnt  versus photon/particle energy.
!mnt  It would be simpler to have a separate tdnpaplt [BH100815].
!..................................................................

      dimension en(nenaa),eflux(nenaa,*),efluxt(nv),inegsxr(nv)
      character*1024 t_
      character*8 softxry

!      write(*,*)'tdsxrplt: en(1:nen)',
!     +     en(1:nen)
!      write(*,*)'tdsxrplt: eflux(1:nen,1:nv)',
!     +    (eflux(1:nen,i),i=1,nv)
!      write(*,*)'tdsxrplt: efluxt(1:nv),inegsxr(1:nv)',
!     +     efluxt(1:nv),inegsxr(1:nv)
!      write(*,*)'tdsxrplt: softxry=',softxry

#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      em100=1.d-100

      fmin=ep100
      fmax=-ep100 ! initialize to negative, will be found below
      do 100  nn=1,nv
        if (inegsxr(nn).le.1) go to 100 ! and then fmax remains negative
        call aminmx(eflux(1:inegsxr(nn),nn),1,inegsxr(nn),1 &
          ,fmin1,fmax1,kmin,kmax)
        fmin=min(fmin,fmin1)
        fmax=max(fmax,fmax1)
 100  continue
      ! Note: if inegsxr(nn).le.1 for ALL nn, then fmax remains equal
      ! to -ep100. This can happen when all sightlines missed plasma.
      !YuP[2018-02-08] Added:
      if(fmax.lt.0.d0)then
        WRITE(*,*)'tdsxrplt: All sightlines missed plasma. Skip plots.'
        return
      endif


      if (fmin .eq. fmax) fmin=.9*fmax-1.e-20
      decades=6.1
      if(decades.ge.0.)  fmin=fmax/10**decades

      DO J=1,nen
         RTAM1(J)=RBOUND(en(j))
      ENDDO

      !write(*,*)'tdsxrplt: fmin,fmax=',fmin,fmax
      ! YuP[2018-02-08] Sometimes fmin<0. Need to check this.
      ! For now, reset to a small value
      !fmin= max(fmin, .9*fmax-1.e-20) !YuP[2018-02-08]

      REMIN=RBOUND(LOG10(fmin))
      REMAX=RBOUND(LOG10(fmax))
      !write(*,*)'tdsxrplt: LOG10(fmin),LOG10(fmax)=',REMIN,REMAX


#ifndef NOPGPLOT
      CALL PGPAGE
#endif
#ifndef NOPGPLOT
!     CALL PGENV(Rtam1(1),Rtam1(nen),Remin,Remax,0,20)
#endif
#ifndef NOPGPLOT
      CALL PGSVP(.2,.8,.45,.9)
#endif
#ifndef NOPGPLOT
      CALL PGSWIN(Rtam1(1),Rtam1(nen),Remin,Remax)
#endif
#ifndef NOPGPLOT
      CALL PGBOX('BCNST',0.,0,'BCNSTL',0.,0)
#endif
#ifndef NOPGPLOT
      CALL PGSAVE
#endif
#ifndef NOPGPLOT
      CALL PGSCH(1.44)
#endif
      if (softxry.eq.'enabled') then
#ifndef NOPGPLOT
      CALL PGLAB('Photon Energy k (keV)', &
           'd\u2\d\(0555)/dtdk (ergs/cm\u2\d-sec-ster-eV)', &
           'SXR Energy Flux versus Photon Energy')
#endif
      else
#ifndef NOPGPLOT
      CALL PGLAB('Particle Energy k (keV)', &
           'd\u2\dN/dtdk (#/cm\u2\d-sec-ster-eV)', &
           'NPA Flux versus Energy')
#endif
      endif

#ifndef NOPGPLOT
      CALL PGUNSA
#endif

      do 200 nn=1,nv ! view lines
        if (inegsxr(nn) .le. 1) go to 200
#ifndef NOPGPLOT
        CALL PGSLW(lnwidth) ! line thickness/width
#endif
        DO J=1,inegsxr(nn)
           rtam2(J)=rbound(eflux(j,nn))
!           write(*,*)'tdsxrplt.f: nn,j,eflux(j,nn),rtam2(j):',
!     +                            nn,j,eflux(j,nn),rtam2(j)
           RTAM2(J)=ABS(RTAM2(J))
           if (RTAM2(J).gt.em100) then
              RTAM2(J)=LOG10(RTAM2(J))
           else
              RTAM2(j)=1.0
           endif
        ENDDO ! J=1,inegsxr(nn)
#ifndef NOPGPLOT
        CALL PGLINE(inegsxr(nn),RTAM1,RTAM2)
#endif
 200  continue ! nn

#ifndef NOPGPLOT
       CALL PGSLW(lnwidth) ! restore line thickness/width
#endif

      if (softxry.eq.'enabled') then
         write(t_,610)
      else
         write(t_,613)
      endif
 610  format("total flux, enmin to enmax (ergs/cm**2-sec-ster):")
 613  format("total flux, enmin_npa to enmax_npa (#/cm**2-sec-ster):")

#ifndef NOPGPLOT
      CALL PGMTXT('B',7.,-0.1,0.,t_)
#endif

      do 300  nn=1,nv,2
         if (nn.eq.nv .and. ((nv/2)*2 .ne. nv)) then
            write(t_,612) efluxt(nv)
         else
            write(t_,611)  (efluxt(in),in=nn,nn+1)
         endif
#ifndef NOPGPLOT
         CALL PGMTXT('B',7.5+0.5*nn,0.,0.,t_)
#endif
 300  continue
!$$$      do 300  nn=1,nv,4
!$$$         if (nn.eq.nv .and. ((nv/4)*4 .ne. nv)) then
!$$$            write(t_,612) efluxt(nv)
!$$$         else
!$$$            write(t_,611)  (efluxt(in),in=nn,nn+3)
!$$$         endif
#ifndef NOPGPLOT
!$$$         CALL PGMTXT('B',7.+nn,0.,0.,t_)
#endif
!$$$ 300  continue
 611  format(1p4e12.2)
 612  format(1pe12.2,12x)

      return
      end subroutine tdsxrplt


!
!
      subroutine tdsxrvw(tempp4,tempp5,tempp6)
        use cqlconf_mod, only : setup0
        use param_mod
        use cqlcomm_mod, only : eqsym, softxry, ez, er, iyjx
        use cqlcomm_mod, only : lorbit, solr, solz, rcontr, zcontr, lensxr
      implicit integer (i-n), real(c_double) (a-h,o-z)
#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

      character*8  pltsxrvw

      real(c_float) ZTOP
      real(c_float) RTAB1(LFIELDA),RTAB2(LFIELDA)
      real(c_float) PGER1,PGERNNR,PGEZNNZ

      real(c_float) RRTAB1,RRTAB2
      dimension rrtab1(:), rrtab2(:)
      dimension tempp4(*),tempp5(*),tempp6(*)
      pointer rrtab1, rrtab2
      allocate(rrtab1(iyjx),STAT=istat)
      allocate(rrtab2(iyjx),STAT=istat)


!..................................................................
!     This routine plots out the contours (flux surfaces) in a
!     poloidal plane, and overplots the SXR view cords.
!..................................................................

#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      pltsxrvw="enabled"
      if (pltsxrvw.eq."disabled") return


#ifndef NOPGPLOT
      CALL PGPAGE
#endif

      if (eqsym.ne."none") then
         ztop=2.*.5*ez(nnz)/(er(nnr)-er(1))+.05
      else
         ztop=.95
      endif

#ifndef NOPGPLOT
      CALL PGSVP(.15,.85,.15,ztop)
#endif

      PGER1=er(1)
      PGERNNR=er(nnr)
!      write(*,*)'tdsxrvw: PGER1,PGERNNR=',PGER1,PGERNNR
      PGEZNNZ=ez(nnz)
#ifndef NOPGPLOT
      CALL PGSWIN(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
#endif
#ifndef NOPGPLOT
      CALL PGWNAD(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
#endif
#ifndef NOPGPLOT
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
#endif
#ifndef NOPGPLOT
      CALL PGSLW(lnwidth) ! line thickness/width
#endif

      if (softxry.eq."enabled") then
#ifndef NOPGPLOT
         CALL PGLAB('Major radius (cms)','Vert height (cms)', &
              'Flux Surfaces and SXR Chords')
#endif
      else  ! NPA
#ifndef NOPGPLOT
         CALL PGLAB('Major radius (cms)','Vert height (cms)', &
              'Flux Surfaces and NPA Chords')
#endif
      endif

      do 10 l=1,setup0%lrzmax
         IF (LORBIT(L).GT.LFIELDA) STOP 'tdsxrvw: CHECK DIM OF RTAB1/2'

        if(eqsym.ne."none") then
!     Up-Down symmetric flux surfaces are plotted:
           do 20 j=1,lorbit(l)
              RTAB1(j)=solr(lorbit(l)+1-j,l)
              RTAB2(j)=abs(solz(lorbit(l)+1-j,l))
 20        continue
#ifndef NOPGPLOT
           CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
#endif
           do 21 j=1,lorbit(l)
              RTAB1(j)=solr(lorbit(l)+1-j,l)
              RTAB2(j)=-abs(solz(lorbit(l)+1-j,l))
 21        continue
#ifndef NOPGPLOT
           CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
#endif

        else
!     Full flux surface contours are plotted:
           write(*,*)
           write(*,*)'tdsxrplt: l,lorbit(l)=',l,lorbit(l)
           do j=1,lorbit(l)
!BH091031             write(*,*)'j,solr(j,l),solz(j,l)=',
!BH091031    +             j,solr(j,l),solz(j,l)
           enddo
           do j=1,lorbit(l)
              RTAB1(j)=solr(lorbit(l)+1-j,l)
              RTAB2(j)=solz(lorbit(l)+1-j,l)
           enddo
#ifndef NOPGPLOT
           CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
#endif
        endif
 10   continue

      if(eqsym.eq.'none' .and. ncontr.gt.1) then
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

      iistep=0
      do 30 nn=1,nv

        do 40 j=1,lensxr(nn)
           if(j.le.iyjx)  then
             RRTAB1(j)=sqrt(tempp4(iistep+j)**2+tempp5(iistep+j)**2)
             RRTAB2(j)=tempp6(iistep+j)
           else
             write(*,*)'tdsxrvw: Increase dimension of rrtab1,2'
             STOP
           endif
 40     continue

!        write(*,*) 'tdsxrvw: RRTAB1',(RRTAB1(jj),jj=1,lensxr(nn))
!        write(*,*) 'tdsxrvw: RRTAB2',(RRTAB2(jj),jj=1,lensxr(nn))

#ifndef NOPGPLOT
        CALL PGLINE(lensxr(nn),RRTAB1,RRTAB2)
#endif

        iistep=iistep+lensxr(nn)

 30   continue

      deallocate(rrtab1,STAT=istat)
      deallocate(rrtab2,STAT=istat)

      return
      end subroutine tdsxrvw


end module tdsxrplt_mod
