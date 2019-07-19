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

module prpprctr_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use pltmain_mod, only : gslnst
  use pltmain_mod, only : gslnsz
  use pltmain_mod, only : gsvp2d
  use pltmain_mod, only : gswd2d
  use pltmain_mod, only : gxglfr

  !---END USE

!
!

contains

      subroutine prpprctr
      use param_mod
      use cqlcomm_mod
      use pltdf_mod, only : cont, tempcntr, nconta
      use pltmain_mod, only : gslnst, gslnsz, gsvp2d, gswd2d, gxglfr
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!
!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!     Still need to finish contour modifications of GPCTV3.
!

      REAL RILIN !-> For PGPLOT (text output positioning)

!...............................................................
!     This routine performs contour plots of data
!     set in tempcntr by subroutine pltcont
!...............................................................
      ipjpxy=ipxy*jpxy
      dmin=0.
      dmax=0.
      do 10 jp=1,jpxy
        do 15 ip=1,ipxy
          xllji(jp,ip)=xpar(jp)
          xppji(jp,ip)=xperp(ip)
 15     continue
 10   continue
      call aminmx(fpn,1,ipjpxy,1,dmin,dmax,kmin,kmax)
      admin=abs(dmin)
      if (dmin.ge.0. .or. admin.lt.contrmin*dmax) then
        k2=1
!990131        smin=alog(contrmin*dmax)
!990131        if (admin/dmax .gt. contrmin) smin=alog(admin)
!990131        smax=alog(dmax)
        smin=log(contrmin*dmax)
        if (admin/dmax .gt. contrmin) smin=log(admin)
        smax=log(dmax)
        dcont=(smax-smin)/float(ncont)
        cont(1)=smin+.5*dcont
        do 20 kc=2,ncont
          cont(kc)=cont(kc-1)+dcont
 20     continue
        do 30 kc=1,ncont
          cont(kc)=exp(cont(kc))
 30     continue
      else
        if (dmax .gt. 0.) then
          k2=ncont/2+1
          ncontp=ncont-k2+1
          ncontm=k2-1
!990131          smaxp=alog(dmax)
!990131          sminp=alog(contrmin*dmax)
          smaxp=log(dmax)
          sminp=log(contrmin*dmax)
        else
          ncontm=ncont
          ncontp=1
          k2=1
        endif
!990131        sminm=alog(-contrmin*dmin)
!990131        if (dmax/dmin.gt.contrmin) sminm=alog(-dmax)
!990131        smaxm=alog(-dmin)
        sminm=log(-contrmin*dmin)
        if (dmax/dmin.gt.contrmin) sminm=log(-dmax)
        smaxm=log(-dmin)
        dcontp=(smaxp-sminp)/float(ncontp)
        dcontm=(smaxm-sminm)/float(ncontm)
        cont(1)=smaxm-.5*dcontm
        do 40 kc=2,ncontm
          cont(kc)=cont(kc-1)-dcontm
 40     continue
        do 50 kc=1,ncontm
          cont(kc)=-exp(cont(kc))
 50     continue
        if (dmax .gt. 0.) then
          cont(k2)=sminp-.5*dcontp
          do 60 kc=k2+1,ncont
            cont(kc)=cont(kc-1)+dcontp
 60       continue
          do 70 kc=k2,ncont
            cont(kc)=exp(cont(kc))
 70       continue
        endif
      endif

      call GXGLFR(0) ! new page
      call GSVP2D(.2,.8,.6,.9)
#ifndef NOPGPLOT
      CALL PGSCH(1.) ! set character size; default is 1.
#endif
      call GSWD2D("linlin$",xpar(1),xpar(jpxy),xperp(1),xperp(ipxy))
#ifndef NOPGPLOT
      call PGLAB('u/unorm_par','u/unorm_perp',' ')
#endif
      if (k2.gt.1) then
        call GSLNST(2)
        call GSLNSZ(.2)
        call GSLNSZ(0.)
      endif
      call GSLNST(1)
      call GSLNSZ(.2)
      call GSLNSZ(0.)

#ifndef NOPGPLOT
      CALL PGSLS(1) ! restore: solid line
#endif
#ifndef NOPGPLOT
      CALL PGSLW(setup0%lnwidth) ! restore linewidth
#endif

      return
      end subroutine prpprctr


end module prpprctr_mod
