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

module dskout_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine dskout(ll)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     At the end of the run, this routine at the option of the user
!     writes out to disk file 'idskf' the namelist input deck
!     and various computed quantities, including the distn functions.
!     Also, write out to disk file 'idskrf' rf diffusion coefficients,
!     and related quantities.
!..................................................................

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

      character*80 line

#ifdef __MPI
      if(mpirank.ne.0) return
#endif

      if (setup0%lrzmax.le.1) then
        if((idskf.eq. "disabled".and.idskrf.eq."disabled") &
          .or. n.ne.nstop+1)  return
      else
        if((idskf.eq. "disabled".and.idskrf.eq."disabled") &
          .or. n.ne.nstop)  return
      endif
      if(ll.gt.1)  go to 4
      if(idskf.ne."disabled") &
        open(unit=4,file=idskf,delim='apostrophe',status='unknown')
      if(idskrf.ne."disabled") &
        open(unit=5,file=idskrf,delim='apostrophe',status='unknown')
!cc      close(unit=2) ! YuP: Why here?
      ilen=0

!..................................................................
!     The input namelist file is transcribed onto the beginning
!     of file idskf and/or idskrf
!..................................................................

      open(unit=2,file='cqlinput',delim='apostrophe',status='old')
 1    read(2,1003) line
      if (line(1:3).eq."end") go to 3
      if(idskf.ne."disabled") write(4,1003) line
      if(idskrf.ne."disabled") write(5,1003) line
      go to 1
 3    close(unit=2)
      if(idskf.ne."disabled") &
           write(4,1006) '***Begin computed output***'
      if(idskrf.ne."disabled") &
           write(5,1006) '***Begin computed output***'
 4    continue

!..................................................................
!     In the following disk write to file named idskf:
!     This subroutine is called to write data for each FP'd
!          flux surface.
!     ll=  FP flux surface number
!          (ll=1:lrors, lrors.le.setup0%lrzmax, see cqlinput_help))
!          (setup0%lrindx(ll) gives flux surface number on the full
!                     radial mesh.(setup0%lrindx(ll)=ll if lrors=setup0%lrzmax,
!                     and using cql3d mode(setup0%cqlpmod="disabled"))).
!     iy,jx= dimensions in theta and u(momentum/mass)
!           (In the case where iy varies with ll, iy(1) will be greatest.)
!     lrors= number of flux surfaces FP'd.
!     setup0%lrzmax= number of flux surfaces, including any not FP'd.
!     x = momentum-per-mass(nomalized to maximum 1.)
!         at  each flux surface setup0%lrindx(ll)
!     y = theta(radians) mesh at  each flux surface setup0%lrindx(ll)
!     rovera= normalized radius (ll)
!             (rho, see Hinton and Haseltine for non-circ).
!             (generally ~sqrt(tor. flux), other coords available.)
!     elecfld = toroidal electric field (volts/cm)
!     bthr(lr_) - the poloidal magnetic field at theta-poloidal = pi/2.
!     btoru(lr_) - the toroidal magnetic field at the same position.
!     bthr0(lr_), btor0(lr_), are the poloidal and  toroidal
!         magnetic fields at the outer midplane of the flux surface.
!     reden= electron density at minimum B point on flux surface.
!     temp= initial electron temperature (keV)
!     radmin= plasma minor radius (cms).
!     vnorm= normalization momentum-per-mass (maximum on grid) (cm/sec)
!     vmaxdvt= vnorm/(temp/mass)**0.5
!     eovedd= electric  field, normalized to Driecer field
!             (calc'd in sub restvty).
!
!     distribution function normalised so that
!         integral( (dx)**3 f) = density at minumum B point.
!
!..................................................................

      if(idskf.ne."disabled") then
        write(4,1004)  ll, iy,jx,lrors,setup0%lrzmax,ngen
        write(4,1004)  itl,itu
        write(4,1005)  (x(j),j=1,jx)
        write(4,1005)  (y(i,ll),i=1,iy)
        do 1000 k=1,ngen
          write(4,1005)  bnumb(k),fmass(k)
          write(4,1005)  rovera(setup0%lrindx(ll)),elecfld(setup0%lrindx(ll)), &
                         bthr(setup0%lrindx(ll)),btoru(setup0%lrindx(ll))
          write(4,1005)  bthr0(setup0%lrindx(ll)),btor0(setup0%lrindx(ll)), &
                         reden(k,setup0%lrindx(ll)),temp(k,setup0%lrindx(ll))
          vmaxdvt=vnorm/(4.19e7*sqrt(2.*1000.*temp(k,setup0%lrindx(ll))))
          write(4,1005)  radmin,vnorm,vmaxdvt,eovedd
          write(4,1005)  ((f(i,j,k,ll),i=1,iy),j=1,jx)
 1000   continue
      endif


!..................................................................
!     In the following disk write:
!     vptb is lambda=u_parallel_o * tau_bounce.
!     (tau_bounce= portion of bounce time from the outer equatorial
!     plane to the bounce point or the inner equatorial plane).
!     temp1 is the bounce averaged uu-diffusion coefficient (cgs units).
!..................................................................
      if(idskrf.ne."disabled") then
        write(5,1004)  ll,iy,jx,lrors,setup0%lrzmax
        write(5,1005)  (y(i,ll),i=1,iy)
        write(5,1005)  (x(j),j=1,jx)
        write(5,1005)  rovera(setup0%lrindx(ll)),vnorm
        write(5,1005)  (vptb(i,setup0%lrindx(ll)),i=1,iy)
        vn2=vnorm*vnorm
        do 10 k=1,mrfn
          vmaxdvt=vnorm/(4.19e7*sqrt(2.*1000.*temp(k,setup0%lrindx(ll))))
          write(5,1005) reden(k,setup0%lrindx(ll)),temp(k,setup0%lrindx(ll)),vmaxdvt
          do 11 i=1,iy
 11       temp1(i,1)=0.0
          do 12 j=2,jx
            x2=x(j)**2
            do 13 i=1,iy
 13         temp1(i,j)=urfb(i,j,indxlr(setup0%lrindx(ll)),k)*vn2/ &
                (vptb(i,setup0%lrindx(ll))*x2)
 12       continue
          write(5,1005) ((temp1(i,j),i=1,iy),j=1,jx)
 10     continue
      endif
 1003 format(a80)
 1004 format(16i5)
 1005 format(5e16.8)
 1006 format(a27)
      if(idskf.ne."disabled".and.ll.eq.lrors)  close(unit=4)
      if(idskrf.ne."disabled".and.ll.eq.lrors)  close(unit=5)
      return
      end subroutine dskout

end module dskout_mod
