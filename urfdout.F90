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

module urfdout_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfdout(ll)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This subroutine writes a disk file which is read by the Ray Tracin
!     code (Brambilla).
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

      character*1 line(132)
! jakub urban 110708: commented out for g95 compiler
! the above blanket save statement should do the job
!      save ifirst
      data ifirst /0/

#ifdef __MPI
      if(mpirank.ne.0) return
#endif

!..................................................................
!     Create or open a disk file.
!..................................................................

      isizeg=iy*jx*(lrors+1)*ngen+20000
      if(ll.eq.1.and.ifirst.eq.0)  then
        open(unit=31,file=idskf,delim='apostrophe',status='new')
      else if (ll.eq.1) then
        open(unit=31,file=idskf,delim='apostrophe',status='old')
      endif
      ifirst=1

!..................................................................
!     Write out ll (flux surf label) iy jx lrors (# of fl. surf)
!     Write out y (theta mesh) and x (momentum mesh).
!     Write out rovera(lr_) (r/radmin), vnorm (x mesh normal. factor).
!..................................................................


      write(31,555)  ll, iy,jx,lrors
      write(31,560)  (y(i,l_),i=1,iy)
      write(31,560)  (x(j),j=1,jx)
      do 1000 k=1,ngen
        write(31,560) rovera(lr_),vnorm,elecfld(lr_),eovedd,reden(k,lr_)
        vmaxdvt=vnorm/(4.19e7*sqrt(2.*1000.*temp(k,lr_)))
        write(31,560)  temp(k,lr_),hnis(k,lr_),vmaxdvt
        write(31,560)  ((f(i,j,k,l_),i=1,iy),j=1,jx)
 1000 continue
 555  format(5i16)
 560  format(5e16.8)
      if(ll.eq.setup0%lrz)  close(unit=31)

      return
      end subroutine urfdout

end module urfdout_mod
