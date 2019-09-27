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

module cqldiag_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine cqldiag
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine controls the calculation of various
!     density and power diagnostics.
!     In particular, it determines density gain and loss due to
!     various physical and numerical processes.
!..................................................................


      if (n .gt. 0) go to 20

!..................................................................
!     compute the line averaged density at time=0 summed over all
!     general species
!..................................................................

      xlndn0(lr_)=0.
      do 10 k=1,ngen
        xlndn00(k,lr_)=xlndn(k,lr_)
        xlndn0(lr_)=xlndn0(lr_)+xlndn(k,lr_)
        do 31 iq=1,8
          sgaint(iq,k,l_)=0.
 31     continue
 10   continue
      sgaint1(l_)=0.
      consn0(l_)=xlndn0(lr_)
      return


 20   continue ! n>0
      dentot=0.

!..................................................................
!     compute cumulative in time particle density diagnostics.
!     see subs diagxswx and diagimpd.
!..................................................................

      do 40 k=1,ngen
        do 41 iq=1,8
          sgaint(iq,k,l_)=sgaint(iq,k,l_)+sgain(iq,k)
          sgaint1(l_)=sgaint1(l_)+sgain(iq,k)
 41     continue
        dentot=dentot+xlndn(k,lr_)
 40   continue

!..................................................................
!     compute density conservation diagnostic - the closer to
!     1.e-14 the better
!..................................................................

      consn(l_)=(dentot-xlndn0(lr_)-sgaint1(l_)) / &
        (.5*(xlndn0(lr_)+dentot))

!cc      if(l_.eq.1)then
!cc       write(*,*)'diag ',
!cc     + dentot-xlndn0(lr_),dentot-xlndn0(lr_)-sgaint1(l_),consn(l_)
!cc       write(*,'(i5,e13.5)') l_,sgain(4,1)
!cc      endif
!..................................................................
!     compute power transfer diagnostics (see sub diagentr)
!..................................................................

      do 60 k=1,ngen
        pelec=pelec+entr(k,2,l_)
        psou=psou+entr(k,5,l_)
        pwrf=pwrf+entr(k,3,l_)
 60   continue

!..................................................................
!     compute the total power being added to the system through
!     d.c. electric fields, RF power, and beam sources. (this does
!     not include power required to sustain backbround species
!     at specified fixed temperatures).
!..................................................................

      pinput=pelec+pwrf+psou
      return
    end subroutine cqldiag


end module cqldiag_mod
