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

module urfrays_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfrays(initrf,nrayptlh,nrayptec,nrayptfw)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.................................................................
!     This routine calls ray tracing code(s), if namelist
!     variables, call_lh, call_ech, or call_fw are "enabled".
!.................................................................


      save icall
      data icall /0/

!..................................................................
!     Determine whether to call the ray tracing code or just to assume
!     that a disk file with r.t. output already exists
!     xbr holds the name of the r.t. controllee.
!..................................................................

      if (call_lh .eq."enabled".and.lh.eq."enabled") then

!..................................................................
!     Set up a disk file to control LH ray tracing code.
!     initrf=0 for first call..
!     nraypts is the maximum number of ray elements allowed per ray.
!..................................................................

        if (nrayptlh.gt.nrayelts) nraypts=nrayelts
        open(unit=21,file='iraylh',status='new')
        write(21,102) initrf,nraypts,ieqbrurf
        close(unit=21)
        write (*,101) n,1
!990131        istat=ishell('xbr')
        write (*,100) n,1
      endif
      if (call_ech.eq."enabled".and.ech.eq."enabled") then

!..................................................................
!     Set up a disk file to control ECH ray tracing code.
!     initrf=0 for first call..
!     nraypts is the maximum number of ray elements allowed per ray.
!..................................................................

        if (nrayptec.gt.nrayelts) nraypts=nrayelts
        open(unit=22,file='irayec',status='new')
        write(22,102) initrf,nraypts,ieqbrurf
        close(unit=22)
        write (*,101) n,2
!990131        istat=ishell('toray')
        write (*,100) n,2
      endif
      if (call_fw.eq."enabled".and.fw.eq."enabled") then

!..................................................................
!     Set up a disk file to control FW ray tracing code.
!     initrf=0 for first call..
!     nraypts is the maximum number of ray elements allowed per ray.
!..................................................................

        if (nrayptfw.gt.nrayelts) nraypts=nrayelts
        open(unit=25,file='irayec',status='new')
        write(25,102) initrf,nraypts,ieqbrurf
        close(unit=25)
        write (*,101) n,2
!990131        istat=ishell('xbr')
        write (*,100) n,2
      endif
 100  format("Call for n=",i5, " mode=",i5,"complete for ray tracing")
 101  format("Call for n=",i5," mode=",i5, "begins for ray tracing")
 102  format(3i5)
      if (icall.eq.0) imprf=1
      icall=1
      return
      end subroutine urfrays

end module urfrays_mod
