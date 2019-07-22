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

module eqinitl_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

  subroutine eqinitl
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
#ifdef __MPI
      include 'mpilib.h'
#endif

      REAL RILIN

!..................................................................
!     This routine does some minor initialization for
!     the "eq" module. Called after the namelist read.
!..................................................................

      if (lfield.gt.lfielda) lfield=lfielda
      nrc=(nnr-1)/2+1
      nzc=(nnz-1)/2+1
      zshift=0.0


#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      if (setup0%noplots.ne."enabled1") then
#ifndef NOPGPLOT
      CALL PGPAGE
#endif
      RILIN=0.
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,"PARAMETER VALUES")
#endif

      write(t_,1000)
 1000 format("EQUILIBRIUM model parameters:")
      RILIN=2.
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif

      write(t_,1001)
 1001 format("nnra,nnza give the Maximum size the eqdsk")
      RILIN=3.
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif

      write(t_,1002) nnra,nnza
 1002 format("====>NNRA = ",i5,"        ====>NNZA = ",i5)
      RILIN=4.
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif

      write(t_,1003) nconteqa
 1003 format("====>NCONTEQA = ",i5)
      RILIN=5.
#ifndef NOPGPLOT
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
#endif

      endif


      return
      end subroutine eqinitl

end module eqinitl_mod
