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

module abchief_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aclear_mod, only : aclear
  use tdchief_mod, only : tdchief

  !---END USE
  !
  !

contains

  subroutine abchief(nml_file)
    use param_mod
    use cqlcomm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save

    character(len=*), intent(in), optional :: nml_file

    !
    !      cputime=second()

    !.......................................................................
    !     2d code controlled by achief1, called at beginning of tdchief
    !.......................................................................
    call tdchief(nml_file)
    call aclear ! reset for STEP mode
    return
  end subroutine abchief

end module abchief_mod
