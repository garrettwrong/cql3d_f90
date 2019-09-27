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

!     wpadvnc.h
!***********************************************************************
!***********************************************************************

!..................................................................
!     wpadvnc contains the statement functions utilized for the parallel
!     transport implicit time advancement routines and its diagnostics
!..................................................................

!     coefficient u_parallel/ds
      wpweigt(i,j,l)=-vnorm*x(j)*coss(i,l)/dsz(l)
      wpweigp(i,j,l)=-vnorm*x(j)*eszp5(l)*coss(i,l)
!     f_s+1/2 (at cst i)
      fsmid(i,j,k,l)=fnp1(i,j,k,l)*dls(i,j,k,l) + &
      fnp1(ilpm1ef(i,l,+1),j,k,lpm1eff(l,+1))*(1-dls(i,j,k,l))

!     contribution to matrix for f(s-1), f(s) and f(s+1)
      wprhsmd(i,j,k,l)=fnhalf(i,j,k,l)/dtreff+velsou(i,j,k,l)
      wprhs(i,j,k,l)=fnhalf(i,j,k,l)/dtreff + &
	0.5*(velsou(i,j,k,l)+velsou(ilpm1ef(i,l,+1),j,k,lpm1eff(l,+1)))

!..................................................................
!     End of statement functions used for parallel transport
!..................................................................
