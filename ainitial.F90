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

module ainitial_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bavgmax_mod, only : bavgmax
  use coefefld_mod, only : coefefld
  use cqldiag_mod, only : cqldiag
  use diagentr_mod, only : diagentr
  use diaggnde_mod, only : diaggnde
  use efield_mod, only : efield
  use finit_mod, only : finit
  use losscone_mod, only : losscone
  use lossegy_mod, only : lossegy
  use micxinim_mod, only : micxinim
  use ntdstore_mod, only : ntdstore
  use pltmain_mod, only : pltmain
  use restvty_mod, only : restvty
  use sourcee_mod, only : sourcee
  use synchrad_mod, only : synchrad
  use tdoutput_mod, only : tdoutput
  use vlf_mod, only : vlf
  use vlh_mod, only : vlh
  use wpvptb_mod, only : wpvptb

  !---END USE

contains

      subroutine ainitial
      use param_mod
      use cqlcomm_mod
      use pltmain_mod, only : pltmain
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!............................................................
!     initializes the driver.
!............................................................

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif
!.......................................................................

      if (n.gt.0) return

!...........................................................
!     call a routine which calculates mesh integration
!     coefficients
! YuP[Sept-2014] moved this part outside of ainitial:
! needed for proper work of FOW-logic.
! See notes in tdinitl, after if(fow.ne.'disabled')...
!yup      call micxinil
!...........................................................

!.....................................................................
!     call orbit bounce time routine.  Also orbit shift calc.
! YuP[Sept-2014] moved this part outside of ainitial:
! needed for proper work of FOW-logic.
! See notes in tdinitl, after if(fow.ne.'disabled')...
!      if (l_ .eq. lmdpln_) then
!         if(taunew.eq."enabled") then
!yup            call baviorbt
!         else
!yup            call baviorbto
!         endif
!      endif
!.....................................................................

!.....................................................................
!     bounce average various geometrical quantities..
!.....................................................................

      if (l_ .eq. lmdpln_) call bavgmax

!.......................................................................
!     Redefine vptb=1.0 if CQLP
!.......................................................................

      if (setup0%cqlpmod.eq."enabled" .and. l_.eq.lmdpln_) call wpvptb

!..............................................................
!     Evaluate an additional integration coefficient.
!..............................................................

      if (l_ .eq. lmdpln_) call micxinim

!...........................................................
!     call routine which initializes distribution functions.
!...........................................................

      call finit
!.....................................................................
!     compute the orbit loss krook operator and the toroidal
!     loss operator
!.....................................................................
!yup [07-14-2014] Moved losscone in front of diaggnde, so that
!                 gone() array is available in diaggnde, at n=0.
      if (l_ .eq. lmdpln_) then
         call losscone
      endif

!...........................................................
!     call routine to calculate energies and densities.
!...........................................................

      call diaggnde

!.....................................................................
!     compute the orbit loss krook operator and the toroidal
!     loss operator
!.....................................................................
!yup      if (l_ .eq. lmdpln_) call losscone

!.....................................................................
!     Compute coefficients for synchrotron radiation loss
!.....................................................................

      if (l_ .eq. lmdpln_) call synchrad

!...................................................................
!     Compute coefficients for phenomenological energy loss
!...................................................................

      if (l_ .eq. lmdpln_) call lossegy

!.....................................................................
!     d.c. ohmic electric field coefficient set-up
!.....................................................................

      call coefefld

!.....................................................................
!     Call plasma resistivity diagnostic...
!.....................................................................

      call efield
      call restvty

!...........................................................
!     call driver routine which initializes analytic sources.
!...........................................................

      if (l_ .eq. lmdpln_) call sourcee

!.......................................................................
!     Call RF module if this is a 2-D calculation and if the RF
!     module is appended.
!.......................................................................

      if (setup0%lrzmax.eq.1 .or.(setup0%cqlpmod.eq."enabled".and.l_.eq.lmdpln_)) then
         if (nrf.ge.1) then
            if (vlfmod.eq."enabled") then
               call vlf("setup")
            elseif (vlhmod.eq."enabled") then
               call vlh("setup")
            else
!     call rf("setup")
               do 850 ku=1,ngen
                  call diagentr(3,ku)
 850           continue
            endif
         endif
      endif

      if (setup0%lrzmax.gt.1 .and. nrf.ge.1 .and. vlhmod.eq."enabled") &
           call vlh("setup")

      if (setup0%lrzmax.gt.1 .and. nrf.ge.1 .and. vlfmod.eq."enabled") &
           call vlf("setup") ! Added YuP[03-2016] Calling vlf for setup0%lrz>1

!.....................................................................
!     Call conservation diagnostic routine..
!.....................................................................

      call cqldiag

!.....................................................................
!     Call time plot storage routine...
!.....................................................................

      call ntdstore

!.....................................................................
!     Call plotting routine...
!.....................................................................

      call pltmain

!.....................................................................
!     print out some parameters
!.......................................................................

      if (setup0%lrzmax .eq. 1) call tdoutput(1)
!

      return
      end subroutine ainitial


end module ainitial_mod
