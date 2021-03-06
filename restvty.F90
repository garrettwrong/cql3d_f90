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

module restvty_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine restvty
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)


!..................................................................
!     Routine computes resistivity and related quantities
!     for electron resistivity calculations. A minimum of
!     jx=125 should be used or the electron-electron fokker-planck
!     energy transfer term will be large (should be 0).
!..................................................................

      if (kelec .gt. ngen) return

!..................................................................
!     [ old: resist is the resistivity - the axis toroidal electric
!     field divided by the flux surface averaged parallel current. ]
!
!
!     If efflag="toroidal", then code electric field is
!                           assumed to be purely toroidal,
!                           varying like elecfld*(Rmag/R).
!     If ifflag="parallel", then code electric field is
!                           assumed to be purely parallel,
!                           varying like elecfld*(Rmag/R).
!     Code current density is purely parallel, and varies
!       on a flux surface as
!       j_par = j_midplane * B(z)/B(0)(=bpsi(l))
!     resist=Resistivity, calc'd from distn fnctn results
!                  =<E_phi/R>/<j_phi/R>, is a "toroidal" resistivity.
!                  Except, if efswtchn.eq."neo_hh" .and.
!                    setup0%cqlpmod.ne."enabled" ==>
!                    resist=(pol cross-section-area avg of E)/currpar
!                    and currpar is sum of Hinton-Hazeltine neoclassical
!                    current + runaway current.

!     resistn=<E_parall*B>/<j_parall*B>/sptzr
!
!
!
!     rovsloc is the local resistivity E_par/j_par(no flux averaged),
!     meaningful in setup0%cqlpmod only
!
!     If efswtchn="neo_hh" and setup0%cqlpmod.ne."enabled":
!                           resist= (area avg of E)/currpar
!                           rovs = resist/(H&H neoclassical value).
!.......................................................................

      if (l_ .eq. lmdpln_) then
        resist=0.0
        resistn=0.0
      endif
      resistlo=0.0
      if (abs(currm(kelec,l_)) .gt. 1.e-10) then

         if (l_ .eq. lmdpln_) then
!%OS  old:   resist=elecfld(lr_)/(300.*curr(kelec,lr_))


            if (efflag.ne."parallel") then

            if (efswtchn.eq."neo_hh" .and. setup0%cqlpmod.ne."enabled") then
            resist=elecfld(lr_) / 300. /(currpar(lr_)*3.e9)  * rmag * &
                 fpsi(lr_)  / bmod0(lr_) * onovrp(2,lr_)/ &
                 psiavg(2,lr_)
            else
               resist=elecfld(lr_) / 300. / currm(kelec,l_) * rmag * &
                    bmod0(lr_) / fpsi(lr_)
            endif
            resistn=elecfld(lr_) / 300. / currm(kelec,l_) * rmag * &
                 fpsi(lr_)  / bmod0(lr_) * onovrp(2,lr_)/ &
                 psiavg(2,lr_)

            else                          !I.E., efflag.eq."parallel"

            if (efswtchn.eq."neo_hh" .and. setup0%cqlpmod.ne."enabled") then
               resist=elecfld(lr_)/300.*rmag*onovrp(2,lr_)/onovrp(1,lr_) &
                    /(currpar(lr_)*3.e9)
            else
               resist=elecfld(lr_) / 300. / currm(kelec,l_) * rmag * &
                    onovpsir3(lr_) / onovrp(2,lr_)
            endif
            resistn=elecfld(lr_) / 300. / currm(kelec,l_) * rmag * &
                 psiovr(lr_)/psiavg(2,lr_)

            endif

         endif

         if (setup0%cqlpmod .eq. "enabled") then   !Local along fld line resist
!%OS  if (mod(nummods,10).le.4 .or. n.ge.nontran) then
            resistlo=elecfld(lr_)/300.*rmag*fpsi(lr_)/bmidplne(lr_)/ &
                 (psis(l_)*solrs(l_)**2) / currm(kelec,l_)
!%OS  else
!%OS  resistlo=elecfld(lr_)/300.*rmag*fpsi(lr_)/bmidplne(lr_)/
!%OS  /           0.125/(psis(l_)+psis(l_+1))/(solrs(l_)+solrs(l_+1))**2
!%OS  /                                                  /currm(kelec,l_)
!%OS  endif
         endif
      endif


      rovsloc(l_)=resistlo/(sptzr(l_)+em90)

!..................................................................
!     rovs(lr_)=computed resistivity / Spitzer resistivity
!     Spitzer resistivity includes Zeff variation.
!..................................................................

      if (l_ .eq. lmdpln_) then
         if (efswtchn.eq."neo_hh" .and. setup0%cqlpmod.ne."enabled") then
            rovs(lr_)=resist/(zreshin(lr_)*sptzr(lmdpln_)+em90)
         else
            rovs(lr_)=resist/(sptzr(lmdpln_)+em90)
         endif
         rovsn(lr_)=resistn/(sptzr(lmdpln_)+em90)
         elect=clight*btor*300.*resist/(radmaj*2.*pi)
         eratio=elect/(elecr(lr_)+em90)

!..................................................................
!     vparl=average parallel velocity of distribution
!..................................................................

         vparl=curr(kelec,lr_)/(charge*reden(kelec,lr_))
         vpovth=vparl/vth(kelec,lr_)

!..................................................................c
!     eovedd is E-toroidal/E-Dreicer at the magnetic field
!..................................................................

         eovedd=elecfld(lr_)/(elecr(lr_)+em90)

      endif
!
      return
      end subroutine restvty

end module restvty_mod
