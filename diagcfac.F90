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

module diagcfac_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use zcunix_mod, only : coeff2
  use zcunix_mod, only : terp1
  use zcunix_mod, only : terp2

  !---END USE

!
!

contains

      real(c_double) function diagcfac(epsil,zeff)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Cordey-Start estimates for current reduction due to electron
!     entrainment. A function of epsil and zeff.
!..................................................................

      parameter(n8=8,n12=12)
      dimension ae(n12),az(n8),af(n8,n12),afx(n8,n12),afy(n8,n12), &
        afxy(n8,n12),iabd(4),wke(100)
      data itimes / 0 /
      data ae / 0.,.01,.02,.04,.07,.1,.2,.3,.4,.5,.6,.9 /
      data az / 1.,1.1,1.2,1.5,2.,4.,8.,16. /
      data iabd / 2,2,2,2 /
      data af / 0.,.0909,.1667,.3333,.5,.75,.875,.9375, &
        .2242,.2892,.344,.4673,.5942,.7921,.8947,.947, &
        .3052,.3615,.4094,.5176,.6303,.8087,.9026,.9509, &
        .4093,.455,.4942,.5838,.6784,.8314,.9135,.9562, &
        .5106,.5468,.5779,.65,.7274,.8551,.9251,.9619, &
        .5824,.6122,.6379,.698,.7633,.8729,.9339,.9663, &
        .7315,.7492,.7646,.8012,.8421,.9132,.9543,.9765, &
        .8185,.8298,.8398,.8637,.8907,.939,.9676,.9833, &
        .8764,.8838,.8903,.9062,.9243,.9573,.9772,.9882, &
        .9173,.9221,.9263,.9367,.9488,.9709,.9843,.9919, &
        .947,.95,.9527,.9592,.9669,.981,.9898,.9947, &
        .9949,.9952,.9954,.996,.9968,.9981,.999,.9995 /
      if (itimes .eq. 1) go to 1
      itimes=1
      call bcast(afx,one,n8*n12)
      call bcast(afy,one,n8*n12)
      call bcast(afxy,zero,n8*n12)
      call coeff2(n8,az,n12,ae,af,afx,afy,afxy,n8,iabd,wke)
 1    continue
      s=terp2(zeff,epsil,n8,az,n12,ae,af,afx,afy,afxy,n8,0,0)
      diagcfac=s
      return
      end function diagcfac
end module diagcfac_mod
