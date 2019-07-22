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

module micgmbnd_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use psif_mod, only : psif

  !---END USE

!
!

contains

      subroutine micgmbnd(xlbmd2)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!......................................................................
!     This routine computes the pass/trapped boundary volume element.
!     See the user manual for mathematical details (pp. 101 and 130-131
!      of Killeen, Kerbel, McCoy, Mirin book).
!     It is called from subroutine micgnbnd.
!......................................................................

      save

      uupls=((y(itl-1,l_)+y(itl,l_))*.5)
      uumns=((y(itl,l_)+y(itl+1,l_))*.5)
      upls=cos(uupls)
      umns=cos(uumns)
      delmu=(upls-umns)*0.005
      hh=zmax(lr_)-zstar
!      write(*,*)'micgmbnd:uupls,uumns,upls,umns,delmu,hh,zmax(),zstar=',
!     +                    uupls,uumns,upls,umns,delmu,hh,zmax(lr_),zstar
      psi1=psif(zstar)
      psi2=psif(zmax(lr_))
!      Write(*,*)'micgmbnd: psi1,psi2=',psi1,psi2
      xlbmd2=0.
!      write(*,*)'micgmbnd: l,ul,umx2,ustar2,alph,ab,xlbmd2='
      do 10 l=1,200
        uz=l
        ul=(uz-.5)*delmu+umns
        ustar2=(1.-psif(zstar)*(1.-ul**2))
        ustar=sqrt(ustar2)
        umx2=(1.-psif(zmax(lr_))*(1.-ul**2))
        alph=1./sqrt(1.-umx2/ustar2)
        ab=abs((alph+1)/(alph-1))
!990131        xlbmd2=xlbmd2+twopi*hh*delmu*ul/ustar*alph*alog(ab)
        xlbmd2=xlbmd2+twopi*hh*delmu*ul/ustar*alph*log(ab)
!BH091031        write(*,100) l,ul,umx2,ustar2,alph,ab,xlbmd
 100    format(i4,3(1pe16.8),3(1pe10.2))
 10   continue
      return
      end subroutine micgmbnd


end module micgmbnd_mod
