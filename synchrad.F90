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

module synchrad_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

  !
  !

contains

  subroutine synchrad
    use param_mod
    use cqlcomm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save

    !..................................................................
    !     This routine computes  synchrotron radiation loss coefficients
    !     for electrons.
    !     syncrad="disabled","gyro","geom", or "gyrogeom".
    !     The relativistic coefficients have been calculated
    !     using the hueristic approach of Bernstein and Baxter.
    !     (I.B. Bernstein and D.C. Baxter, Phys. Fluids 24,
    !      p. 108 (1981), Appendix A (not equation given in main text).)
    !     Plus
    !     synchrotron radiation due to toroidal and poloidal rotation
    !     (BH + Paul Parks, March 2000).
    !..................................................................



    !     Temporary, for writing:
    dimension synca_gyro(jx),syncd_gyro(jx)
    dimension synca_geom(jx),syncd_geom(jx)

    if (syncrad .eq. "disabled" .or. kelecg .eq. 0) return
    call bcast(synca,zero,iyjx*setup0%lrz)
    call bcast(syncd,zero,iyjx*setup0%lrz)

    if (syncrad.eq."gyro" .or. syncrad.eq."gyrogeom") then

       coefnt=(2.*(charge*bnumb(kelecg))**4*bmod0(lr_)**2) / (3.*clight**5*fmass(kelecg)**3)
       do  j=1,jx
          do  i=1,iy
             synca(i,j,indxlr_)=coefnt*vptb(i,lr_)*sinn(i,l_)**2*xcu(j) &
                  * psicu(i,lr_)*gamma(j)
             syncd(i,j,indxlr_)=coefnt*vptb(i,lr_)*xsq(j)*gammi(j) &
                  * sinn(i,l_)**2/coss(i,l_) * (psisq(i,lr_)-sinn(i,l_)**2*psicu(i,lr_))
          enddo
       enddo

       call bcast(synca_gyro,zero,jx)
       call bcast(syncd_gyro,zero,jx)
       do j=1,jx
          !         do i=1,iy
          i=11 !   comment from YuP[2018-01-16] Why i=11 only ?
          synca_gyro(j)=synca_gyro(j) + synca(i,j,indxlr_)
          !     *                    sinn(i,l_)*dy(i,l_)*synca(i,j,indxlr_) &

          syncd_gyro(j)=syncd_gyro(j) + syncd(i,j,indxlr_)
          !     *                    sinn(i,l_)*dy(i,l_)*syncd(i,j,indxlr_) &

          !         enddo
       enddo

    endif

    if (syncrad.eq."geom" .or. syncrad.eq."gyrogeom") then

       coefnt1=2.*(charge*bnumb(kelecg))**2/(3.*clight**3*fmass(kelecg))
       bsq=bthr(lr_)**2+btoru(lr_)**2
       bphiob4=(btoru(lr_)**2/bsq)**2
       bthob4=(bthr(lr_)**2/bsq)**2
       reffm2=onovrp(2,lr_)*bphiob4 + pi/areacon(lr_)*bthob4

       call bcast(synca_geom,zero,jx)
       call bcast(syncd_geom,zero,jx)
       do  j=1,jx
          do  i=1,iy

             sync1=  +coefnt1*vptb(i,lr_) * &
                  xcu(j)*xsq(j)*vnorm2*gammi(j)*reffm2* &
                  (1.-2.*sinn(i,l_)**2*psiba(i,lr_)+ &
                  sinn(i,l_)**4*psisq(i,lr_))
             synca(i,j,indxlr_)=synca(i,j,indxlr_)+sync1
             sync2=  -coefnt1*vptb(i,lr_)* &
                  xcu(j)*x(j)*vnorm2*gammi(j)*reffm2* &
                  sinn(i,l_)**2/coss(i,l_)* &
                  (1.-2.*sinn(i,l_)**2*psiba(i,lr_)+ &
                  sinn(i,l_)**4*psisq(i,lr_))
             syncd(i,j,indxlr_)=syncd(i,j,indxlr_)+sync2

             if (i.eq.1) then
                synca_geom(j)=synca_geom(j) + sync1
                !     *                    sinn(i,l_)*dy(i,l_)*sync1 &

                syncd_geom(j)=syncd_geom(j) + sync2
                !     *                    sinn(i,l_)*dy(i,l_)*sync1 &

             endif
          enddo
       enddo

    endif

      write(*,*) 'syncrad: synca_gyro: ',synca_gyro
      write(*,*) 'syncrad: syncd_gyro: ',syncd_gyro
      write(*,*) 'syncrad: synca_geom: ',synca_geom
      write(*,*) 'syncrad: syncd_geom: ',syncd_geom


    return
  end subroutine synchrad

end module synchrad_mod
