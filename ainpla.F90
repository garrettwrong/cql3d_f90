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

module ainpla_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use cqlcomm_mod
  use cfpmodbe_mod, only : cfpmodbe
  use param_mod

  !---END USE

contains

  subroutine ainpla
    !implicit integer (i-n), real(c_double) (a-h,o-z)
    implicit none
    integer :: k
    integer :: k1
    integer :: l
    real(c_double) :: bk1
    real(c_double) :: bk2
    real(c_double) :: rstmss
    real(c_double) :: thta
    real(c_double) :: xq
    real(c_double) :: zeff1
    save

    !..........................................................
    !     This routine initialize some plasma parameter profiles
    !     on the whole setup0%lrzmax radial mesh
    !     (some where defined in diaggnde before)
    !.............................................................


    !.......................................................................
    !l    1. Energy, v-thermal
    !.......................................................................

    !l    1.1 radial mesh

    do 110 k=1,ntotal
       rstmss=fmass(k)*clite2/ergtkev
       do 111 l=1,setup0%lrzmax
          thta=rstmss/temp(k,l)
          if (thta.gt.100. .or. relativ.eq."disabled") then
             energy(k,l)=1.5*temp(k,l)
             !write(*,*)'ainpla.1:  energy(k,l)/1.5=',energy(k,l)/1.5
          else
             call cfpmodbe(thta,bk1,bk2)
             energy(k,l)=rstmss*(bk1/bk2-1.+3./thta)
             !write(*,*)'ainpla.2:  energy(k,l)/1.5=',energy(k,l)/1.5
          endif
          !          write(*,*)'ainpla: k,l,energy(k,l):',k,l,energy(k,l)
          vth(k,l)=((temp(k,l)*ergtkev)/fmass(k))**.5
          if (k .eq. kelec) vthe(l)=vth(kelec,l)
111    end do
110 end do

    !.......................................................................
    !l    1.2 parallel mesh
    !.......................................................................

    if (setup0%cqlpmod .eq. "enabled") then

       do 120 k=1,ntotal
          rstmss=fmass(k)*clite2/ergtkev
          do 121 l=1,setup0%lsmax
             thta=rstmss/temppar(k,l)
             if (thta.gt.100. .or. relativ.eq."disabled") then
                enrgypa(k,l)=1.5*temppar(k,l)
             else
                call cfpmodbe(thta,bk1,bk2)
                enrgypa(k,l)=rstmss*(bk1/bk2-1.+3./thta)
             endif
             vthpar(k,l)=((temppar(k,l)*ergtkev)/fmass(k))**.5
121       end do
          if (sbdry.eq."periodic" .and. transp.eq."enabled") then
             enrgypa(k,0)=enrgypa(k,setup0%lsmax)
             enrgypa(k,setup0%lsmax+1)=enrgypa(k,1)
             vthpar(k,0)=vthpar(k,setup0%lsmax)
             vthpar(k,setup0%lsmax+1)=vthpar(k,1)
          endif
120    end do

    endif
    !.......................................................................
    !     2. Compute radial Z-effective
    !.......................................................................

    if (izeff.eq."ion") then
        k1=ngen+1
     else
        k1=1
     endif
     do 200 l=1,setup0%lrzmax
        zeff(l)=0.
        zeff1=0.
        zeff4(l)=0.d0 !Yup[2014-05-27] Initialize to 0.
        xq=0.
        do 210 k=k1,ntotal
           if (k.eq.kelecg .or. k.eq.kelecm) goto 210
           !BobH990128          if (k.eq.izeff) goto 210
           xq=xq+1.
           zeff(l)=zeff(l)+bnumb(k)**2*reden(k,l)
           zeff4(l)=bnumb(k)**4*reden(k,l)+zeff4(l)
           zeff1=zeff1+bnumb(k)*reden(k,l)
210     end do
        zeff4(l)=zeff4(l)/xq
        zeff(l)=zeff(l)/zeff1
200  end do

     return
   end subroutine ainpla

end module ainpla_mod
