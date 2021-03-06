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

module diagdenz_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use cfpleg_mod, only : cfpleg
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine diagdenz
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Computes density as a function of angle= twopi*z(l,lr_)/zmax(lr_)
!     and as a function of energy range. The results are used as
!     diagnostics in plots, on in the event that locquas="enabled"
!     to establish the background electron density as a function of
!     poloidal angle to maintain LOCAL (in z) quasi-neutrality.
!     negyrg is the number of energy bins over which we will be
!     computing density. For each energy range [eegy(ny,1,k,lr_)
!     ,eegy(ny,2,k,lr_)]
!     the velocity range associated with it must be determined.
!     The indicator will be jegy(ny,1,k,lr_) (lower) and
!     jegy(ny,2,k,lr_) (upper). Note that the
!     final energy bin or range is forced to be the entire energy range.
!..................................................................


      if (negyrg.lt.1) return
      do 103 k=1,ngen

!..................................................................
!     Determine the lower and upper momentum indices for the given
!     energy range.
!..................................................................

        do 100 ny=1,negyrg

!..................................................................
!     ny=negyrg forces consideration of complete energy range.
!..................................................................

          if (ny .eq. negyrg) then
            jegy(ny,1,k,lr_)=1
            jegy(ny,2,k,lr_)=jx
            eegy(negyrg,1,k,lr_)=0.
            eegy(negyrg,2,k,lr_)=fions(k)*tcsgm1(jx)
            go to 100
          endif
          do 101 j=2,jx
            if (fions(k)*tcsgm1(j) .gt. eegy(ny,1,k,lr_)) go to 104
 101      continue
          jegy(ny,1,k,lr_)=0
          go to 108
 104      continue
          jegy(ny,1,k,lr_)=j-1
          do 106 j=2,jx
            if (fions(k)*tcsgm1(j) .gt. eegy(ny,2,k,lr_)) go to 107
 106      continue
          jegy(ny,2,k,lr_)=jx
          go to 108
 107      jegy(ny,2,k,lr_)=j
 108      continue
          if (jegy(ny,1,k,lr_) .ge. jegy(ny,2,k,lr_)) jegy(ny,1,k,lr_)=0
 100    continue
 103  continue
      call bcast(xlndnz(1:(ngen+1),1:negyrga),zero,(ngen+1)*negyrga)

!..................................................................
!     xlndnz(k,ny) will eventually hold the line integrated density
!
!     densz(l,k,ny,lr_) will contain the local density at z(l,lr_) (partic
!..................................................................

      if (setup0%cqlpmod .ne. "enabled") then
        iorbstr=1
        iorbend=lz
      else
        iorbstr=l_
        iorbend=l_
      endif
      do 10 k=1,ngen
         ! xxx just assign it, sigh
        call dcopy(iyjx2,f(0:iy+1,0:jx+1,k,l_),1,temp3(0:iy+1,0:jx+1),1)
        do 11 l=iorbstr,iorbend
          ileff=l
          if (setup0%cqlpmod .eq. "enabled") ileff=ls_
          call cfpleg(0,ileff,1)
          do 40 ny=1,negyrg
            if (jegy(ny,1,k,lr_).eq.0 .or. eegy(ny,2,k,lr_).lt.1.e-15) &
              go to 40
            densz(l,k,ny,lr_)=0.
            do 6 j=jegy(ny,1,k,lr_),jegy(ny,2,k,lr_)
 6          densz(l,k,ny,lr_)=densz(l,k,ny,lr_)+cint2(j) &
                *tam1(j)*one_*twopi*2.
            xlndnz(k,ny)=xlndnz(k,ny)+densz(l,k,ny,lr_) &
              *dz(ileff,lr_)/bbpsi(ileff,lr_)
 40       continue
 11     continue
 10   continue


!      write(*,*) 'DENSZ(1,1,1:3,LR_),LR_',
!     +                 (DENSZ(1,1,NY,LR_),NY=1,NEGYRG),LR_

!..................................................................
!     Now determine local electron density for locquas="enabled"
!     calculations
!..................................................................

      if (kelecg.gt.0 .or. locquas .eq. "disabled") return
      nw=negyrg
      xlndnz(ngen+1,nw)=0.
      do 20 l=iorbstr,iorbend
        densz(l,ngen+1,nw,lr_)=0.0
        ileff=l
        if (setup0%cqlpmod .eq. "enabled") ileff=ls_
        do 30 k=1,ngen
          densz(l,ngen+1,nw,lr_)=densz(l,ngen+1,nw,lr_) &
            +densz(l,k,nw,lr_)*bnumb(k)
 30     continue
        xlndnz(ngen+1,nw)=xlndnz(ngen+1,nw) &
          +densz(l,ngen+1,nw,lr_)*dz(ileff,lr_) &
          /bbpsi(ileff,lr_)
 20   continue
      denn=0.
      do 200 ku=ngen+1,ntotal
        if (ku .eq. kelecm) go to 200
        denn=denn+reden(ku,lr_)
 200  continue
      xlndnz(ngen+1,nw)=xlndnz(ngen+1,nw)+denn*zmaxpsi(lr_)
      do 201 l=iorbstr,iorbend
        densz(l,ngen+1,nw,lr_)=densz(l,ngen+1,nw,lr_)+denn
 201  continue
      return
      end subroutine diagdenz


end module diagdenz_mod
