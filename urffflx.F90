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

module urffflx_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : ibcast
  use zcunix_mod, only : terp2

  !---END USE
#ifdef __MPI
  include 'cql3d_mpilib.h'
#endif

!
!

contains

      subroutine urffflx
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : luf
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!..................................................................
!     Takes ray data from ray tracing codes and segregates
!     information by flux surface.
!     Three arrays are calculated:
!     -ncontrib(l) will hold the number of ray elements that contribute
!      to the diffusion coefficient for flux surface l, i.e., in vol(l).
!      lloc(is,iray,krf) is the index of the flux surface to which ray
!      element (is,iray,krf) contributes.  The counting is restricted
!      to wave types, and does not include harmonics.
!     -psiloc(,,) is the value of pol. flux psi associated with the
!      ray element;     krf is the index of the excitation mode.
!     -lloc(,,) is the radial mesh rho bin.
!     NOTE: psivalm(l) demarks the outer edge of the flux volume associated
!      with flux surface l (vol(l)).
!..................................................................



      do 1 l=1,setup0%lrzmax
 1    tr2(l)=psimag-psivalm(l)

!..................................................................
!     According to the psi value determined for a particular ray element
!     assign it to a given flux surface.
!..................................................................

      call ibcast(ncontrib(1),0,setup0%lrzmax)
      icount_outside_lim=0 ! only for printout
      icount_outside_ez=0  ! only for printout
      icount_outside_er=0  ! only for printout
      do 100 krf=1,mrf
        do 20 iray=1,nray(irfn(krf))
          do 30 is=lrayelt(iray,irfn(krf))+1,nrayelt(iray,irfn(krf))
            !-------------------------
            zray=wz(is,iray,irfn(krf))
            if(zray.gt.ez(nnz))then
              !For a mirror machine: ray can get outside of zbox
              !which defines the border of ez() equilibrium grid
              ![so that zbox=ez(nnz)-ez(1)]
              if (icount_outside_ez.eq.0) &
              write(*,*)'urffflx: Ray elements outside of ez grid'
              icount_outside_ez=icount_outside_ez+1 !for a printout
              write(*,'(a,i4,2i7)') &
                  'urffflx: zray>ez; iray,is,icount_outside_ez', &
                                     iray,is,icount_outside_ez
              ! Make an adjustment:
              zray=ez(nnz)
              !This correction is ok for a tokamak, too,
              !although not likely to happen.
            endif
            if(zray.lt.ez(1))then
              if (icount_outside_ez.eq.0) &
              write(*,*)'urffflx: Ray elements outside of ez grid'
              icount_outside_ez=icount_outside_ez+1 !for a printout
              write(*,'(a,i4,2i7)') &
                  'urffflx: zray<ez; iray,is,icount_outside_ez', &
                                     iray,is,icount_outside_ez
              ! Similarly, Make an adjustment:
              zray=ez(1)
            endif
            !-------------------------
            rray=wr(is,iray,irfn(krf))
            if(rray.gt.er(nnr))then
              !For a mirror machine: ray can get outside of
              !er() equilibrium grid
              ![so that zbox=ez(nnz)-ez(1)]
              ! Make an adjustment:
              rray=er(nnr)
              !This correction is ok for a tokamak, too,
              !although not likely to happen.
            endif
            if(rray.lt.er(1))then ! this cannot happen, but ok to add.
              ! Similarly, Make an adjustment:
              rray=er(1)
            endif
            !-------------------------
            psiloc(is,iray,irfn(krf))= terp2(rray,zray,nnr,er,nnz,ez, &
                                     epsi,epsirr,epsizz,epsirz,nnra,0,0)
            apsi=psimag-psiloc(is,iray,irfn(krf))
            l=luf(apsi,tr2(1),setup0%lrzmax)
!BH090602   Ray elements outside LCFS (rho=1 surface) will be attributed to setup0%lrzmax
            if (l.gt.setup0%lrzmax) then
               if (icount_outside_lim.eq.0) &
               write(*,*)'urffflx: Ray elements outside of rho=1'
               icount_outside_lim=icount_outside_lim+1 !for a printout
#ifdef __MPI
               ! for MPI, only print from master
               if(mpirank.eq.0) then
#endif
                  if(icount_outside_lim.LT.10) then
                     write(*,'(a,i4,2i7)') &
                          'urffflx:l>setup0%lrzmax; iray,is,icount_outside_lim', &
                          iray,is,icount_outside_lim
                  end if
                  if(icount_outside_lim.EQ.10) then
                     write(*,'(a)') &
                          'urffflx:l>setup0%lrzmax; GT 10 elem outside, quieting...'
                  end if
#ifdef __MPI
         endif !mpirank 0
#endif
               l=setup0%lrzmax ! Adjusted
            endif
!$$$            if (l.le.0) then
!$$$               write(*,*)'urffflx:l,setup0%lrzmax,k,iray,is',l,setup0%lrzmax,k,iray,is
!$$$               go to 30
!$$$            endif
            lloc(is,iray,irfn(krf))=l
!BH090602            if (l.gt.setup0%lrzmax) go to 30
            ncontrib(l)=ncontrib(l)+1
30       end do
20    end do
100 end do !  krf=1,mrf

! if we "quietted" log, report totals now
#ifdef __MPI
         ! for MPI, only print from master
         if(mpirank.eq.0) then
#endif
            if(icount_outside_lim.GT.0) then
               write(*,'(a,i4)')  &
                    'urffflx:l>setup0%lrzmax; total icount_outside_lim: ', &
                    icount_outside_lim
            end if
#ifdef __MPI
         endif !mpirank 0
#endif

!     Duplicate data for psiloc and lloc into multi-harmonic
!     related arrays.
      do 120 krf=1,mrf
      if (nharms(krf).gt.1) then
        do 110  i=1,nharms(krf)-1
          do 21  iray=1,nray(irfn(krf))
            do 31  is=lrayelt(iray,irfn(krf))+1,nrayelt(iray,irfn(krf))
              psiloc(is,iray,irfn(krf)+i)=psiloc(is,iray,irfn(krf))
              lloc(is,iray,irfn(krf)+i)=lloc(is,iray,irfn(krf))
 31         continue
 21       continue
 110    continue
      endif
 120  continue !  krf=1,mrf

!     Temporary print out checking number of elements at each rad bin:
!      write(*,*)'urffflx:ncontrib(1:setup0%lrzmax):',(ncontrib(l),l=1,setup0%lrzmax)

      return
      end subroutine urffflx

end module urffflx_mod
