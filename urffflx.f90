module urffflx_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use bcast_mod, only : ibcast
  use zcunix_mod, only : terp2

  !---END USE

!
!

contains

      subroutine urffflx
      use param_mod
      use comm_mod
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



      do 1 l=1,lrzmax
 1    tr2(l)=psimag-psivalm(l)

!..................................................................
!     According to the psi value determined for a particular ray element
!     assign it to a given flux surface.
!..................................................................

      call ibcast(ncontrib(1),0,lrzmax)
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
            l=luf(apsi,tr2(1),lrzmax)
!BH090602   Ray elements outside LCFS (rho=1 surface) will be attributed to lrzmax
            if (l.gt.lrzmax) then
               if (icount_outside_lim.eq.0) &
               write(*,*)'urffflx: Ray elements outside of rho=1'
               icount_outside_lim=icount_outside_lim+1 !for a printout
               write(*,'(a,i4,2i7)') &
                   'urffflx:l>lrzmax; iray,is,icount_outside_lim', &
                                      iray,is,icount_outside_lim
               l=lrzmax ! Adjusted
            endif
!$$$            if (l.le.0) then
!$$$               write(*,*)'urffflx:l,lrzmax,k,iray,is',l,lrzmax,k,iray,is
!$$$               go to 30
!$$$            endif
            lloc(is,iray,irfn(krf))=l
!BH090602            if (l.gt.lrzmax) go to 30
            ncontrib(l)=ncontrib(l)+1
 30       continue
 20     continue
 100  continue !  krf=1,mrf

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
!      write(*,*)'urffflx:ncontrib(1:lrzmax):',(ncontrib(l),l=1,lrzmax)

      return
      end
end module urffflx_mod
