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

module urfdamp0_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use r8subs_mod, only : dcopy
  use urfdamp1_mod, only : urfdamp1
  use urfdamp2_mod, only : urfdamp2
  use urfdampa_mod, only : urfdampa

  !---END USE

!
!

contains

      subroutine urfdamp0(irfpwr,kopt)
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
!yup      save

!..................................................................
!     This routine calculates the power flowing along each ray
!     accounting for the damping, starting from an edge value
!     equal to (1/2)**irfpwr*(initialized input value).
!
!     kopt  .ne. 2: Modify local power according to power already
!     damped along the ray
!     =   2: Keep delpwr as initially given. Used for diagnostic rerun
!
!     Reminder:
!        mrf= number of rf "types"
!        mrfn= number of rf "modes" (sum over mrf of the nharms())
!        irfn(1:mrf)=  rf mode index (in 1:mrfn) of the lowest
!                      harmonic for each wave type
!        irfm(1:mrf)=  number of modes for each distinct wave
!                      type (where there may be more than one
!                      identical wave file applied to multiple
!                      species.
!        irffile(1:mrf)= "separate", "combine1", or "combine2",
!                      indicating separate, first of series of
!                      identical files applied to different
!                      species, continuation of the identical files
!                      applied to different species.
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

!yup      character*8 ifirst
!yup      data ifirst/"first"/
      call bcast(powrft(1:setup0%lrzmax),zero,setup0%lrzmax)
      call bcast(powurf(0:nmodsa),zero,nmodsa+1)
      call bcast(powurfc(0:nmodsa),zero,nmodsa+1)
      call bcast(powurfl(0:nmodsa),zero,nmodsa+1)
      call bcast(powrf( 1:lrza, 1:nmodsa),zero,lrza*nmodsa)
      call bcast(powrfc(1:lrza, 1:nmodsa),zero,lrza*nmodsa)
      call bcast(powrfl(1:lrza, 1:nmodsa),zero,lrza*nmodsa)
      call bcast(powurfi(0:lrza,0:nmodsa),zero,(lrza+1)*(nmodsa+1))

      call bcast(urfpwr,zero,nrayelts*nrayn*mrfn)
!yup      if (ifirst.eq."first") then
         call bcast(urfpwrc,zero,nrayelts*nrayn*mrfn)
         call bcast(urfpwrl,zero,nrayelts*nrayn*mrfn)
         ! urfpwrc, urfpwrl, scalurf are dimensioned as (1:nrayelts,1:nrayn,1:mrfn)
!yup      endif
      call bcast(scalurf,one, nrayelts*nrayn*mrfn)
!BH120815      call bcast(salphac,zero,nrayelts*nrayn*mrfn)

      do krfmode=1,mrfn ! loop over wave modes.
!..................................................................
!     Apply collisional or added linear damping to first harmonic
!     of separate wave types.
!     Indicate whether this wave mode, krfmode, corresponds to the
!     lowest harmonic for a ray type (lk=1, if yes).
!..................................................................
        lk=0
        do ktype=1,mrf ! loop over wave types.
           if (krfmode.eq.irfn(ktype)) lk=1
        enddo
!..................................................................
!     Get damping coefficients for all the ray elements
!..................................................................
        if (lk.eq.1) call urfdampa(krfmode) !collisional or added linear
             ! damping: urfpwrc(is,iray,krfmode) urfpwrl(is,iray,krfmode)
      enddo

!cc        call bcast(da(1,0),zero,iyjxp1) ! YuP: Why needed? No effect.
      if (n.ge.nondamp) then
!cc      do 6 ll=1,setup0%lrz ! Loop over flux surfaces
      do 5 krfmode=1,mrfn ! loop over wave modes.
!cc          call tdnflxs(lmdpln(ll)) !-> Determine l_, lr_, etc.
          if (n.gt.0 .and. urfdmp.eq."secondd") then
             call urfdamp2(krfmode)
          else
             call urfdamp1(krfmode)
          endif
 5    continue ! krfmode=1,mrfn ! loop over wave modes.
!cc 6    continue ! ll
      endif

#ifdef __MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfpwr, nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfpwrc,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfpwrl,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(scalurf,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(salphac,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
#endif


!yup      ifirst="notfirst"

!..................................................................
!
!     If nharms(1:mrf).gt.0, then power flow along delpwr for
!     the given wave type is calculated only for the lowest harmonic.
!     The resulting delpwr is copied to the harmonics
!     at end of subroutine.
!
!     Explanation of some of following variables:
!
!     delpwr(is=1:nrayelt(iray,krf),iray=1:nray(krf),krf=1:mrfn), where
!     mrfn=sum_over_wave_types(nharms), is the power flowing
!     along each ray (ergs/sec), for each mode krf.
!
!     urfpwr(is,iray,krf) is the fractional power absorbed by each
!     ray element, i.e., fraction of power flowing into the ray
!     element which is absorbed in the ray element.
!     In nharms=0,1 case, this is for each mode krf.
!     In nharms.gt.1 case, this is for each harmonic, krf=1:nharms.
!
!     urfpwrc(is,iray,krf) fractional collisional power absorption,
!     passed from ray tracing code.
!
!     urfpwrl(is,iray,krf) additional linear power absorption passed
!     from ray tracing code.
!
!     powrf(lr,krf) is the sum of urf power deposited from each mode
!     (or harmonic for nharms.gt.1) in each radial bin, divided
!     by bin volume. (watts/cm**3).
!
!     powrfc(lr,krf) collisional power/volume (watts/cm**3).
!
!     powrfl(lr,krf) additional linear power/volume (watts/cm**3).
!
!     powrft(lr) power per volume (watts/cm**3) summed over modes
!     or harmonics, dur to urf, collisional and add. linear abs.
!
!     powurf(0) total rf power summed over modes and radius
!     powurf(1:mrfn) total power summed over radius, for each
!     mode or harmonic.  (watts)
!
!     powurfi(0:lr,0:mrfn) partial powers integrated up to rho(lr).
!     powurf(krf)=powurfi(setup0%lrzmax,krf).    (watts)
!
!     powurfc(0) total coll. urf power summed over modes and radius
!     powurfc(1:mrfn) total coll. urf power summed over radius, for
!     each mode or harmonic. (watts)
!
!     powurfl(0:mrfn) ditto for urf linear power.
!..................................................................


!..................................................................
!     Calculate time-dependent pwrscalet, based on
!     nurftime,urftime(),pwrscale1()
!..................................................................

         if (nurftime.le.0) then
            pwrscalet=one
         else
            itme=1
            do jtm=1,nurftime
               if (timet.ge.urftime(jtm)) itme=jtm
            enddo
            itme1=itme+1
            if (itme.lt.nurftime) then
               pwrscalet=pwrscale1(itme) &
                    +(pwrscale1(itme1)-pwrscale1(itme)) &
                    /(urftime(itme1)-urftime(itme)) &
                    *(timet-urftime(itme))
            else
               pwrscalet=pwrscale1(nurftime)
            endif
         endif



      do 500 ktype=1,mrf ! loop over wave types.
         krf=irfn(ktype) !index of lowest harmonic in 1:mrfn (all modes)
                         !for each ktype
!         write(*,*)'urfdamp0: ktype, krf, pwrscale(ktype)',
!     +                        ktype, krf, pwrscale(ktype)
!.......................................................................
!     Update power flowing along the rays using run indicator irffile,
!     set in urfread.
!
!     If irffile()="separate" and nharms(krf).gt.0,
!     then sum power absorption over harmonics
!     for the given wave type.
!
!     If irffile()="combine1"
!     then sum power absorption over all modes,
!     and set pwr absorption data for irffile()="combine2" files
!     equal to the "combined1" data.
!.......................................................................

!        maxmodes=nharms(ktype)
!        if (irffile.eq.'combined') then
!           maxmodes=mrfn
!           if (ktype.ge.2) goto 501 !I.E., finished with the loop
!        endif ! if (irffile.eq.'combined')

        if (irffile(ktype).eq."separate") then
           maxmodes=nharms(ktype)
        elseif (irffile(ktype).eq."combine1") then
              maxmodes=irfm(ktype)
        elseif (irffile(ktype).eq."combine2") then
           goto 500  !i.e., finished with this distinct wave type
        endif


        if(kopt.ne.2)then ! modify local power, and determine nrayelt0
        pwrfrac=0.5**irfpwr*pwrscale(ktype)*pwrscalet
        powray=0.d0
        do iray=1,nray(krf)
          nrayelt0(iray,krf)=nrayelt(iray,krf)
          delpwr(1,iray,krf)=pwrfrac*delpwr0(iray,krf) ! is=1
          powray=powray+delpwr(1,iray,krf) ! Not used???
          if(nrayelt(iray,krf).ge.2)  then ! A normal long ray
            do is=2,nrayelt(iray,krf)
               ism=is-1
               urfpwrt=0.d0
               do kk=1,MAX(1,maxmodes) ! if maxmodes<1, then kk=1
                  kkk=krf+kk-1 !mode index
                  urfpwrt= urfpwrt + urfpwr(ism,iray,kkk) &
                          +urfpwrc(ism,iray,kkk)+urfpwrl(ism,iray,kkk)
               enddo ! kk=1,MAX(1,maxmodes)
               ! Redefine delpwr in this element through the previous:
               delpwr(is,iray,krf)=delpwr(ism,iray,krf)*(1.d0-urfpwrt)
               !!!!!!! here some of delpwr may become <0, because urfpwr>1.
               !(but normally, urfpwr<<1.)
               !write(*,*)'kopt=1: iray,is, urfpwrt=', iray,is, urfpwrt
               !if(urfpwrt.ge.1.d0) pause
               if(delpwr(is,iray,krf).le.1.d-6*delpwr(1,iray,krf))then
                 ! The left-over power in the ray element is negligible
                 delpwr(is:nrayelt(iray,krf),iray,krf)=em100
                 nrayelt0(iray,krf)=ism ! The last element where power is still not small
                 if(istarts(iray,krf).eq.1)  istarts(iray,krf)=-1
                 ! Finished with this ray, go to next ray:
                 go to 22  !!!!!!! this jump also happens when delpwr<0
               endif
            enddo ! is
 22       continue   ! for skipping all the rest ray-elements
          else ! A short(one point) ray :
            delpwr(1,iray,krf)=em100
          endif
        enddo ! iray
        endif ! kopt=1


        if(kopt.eq.2)then !do NOT modify delpwr, only determine nrayelt0
        do iray=1,nray(krf)
          nrayelt0(iray,krf)=nrayelt(iray,krf)
          if(nrayelt(iray,krf).ge.2)  then ! A normal long ray
            do is=2,nrayelt(iray,krf)
               ism=is-1
               if(delpwr(is,iray,krf).le.1.d-6*delpwr(1,iray,krf))then
                 ! The left-over power in the ray element is negligible
                 nrayelt0(iray,krf)=ism ! The last element where power is still not small
                 if(istarts(iray,krf).eq.1)  istarts(iray,krf)=-1
                 go to 23  ! Finished with this ray, go to next ray:
               endif
            enddo ! is
 23         continue   ! for skipping all the rest ray-elements
          endif ! if nrayelt(iray,krf)=1, do nothing
        enddo ! iray
        endif ! kopt=2


        ! Now sum-up deposited RF powers, sum over rays and ray-elements
        do 20  iray=1,nray(krf)
            do 15  is=1,nrayelt0(iray,krf) ! up to the element where delpwr is still not too small
              lrl=lloc(is,iray,krf) ! Only for ZOW:
              !ZOW: one ray-element contributes to only one radial bin.
              do 17 kk=1,MAX(1,maxmodes) ! if maxmodes<1, then kk=1
                kkk=krf+kk-1
                powrf(lrl,kkk)=powrf(lrl,kkk) &
                  +delpwr(is,iray,krf)*urfpwr(is,iray,kkk)*1.d-7
                !YuP[02/08/2015] Correction for powrfc and powrfl: is->is-1
                !Should there be a factor of scalurf present?
                powrfc(lrl,kkk)=powrfc(lrl,kkk) &
                  +delpwr(is,iray,krf)*urfpwrc(is,iray,kkk)*1.d-7
                powrfl(lrl,kkk)=powrfl(lrl,kkk) &
                  +delpwr(is,iray,krf)*urfpwrl(is,iray,kkk)*1.d-7
 17           continue
 15         continue ! is=1,nrayelt0(iray,krf)
 20     continue ! do iray=1,nray(krf)

        ! Total powers (sum over radial index)
        do 30 l=1,setup0%lrzmax
            do 19 kk=1,MAX(1,maxmodes) ! if maxmodes<1, then kk=1
              kkk=krf+kk-1
              powurfi(l,kkk)=powurfi(l-1,kkk)+powrf(l,kkk)
              powurfi(l,0)=powurfi(l,0)+powurfi(l,kkk)
              powrf(l,kkk)=powrf(l,kkk)/dvol(l)  !Converted to pwr dens
              powrft(l)=powrft(l)+powrf(l,kkk) !Pwr dens summed on modes
              !-> Coll.damping, if any:
              ! Total coll.power (integral over rho and sum over modes):
              powurfc(0)=powurfc(0)+powrfc(l,kkk)
              ! Integral over rho, for each mode (Watts):
              powurfc(kkk)=powurfc(kkk)+powrfc(l,kkk)
              powrfc(l,kkk)=powrfc(l,kkk)/dvol(l) ! watt/cm^3
              !-> Linear damping, if any:
              powurfl(0)=powurfl(0)+powrfl(l,kkk)
              powurfl(kkk)=powurfl(kkk)+powrfl(l,kkk)
              powrfl(l,kkk)=powrfl(l,kkk)/dvol(l)
 19         continue ! kk=1,maxmodes
 30     continue ! l=1,setup0%lrzmax
        do 29 kk=1,MAX(1,maxmodes) ! if maxmodes<1, then kk=1
          kkk=krf+kk-1
          powurf(kkk)=powurfi(setup0%lrzmax,kkk)
 29     continue ! kk=1,maxmodes
      powurf(0)=powurfi(setup0%lrzmax,0)

      if (maxmodes.gt.1) then
        do 600 kk=2,maxmodes
          kkk=krf+kk-1
          call dcopy(nrayn*nrayelts,delpwr(1:nrayelts,1:nrayn,krf),1, &
                                    delpwr(1:nrayelts,1:nrayn,kkk),1)
          do iray=1,nrayn
             nrayelt0(iray,kkk)=nrayelt0(iray,krf)
          enddo
 600    continue ! kk=1,maxmodes
      endif


 500  continue ! do ktype=1,mrf




 501  continue ! After end of do 500 ! (Not used)

!..................................................................
!     Print out min/max , for diagnostic purposes
!      do kkk=1,5
!      urfb_max=zero
!      urfb_min=ep100
!      do iray=1,nrayn
!         do is=1,nrayelts
!            yuri= urfpwr(is,iray,kkk)
!            urfb_max1=max(urfb_max,yuri)
!            urfb_min1=min(urfb_min,yuri)
!            if (urfb_max1.gt.urfb_max) then
!               urfb_max=urfb_max1
!            endif
!            if (urfb_min1.lt.urfb_min) then
!               urfb_min=urfb_min1
!            endif
!         enddo
!      enddo
!      write(*,*)'powurf-e->:',kkk,powurf(kkk),powurfl(kkk)
!      enddo
!...................................................................

      return
      end subroutine urfdamp0


end module urfdamp0_mod
