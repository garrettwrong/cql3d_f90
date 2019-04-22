module urfpack_mod


contains

      subroutine urfpack
      use param_mod
      use comm_mod
      use r8subs_mod, only : luf
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!.......................................................................
!     This routine determines arrays which contain data
!     distilled from the ray tracing routines. Each time the subroutine
!     urfb0 is called it uses and reuses the information computed
!     and stored in this routine.
!.......................................................................

      data nray0 /1/

!     Temporary printout, checking ray data is on radial mesh.
!      krf=1
!      do iray=nray0,nray(krf)
!      write(*,*)'urfpack:lloc(*,iray,1):',
!     1        (lloc(is,iray,1),is=1,nrayelt(iray,krf))
!      enddo

!.......................................................................
!     Uses library routines:
!     luf(x,table,n) (MATHLIB) which is a function returning the index
!        of the first element in the table that is greater than x.
!        Elements must be strictly increasing. x.gt.table(n)==>n+1.
!     lug(x,table,n,iguess) (MATHLIB) same as luf,
!        but with guess index iguess.
!     pack/unpack, as below.
!.......................................................................


      wx=.5
!**bh930729bcnst=charge/(fmass(kelecg)*clight)

!..................................................................
!     Loop over modes (or harmonics) and flux surfaces..
!..................................................................

      do 500 krf=1,mrfn
        write(*,*)'urfpack: krf(1:mrfn) = ', krf

!       k is general species to which this krf mode is to be applied:
        k=nrfspecies(krfn(krf))
        bcnst=abs(bnumb(k))*charge/(fmass(k)*clight)

        do 9 l=1,lrz
          call tdnflxs(lmdpln(l))

!..................................................................
!     Loop over rays
!..................................................................

          icount_outside_ez=0  ! only for printout
          icount_outside_er=0  ! only for printout

          do 10 iray=nray0,nray(krf)

!.......................................................................
!     Loop of ray elements - jump out if ray element does not contribute
!     to current flux surface or if data has been stored previously.
!.......................................................................

            do 20 is=1,nrayelt(iray,krf)

              if(lr_.ne.lloc(is,iray,krf)) go to 20

!..................................................................
!     lrayelt is 0 during the first call to urfpack.  If rays are
!     extended, lrayelt will contain the previous values of nrayelt.
!..................................................................

              if (lrayelt(iray,krf).ge.is) go to 20

!.......................................................................
!     The arrays ilim1 and ilim2 defined in this routine will be stored
!     in packed form in array ilowp and iupp. The number of bits
!     utilized from each number in ilim1 or ilim2 will be the rightmost
!     8 bits of each word. The routine pack accomplishes
!     compression. The variable "locatn" pinpoints the address
!     in the packed array where the first of the jx contributions
!     from a given ray element will be stored.
!     jjx is the first multiple of 8 greater than jx .
!     Similary, ifct1 and ifct2 are stored in packed form in ifct1_
!     and ifct2_, but in 16 bit (i.e., 2 byte) form. locatn16 specifies
!     the storage location in the compressed arrays.
!
!     From man pack on Cray C90 and J90 (after one or more "q"):
!      NAME
!           PACK - Compresses stored data
!      SYNOPSIS
!           CALL PACK(p, nbits, u, nw)
!      IMPLEMENTATION
!           Cray PVP systems
!      DESCRIPTION
!           PACK compresses stored data.  The following is a list
!           of valid arguments for this routine.
!           p         On exit, vector of packed data.
!           nbits     Number of rightmost bits of data in each
!                     partial word; must be 1, 2, 4, 8, 16, or 32.
!           u         Vector of partial words to be compressed.
!           nw        Number of partial words to be compressed.
!
!           PACK takes the 1, 2, 4, 8, 16, or 32 rightmost bits
!           of several partial words and concatenates them
!           into full 64-bit words.
!     The code has been modified to also treat the 32-bit interger
!        (machinea=2) case.  (BobH, 990308).
!     subroutines pack/unpack are for 8-bit words, pack16/unpack16
!     are for 16-bit words.
!.......................................................................

              locatn=  (jjx*(is-1)+jjx*nrayelts*(iray-1))/ibytes +1
              locatn16=(jjx*(is-1)+jjx*nrayelts*(iray-1))/ibytes16 +1

!.......................................................................
!     Determine the poloidal angle on the flux surface where the ray
!     element contributes.
!     rray is the major radius R location where it contributes and
!     zray is the Z location (distance above or below the midplane).
!..................................................................

            !-------------------------
            !YuP[04-2016]: adjust Zray,Rray if got outside of equil grid
            zray=wz(is,iray,krf)
            if(zray.gt.ez(nnz))then
              !For a mirror machine: ray can get outside of zbox
              !which defines the border of ez() equilibrium grid
              ![so that zbox=ez(nnz)-ez(1)]
              if (icount_outside_ez.eq.0) &
              write(*,*)'urfpack: Ray elements outside of ez grid'
              icount_outside_ez=icount_outside_ez+1 !for a printout
              write(*,*)'urfpack: zray>ez; iray,is,icount_outside_ez', &
                                           iray,is,icount_outside_ez
              ! Make an adjustment:
              zray=ez(nnz)
              !This correction is ok for a tokamak, too,
              !although not likely to happen.
            endif
            if(zray.lt.ez(1))then
              if (icount_outside_ez.eq.0) &
              write(*,*)'urfpack: Ray elements outside of ez grid'
              icount_outside_ez=icount_outside_ez+1 !for a printout
              write(*,*)'urfpack: zray<ez; iray,is,icount_outside_ez', &
                                           iray,is,icount_outside_ez
              ! Similarly, Make an adjustment:
              zray=ez(1)
            endif
            !-------------------------
            rray=wr(is,iray,krf)
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

!..................................................................
!     rmag is the location of the magnetic axis
!..................................................................

!BH091012:              rnew=rray-rmag
!BH091012:  For consistency with other thetapol defns
!BH091012:  (and non-updown symmetry).
!$$$             thetapol=atan(abs(zray)/abs(rnew))
!$$$              if (rnew.lt.0.) then
!$$$                thetapol=pi-thetapol
!$$$              endif
            if (eqsym.ne."none") then  !up-down symmetrize: abs(zz-zmag)
               thetapol=atan2(abs(zray-zmag),rray-rmag) ! in [0,pi]
               !but can be only [0; +pi/2) range in mirror machine
            else  !i.e., no up-down symmetrization
               thetapol=atan2(zray-zmag,rray-rmag)
               if (thetapol.lt.zero) thetapol=thetapol+twopi !in [0,pi2)
            endif

!..................................................................
!     Determine the relevant index "ll" on the computational poloidal
!     mesh - this will differ slightly from the exact value thetapol.
!..................................................................

              ll=luf(thetapol,pol(1:lz,lr_),lz)
              if (ll.gt.lz) then
!      write(*,*)
!     +'urfpack: thetapol,lr_,lz,pol(lz,lr_),pol(lz,lr_-1)=',
!     +          thetapol,lr_,lz,pol(lz,lr_),pol(lz,lr_-1)
                 if (abs(thetapol-pol(lz,lr_)).lt.em12*thetapol) then
                    ! just a small over-shooting
                    ll=lz-1 ! YuP: or maybe ll=lz ?
                    thetapol= pol(lz,lr_) !YuP: adjusted
                 else
                    if(machine.eq."mirror")then
                      ! Field lines can only go from 0 to pi/2,
                      ! or (normally) less. Ray element can be between
                      ! this given lr surface and the lower lr-1.
                      ! Try lr-1
                      if(lr_.eq.1)then !ray is inside 1st surface, near Z=Zmax
                         ll=lz
                         thetapol= pol(lz,lr_) !adjusted
                      else !(lr_.ge.2) then
                         !ll=luf(thetapol,pol(1,lr_-1),lz) ! try the lower lr
                         ! Or do same as for lr_=1 :
                         ll=lz
                         thetapol= pol(lz,lr_) !adjusted
                      endif
                      if(ll.gt.lz) then ! still could not find: stop the run
                         WRITE(*,*)'urfpack/urfwrong: lr_,ll,Rray,Zray', &
                             lr_,ll, rray,zray, &
                             solrz(lz,lr_),  solzz(lz,lr_), &
                             solrz(lz,lr_-1),solzz(lz,lr_-1)
                         call urfwrong(3)
                      endif
                    else ! tokamak
                      call urfwrong(3)
                    endif
                 endif
              endif ! ll>lz problem

              if (ll.gt.1) then
                if (pol(ll,lr_)-thetapol .gt. thetapol-pol(ll-1,lr_)) &
                  ll=ll-1
              endif
              llray(is,iray,krf)=ll

              signn=1.
              if (wnpar(is,iray,krf).lt.0) signn=-1.

              omomn=1.0
              signom=1.0

!..................................................................
!BH080920              if (nharm(krf).eq.0 .or. kiong(k).eq.1)  then
              if (nharm(krf).eq.0 .or. kiongg(k).eq.k)  then
!..................................................................

!..................................................................
!     Lower Hybrid/Fast Wave, or non-relativistic ions
!..................................................................


!..................................................................
!     Determine the local minimum and maximum parallel speed for
!     which we assume a resonance.
!     For nharm(krf)=0, if vparl has opposite sign from vparu, damping
!     will be be omitted.  Generally, this occurs when npar goes
!     through zero and hence the wave resonance is beyond the grid.
!     Sign reversal on vl,vu (and cosmz) is to accomodate luf.
!     Velocity integration limits are worked out using the
!     abs(wnpar).  Limits for negative wnpar can be obtained from
!     the abs(wnpar)-results by reflection about theta=pi/2.
!
!     The ion damping case is treated non-relativistically,
!     including the cyclotron damping.  The NR resonances are similar
!     to the electron case, in that the resonant velocity is
!     independent of perpendicular velocity.
!..................................................................

!BH080920                if (kiong(k).eq.1) then
                if (kiongg(k).eq.k) then
                  omn=nharm(krf)*bcnst*sbtot(is,iray,krf)/omega(krf)
                  omomn=1.0-omn
!990131                  signom=sign(1.,omomn)
                  signom=sign(one,omomn)
                endif

                vparu=clight*signom*omomn/ &
                  (signn*wnpar(is,iray,krf)-wdnpar(is,iray,krf)*wx)
                vparl=clight*signom*omomn/ &
                  (signn*wnpar(is,iray,krf)+wdnpar(is,iray,krf)*wx)

                if(vparl*vparu .gt. 0.)  then
                  vl=-abs(vparl/vnorm)
                  vu=-abs(vparu/vnorm)
                else
                  vl=-1.1
                  vu=-1.2
                endif

!..................................................................
!     Determine the lowest index in the velocity mesh that experiences
!     a resonance. sx is the normalized speed mesh, whereas x in the
!     momentum mesh.
!     jmin will be .gt.1.  If abs(vl) is off the mesh, jmin=jx+1
!..................................................................

                jmin=luf(-vl,sx,jx)
                jminray(is,iray,krf)=jmin
                jmaxray(is,iray,krf)=jx
                if (jmin.gt.jx) go to 20
                if (jmin.eq.1) call urfwrong(4)
                do 25 j=1,jmin-1
                  ilim2(j)=0
                  ilim1(j)=0
 25             continue

!..................................................................
!     Determine the lowest(highest) theta index for which there is a
!     resonance at speed mesh sx(j)
!..................................................................

                ilim2(jmin)=luf(vu/sx(jmin),cosmz(1:iyh,ll,lr_),iyh)
                ilim1(jmin)=luf(vl/sx(jmin),cosmz(1:iyh,ll,lr_),iyh)

                !YuP[04-2016] adjusted the range in theta around resonance
                !to have at least 3 points
!                i2=ilim2(jmin)
!                i1=ilim1(jmin)  ! i2<=i1
!                if( i2-i1 .eq. 0) then
!                   ! One point only i2=i1. Add two more points:
!                   if(i2.ge.2 .and. i1.le.iy-1) then
!                     i2=i2-1
!                     i1=i1+1
!                   endif
!                   if(i2.eq.1) then
!                     i2=i2
!                     i1=i1+2
!                   endif
!                   if(i1.eq.iy) then
!                     i2=i2-2
!                     i1=i1
!                   endif
!                endif
!                ilim2(jmin)=i2
!                ilim1(jmin)=i1
                !YuP[04-2016] done: adjusted
                call urfedge(ilim1(jmin),ilim2(jmin),vl/sx(jmin), &
                  vu/sx(jmin),ll,lr_,jmin)

                ilim2(jx)=luf(vu/sx(jx),cosmz(1:iyh,ll,lr_),iyh)
                ilim1(jx)=luf(vl/sx(jx),cosmz(1:iyh,ll,lr_),iyh)
                !YuP[04-2016] adjusted the range in theta around resonance
                !to have at least 3 points
!                i2=ilim2(jx)
!                i1=ilim1(jx)  ! i2<=i1
!                if( i2-i1 .eq. 0) then
!                   ! One point only i2=i1. Add two more points:
!                   if(i2.ge.2 .and. i1.le.iy-1) then
!                     i2=i2-1
!                     i1=i1+1
!                   endif
!                   if(i2.eq.1) then
!                     i2=i2
!                     i1=i1+2
!                   endif
!                   if(i1.eq.iy) then
!                     i2=i2-2
!                     i1=i1
!                   endif
!                endif
!                ilim2(jx)=i2
!                ilim1(jx)=i1
                !YuP[04-2016] done: adjusted
                call urfedge(ilim1(jx),ilim2(jx),vl/sx(jx),vu/sx(jx), &
                  ll,lr_,jx)

                iupjx=ilim2(jx)
                ilwjx=ilim1(jx)
                do 30 j=jmin+1,jx-1
                  ilim2(j)=luf(vu/sx(j),cosmz(ilim2(j-1): &
                    ilim2(j-1)+iupjx-ilim2(j-1),ll,lr_), &
                    iupjx-ilim2(j-1))-1+ilim2(j-1)
                  ilim1(j)=luf(vl/sx(j),cosmz(ilim1(j-1): &
                    ilim1(j-1)+ilwjx-ilim1(j-1),ll,lr_), &
                    ilwjx-ilim1(j-1))-1+ilim1(j-1)
                    !YuP[04-2016] adjusted the range in theta around resonance
                    !to have at least 3 points
!                    i2=ilim2(j)
!                    i1=ilim1(j)  ! i2<=i1
!                    if( i2-i1 .eq. 0) then
!                       ! One point only i2=i1. Add two more points:
!                       if(i2.ge.2 .and. i1.le.iy-1) then
!                         i2=i2-1
!                         i1=i1+1
!                       endif
!                       if(i2.eq.1) then
!                         i2=i2
!                         i1=i1+2
!                       endif
!                       if(i1.eq.iy) then
!                         i2=i2-2
!                         i1=i1
!                       endif
!                    endif
!                    ilim2(j)=i2
!                    ilim1(j)=i1
                    !YuP[04-2016] done: adjusted
                  call urfedge(ilim1(j),ilim2(j),vl/sx(j),vu/sx(j), &
                    ll,lr_,j)
 30             continue

!..................................................................
              elseif (nharm(krf).gt.0)  then
!..................................................................

!..................................................................
!     Relativistic cyclotron damping case.  rnpar1 is the lower n_parallel
!     and rnpar2 is the upper one.  The associated resonance velocities
!     (not momentum-per-mass) are vpar1 and vpar2.  ilim1(j) and ilim2(j)
!     are the lower and upper limits of the pitch angle index.
!
!     Also, if rnpar1 and rnpar2 have opposite signs, damping is again
!       omitted.  This can give inaccuracies for very close to
!       perpendicular propagation.  It can be mitigated by choosing
!       small wdnpar.
!..................................................................

                rnpar1=signn*wnpar(is,iray,krf)-wx*wdnpar(is,iray,krf)
                rnpar2=signn*wnpar(is,iray,krf)+wx*wdnpar(is,iray,krf)
                r1mn12=1.-rnpar1*rnpar1
                r1mn22=1.-rnpar2*rnpar2

                omn=nharm(krf)*bcnst*sbtot(is,iray,krf)/omega(krf)
                omn2=omn*omn
                rad1=omn2-r1mn12
                rad2=omn2-r1mn22
!
                if (rnpar1*rnpar2.le.0.) then
                  jminray(is,iray,krf)=jx+1
                  go to 20
                endif

!..................................................................
!     Ray outside of pinch point. No reson. (omega(krf) will be
!     .gt. omega_ce). Indicate skipping of the ray element by taking
!     jminray(is,iray,krf)=jx+1.
!..................................................................

                if (rad1.le.0.0.and.rad2.le.0.0)  then

!..................................................................
!     set flag t skip this ray element
!..................................................................

                  jminray(is,iray,krf)=jx+1
                  go to 20
                endif

!..................................................................
!     The following case should not be possible:
!..................................................................

                if (rad1.gt.0.0.and.rad2.lt.0.0)  call urfwrong(6)

!..................................................................
!     Note: We let uu01 and uu02 become negative for npar**2 > 1.
!           This helps the logic go through for this case,
!           with minimum changes from the npar**2<1 case.
!..................................................................

                if (rad2.gt.0.0.and.rad1.le.0.0)  then
                  upar01=0.0
                  uu01=0.0
                else
                  upar01=clight*omn*rnpar1/(r1mn12*vnorm)
                  uu01=clight*sqrt(rad1)/(r1mn12*vnorm)
                endif
                upar02=clight*omn*rnpar2/(r1mn22*vnorm)
                uu02=clight*sqrt(rad2)/(r1mn22*vnorm)

!..................................................................
!     upar11 and upar12 are parallel momentum-per-mass limits of the
!     resonance ellipse associated with rnpar1..., upar21 and upar22,
!     associated with rnpar2.
!..................................................................

                upar11=upar01-uu01
                upar12=upar01+uu01
                upar21=upar02-uu02
                upar22=upar02+uu02

!               Set rnpar1, rnpar2 .gt.1.0, for npar**2>1 case:
                if (rnpar1.gt.1.) upar12=1.1
                if (rnpar2.gt.1.) upar22=1.1

                jmin=luf(abs(upar21),x,jx)
                jminray(is,iray,krf)=jmin
                jmax=luf(upar22,x,jx)-1
                jmaxray(is,iray,krf)=jmax
!***  if(jmin.gt.jx.or.jmin.eq.jmax)  go to 20
!***  Try a less restricted condition, compatible with urfb0,...
                if(jmin.gt.jx)  go to 20

                if(jmin.eq.1)  call urfwrong(4)

                do 41  j=1,jmin-1
                  ilim1(j)=0
                  ilim2(j)=0
 41             continue
                do 42  j=jmax+1,jx
                  ilim1(j)=0
                  ilim2(j)=0
 42             continue

!..................................................................
!     Treat cases of (nharm(krf)*omega_ce/omega(krf)).le.1,
!     and .gt.1, separately
!..................................................................

                if (omn.le.1.)  then
!*bh*940227if(x(jmin).le.upar11)then
                  if(x(jmin).le.upar11 .or. rad1.le.0.0)  then
                    vres1=clight/vnorm ! c*bh*940227
                    ilim2(jmin)=1
                  else
                    vres1=clight*(1.-omn/gamma(jmin))/(rnpar1*vnorm)
                    ilim2(jmin)=luf(-vres1/sx(jmin), &
                        cosmz(1:iyh,ll,lr_),iyh)
                  endif
                  vres2=clight*(1.-omn/gamma(jmin))/(rnpar2*vnorm)
                  ilim1(jmin)=luf(-vres2/sx(jmin), &
                      cosmz(1:iyh,ll,lr_),iyh)
                  ii1=ilim1(jmin)
                  call urfedge(ii1,ilim2(jmin),-vres2/sx(jmin), &
                    -vres1/sx(jmin),ll,lr_,jmin)
                  do 45  j=jmin+1,jmax
                    if(x(j).le.upar11)  then
                      ilim2(j)=1
                      vres2=clight*(1.-omn/gamma(j))/(rnpar2*vnorm)
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1:iyh,ll,lr_),iyh)
                    elseif (x(j).lt.upar12)  then
                      vres1=clight*(1.-omn/gamma(j))/(rnpar1*vnorm)
                      vres2=clight*(1.-omn/gamma(j))/(rnpar2*vnorm)
                      ilim2(j)=luf(-vres1/sx(j),cosmz(1:iyh,ll,lr_),iyh)
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1:iyh,ll,lr_),iyh)
                    else
                      ilim2(j)=1
                      vres2=clight*(1.-omn/gamma(j))/(rnpar2*vnorm)
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1:iyh,ll,lr_),iyh)
                    endif
                    ii1=ilim1(j)
                    call urfedge(ii1,ilim2(j),-vres2/sx(j),-vres1/sx(j), &
                      ll,lr_,j)
 45               continue


                else

!...............................................................
!     nharm(krf)*omega_ce/omega(krf) .gt.1 - case:
!...............................................................

                  uperpstr=sqrt(omn*omn-1.)*clight/vnorm
                  if(x(jmin).le.abs(upar11))  then
!bh940227
                    vres1=clight*(1.-omn/gamma(jmin))/(rnpar1*vnorm)
                    ilim1(jmin)=iy
                  else
                    vres1=clight*(1.-omn/gamma(jmin))/(rnpar1*vnorm)
                    ilim1(jmin)=luf(-vres1/sx(jmin), &
                        cosmz(1:iy,ll,lr_),iy)
                  endif
                  vres2=clight*(1.-omn/gamma(jmin))/(rnpar2*vnorm)
                  ilim2(jmin)=luf(-vres2/sx(jmin), &
                    cosmz(1:iy,ll,lr_),iy)
                  ii1=ilim1(jmin)
                  call urfedge(ii1,ilim2(jmin),-vres1/sx(jmin), &
                    -vres2/sx(jmin),ll,lr_,jmin)
                  do 55  j=jmin+1,jmax
                    vres1=clight*(1.-omn/gamma(j))/(rnpar1*vnorm)
                    vres2=clight*(1.-omn/gamma(j))/(rnpar2*vnorm)
                    if(x(j).le.abs(upar11))  then
                      ilim1(j)=iy
                      ilim2(j)=luf(-vres2/sx(j), &
                       cosmz(1:iy,ll,lr_),iy)
                      ii1=ilim1(j)
                      call urfedge(ii1,ilim2(j),-vres1/sx(j), &
                        -vres2/sx(j),ll,lr_,j)
                    elseif (x(j).lt.uperpstr)  then
                      ilim2(j)=luf(-vres2/sx(j), &
                       cosmz(1:iy,ll,lr_),iy)
                      ilim1(j)=luf(-vres1/sx(j),cosmz(1:iy,ll,lr_),iy)
                      ii1=ilim1(j)
                      call urfedge(ii1,ilim2(j),-vres1/sx(j), &
                        -vres2/sx(j),ll,lr_,j)
                    elseif(x(j).lt.upar12)  then
                      ilim2(j)=luf(-vres1/sx(j),cosmz(1:iy,ll,lr_),iy)
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1:iy,ll,lr_),iy)
                      ii1=ilim1(j)
                      call urfedge(ii1,ilim2(j),-vres2/sx(j), &
                        -vres1/sx(j),ll,lr_,j)
                    elseif (x(j).lt.upar22)  then
                      ilim2(j)=1
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1:iy,ll,lr_),iy)
                      ii1=ilim1(j)
                      call urfedge(ii1,ilim2(j),-vres2/sx(j), &
                        -vres1/sx(j),ll,lr_,j)
                    endif


 55               continue

                endif

!..................................................................
!     End of nharm(krf) if construction:
!..................................................................

              endif

!..................................................................
!     Check for errors
!..................................................................

              do 70 j=1,jx
                if(ilim1(j).gt.iy) call urfwrong(7)
                if(ilim2(j).gt.iy) call urfwrong(7)
 70           continue

!..................................................................
!     Reflect about pi/2 for negative k_par
!..................................................................

              if (signn*signom.lt.0.) then
                do 80 j=1,jx
                  iota(j)=ilim2(j)
 80             continue
                do 90 j=1,jx
                  ilim2(j)=iy+1-ilim1(j)
                  ilim1(j)=iy+1-iota(j)
 90             continue
                do 81 j=1,jx
                  iota(j)=ifct1(j)
 81             continue
                do 91 j=1,jx
                  ifct1(j)=ifct2(j)
                  ifct2(j)=iota(j)
 91             continue
              endif

!..................................................................
!     Pack results in 1 byte chunks (ilowp and iupp) and
!     2 byte chunks (ifct1_ and ifct2_) to save space.
!..................................................................

              !if(urfb_version.eq.1)then ! 2 is the new version developed by YuP
                ! if 1, it will use the original version
                call pack(ilowp(locatn,krf),8,ilim1,jjx)
                call pack(iupp(locatn,krf),8,ilim2,jjx)
                call pack16(ifct1_(locatn16,krf),8,ifct1,jjx)
                call pack16(ifct2_(locatn16,krf),8,ifct2,jjx)
              !endif

!     Temporary printout, checking elements fall on vel grid:
!              write(*,*)'urfpack: is,iray,krf,jminray,jmaxray',
!     1             is,iray,krf,jminray(is,iray,krf),jmaxray(is,iray,krf)


 20         continue ! is=1,nrayelt(iray,krf)
 10       continue ! iray=nray0,nray(krf)

!.......................................................................
!     end of loop over flux surfaces
!.......................................................................
 9      continue !  l=1,lrz

        do 50 iray=1,nray(krf)
          lrayelt(iray,krf)=nrayelt(iray,krf)
 50     continue

 500  continue     !end of loop over modes krf=1,mrfn


!$$$      do is=1,nrayelt(1,1)
!$$$         jmin=jminray(is,1,1)
!$$$         jmax=jmaxray(is,1,1)
!$$$         write(*,1000) is,lloc(is,1,1),jmin,jmax
!$$$      enddo
!$$$ 1000 format(10i5)
!$$$      write(*,*)'urfpack: STOP'
!$$$      STOP


      !if(urfb_version.eq.1)then ! 2 is the new version developed by YuP
         ! if 1, it will use the original version
          write(*,*)'URFpacked ifct1 to ifct1_; Sizes=', &
                     size(ifct1),size(ifct1_)
      !endif

      write(*,*)'urfpack:  END'
      return
      end
end module urfpack_mod
