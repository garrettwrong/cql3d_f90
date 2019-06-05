module urfbes_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use bcast_mod, only : bcast
  use zcunix_mod, only : zzbeslri

  !---END USE


!
!

contains

      subroutine urfbes
      use param_mod
      use cqlcomm_mod
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     Routine computes three Bessel tables, J_nharm-1,J_nharm, and
!     J_nharm+1. Space is allocated. Once tables are created, subsequent
!     calls check to see if the maximum argument has increased beyond
!     table limits. If so, table is augmented.
!..................................................................


      pointer besl
      dimension besl(:) ! Now it's a local working array

      nharmx = 1
      do krf=1,mrfn
         nharmx=MAX(nharmx,nharm(krf))
      enddo
      write(*,*)'urfbes:  Before allocation, nharmx=',nharmx
      allocate(besl(nharmx+2),STAT=istat) !-YuP->added
      !XXXcall bcast(besl,0.0,SIZE(besl))     !-YuP->added
      besl = 0.d0
      allocate(jbm1(nbssltbl,mrfn),STAT=istat) !-YuP->added
      !XXXcall bcast(jbm1,0.0,SIZE(jbm1))          !-YuP->added
      jbm1 = 0.d0 ! for all (1:nbssltbl,1:mrfn)
      allocate(jb0(nbssltbl,mrfn),STAT=istat)  !-YuP->added
      !XXXcall bcast(jb0,0.0,SIZE(jbm1))           !-YuP->added
      jb0 = 0.d0 ! for all (1:nbssltbl,1:mrfn)
      allocate(jbp1(nbssltbl,mrfn),STAT=istat) !-YuP->added
      !XXXcall bcast(jbp1,0.0,SIZE(jbm1))          !-YuP->added
      jbp1 = 0.d0 ! for all (1:nbssltbl,1:mrfn)
      write(*,*)'urfbes:  After allocation'

!..................................................................
!     Loop over excitation modes
!..................................................................

      do 200 krf=1,mrfn
!         write(*,*)'urfbes: krf(1:mrfn) = ',krf

!..................................................................
!     Find maximum k-perp
!..................................................................

        call aminmx(wnper(1:nrayelts,1:nrayn,krf), &
              1,nrayelts*nrayn,1,emin,emax, &
              kmin,kmax)
        em=emax
        if (abs(emin).gt.emax) em=abs(emin)
        xkp=em*freqcy(krf)*twopi/clight

!..................................................................
!     Guess at minimum cyclotron frequency - leave bpol out of estimate
!     this will decrease guess.
!..................................................................

!        bvalmin=fpsiz(lrindx(lrors))/er(nnr)
!c**bh930729wcemin=charge*bvalmin/(fmass(kelecg)*clight)
!        wcemin=abs(bnumb(1))*charge*bvalmin/(fmass(1)*clight)
!**bh970620  Using ray data to get bvalmin:
!**bh970620  (This was necessary, to avoid neg. tor. B problem with
!**bh970620   the RFP).

        bvalmin=ep100
        do iray=1,nray(krf)
           do is=1,nrayelt(iray,krf)
!990131              bvalmin=amin1(bvalmin,abs(sbtot(is,iray,krf)))
              bvalmin=min(bvalmin,abs(sbtot(is,iray,krf)))
           enddo
!           write(*,*)'urfbes: iray,krf,sbtot(is,iray,krf) = ',
!     1             iray,krf,(sbtot(is,iray,krf),is=1,nrayelt(iray,krf))
        enddo
!       k is species to which this krf mode is to be applied:
        k=nrfspecies(krfn(krf))
        wcemin=abs(bnumb(k))*charge*bvalmin/(fmass(k)*clight)

!..................................................................
!     Maximum argument in table will be argmax - give it a 40% boost
!..................................................................

        argmax=xkp*vnorm/wcemin*1.4
        write(*,*)'urfbes: argmax,nharm(krf),bvalmin = ', &
                           argmax,nharm(krf),bvalmin
        bsslstp(krf)=argmax/(nbssltbl-1)
        arg=0.
        do 10 i=1,nbssltbl
          call zzbeslri(arg,nharm(krf)+2,0,besl,ncalc)
!BH150620: Having problem with bessel calc.
!BH150620: Try letting simple accuracy loss get by ok
!BH150620:
           if (ncalc.ne.nharm(krf)+2) then
             write(*,*)'zzbeslri accuracy loss: arg,nharm(krf)+2,ncalc=' &
                  ,arg,nharm(krf)+2,ncalc
             if (ncalc.lt.0) then
                stop ' in urfbes; ncalc.lt.0'
             endif
          endif
          if(nharm(krf).eq.0)  then
            jbm1(i,krf)=-besl(2)
            jb0(i,krf)=besl(1)
            jbp1(i,krf)=besl(2)
          elseif (nharm(krf).gt.0)  then
            jbm1(i,krf)=besl(nharm(krf))
            jb0(i,krf)=besl(nharm(krf)+1)
            jbp1(i,krf)=besl(nharm(krf)+2)
          endif
          arg=arg+bsslstp(krf)
 10     continue
        argmax=arg-bsslstp(krf)

 200  continue ! krf=1,mrfn

      write(*,*)'End urfbes'

      return
      end subroutine urfbes
      
end module urfbes_mod
