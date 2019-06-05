module bavgmax_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use psif_mod, only : psif

  !---END USE

!
!

contains

      subroutine bavgmax
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!.................................................................
!     This routine computes the bounce average of several bbpsi
!     related quantities.
!.................................................................
      save


      if (n .ge. 2) return
      call bcast(psiiv(1:iy,lr_),zero,iy)
      call bcast(psiba(1:iy,lr_),zero,iy)
      call bcast(psisq(1:iy,lr_),zero,iy)
      call bcast(psiir(1:iy,lr_),zero,iy)
      call bcast(vderb(1:iy,lr_),zero,iy)
      call bcast(psicu(1:iy,lr_),zero,iy)
      call bcast(psiqu(1:iy,lr_),zero,iy)
      call bcast(sincosba(1:iy,1:ngen,lr_),zero,iy*ngen)


      lrange=lz
      if (numclas .eq. 1) lrange=lz/2+1 !!! YuP: needs to check ?


      do 30 l=1,lrange
        do 31 i=1,imax(l,lr_)
          iii=iy+1-i

          ! 1. All orbits that reach/pass a given poloidal node l:
          ! passing (i<itl), or trapped that could reach/pass this l;
          ! also includes Last Trapped, with tip at l=lz_bmax
          !(for such particle, lmax(itl)=lz_bmax; see micxinil)
          ztt=bbpsi(l,lr_)
          dtau_tau=dtau(i,l,lr_)/tau(i,lr_)
          vddert=sin(pol(l,lr_))*(2.+ztt*sinn(i,l_)**2)
          !psiir(i,lr_)=psiir(i,lr_)+sqrt(ztt)*dtau_tau    !-YuP Not used; ztt or 1/ztt?
          psiiv(i,lr_)=psiiv(i,lr_)+(1./ztt)*dtau_tau
          psiba(i,lr_)=psiba(i,lr_)+ztt*dtau_tau
          psisq(i,lr_)=psisq(i,lr_)+ztt**2*dtau_tau
          psicu(i,lr_)=psicu(i,lr_)+ztt**3*dtau_tau
          psiqu(i,lr_)=psiqu(i,lr_)+ztt**4*dtau_tau
          vderb(i,lr_)=vderb(i,lr_)+vddert*dtau_tau
          !write(*,*)'bavg oo',lr_,l,i,dtau(i,l,lr_)

          ! 2a. Trapped, with tips between l and l+1 (ABOVE midplane):
          if (l.eq.lmax(i,lr_) .and. l.lt.lz_bmax(lr_)) then
            ! Add contribution from orbit tip:
            ztt=psif(zboun(i,lr_))
            dtau_tau=dtau(i,l+1,lr_)/tau(i,lr_)
            vddert=sin(pi*zboun(i,lr_)/z_bmax(lr_)) !-YuP  Was /zmax()
            !psiir(i,lr_)=psiir(i,lr_)+sqrt(ztt)*dtau_tau    !-YuP Not used; ztt or 1/ztt?
            psiiv(i,lr_)=psiiv(i,lr_)+(1./ztt)*dtau_tau
            psiba(i,lr_)=psiba(i,lr_)+ztt*dtau_tau
            psisq(i,lr_)=psisq(i,lr_)+ztt**2*dtau_tau
            psicu(i,lr_)=psicu(i,lr_)+ztt**3*dtau_tau
            psiqu(i,lr_)=psiqu(i,lr_)+ztt**4*dtau_tau
            vderb(i,lr_)=vderb(i,lr_)+vddert*dtau_tau
            do  k=1,ngen
              sincost=0.0
              if(paregy(k).le.0.) sincost=1.
              sincosba(i,k,lr_)=sincosba(i,k,lr_)+sincost*dtau_tau
              sincosba(iii,k,lr_)=sincosba(i,k,lr_)
            enddo
            !write(*,*)'bavg <<',lr_,l,i,dtau(i,l+1,lr_)
          endif

          ! 2b. Trapped, with tips between l and l-1 (BELOW midplane):
          if (l.eq.lmax(i+iyh,lr_) .and. l.gt.lz_bmax(lr_)) then
            ! Add contribution from orbit tip:
            ztt=psif(zboun(i,lr_))
            dtau_tau=dtau(i,l-1,lr_)/tau(i,lr_)
            vddert=sin(pi*zboun(i,lr_)/z_bmax(lr_)) !-YuP  Was /zmax()
            !psiir(i,lr_)=psiir(i,lr_)+sqrt(ztt)*dtau_tau    !-YuP Not used; ztt or 1/ztt?
            psiiv(i,lr_)=psiiv(i,lr_)+(1./ztt)*dtau_tau
            psiba(i,lr_)=psiba(i,lr_)+ztt*dtau_tau
            psisq(i,lr_)=psisq(i,lr_)+ztt**2*dtau_tau
            psicu(i,lr_)=psicu(i,lr_)+ztt**3*dtau_tau
            psiqu(i,lr_)=psiqu(i,lr_)+ztt**4*dtau_tau
            vderb(i,lr_)=vderb(i,lr_)+vddert*dtau_tau
            do  k=1,ngen
              sincost=0.0
              if(paregy(k).le.0.) sincost=1.
              sincosba(i,k,lr_)=sincosba(i,k,lr_)+sincost*dtau_tau
              sincosba(iii,k,lr_)=sincosba(i,k,lr_)
            enddo
            !write(*,*)'bavg <<',lr_,l,i,dtau(i,l-1,lr_)
          endif

          ! Symmetrize around pitch=pi/2
!          psiir(iii,lr_)=psiir(i,lr_)    !-YuP Not used
          psiiv(iii,lr_)=psiiv(i,lr_)
          psiba(iii,lr_)=psiba(i,lr_)
          psisq(iii,lr_)=psisq(i,lr_)
          psicu(iii,lr_)=psicu(i,lr_)
          psiqu(iii,lr_)=psiqu(i,lr_)
          vderb(iii,lr_)=vderb(i,lr_)

 31     continue ! i
 30   continue ! l



      do 70 i=2,iyh
        iii=iy+1-i
        swwswwsw=psiiv(i,lr_)/sinn(i,l_)**2-1.
        batot(i,lr_)=tann(i,l_)**2*swwswwsw
        batot(iii,lr_)=batot(i,lr_)
 70   continue
      batot(1,lr_)=psiiv(1,lr_)
      batot(iy,lr_)=batot(1,lr_)

!      do i=1,iy
!         write(*,*)lr_,i,tau(i,lr_)
!      enddo
!      pause
!         write(*,*)'bavgmax psimx=',lr_,psimx(lr_)
!         if(lr_.eq.1) pause
      return
      end subroutine bavgmax
      
      
end module bavgmax_mod
