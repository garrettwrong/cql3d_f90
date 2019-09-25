!
!
      subroutine frstup(mf,mfm1,mi,mj,nion,potsid,codeid,rin, &
        rmax,zax,zminn,zmaxx,zni,zne,zte,zti,zzi,frpsi,psivol, &
        xxx,yyy,nprim,nimp,zeffctv,zshift1)

      use bcast_mod, only : bcast
      use cqlcomm_mod
      use cqlconf_mod, only : setup0
      use param_mod

      implicit none
      integer :: mf,mfm1,mi,mj,nion,nprim,nimp ! input and output
      character*8 codeid ! output
      real(c_double) :: psivol(*),potsid(*),frpsi(nnra,*) ! output
      real(c_double) :: zne(kz,*),zni(kz,*),zte(*),zzi(kz,*) ! input
      real(c_double) :: xxx(*),yyy(*) ! output
      real(c_double) :: zti(*),zeffctv(*) ! input
      real(c_double) :: rin,rmax,zminn,zmaxx,zax,zshift1 ! output
      integer :: i,j,l,ll,k,kk,kimp_ ! local
      real(c_double) :: tp,dn ! local

      ! YuP[2019-06-12] Getting lrz and lrzmax from setup0 type:
      integer :: lrz
      integer :: lrzmax
      lrz = setup0%lrz
      lrzmax =  setup0%lrzmax

!..................................................................
!     This routine should be called just before the call to FREYA.
!     It defines a number of time (or iteration) dependent input
!     variables required in the call to FREYA.
!..................................................................


!..................................................................
!     Set the frpsi (epsi) array as required by FREYA
!     Maximum will be at limiter not the axis.
!..................................................................

      do 10 i=1,nnr
        xxx(i)=er(i)
        do 11 j=1,nnz
          frpsi(i,j)=-epsi(i,j)
 11     continue
 10   continue
      do 15 j=1,nnz
        yyy(j)=ez(j)
 15   continue
!.....................................................................
!     eqpsi(1:nconteqn) is set up in cql3d routine eqrhopsi.f.
!     It is used to give the flux zones for freya returned in the
!     izone argument of subroutine inject (now inject_old and inject1),
!     called in subroutine freya.
!     Change the sign of (positive) eqpsi array to get it in ascending
!     order for splines.  (Changed back at end of subroutine).
!.....................................................................

      do 2 l=1,nconteqn
        eqpsi(l)=-eqpsi(l)
 2    continue
!..................................................................
!     The number of flux zones (mfm1) used in FREYA
!..................................................................

      mfm1=nconteqn-1
      mf=mfm1+1

      write(*,*)
      write(*,*)'Number of flux zones (mfm1) used in FREYA=',mfm1
      write(*,*)'Adjust mfm1=nconteqn-1 throught NL input of nconteqn'
      write(*,*)

!..................................................................
!     Interpolate densities, temperatures etc over to the
!     (high resolution?) mesh utilized by FREYA.
!
!..................................................................

      psivol(1)=.5*(eqvol(2)+eqvol(1))
      psivol(mfm1)=eqvol(mf)-.5*(eqvol(mfm1)+eqvol(mfm1-1))
      do 20 ll=2,mfm1-1
        psivol(ll)=.5*(eqvol(ll+1)-eqvol(ll-1))
 20   continue

!..................................................................
!     Electron density, energy
!..................................................................

      do 30 ll=1,lrzmax
        tr1(ll)=reden(kelec,ll)
        tr2(ll)=energy(kelec,ll)*2./3.
 30   continue

!..................................................................
!     Call spline interpolation routine to get values on FREYA
!     eqspi mesh
!     (which is same as intemediary radial mesh used in CQL3D proper).
!..................................................................

      call frsplft(lrzmax,equilpsp(1),tr1(1),mfm1,eqpsi,zne)
      call frsplft(lrzmax,equilpsp(1),tr2(1),mfm1,eqpsi,zte)

!..................................................................
!     Average ion temperature..
!..................................................................

      kimp_=ntotal-nimp+1
      tr1=zero !YuP[2019-06-08]was call bcast(tr1(1),zero,lrzmax)
      do 50 ll=1,lrzmax
        tp=0.
        dn=0.
        kk=0
        do 40 k=1,ntotal
          if (k.eq.kelecg .or. k.eq.kelecm) go to 40
          kk=kk+1
          if (kk.gt.nprim .and. k.lt.kimp_) then
            kk=kk-1
            go to 40
          endif
          tp=tp+2./3.*energy(k,ll)*reden(k,ll)
          dn=dn+reden(k,ll)
 40     continue
        tr1(ll)=tp/dn
 50   continue
      call frsplft(lrzmax,equilpsp(1),tr1(1),mfm1,eqpsi,zti(1))

!..................................................................
!     Ion densities..
!..................................................................
!
      kk=0
      do 70 k=1,ntotal
        tr1=zero !YuP[2019-06-08]was call bcast(tr1(1),zero,lrzmax)
        if (k.eq.kelecm .or. k.eq.kelecg) go to 70
        kk=kk+1
!        if (kk.gt.nprim .and. k.lt.kimp_) then
        if (kk.gt.nprim+nimp+1 .and. k.lt.kimp_) then
          kk=kk-1
          go to 70
        endif
        do 60 ll=1,lrzmax
          tr1(ll)=reden(k,ll)
!         write(*,*) k, ll, tr1(ll), ' ion density'
 60     continue
        call frsplft(lrzmax,equilpsp(1),tr1(1),mfm1,eqpsi,zni(1,kk))
 70   continue
!      stop
!..................................................................
!     Average charge state..
!..................................................................
      kk=0
      do 90 k=1,ntotal
        if (k.eq.kelecm .or. k.eq.kelecg) go to 90
        kk=kk+1
!        if (kk.gt.nprim .and. k.lt.kimp_) then
        if (kk.gt.nprim+nimp+1 .and. k.lt.kimp_) then
          kk=kk-1
          go to 90
        endif
        do 80 ll=1,lrzmax
          tr1(ll)=bnumb(k)
!         write(*,*) k, ll, tr1(ll), ' ion charge'
 80     continue
        call frsplft(lrzmax,equilpsp(1),tr1(1),mfm1,eqpsi,zzi(1,kk))
 90   continue
!      stop
!..................................................................
!     Compute zeff(lr_) on Freya mesh.
!..................................................................

      call frsplft(lrzmax,equilpsp(1),zeff(1),mfm1,eqpsi,zeffctv(1))
!
       if(nprim.eq.1 .and. nimp.eq.1) then
         write(*,*) 'copying Zeff into zzi...'
         do k=1,mfm1
           zzi(k,nion+2)=zeff(k)
!           write(*,*) 'zzi(k,nion+2) = ',zzi(k,nion+2)
         enddo
       endif
!
!..................................................................
!     Change eqpsi back to convention at call to this subroutine.
!..................................................................

      do 101 l=1,nconteqn
        eqpsi(l)=-eqpsi(l)
 101  continue

!..................................................................
!     For running NFREYA with CQL3D, we assume iborb=0; so
!     it is necessary to set only the endpoints of potsid
!..................................................................

      potsid(1)=-psimag
      potsid(mf)=-psilim
      !elong=0 ! YuP: not used here?
      codeid="twodee"
      rin=rmincon
      write(*,*)'frstup: rmincon,rmaxcon=',rmincon,rmaxcon
      rmax=rmaxcon
      zminn=zmincon
      zmaxx=zmaxcon
      zax=0.
      mi=nnr
      mj=nnz
      zshift1=zshift

      return
      end subroutine frstup
