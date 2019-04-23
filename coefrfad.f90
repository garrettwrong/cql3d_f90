module coefrfad_mod

  !---BEGIN USE

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine coefrfad(k,xrf)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!...............................................................
!     This routine adds the RF contributions to the Fokker-Planck
!     coefficients for species k.
!...............................................................

      save

!..................................................................
!     rf module contributions first.
!..................................................................


      call bcast(dbb(1:iyjxp1,0),zero,iyjxp1)
      call bcast(dff(0:iyp1jx-1,1),zero,iyp1jx)

      if (cqlpmod.ne."enabled") then
      if (vlhmod.eq."enabled" .or. vlfmod.eq."enabled") then
!     vlhmod and vlfmod are only set up for nrf=1
      if (nrf.gt.1) then
         write(*,*)'coefrfad: STOP, Problem with nrf'
         STOP
      endif
      do 10 kk=1,nrf
        if (k .eq. kk) then
          do krf=1,mrfn
          do 2 j=1,jx
            do 3 i=1,iy
              dbb(i,j)=dbb(i,j)+cqlb(i,j,indxlr_,krf)
              dff(i,j)=dff(i,j)+cqlf(i,j,indxlr_,krf)
              db(i,j)=  db(i,j)+cqlb(i,j,indxlr_,krf)
              dc(i,j)=  dc(i,j)+cqlc(i,j,indxlr_,krf)
              de(i,j)=  de(i,j)+cqle(i,j,indxlr_,krf)
              df(i,j)=  df(i,j)+cqlf(i,j,indxlr_,krf)
 3          continue
 2        continue
          enddo  !On krf
          xrf=1.
        endif
 10   continue
      endif   ! On vlhmod.or.vlfmod enabled

      elseif (cqlpmod.eq."enabled" .and. vlfmod.eq."enabled") then

      do 20 kk=1,nrf
!cc YuP-101220: Commenting next line, because wcqlb-wcqlf
!cc are defined for vlfmod="enabled" (but not for vlhmod="enabled")
!cc      if (vlhmod.eq."enabled" .or. vlfmod.eq."enabled") then
!     vlhmod and vlfmod are only set up for nrf=1
        if (nrf.gt.1) then
           write(*,*)'coefrfad: STOP, Problem with nrf'
           STOP
        endif
        if (k .eq. kk) then
          do krf=1,mrfn
          do 12 j=1,jx
            do 13 i=1,iy
              dbb(i,j)=dbb(i,j)+wcqlb(i,j,krf,l_)
              dff(i,j)=dff(i,j)+wcqlf(i,j,krf,l_)
              db(i,j)=db(i,j)+wcqlb(i,j,krf,l_)
              dc(i,j)=dc(i,j)+scatfrac*wcqlc(i,j,krf,l_)
              de(i,j)=de(i,j)+scatfrac*wcqle(i,j,krf,l_)
              df(i,j)=df(i,j)+scatfrac*wcqlf(i,j,krf,l_)
 13         continue
 12       continue
          enddo  !On krf
          xrf=1.
        endif
!cc      endif   ! On vlhmod,vlfmod enabled
 20   continue

      endif   ! On cqlpmod,vlfmod enabled

!..................................................................
!     urf module contribution next
!..................................................................

!**bh930729if(urfmod.eq."enabled".and.k.eq.kelecg) then
      if (urfmod.eq."enabled") then
        l=indxlr_
        do 30 krf=1,mrfn
!         Check if this rf mode is to applied to the current species k:
          if (nrfspecies(krfn(krf)) .eq. k) then
          do 22 j=1,jx
            do 23 i=1,iy
             delb0= urfb(i,j,l,krf)
             if(delb0.ne.0.d0)then
                !YuP[03/18/2015] urfe,urff are expr.through urfb,urfc
                !urfR23,urfR33 are expr.through urfb,urfc,urfR13
                !It can be shown that if the local delB=0, then also
                ! delC=0, delE=0, delF=0, and then all urfR** coeffs are zero.
                delc0=  urfc(i,j,l,krf)
                dele0=  delc0*sinn(i,l_)   ! urfe ! valid for ZOW&FOW
                delf0= (delc0*dele0)/delb0 ! urff ! valid for ZOW&FOW
                dbb(i,j)= dbb(i,j)+delb0
                dff(i,j)= dff(i,j)+delf0
                db(i,j)=  db(i,j) +delb0
                dc(i,j)=  dc(i,j) +scatfrac*delc0
                de(i,j)=  de(i,j) +scatfrac*dele0
                df(i,j)=  df(i,j) +scatfrac*delf0
             endif
 23        continue
 22     continue
          endif  ! on nrfspecies
 30     continue ! krf
      xrf=1.
      endif  ! on urfmod

!.......................................................................

!..................................................................
!     rdc module contribution next
!..................................................................

      if (rdcmod.ne."disabled") then
!         rdcbmax=0.d0
!         do j=1,jx
!            do i=1,iy
!               rdcbmax=max(rdcbmax,rdcb(i,j,indxlr_))
!            enddo
!         enddo
!         write(*,*)'coefrfad: Before rdcmod add, rdcbmax ,indxlr_=',
!     &        rdcbmax,indxlr_

!$$$         do j=1,jx
!$$$            do i=1,iy
!$$$              dbb(i,j)=dbb(i,j)+rdcb(i,j,indxlr_)
!$$$              dff(i,j)=dff(i,j)+rdcf(i,j,indxlr_)
!$$$              db(i,j)=db(i,j)+rdcb(i,j,indxlr_)
!$$$              dc(i,j)=dc(i,j)+scatfrac*rdcc(i,j,indxlr_)
!$$$              de(i,j)=de(i,j)+scatfrac*rdce(i,j,indxlr_)
!$$$              df(i,j)=df(i,j)+scatfrac*rdcf(i,j,indxlr_)
!$$$           enddo
!$$$        enddo
!$$$        xrf=1.

        l=indxlr_
        do 40 krf=1,nrdc
!         Check if this rf mode is to applied to the current species k:
          if (nrdcspecies(krf) .eq. k) then
          do 42 j=1,jx
            do 43 i=1,iy
             delb0= rdcb(i,j,l,krf)
             if(delb0.ne.0.d0)then
!$$$                !YuP[03/18/2015] urfe,urff are expr.through urfb,urfc
!$$$                !urfR23,urfR33 are expr.through urfb,urfc,urfR13
!$$$                !It can be shown that if the local delB=0, then also
!$$$                ! delC=0, delE=0, delF=0, and then all urfR** coeffs are zero.
                delc0=  rdcc(i,j,l,krf)
!$$$            BH180511: Could try following, as an alternative to direct calc.
!$$$                dele0=  delc0*sinn(i,l_)   ! urfe ! valid for ZOW&FOW
!$$$                delf0= (delc0*dele0)/delb0 ! urff ! valid for ZOW&FOW
                dele0=rdce(i,j,l,krf)
                delf0=rdcf(i,j,l,krf)
                dbb(i,j)= dbb(i,j)+delb0
                dff(i,j)= dff(i,j)+delf0
                db(i,j)=  db(i,j) +delb0
                dc(i,j)=  dc(i,j) +scatfrac*delc0
                de(i,j)=  de(i,j) +scatfrac*dele0
                df(i,j)=  df(i,j) +scatfrac*delf0
             endif
 43        continue
 42     continue
          endif  ! on nrfspecies
 40     continue ! krf
      xrf=1.

      endif  !On rdcmod

      return
      end
end module coefrfad_mod
