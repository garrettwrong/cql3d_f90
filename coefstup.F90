module coefstup_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use coefefad_mod, only : coefefad
  use coefegad_mod, only : coefegad
  use coeffpad_mod, only : coeffpad
  use coefload_mod, only : coefload
  use coefrfad_mod, only : coefrfad
  use coefsyad_mod, only : coefsyad
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine coefstup(k)
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine stores all the relevant coefficients into the
!     proper arrays for time advancement
!..................................................................


      call bcast(da,zero,iyjxp1)
      call bcast(db,zero,iyjxp1)
      call bcast(dc,zero,iyjxp1)
      call bcast(dd,zero,iyp1jx)
      call bcast(de,zero,iyp1jx)
      call bcast(df,zero,iyp1jx)
      call bcast(dbb,zero,iyjxp1)
      call bcast(dff,zero,iyp1jx)

!..................................................................
!     particle source
!..................................................................
      ifag=0
      do 10 is=1,nso
        if (n.ge.nonso(k,is) .and. n.lt.noffso(k,is)) then
          ifag=1
        endif
 10   continue
      if (ifag.eq.0) then
        call bcast(so,zero,iyjx2)
      else
        call dcopy(iyjx2,source(0:iy+1,0:jx+1,k,indxlr_),1,so,1)
      endif
      xrf=0.

!..................................................................
!     rf coefficients
!..................................................................

      if (n.eq.nonrf(k).and.(nrf.gt.0)) imprf=1
      if (urfmod.eq."enabled" .and. n.eq.nonrf(k)) imprf=1
      if (vlhplse.ne."disabled" .and. (vlhmod.eq."enabled".or. &
        vlfmod.eq."enabled")) then
        if (vlhplse.eq."tauee") then
          timespn=tauee(lr_)*(vlhpon+vlhpoff)
          timeon=tauee(lr_)*vlhpon
        else
          timespn=vlhpon+vlhpoff
          timeon=vlhpon
        endif
        iperiod=timet/timespn
        tim1=iperiod*timespn
        timedf=timet-tim1
        if (timedf.le.timeon) then
          call coefrfad(k,xrf)
        endif
      else if (n .ge. nonrf(k) .and. n .lt. noffrf(k)) then
        call coefrfad(k,xrf)
      endif

!..................................................................
!     parallel electric field
!..................................................................

      impelec=0
      if (elecfld(lr_) .ne. 0. .and.(nonel.le.n) .and.(noffel.gt.n))then
        call coefefad(k)
      endif
      if (n.eq.nonel .or. n.eq.noffel) impelec=1

!..................................................................
!     Add in the krook operator loss term..
!..................................................................

      call coefload(k)

!..................................................................
!     Add in phenomenological energy loss operator
!..................................................................

      call coefegad(k)

!..................................................................
!     Add in collisional coefficients
!..................................................................

      call coeffpad(k)

!..................................................................
!     add in synchrotron radiation losses
!..................................................................

      call coefsyad(k)

!.......................................................................
!     Add contribution of theta operator due to (mu*grad_parallel B) force
!.......................................................................

!%OS  if (setup0%cqlpmod.eq."enabled" .and. transp.eq."enabled" .and.
!%OS  +                                             n.ge.nontran+1)  then
!%OS  call wpcthta
!%OS  else
      call bcast(cthta,zero,iyjx)
!%OS  endif

      return
      end subroutine coefstup


end module coefstup_mod
