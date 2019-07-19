module urfdampa_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfdampa(krf)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine computes the collisional or additional linear
!     damping, depending on iurfcoll and iurfl, for all
!     ray elements on all rays for UNIT power delpwr.
!     It only needs to be done once.
!..................................................................


      character*8 ifirst
      save ifirst,nray0
      data ifirst/"first"/
      data nray0 /1/


!-YuP      if (ifirst.ne."first") go to 110 ! urfpwrl and urfpwrc are defined
!-YuP      ifirst="notfirst"  ! only for the first called krf (usually krf=1)
!      if (nharms.gt.0 .and. krf.gt.1) go to 110

!..................................................................
!     Loop over rays
!..................................................................

      do 10 iray=nray0,nray(krf)

!..................................................................
!     Loop of ray elements
!..................................................................

        do 20 is=1,nrayelt(iray,krf)

!..................................................................
!     Determine the length of the ray element, slngth
!..................................................................

          slngth=0.0
          if (is.eq.1) then
            slngth=ws(1,iray,krf)
          elseif (is.lt.nrayelt(iray,krf)) then
            slngth=0.5*(ws(is+1,iray,krf)-ws(is-1,iray,krf))
          elseif (is.eq.nrayelt(iray,krf)) then
            slngth=0.5*(ws(is,iray,krf)-ws(is-1,iray,krf))
          endif

!.......................................................................
!     Damping due to collisional absorption of the wave:
!.......................................................................

          if (iurfcoll(krfn(krf)).eq."enabled")  then
            urfpwrc(is,iray,krf)=1.d0-exp(-salphac(is,iray,krf)*slngth)
          endif

!.......................................................................
!     damping due to additional linear absorption (e.g., ions)
!     passed from rayop file:
!.......................................................................

          if (iurfl(krfn(krf)).ne."disabled") then !"enabled" or "split"
            urfpwrl(is,iray,krf)=1.d0-exp(-salphal(is,iray,krf)*slngth)
          endif

 20     continue
 10   continue
      iray=1

 110  continue

      return
      end subroutine urfdampa

end module urfdampa_mod
